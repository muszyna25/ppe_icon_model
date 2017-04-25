!>
!! Computation of vertical tracer flux
!!
!! Vertical upwind fluxes are calculated at triangle centers on half-levels.
!! Possible options for vertical flux calculation include
!! - first order Godunov method (UP1)
!! - third order PPM method with or without CFL restriction
!!
!! Semi-monotone and monotone limiters are available for MUSCL and PPM
!!
!! These routines compute only the correct half level value of
!! 'c*weta'. The vertical divergence is computed in step_advection.
!!
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Jochen Foerstner, DWD (2008-05-15)
!! Modification by Daniel Reinert, DWD (2009-08-06)
!! - code restructured
!! Modification by Daniel Reinert, DWD (2009-08-12)
!! - included piecewise parabolic method (PPM)
!! Modification by Daniel Reinert, DWD (2009-08-12)
!! - recoding of MUSCL in order to account for time and space
!!   dependent surface pressure and layer thickness
!! Modification by Daniel Reinert, DWD (2010-01-22)
!! - modified MUSCL scheme which handles CFL>1 (see Lin and Rood (1996))
!!   added optional semi-monotone limiter.
!! Modification by Daniel Reinert, DWD (2010-01-09)
!! - transferred vertical flux calculation (UP, MUSCL, PPM) to this
!!   new module
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - moved slope limiter to new module mo_advection_limiter
!! Modification by Daniel Reinert, DWD (2013-05-07)
!! - removed unused second order MUSCL scheme
!! Modification by Daniel Reinert, DWD (2016-03 ?)
!! - refactoring in upwind_vflux_ppm_cfl
!! Modification by Will Sawyer, CSCS (2016-07-15)
!! - added OpenACC support
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
MODULE mo_advection_vflux

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, SUCCESS, min_rlcell_int,   &
    &                               iup_v, ippm_v, ippm_vcfl, islopel_vsm,      &
    &                               islopel_vm, ifluxl_vpd, ino_flx, izero_grad,&
    &                               iparent_flx
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_dynamics_config,     ONLY: iequations 
  USE mo_run_config,          ONLY: msg_level, lvert_nest, timers_level, iqtke
  USE mo_advection_config,    ONLY: advection_config, lcompute, lcleanup, t_trList 
  USE mo_advection_utils,     ONLY: laxfr_upflux_v
  USE mo_advection_limiter,   ONLY: v_ppm_slimiter_mo, v_ppm_slimiter_sm,     &
   &                                vflx_limiter_pd, vflx_limiter_pd_ha
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: global_max
  USE mo_mpi,                 ONLY: process_mpi_stdio_id, my_process_is_stdio, get_my_mpi_work_id, &
                                    get_glob_proc0, comm_lev
#ifdef _OPENACC
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, check_patch_array
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif
  USE mo_timer,               ONLY: timer_adv_vert, timer_start, timer_stop


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: vert_upwind_flux
  PUBLIC :: upwind_vflux_up
  PUBLIC :: upwind_vflux_ppm
  PUBLIC :: upwind_vflux_ppm_cfl

#if defined( _OPENACC )
#define ACC_DEBUG NOACC
#if defined(__ADVECTION_VFLUX_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
#endif

  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !
  !

  !>
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !!
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !! using either
  !! - the first order Godunov method (UP1)
  !! - the third order PPM method
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-02-10)
  !! Modification by Daniel Reinert (2010-11-05)
  !! - tracer loop moved from step_advection to vert_upwind_flux
  !! Modification by Daniel Reinert (2010-11-08)
  !! - lcompute and lcleanup are precomputed in setup_transport.
  !!   This was necessary to allow for different flux-methods
  !!   for different tracers.
  !! Modification by Daniel Reinert, DWD (2011-02-15)
  !! - new field providing the upper margin tracer flux (required)
  !! Modification by Daniel Reinert, DWD (2013-05-07)
  !! - removed unused second order MUSCL scheme
  !!
  !
  ! !LITERATURE
  ! PPM  : Colella and Woodward (1984), JCP, 54, 174-201
  !        Carpenter et al. (1989), MWR, 118, 586-612
  !        Lin and Rood (1996), MWR, 124, 2046-2070 (see also for CFL-
  !                                                  independent versions)
  !
  SUBROUTINE vert_upwind_flux( p_patch, p_cc, p_mflx_contra_v, p_w_contra,    &
    &                      p_dtime, p_cellhgt_mc_now,                         &
    &                      p_cellmass_now, p_ivadv_tracer, p_itype_vlimit,    &
    &                      p_iubc_adv, p_iadv_slev, lprint_cfl, p_upflux,     &
    &                      opt_topflx_tra, opt_q_int, opt_rlstart, opt_rlend  )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is 
      &  p_patch                             !< performed

    REAL(wp), INTENT(IN) ::  &      !< advected cell centered variable
      &  p_cc(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT) ::  &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)     !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &      !< contravariant vertical velocity
      &  p_w_contra(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    REAL(wp), INTENT(IN) ::  &      !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)    !< time step n (either \Delta p or \Delta z)
                                    !< dim: (nproma,nlev,nblks_c)

    REAL(wp), TARGET, INTENT(IN)::& !< NH: density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)      !< at time step n [kg/m**2]
                                    !< HA: pressure thickness for full levels 
                                    !< at time step n [Pa]
                                    !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN) ::   &      !< parameter to select numerical
      &  p_ivadv_tracer(:)          !< scheme for vertical transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< parameter to select the limiter
      &  p_itype_vlimit(:)          !< for vertical transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< vertical start level for transport
      &  p_iadv_slev(:)             !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< selects upper boundary condition
      &  p_iubc_adv

    LOGICAL, INTENT(IN) ::   &      !< determines if vertical CFL number shall be printed
      &  lprint_cfl                 !< in routine upwind_vflux_ppm_cfl

    REAL(wp), INTENT(INOUT) :: &    !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:,:)          !< dim: (nproma,nlevp1,nblks_c,ntracer)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:,:)          !< NH: [kg/m**2/s]
                                        !< HA: [Pa/s]
                                        !< dim: (nproma,nblks_c,ntracer)

    REAL(wp), INTENT(OUT), OPTIONAL :: & !< tracer value at upper boundary of child nest 
      &  opt_q_int(:,:,:)               !< NH: [kg/kg]
                                        !< HA: [kg/kg]
                                        !< dim: (nproma,nblks_c,ntracer)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    TYPE(t_trList), POINTER ::  &      !< pointer to tracer sublist
      trAdvect

    INTEGER :: jt, nt                  !< tracer index and loop index
    INTEGER :: jg                      !< patch ID
    INTEGER :: jc, jb                  !< cell and block loop index
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart_c, i_rlend_c, i_nchdom
    INTEGER :: iadv_min_slev           !< scheme specific minimum slev

    REAL(wp) :: z_mflx_contra_v(nproma) !< auxiliary variable for computing vertical nest interface quantities
    !-----------------------------------------------------------------------

    IF (timers_level > 2) CALL timer_start(timer_adv_vert)

    ! get patch ID
    jg = p_patch%id

    trAdvect => advection_config(jg)%trAdvect

    !
    ! Loop over different tracers
    !
    ! Note (DR): Since we define p_ivadv_tracer as a 1D array of dimension
    ! ntracer, we can select different flux calculation methods for
    ! different tracers. We may also decide not to advect special
    ! tracers. Furthermore, options regarding the desired limiter can be passed 
    ! via the argument list. The same is true for precomputed lcompute/lcleanup
    ! values.
    IF (PRESENT(opt_topflx_tra)) THEN

      DO nt = 1, trAdvect%len

        jt = trAdvect%list(nt)

        IF (.NOT. PRESENT(opt_rlend) .OR. (jt == iqtke .AND. advection_config(jg)%iadv_tke == 1)) THEN
          i_rlend_c = min_rlcell_int
        ELSE
          i_rlend_c = opt_rlend
        ENDIF

        ! Select desired flux calculation method
        SELECT  CASE( p_ivadv_tracer(jt) )

        CASE( iup_v )
          ! CALL first order upwind
          CALL upwind_vflux_up( p_patch, p_cc(:,:,:,jt), p_iubc_adv,   &! in
            &                   p_mflx_contra_v, p_upflux(:,:,:,jt),   &! in,out
            &                   opt_topflx_tra=opt_topflx_tra(:,:,jt), &! in
            &                   opt_slev=p_iadv_slev(jt),              &! in
            &                   opt_rlstart=opt_rlstart,               &! in
            &                   opt_rlend=i_rlend_c                    )! in


        CASE( ippm_vcfl )

          iadv_min_slev = advection_config(jg)%ppm_v%iadv_min_slev

          ! CALL third order PPM (unrestricted timestep-version) (i.e. CFL>1)
          CALL upwind_vflux_ppm_cfl( p_patch, p_cc(:,:,:,jt), p_iubc_adv,    &! in
            &                  p_mflx_contra_v, p_dtime, lcompute%ppm_v(jt), &! in
            &                  lcleanup%ppm_v(jt), p_itype_vlimit(jt),       &! in
            &                  p_cellhgt_mc_now, p_cellmass_now, lprint_cfl, &! in
            &                  p_upflux(:,:,:,jt),                           &! out
            &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),        &! in
            &                  opt_slev=p_iadv_slev(jt),                     &! in
            &                  opt_ti_slev=iadv_min_slev,                    &! in
            &                  opt_rlstart=opt_rlstart,                      &! in
            &                  opt_rlend=i_rlend_c                           )! in
        CASE( ippm_v )
          ! CALL third order PPM
          CALL upwind_vflux_ppm( p_patch, p_cc(:,:,:,jt), p_iubc_adv,  &! in
            &                  p_mflx_contra_v, p_w_contra, p_dtime,   &! in
            &                  p_itype_vlimit(jt), p_cellhgt_mc_now,   &! in
            &                  p_upflux(:,:,:,jt),                     &! out
            &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),  &! in
            &                  opt_slev=p_iadv_slev(jt),               &! in
            &                  opt_rlstart=opt_rlstart,                &! in
            &                  opt_rlend=i_rlend_c                     )! in
        END SELECT
      END DO  ! Tracer loop


      !
      ! get face value "q_int" at vertical boundary of child nest
      !
      ! determine if upper boundary values are needed
      IF (lvert_nest .AND. (p_patch%nshift_child > 0)) THEN 

        ! refinement control start/end level for cells
        i_rlstart_c = 1
        i_rlend_c   = min_rlcell_int

        ! number of child domains
        i_nchdom = MAX(1,p_patch%n_childdom)

        i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

#ifndef _OPENACC
!$OMP PARALLEL DO PRIVATE(jb,jt,jc,nt,i_startidx,i_endidx,z_mflx_contra_v) ICON_OMP_DEFAULT_SCHEDULE
#endif
        DO jb = i_startblk, i_endblk
          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
            &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

          ! Be sure to avoid division by zero
          DO jc = i_startidx, i_endidx
            z_mflx_contra_v(jc) = SIGN( MAX(ABS(p_mflx_contra_v(jc,p_patch%nshift_child,jb)),dbl_eps), &
              &                                 p_mflx_contra_v(jc,p_patch%nshift_child,jb) )
          ENDDO

          DO nt = 1, trAdvect%len
            jt = trAdvect%list(nt)
            DO jc = i_startidx, i_endidx
              opt_q_int(jc,jb,jt) = p_upflux(jc,p_patch%nshift_child,jb,jt) / z_mflx_contra_v(jc)
            ENDDO
          ENDDO
        ENDDO
#ifndef _OPENACC
!$OMP END PARALLEL DO
#endif

      ENDIF


    ELSE ! opt_topflx_tra not present (i.e. for hydrostatic model)

      DO nt = 1, trAdvect%len

        jt = trAdvect%list(nt)

        ! Select desired flux calculation method
        SELECT  CASE( p_ivadv_tracer(jt) )

        CASE( iup_v )
          ! CALL first order upwind
          CALL upwind_vflux_up( p_patch, p_cc(:,:,:,jt), p_iubc_adv,   &! in
            &                   p_mflx_contra_v, p_upflux(:,:,:,jt),   &! in,out
            &                   opt_slev=p_iadv_slev(jt),              &! in
            &                   opt_rlstart=opt_rlstart,               &! in
            &                   opt_rlend=opt_rlend                    )! in


        CASE( ippm_vcfl )

          iadv_min_slev = advection_config(jg)%ppm_v%iadv_min_slev

          ! CALL third order PPM which handles long time steps (i.e. CFL>1)
          CALL upwind_vflux_ppm_cfl( p_patch, p_cc(:,:,:,jt), p_iubc_adv,    &! in
            &                  p_mflx_contra_v, p_dtime, lcompute%ppm_v(jt), &! in
            &                  lcleanup%ppm_v(jt), p_itype_vlimit(jt),       &! in
            &                  p_cellhgt_mc_now, p_cellmass_now, lprint_cfl, &! in
            &                  p_upflux(:,:,:,jt),                           &! out
            &                  opt_slev=p_iadv_slev(jt),                     &! in
            &                  opt_ti_slev=iadv_min_slev,                    &! in
            &                  opt_rlstart=opt_rlstart,                      &! in
            &                  opt_rlend=opt_rlend                           )! in

        CASE( ippm_v )
          ! CALL third order PPM
          CALL upwind_vflux_ppm( p_patch, p_cc(:,:,:,jt), p_iubc_adv,  &! in
            &                  p_mflx_contra_v, p_w_contra, p_dtime,   &! in
            &                  p_itype_vlimit(jt), p_cellhgt_mc_now,   &! in
            &                  p_upflux(:,:,:,jt),                     &! out
            &                  opt_slev=p_iadv_slev(jt),               &! in
            &                  opt_rlstart=opt_rlstart,                &! in
            &                  opt_rlend=opt_rlend                     )! in

        END SELECT
      END DO  ! Tracer loop

    END IF ! PRESENT(opt_topflx_tra)

    IF (timers_level > 2) CALL timer_stop(timer_adv_vert)

  END SUBROUTINE vert_upwind_flux




  !-------------------------------------------------------------------------
  !>
  !! The first order Godunov method
  !!
  !! Calculation of time averaged vertical tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - transferred to separate subroutine
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized to height based vertical coordinate systems
  !!
  SUBROUTINE upwind_vflux_up( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v, &
    &                         p_upflux, opt_topflx_tra, opt_slev,         &
    &                         opt_rlstart, opt_rlend )

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_advection_vflux: upwind_vflux_up'

    TYPE(t_patch), TARGET, INTENT(IN) ::  & !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::   &   !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(IN) ::   &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(INOUT) ::  & !< vertical tracer flux at half levels
      &  p_upflux(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    REAL(wp) ::  &                             !< necessary, to make this routine
     &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 
                                       
    INTEGER  :: slev                   !< vertical start level
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: jc, jk, jb             !< index of cell, vertical level and block
    INTEGER  :: jg                     !< patch ID
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    !-------------------------------------------------------------------------

#ifdef _OPENACC
!$ACC DATA CREATE( zparent_topflx ), PCOPYIN( p_cc, p_mflx_contra_v ), PCOPYOUT( p_upflux ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cc, p_mflx_contra_v ), IF( i_am_accel_node .AND. acc_on )
#endif

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_topflx_tra) ) THEN
!$ACC KERNELS PRESENT( opt_topflx_tra, zparent_topflx ), IF ( i_am_accel_node .AND. acc_on )
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
!$ACC END KERNELS
    ELSE
!$ACC KERNELS PRESENT( zparent_topflx ), IF ( i_am_accel_node .AND. acc_on )
      zparent_topflx(:,:) = 0._wp
!$ACC END KERNELS
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! get patch ID
    jg = p_patch%id

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( p_patch, advection_config, p_cc, p_mflx_contra_v, p_upflux ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

!$ACC LOOP VECTOR COLLAPSE(2)
      DO jk = slev+1, nlev
        DO jc = i_startidx, i_endidx
          ! calculate vertical tracer flux
          p_upflux(jc,jk,jb) =                                  &
            &  laxfr_upflux_v( p_mflx_contra_v(jc,jk,jb),       &
            &                p_cc(jc,jk-1,jb), p_cc(jc,jk,jb),  &
            &                advection_config(jg)%coeff_grid )

        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels


      !
      ! set upper and lower boundary condition
      !
!
! With OpenACC this sometimes causes the compiler to crash, apparently due to the array syntax of one argument.  
!
      CALL set_bc_vadv(nproma, p_upflux(:,slev+1,jb),    &! in
        &              p_mflx_contra_v(:,slev+1,jb),     &! in
        &              p_mflx_contra_v(:,slev  ,jb),     &! in
        &              p_iubc_adv, i_startidx, i_endidx, &! in
        &              zparent_topflx(:,jb),             &! in
        &              p_upflux(:,slev,jb),              &! out
        &              p_upflux(:,nlevp1,jb), .TRUE.)     ! out

    ENDDO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
!ACC_DEBUG UPDATE HOST( p_upflux ), IF( i_am_accel_node .AND. acc_on )
!$ACC END DATA
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

  END SUBROUTINE upwind_vflux_up



  !-------------------------------------------------------------------------
  !>
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
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin et al (1994), MWR, 122, 1575-1593 (slope limiter)
  ! - Lin and Rood (1996), MWR, 124, 2046-2070
  !
  SUBROUTINE upwind_vflux_ppm( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v,  &
    &                      p_w_contra, p_dtime, p_itype_vlimit,             &
    &                      p_cellhgt_mc_now, p_upflux,                      &
    &                      opt_lout_edge, opt_topflx_tra, opt_slev,         &
    &                      opt_rlstart, opt_rlend )

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_advection_vflux: upwind_vflux_ppm'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(IN) ::  &    !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< contravariant vertical velocity
      &  p_w_contra(:,:,:)        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime  !< time step


    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) :: &  !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    INTEGER, INTENT(IN) :: p_itype_vlimit  !< parameter to select the limiter for
                                           !< vertical transport

    LOGICAL, INTENT(IN), OPTIONAL :: & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                   !< or the flux across the edge 
                                       !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL  :: l_out_edgeval          !< corresponding local variable; default 
                                       !< .FALSE. i.e. output flux across the edge

    REAL(wp) :: &                      !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) :: &                      !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                      !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) :: &                      !< linear extrapolation value 1
      &  z_lext_1
 
    REAL(wp) :: &                      !< linear extrapolation value 2
      &  z_lext_2
 
    REAL(wp) :: &                      !< CFL number (weta>0, w<0)
      &  z_cfl_m(nproma,p_patch%nlevp1,p_patch%nblks_c) !< (nproma,nlevp1,p_patch%nblks_c)
                                                         
    REAL(wp) :: &                      !< CFL number (weta<0, w>0)
      &  z_cfl_p(nproma,p_patch%nlevp1,p_patch%nblks_c) !< (nproma,nlevp1,p_patch%nblks_c)
                                                       
    REAL(wp) :: &                      !< monotonized slope
      &  z_slope(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) ::  &                             !< necessary, to make this routine
     &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core

    REAL(wp) :: p_cc_min, p_cc_max       !< 3-point max/min values

    REAL(wp) :: z_delta_m, z_delta_p   !< difference between lower and upper face value
                                       !< for weta >0 and weta <0
    REAL(wp) :: z_a11, z_a12           !< 1/6 * a6,i (see Colella and Woodward (1984))
    REAL(wp) :: z_weta_dt              !< weta times p_dtime

    INTEGER  :: slev, slevp1           !< vertical start level and start level +1
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: ikm1, ikp1, ikp2       !< vertical level minus and plus one, plus two
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom
    INTEGER  :: jc, jk, jb              !< index of cell, vertical level and block
    INTEGER  :: jg                      !< patch ID
    REAL(wp) :: coeff_grid              !< parameter which is used to make the vertical 
                                        !< advection scheme applicable to a height      
                                        !< based coordinate system (coeff_grid=-1)

    !-----------------------------------------------------------------------

#ifdef _OPENACC
!$ACC DATA CREATE( z_face, z_face_up, z_face_low, z_cfl_m, z_cfl_p, z_slope ), &
!$ACC      PCOPYIN( p_cc, p_mflx_contra_v, p_w_contra, p_cellhgt_mc_now ), &
!$ACC      PCOPYOUT( p_upflux ), IF( i_am_accel_node .AND. acc_on )
!ACC_DEBUG UPDATE DEVICE( p_cc, p_mflx_contra_v, p_w_contra, p_cellhgt_mc_now ), IF( i_am_accel_node .AND. acc_on )
#endif

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
      slevp1= opt_slev + 1
    ELSE
      slev  = 1
      slevp1= 2
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_topflx_tra) ) THEN
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
    ELSE
      zparent_topflx(:,:) = 0._wp
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! get patch ID
    jg = p_patch%id

    ! save some paperwork
    coeff_grid = advection_config(jg)%coeff_grid


    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    !
    ! advection is done with an upwind scheme and a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    !
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    !
    ! 1. Calculate Courant number for weta>0 (w<0) and weta<0 (w>0)
    !
#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( p_patch, advection_config, p_cc, p_mflx_contra_v, p_cellhgt_mc_now ), &
!$ACC PRESENT( z_cfl_p, z_cfl_m, z_slope, z_face, p_upflux ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikm1,z_weta_dt,ikp1, &
!$OMP            p_cc_min,p_cc_max,ikp2) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! Courant number at top
      z_cfl_p(i_startidx:i_endidx,slev,jb) = 0._wp
      z_cfl_m(i_startidx:i_endidx,slev,jb) = 0._wp
      ! Courant number at bottom
      z_cfl_p(i_startidx:i_endidx,nlevp1,jb) = 0._wp
      z_cfl_m(i_startidx:i_endidx,nlevp1,jb) = 0._wp

!$ACC LOOP VECTOR, COLLAPSE(2)
      DO jk = slevp1, nlev
        DO jc = i_startidx, i_endidx

          ! WS 2015-05-07: moved to inner loop to ensure loops can be collapsed
          ! index of top half level
          ikm1 = jk - 1

          ! Calculate local Courant number at half levels
          ! z_cfl_m for weta >0 (w <0)
          ! z_cfl_p for weta <0 (w >0)
          z_weta_dt = ABS(p_w_contra(jc,jk,jb)) * p_dtime

          z_cfl_m(jc,jk,jb) = z_weta_dt / p_cellhgt_mc_now(jc,ikm1,jb)

          z_cfl_p(jc,jk,jb) = z_weta_dt / p_cellhgt_mc_now(jc,jk,jb)

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels


      !
      ! 2. Calculate monotonized slope
      !
      z_slope(i_startidx:i_endidx,slev,jb) = 0._wp

!$ACC LOOP VECTOR, COLLAPSE(2)
      DO jk = slevp1, nlev
        DO jc = i_startidx, i_endidx

          ! WS 2015-05-07: moved to inner loop to ensure loops can be collapsed
          ! index of top half level
          ikm1    = jk - 1
          ! index of bottom half level
          ikp1    = MIN( jk+1, nlev )


          z_slope(jc,jk,jb) = ( p_cellhgt_mc_now(jc,jk,jb)                             &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)            &
              &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                       &
              &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
              &  / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                                   &
              &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb))   &
              &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))           &
              &  * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb)) )

          ! equivalent formulation of Colella and Woodward (1984) slope limiter 
          ! following Lin et al (1994).
          p_cc_min = MIN(p_cc(jc,ikm1,jb),p_cc(jc,jk,jb),p_cc(jc,ikp1,jb))
          p_cc_max = MAX(p_cc(jc,ikm1,jb),p_cc(jc,jk,jb),p_cc(jc,ikp1,jb))
          z_slope(jc,jk,jb) = SIGN(                                            &
            &  MIN( ABS(z_slope(jc,jk,jb)), 2._wp*(p_cc(jc,jk,jb)-p_cc_min),   &
            &                               2._wp*(p_cc_max-p_cc(jc,jk,jb)) ), &
            &    z_slope(jc,jk,jb))

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
!$ACC LOOP VECTOR
      DO jc = i_startidx, i_endidx

        z_face(jc,slevp1,jb) = p_cc(jc,slev,jb)*(1._wp - (p_cellhgt_mc_now(jc,slev,jb)&
          &       / p_cellhgt_mc_now(jc,slevp1,jb))) + (p_cellhgt_mc_now(jc,slev,jb)  &
          &       /(p_cellhgt_mc_now(jc,slev,jb) + p_cellhgt_mc_now(jc,slevp1,jb)))   &
          &       * ((p_cellhgt_mc_now(jc,slev,jb) / p_cellhgt_mc_now(jc,slevp1,jb))  &
          &       * p_cc(jc,slev,jb) + p_cc(jc,slevp1,jb))

        z_face(jc,nlev,jb) = p_cc(jc,nlev-1,jb)*( 1._wp                               &
          &       - (p_cellhgt_mc_now(jc,nlev-1,jb) / p_cellhgt_mc_now(jc,nlev,jb)))  &
          &       + (p_cellhgt_mc_now(jc,nlev-1,jb)/(p_cellhgt_mc_now(jc,nlev-1,jb)   &
          &       + p_cellhgt_mc_now(jc,nlev,jb))) * ((p_cellhgt_mc_now(jc,nlev-1,jb) &
          &       / p_cellhgt_mc_now(jc,nlev,jb)) * p_cc(jc,nlev-1,jb)                &
          &       + p_cc(jc,nlev,jb))

        z_face(jc,slev,jb)   = p_cc(jc,slev,jb)
        z_face(jc,nlevp1,jb) = p_cc(jc,nlev,jb)

      ENDDO


!$ACC LOOP VECTOR, COLLAPSE(2)
      DO jk = slevp1, nlev-2

        DO jc = i_startidx, i_endidx

          ! index of top half level
          ikm1 = jk - 1
          ! index of bottom half level
          ikp1 = jk + 1
          ikp2 = jk + 2

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

        END DO ! end loop over cells

      END DO ! end loop over vertical levels

    END DO
#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO
#endif



    !
    ! 4. Limitation of first guess parabola (which is based on z_face)
    ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
    ! the limitation procedure.
    ! Therefore 2 additional fields z_face_up and z_face_low are
    ! introduced.
    !
#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( p_patch, p_cc, z_face, z_face_up, z_face_low, z_slope ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP DO PRIVATE(jk,ikp1,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      IF (p_itype_vlimit == islopel_vsm) THEN
        ! semi-monotonic (sm) limiter
! TODO (WS): if the compiler cannot inline this, make it an ACC ROUTINE
        CALL v_ppm_slimiter_sm( p_cc(:,:,jb), z_face(:,:,jb),           & !in
          &                   z_face_up(:,:,jb), z_face_low(:,:,jb),    & !inout
          &                   i_startidx, i_endidx, slev, nlev          ) !in
      ELSE IF (p_itype_vlimit == islopel_vm) THEN
        ! monotonic (mo) limiter
! TODO (WS): if the compiler cannot inline this, make it an ACC ROUTINE
        CALL v_ppm_slimiter_mo( p_cc(:,:,jb), z_face(:,:,jb), z_slope(:,:,jb), & !in
          &                   z_face_up(:,:,jb), z_face_low(:,:,jb),           & !inout
          &                   i_startidx, i_endidx, slev, nlev                 ) !in
      ELSE
       ! simply copy face values to 'face_up' and 'face_low' arrays
!$ACC LOOP VECTOR
        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          z_face_up(i_startidx:i_endidx,jk,jb)  = z_face(i_startidx:i_endidx,jk,jb)
          z_face_low(i_startidx:i_endidx,jk,jb) = z_face(i_startidx:i_endidx,ikp1,jb)
        ENDDO
      ENDIF
    ENDDO
#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP ENDDO
#endif


    IF ( l_out_edgeval ) THEN

      ! 5a. Compute edge value of advected quantity by applying a 
      !     piecewise parabolic reconstruction

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( p_patch, p_cc, p_mflx_contra_v, z_face_up, z_face_low, z_cfl_m ), &
!$ACC PRESENT( p_upflux ),                                                       &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,ikm1,z_delta_m, &
!$OMP            z_delta_p,z_a11,z_a12) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        ! for the time being, when computing the edge value, we set 
        ! the top and bottom values to 0
        p_upflux(i_startidx:i_endidx,  slev,jb) = 0._wp
        p_upflux(i_startidx:i_endidx,nlevp1,jb) = 0._wp

!$ACC LOOP VECTOR, COLLAPSE(2)
        DO jk = slevp1, nlev
          DO jc = i_startidx, i_endidx

            ! WS 2015-05-07: moved to inner loop to ensure loops can be collapsed
            ! index of top half level
            ikm1 = jk -1

            ! linear extrapolated values
            ! for the height based coordinate system multiplication by coeff_grid
            ! is not necessary due to compensating (-) signs.
            ! first (of cell above) (case of w < 0; weta > 0)
            z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
            z_a11     = p_cc(jc,ikm1,jb)                                  &
              &       - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))

            z_lext_1 = p_cc(jc,ikm1,jb)                                   &
              &  + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
              &  - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))


            ! second (of cell below) (case of w > 0; weta < 0)
            z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
            z_a12     = p_cc(jc,jk,jb)                                    &
              &       - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))

            z_lext_2 = p_cc(jc,jk,jb)                                     &
              &  - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
              &  - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))

            !
            ! calculate 'edge value' of advected quantity
            !
            p_upflux(jc,jk,jb) =                                         &
              &  laxfr_upflux_v( SIGN(1._wp, p_mflx_contra_v(jc,jk,jb)), &
              &                  z_lext_1, z_lext_2,                     &
              &                  -1.0_wp )  ! for test purposes only in nh model

            ! sign of the edge value
            p_upflux(jc,jk,jb) =  p_upflux(jc,jk,jb)*SIGN(1._wp, p_mflx_contra_v(jc,jk,jb))

          END DO ! end loop over cells

        ENDDO ! end loop over vertical levels

      ENDDO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO NOWAIT
#endif

    ELSE
    !
    ! 5b. extrapolation using piecewise parabolic approx. of the transported
    ! quantity to the edge and finally, calculation of the upwind fluxes
    !

#ifdef _OPENACC
!$ACC PARALLEL &
!$ACC PRESENT( p_patch, p_cc, p_mflx_contra_v, z_face_up, z_face_low, z_cfl_m, z_cfl_p ), &
!$ACC PRESENT( p_upflux ),                                                                &
!$ACC IF( i_am_accel_node .AND. acc_on )

!$ACC LOOP GANG
#else
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,z_lext_1,z_lext_2,ikm1,z_delta_m, &
!$OMP            z_delta_p,z_a11,z_a12) ICON_OMP_DEFAULT_SCHEDULE
#endif
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

!$ACC LOOP VECTOR, COLLAPSE(2)
        DO jk = slevp1, nlev
          DO jc = i_startidx, i_endidx

            ! WS 2015-05-07: moved to inner loop to ensure loops can be collapsed
            ! index of top half level
            ikm1 = jk -1

            ! linear extrapolated values
            ! for the height based coordinate system multiplication by coeff_grid
            ! is not necessary due to compensating (-) signs.
            ! first (of cell above) (case of w < 0; weta > 0)
            z_delta_m = z_face_low(jc,ikm1,jb) - z_face_up(jc,ikm1,jb)
            z_a11     = p_cc(jc,ikm1,jb)                                  &
              &       - 0.5_wp * (z_face_low(jc,ikm1,jb) + z_face_up(jc,ikm1,jb))

            z_lext_1 = p_cc(jc,ikm1,jb)                                   &
              &  + (0.5_wp * z_delta_m * (1._wp - z_cfl_m(jc,jk,jb)))     &
              &  - z_a11*(1._wp - 3._wp*z_cfl_m(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_m(jc,jk,jb)*z_cfl_m(jc,jk,jb))


            ! second (of cell below) (case of w > 0; weta < 0)
            z_delta_p = z_face_low(jc,jk,jb) - z_face_up(jc,jk,jb)
            z_a12     = p_cc(jc,jk,jb)                                    &
              &       - 0.5_wp * (z_face_low(jc,jk,jb) + z_face_up(jc,jk,jb))

            z_lext_2 = p_cc(jc,jk,jb)                                     &
              &  - (0.5_wp * z_delta_p * (1._wp - z_cfl_p(jc,jk,jb)))     &
              &  - z_a12*(1._wp - 3._wp*z_cfl_p(jc,jk,jb)                 &
              &  + 2._wp*z_cfl_p(jc,jk,jb)*z_cfl_p(jc,jk,jb))

            !
            ! calculate vertical tracer flux
            !
            p_upflux(jc,jk,jb) =                                  &
              &  laxfr_upflux_v( p_mflx_contra_v(jc,jk,jb),       &
              &                z_lext_1, z_lext_2,                &
              &                coeff_grid )

          END DO ! end loop over cells

        ENDDO ! end loop over vertical levels


        !
        ! set upper and lower boundary condition
        !
        CALL set_bc_vadv(nproma, p_upflux(:,slev+1,jb),    &! in
          &              p_mflx_contra_v(:,slev+1,jb),     &! in
          &              p_mflx_contra_v(:,slev  ,jb),     &! in
          &              p_iubc_adv, i_startidx, i_endidx, &! in
          &              zparent_topflx(:,jb),             &! in
          &              p_upflux(:,slev,jb),              &! out
          &              p_upflux(:,nlevp1,jb), .TRUE.)     ! out

      ENDDO ! end loop over blocks
#ifdef _OPENACC
!$ACC END PARALLEL
#else
!$OMP END DO NOWAIT
#endif


    ENDIF
#ifndef _OPENACC
!$OMP END PARALLEL
#endif

    !
    ! 6. If desired, apply a flux limiter to limit computed fluxes.
    !    These flux limiters are based on work by Zalesak (1979)
    !
    IF (iequations /= 3) THEN
      IF (p_itype_vlimit == ifluxl_vpd) THEN
        ! positive-definite (pd) limiter
        print *, "Calling VFLX_LIMITER_PD_HA"
        CALL vflx_limiter_pd_ha( p_patch, p_dtime, p_cc, p_upflux,    & !in,inout
          &                 opt_rlstart=i_rlstart, opt_rlend=i_rlend, & !in
          &                 opt_slev=slev                             ) !in
      ENDIF
    ELSE
      IF (p_itype_vlimit == ifluxl_vpd) THEN
        ! positive-definite (pd) limiter
        print *, "Calling VFLX_LIMITER_PD"
        CALL vflx_limiter_pd( p_patch, p_dtime, p_cc, p_upflux,       & !in,inout
          &                opt_rlstart=i_rlstart, opt_rlend=i_rlend,  & !in
          &                opt_slev=slev                              ) !in
      ENDIF
    ENDIF

#ifdef _OPENACC
!ACC_DEBUG UPDATE HOST(p_upflux), IF (i_am_accel_node .AND. acc_on)
!$ACC END DATA
#endif

  END SUBROUTINE upwind_vflux_ppm



  !-------------------------------------------------------------------------
  !>
  !! The third order PPM scheme for large time steps (CFL>1)
  !!
  !! Calculation of time averaged vertical tracer fluxes or tracer edge 
  !! values using the third order PPM scheme. This scheme can handle 
  !! large time steps (i.e. CFL>1)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-01-14)
  !!
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin et al (1994), MWR, 122, 1575-1593 (slope limiter)
  ! - Lin and Rood (1996), MWR, 124, 2046-2070 (CFL-independent version)
  !
  SUBROUTINE upwind_vflux_ppm_cfl( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v, &
    &                      p_dtime,  ld_compute, ld_cleanup, p_itype_vlimit,   &
    &                      p_cellhgt_mc_now, p_cellmass_now, lprint_cfl,       &
    &                      p_upflux, opt_lout_edge, opt_topflx_tra, opt_slev,  &
    &                      opt_ti_slev, opt_rlstart, opt_rlend, opt_elev )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux:upwind_vflux_ppm_cfl'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(INOUT) ::  & !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< time step
      &  p_dtime

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN)  ::  &    !< parameter to select the limiter for
      &  p_itype_vlimit           !< vertical transport

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< NH: density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)    !< at time step n [kg/m**2]
                                  !< HA: pressure thickness for full levels 
                                  !< at time step n [Pa]
                                  !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &    !< determines if vertical CFL number shall be written out
      &  lprint_cfl

    REAL(wp), INTENT(INOUT) :: &  !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    LOGICAL, INTENT(IN), OPTIONAL ::  & !< optional: output edge value (.TRUE.),
     & opt_lout_edge                    !< or the flux across the edge 
                                        !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end   level (for sedimentation)
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)

    LOGICAL  :: l_out_edgeval     !< corresponding local variable; default 
                                  !< .FALSE. i.e. output flux across the edge

    REAL(wp) :: &                 !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< integer fluxes for w>0 (weta<0)
      &  z_iflx_p(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< integer fluxes for w<0 (weta>0)
      &  z_iflx_m(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< monotonized slope
      &  z_slope(nproma,p_patch%nlev)

    REAL(wp) :: p_cc_min, p_cc_max       !< 3-point max/min values

    REAL(wp) :: z_delta_p, z_delta_m     !< difference between upper and lower face value
                                         !< for w>0 and w<0
    REAL(wp) :: z_a11, z_a12             !< 1/6 * a6,i (see Colella and Woodward (1984))

    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: ikm1, ikp1, ikp2         !< vertical level minus and plus one, plus two
    INTEGER  :: slev, slevp1             !< vertical start level and start level +1
    INTEGER  :: slev_ti, slevp1_ti       !< vertical start level (+1)  (tracer independent part)
    INTEGER  :: nlev, nlevp1             !< number of full and half levels

    ! JF: for treatment of sedimentation
    INTEGER  :: elev, elev_lim           !< vertical end level
    LOGICAL  :: llbc_adv                 !< apply lower boundary condition?
    INTEGER  :: ik                       !< = MIN(jk,nlev)

    INTEGER  :: ji_p, ji_m               !< loop variable for index list
    INTEGER  :: ist                      !< status variable
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

    INTEGER  :: nlist_max                !< maximum number of index lists
    INTEGER  :: nlist_p, nlist_m         !< list number
    INTEGER  :: nlist                    !< list loop variable

    REAL(wp), ALLOCATABLE, SAVE  ::   &  !< fractional (mass weighted) Courant number 
      &  z_cflfrac_p(:,:,:)              !< for w>0

    REAL(wp), ALLOCATABLE, SAVE  ::   &  !< fractional (mass weighted) Courant number 
      &  z_cflfrac_m(:,:,:)              !< for w<0

    REAL(wp), ALLOCATABLE, SAVE ::    &  !< maximum vertical Courant number
      &  max_cfl_blk(:)                  !< per block

    INTEGER, ALLOCATABLE, SAVE  ::    &  !< Index and level lists and list dimensions 
      &  i_indlist_p(:,:,:),  &          !< for points with CFL>1,2,3
      &  i_indlist_m(:,:,:),  &
      &  i_levlist_p(:,:,:),  &
      &  i_levlist_m(:,:,:),  &
      &  i_listdim_p(:,:),    &
      &  i_listdim_m(:,:)

    INTEGER, ALLOCATABLE, SAVE  ::    &  !< shifted indices
      &  jk_int_p(:,:,:),             &  ! jk+s
      &  jk_int_m(:,:,:)                 ! jk-s, with shift index s

    INTEGER  :: jk_shift
    INTEGER  :: counter_p, counter_m, &  !< check whether any of the points has 
      &  counter_jip, counter_jim        !< CFL>nlist_p/m

    INTEGER  :: jg                       !< patch ID

    REAL(wp) ::   &                      !< high order flux
      &  z_flx_frac_high

    REAL(wp) ::   &                      !< maximum vertical Courant number
      &  max_cfl(nproma,p_patch%nlevp1)

    REAL(wp) :: &                        !< linear extrapolation value 1
      &  z_lext_1
 
    REAL(wp) :: &                        !< linear extrapolation value 2
      &  z_lext_2

    REAL(wp) ::  &                              !< necessary, to make this routine
      &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 

    REAL(wp) ::   &                      !< dummy variable
      &  z_dummy

    REAL(wp) ::   &                      !< maximum CFL within one layer, and domain-wide maximum
      &  max_cfl_lay(p_patch%nlevp1,p_patch%nblks_c), max_cfl_tot, max_cfl_lay_tot(p_patch%nlevp1)

    REAL(wp) ::   &                      !< auxiliaries for fractional CFL number computation
      &  z_aux_p(nproma), z_aux_m(nproma)

    REAL(wp) ::   &                      !< auxiliaries for optimization
      &   zfac, zfac_n(nproma), zgeo1, zgeo2, zgeo3, zgeo4

    REAL(wp) :: coeff_grid              !< parameter which is used to make the vertical 
                                        !< advection scheme applicable to a height      
                                        !< based coordinate system (coeff_grid=-1)

    REAL(wp) :: rdtime                  !< 1/dt

    !-----------------------------------------------------------------------

    ! inverse of time step for computational efficiency
    rdtime = 1._wp/p_dtime

    ! get patch ID
    jg = p_patch%id

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
      slevp1= opt_slev + 1
    ELSE
      slev  = 1
      slevp1= 2
    END IF

    ! check optional arguments
    IF ( PRESENT(opt_ti_slev) ) THEN
      slev_ti  = opt_ti_slev
      slevp1_ti= opt_ti_slev + 1
    ELSE
      slev_ti  = 1
      slevp1_ti= 2
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF (l_out_edgeval) THEN
      coeff_grid = -1._wp ! needs to be set in case of a call from the hexagonal NH code
    ELSE
      coeff_grid = advection_config(jg)%coeff_grid
    ENDIF

    IF ( PRESENT(opt_topflx_tra) ) THEN
      zparent_topflx(:,:) = opt_topflx_tra(:,:)
    ELSE
      zparent_topflx(:,:) = 0._wp
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! maximum number of lists
    nlist_max = advection_config(jg)%ivcfl_max

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! check optional arguments
    llbc_adv = .TRUE.
    IF ( PRESENT(opt_elev) ) THEN
      IF ( opt_elev == nlevp1 ) THEN
        elev = nlevp1
        llbc_adv = .FALSE.
      ELSE
        elev = nlev
      END IF
    ELSE
      elev = nlev
    END IF
    elev_lim = nlev

    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


    !
    ! advection is done with an upwind scheme where a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    !
    ! 3 options:  standard without limiter
    !             standard with semi-monotone or monotone limiter
    !             special version with limiter which handles CFL >1

    IF ( ld_compute ) THEN
      ! allocate temporary arrays 
      ALLOCATE( i_indlist_p(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_indlist_m(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_levlist_p(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_levlist_m(nproma*nlevp1,nlist_max,p_patch%nblks_c),  &
        &       i_listdim_p(nlist_max,p_patch%nblks_c),                &
        &       i_listdim_m(nlist_max,p_patch%nblks_c),                &
        &       jk_int_p(nproma,nlev,p_patch%nblks_c),                 &
        &       jk_int_m(nproma,nlevp1,p_patch%nblks_c),               &
        &       z_cflfrac_p(nproma,nlevp1,p_patch%nblks_c),            &
        &       z_cflfrac_m(nproma,nlevp1,p_patch%nblks_c),            &
        &       max_cfl_blk(p_patch%nblks_c), STAT=ist                 )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                         &
          &  'allocation for i_indlist_p, i_indlist_m, i_levlist_p, '  //  &
          &  'i_levlist_m, i_listdim_p, i_listdim_m, jk_int_p, '       //  &
          &  'jk_int_m, z_cflfrac_p, z_cflfrac_m, max_cfl_blk  failed '    )
      ENDIF
    END IF

#ifdef _OPENACC
    PRINT *, "Sorry, you are out of luck: the OpenACC version is not yet available"
    IF ( .FALSE. ) THEN
#else
!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jc,ik,ikm1,i_startidx,i_endidx,z_dummy,nlist_p,nlist_m,        &
!$OMP            counter_p,counter_m,counter_jip,counter_jim,max_cfl,                 &
!$OMP            z_aux_p,z_aux_m,ikp1,p_cc_min,p_cc_max,ikp2,nlist,ji_p,    &
!$OMP            ji_m,jk_shift,z_iflx_m,z_iflx_p,z_delta_m,z_delta_p,z_a11,z_a12,     &
!$OMP            zfac, zfac_n, zgeo1, zgeo2, zgeo3, zgeo4,                            &
!$OMP            z_lext_1,z_lext_2,z_slope,z_face,z_face_up,z_face_low,z_flx_frac_high) ICON_OMP_GUIDED_SCHEDULE
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
      &                 i_startidx, i_endidx, i_rlstart, i_rlend )


    ! The contravariant mass flux should never exactly vanish
    !
    IF (l_out_edgeval) THEN
      DO jk = slevp1, elev
        p_mflx_contra_v(i_startidx:i_endidx,jk,jb) =                            &
        &              p_mflx_contra_v(i_startidx:i_endidx,jk,jb)               &
        &              + SIGN(dbl_eps,p_mflx_contra_v(i_startidx:i_endidx,jk,jb))
      ENDDO
    ENDIF

    !
    ! 1. Compute density weighted (fractional) Courant number 
    !    for w<0 (weta>0) and w>0 (weta<0) and integer shift s
    !
    IF (ld_compute) THEN

      ! set start values for index list dimension
      i_listdim_p(1:nlist_max,jb) = 0
      i_listdim_m(1:nlist_max,jb) = 0

      ! (fractional) Courant number at top
      z_cflfrac_p(i_startidx:i_endidx,slev_ti,jb) = 0._wp
      z_cflfrac_m(i_startidx:i_endidx,slev_ti,jb) = 0._wp
      ! (fractional) Courant number for w<0 at bottom
      z_cflfrac_m(i_startidx:i_endidx,nlevp1,jb) = 0._wp


      !
      ! compute (fractional) Courant number
      !
      DO jk = slevp1_ti, elev

        ik   = MIN(jk,nlev)
        ikm1 = jk-1

        DO jc = i_startidx, i_endidx

          ! initialize shift index
          jk_int_p(jc,ik,jb) = ik
          jk_int_m(jc,jk,jb) = ikm1

          z_dummy = p_dtime * ABS(p_mflx_contra_v(jc,jk,jb))

          ! compute Courant number
          z_cflfrac_p(jc,jk,jb) = z_dummy / p_cellmass_now(jc,ik,jb)
          z_cflfrac_m(jc,jk,jb) = z_dummy / p_cellmass_now(jc,ikm1,jb)

          max_cfl(jc,jk) = MAX(z_cflfrac_p(jc,jk,jb),z_cflfrac_m(jc,jk,jb))
        ENDDO

        max_cfl_lay(jk,jb) = MAXVAL(max_cfl(i_startidx:i_endidx,jk))

      ENDDO

      ! (fractional) Courant number  for w>0 at bottom
      z_cflfrac_p(i_startidx:i_endidx,nlevp1,jb) = 0._wp

      max_cfl_blk(jb) = MAXVAL(max_cfl_lay(slevp1_ti:elev,jb))


      ! If CFL>1 then split the CFL number into the fractional CFL number 
      ! and the index shift s.
      IF ( max_cfl_blk(jb) > 1._wp ) THEN

        DO jk = slevp1_ti, nlev

          IF (max_cfl_lay(jk,jb) <= 1._wp) CYCLE

          ! start construction of fractional Courant number
          z_aux_p(i_startidx:i_endidx) = p_dtime*ABS(p_mflx_contra_v(i_startidx:i_endidx,jk,jb))

          ! initialize list number
          nlist_p = 0

          ! checks whether there exists any point with 'large CFL number'
          counter_p = 1

          ! loop until no point has CFL > 1, or nlist_p > nlist_max
          DO WHILE(counter_p > 0 .AND. nlist_p < nlist_max )

            ! get number of current list
            nlist_p     = nlist_p + 1
            ! re-initialize counter for CFL>nlist_p
            counter_p   = 0

            ! copy value from counter in vector form to counter in scalar form, 
            ! since otherwise the following DO-Loop will not vectorize.
            counter_jip = i_listdim_p(nlist_p,jb)

            DO jc = i_startidx, i_endidx

              IF ( z_aux_p(jc) > p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)    &
                &  .AND. jk_int_p(jc,jk,jb) <= nlev-1                        &
                &  .AND. coeff_grid * p_mflx_contra_v(jc,jk,jb) < 0._wp ) THEN

                z_aux_p(jc) = z_aux_p(jc) - p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)

                ! update shift index
                jk_int_p(jc,jk,jb) = jk_int_p(jc,jk,jb) + 1

                ! tests whether we need to loop once again
                counter_p = counter_p + 1

                ! Fill index lists with those points that need index shifts
                ! Note that we have to use a scalar counter instead of a vector, like
                ! i_listdim_p(nlist_p,jb). Otherwise this loop will not vectorize.  
                counter_jip = counter_jip + 1
                i_indlist_p(counter_jip,nlist_p,jb) = jc
                i_levlist_p(counter_jip,nlist_p,jb) = jk

                ! compute fractional Courant number
                z_cflfrac_p(jc,jk,jb) = z_aux_p(jc) / p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)

              ENDIF

            END DO ! end loop over cells

            ! store index of current index list
            ! after the last jk-loop this will be the list dimension
            i_listdim_p(nlist_p,jb) = counter_jip

          ENDDO  ! DO WHILE loop
 
        ENDDO ! end loop over vertical levels

        DO jk = slevp1_ti, elev

          IF (max_cfl_lay(jk,jb) <= 1._wp) CYCLE

          ! start construction of fractional Courant number
          z_aux_m(i_startidx:i_endidx) = p_dtime*ABS(p_mflx_contra_v(i_startidx:i_endidx,jk,jb))

          ! initialize list number
          nlist_m = 0

          ! checks whether there exists any point with 'large CFL number'
          counter_m = 1

          ! loop until no point has CFL>nlist_m, or nlist_m > nlist_max
          DO WHILE(counter_m > 0 .AND. nlist_m < nlist_max )

            ! set number of current list
            nlist_m     = nlist_m + 1
            ! re-initialize counter for CFL>nlist_m
            counter_m   = 0

            ! copy value from counter in vector form to counter in scalar form, 
            ! since otherwise the following DO-Loop will not vectorize.
            counter_jim = i_listdim_m(nlist_m,jb)

            DO jc = i_startidx, i_endidx

              IF ( z_aux_m(jc) > p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)    &
                &  .AND. jk_int_m(jc,jk,jb) >= slevp1_ti                     &
                &  .AND. coeff_grid * p_mflx_contra_v(jc,jk,jb) > 0._wp ) THEN

                z_aux_m(jc) = z_aux_m(jc) - p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)

                ! Index shift
                jk_int_m(jc,jk,jb) = jk_int_m(jc,jk,jb) - 1

                ! checks whether we need to loop once again
                counter_m = counter_m + 1

                ! Fill index lists with those points that need index shifts
                ! Note that we have to use a scalar instead of a vector, like
                ! i_listdim_m(nlist_m,jb). Otherwise it will not vectorize  
                counter_jim = counter_jim + 1
                i_indlist_m(counter_jim,nlist_m,jb) = jc
                i_levlist_m(counter_jim,nlist_m,jb) = jk

                ! compute fractional Courant number
                z_cflfrac_m(jc,jk,jb) = z_aux_m(jc) / p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)
              ENDIF

            END DO ! end loop over cells

            ! store index of current index list.
            ! after the last jk-loop this will be the list dimension
            i_listdim_m(nlist_m,jb) = counter_jim

          ENDDO  ! DO WHILE loop


        ENDDO ! end loop over vertical levels

      END IF  ! IF ( max_cfl_blk(jb) > 1._wp )

    END IF ! ld_compute

    !
    ! 2. Compute monotonized slope
    !

      ! Initialize z_slope and zfac_n for jk=slev
      z_slope(i_startidx:i_endidx,slev) = 0._wp
      zfac_n(i_startidx:i_endidx) = 1._wp/(p_cellhgt_mc_now(i_startidx:i_endidx,slevp1,jb) &
        &                         + p_cellhgt_mc_now(i_startidx:i_endidx,slev,jb))         &
        &                         * (p_cc(i_startidx:i_endidx,slevp1,jb) - p_cc(i_startidx:i_endidx,slev,jb))

      DO jk = slevp1, nlev

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, nlev )

        DO jc = i_startidx, i_endidx
          zfac = 1._wp / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))

          z_slope(jc,jk) = ( p_cellhgt_mc_now(jc,jk,jb)                                          &
            &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)                      &
            &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                                 &
            &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) * zfac    &
            &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb)) * zfac_n(jc))

          zfac_n(jc) = zfac

          ! equivalent formulation of Colella and Woodward (1984) slope limiter 
          ! following Lin et al (1994).
          p_cc_min = MIN(p_cc(jc,ikm1,jb),p_cc(jc,jk,jb),p_cc(jc,ikp1,jb))
          p_cc_max = MAX(p_cc(jc,ikm1,jb),p_cc(jc,jk,jb),p_cc(jc,ikp1,jb))
          z_slope(jc,jk) = SIGN(                                            &
            &  MIN( ABS(z_slope(jc,jk)), 2._wp*(p_cc(jc,jk,jb)-p_cc_min),   &
            &                            2._wp*(p_cc_max-p_cc(jc,jk,jb)) ), &
            &    z_slope(jc,jk))
           
        END DO

      END DO




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

        z_face(jc,slevp1) = p_cc(jc,slev,jb)*(1._wp - (p_cellhgt_mc_now(jc,slev,jb)   &
          &       / p_cellhgt_mc_now(jc,slevp1,jb))) + (p_cellhgt_mc_now(jc,slev,jb)  &
          &       /(p_cellhgt_mc_now(jc,slev,jb) + p_cellhgt_mc_now(jc,slevp1,jb)))   &
          &       * ((p_cellhgt_mc_now(jc,slev,jb) / p_cellhgt_mc_now(jc,slevp1,jb))  &
          &       * p_cc(jc,slev,jb) + p_cc(jc,slevp1,jb))

        z_face(jc,nlev) = p_cc(jc,nlev-1,jb)*( 1._wp                                  &
          &       - (p_cellhgt_mc_now(jc,nlev-1,jb) / p_cellhgt_mc_now(jc,nlev,jb)))  &
          &       + (p_cellhgt_mc_now(jc,nlev-1,jb)/(p_cellhgt_mc_now(jc,nlev-1,jb)   &
          &       + p_cellhgt_mc_now(jc,nlev,jb))) * ((p_cellhgt_mc_now(jc,nlev-1,jb) &
          &       / p_cellhgt_mc_now(jc,nlev,jb)) * p_cc(jc,nlev-1,jb)                &
          &       + p_cc(jc,nlev,jb))

        z_face(jc,slev)   = p_cc(jc,slev,jb)
        z_face(jc,nlevp1) = p_cc(jc,nlev,jb)

      ENDDO


      DO jk = slevp1, nlev-2

        ! index of top half level
        ikm1 = jk - 1
        ! index of bottom half level
        ikp1 = jk + 1
        ikp2 = jk + 2

        DO jc = i_startidx, i_endidx
          zgeo1 = p_cellhgt_mc_now(jc,jk,jb)                                         &
            &   / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))
          zgeo2 = 1._wp / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb) &
            &   + p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,ikp2,jb))
          zgeo3 = (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb))        &
            &   / (2._wp*p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikp1,jb))
          zgeo4 = (p_cellhgt_mc_now(jc,ikp2,jb) + p_cellhgt_mc_now(jc,ikp1,jb))      &
            &   / (2._wp*p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb))


          z_face(jc,ikp1) = p_cc(jc,jk,jb)                                  &
            &  + zgeo1 * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))                &
            &  + zgeo2 * ( (2._wp * p_cellhgt_mc_now(jc,ikp1,jb) * zgeo1)   &
            &  * ( zgeo3 - zgeo4 ) * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))    &
            &  - zgeo3 * p_cellhgt_mc_now(jc,jk,jb)   * z_slope(jc,ikp1)    &
            &  + zgeo4 * p_cellhgt_mc_now(jc,ikp1,jb) * z_slope(jc,jk) )

        END DO

      END DO



      !
      ! 4. Limitation of first guess parabola (which is based on z_face)
      ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
      ! the limitation procedure.
      ! Therefore 2 additional fields z_face_up and z_face_low are
      ! introduced.
      !
      IF (p_itype_vlimit == islopel_vsm) THEN
        ! semi-monotonic (sm) limiter
        CALL v_ppm_slimiter_sm( p_cc(:,:,jb), z_face(:,:),          & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in
      ELSE IF (p_itype_vlimit == islopel_vm) THEN
        ! monotonic (mo) limiter
        CALL v_ppm_slimiter_mo( p_cc(:,:,jb), z_face(:,:), z_slope(:,:), & !in
          &                   z_face_up(:,:), z_face_low(:,:),           & !inout
          &                   i_startidx, i_endidx, slev, elev_lim       ) !in
      ENDIF


      IF (p_itype_vlimit /= islopel_vsm .AND. p_itype_vlimit /= islopel_vm) THEN
        ! simply copy face values to 'face_up' and 'face_low' arrays

        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          z_face_up(i_startidx:i_endidx,jk)  = z_face(i_startidx:i_endidx,jk)
          z_face_low(i_startidx:i_endidx,jk) = z_face(i_startidx:i_endidx,ikp1)
        ENDDO

      ENDIF


      !
      ! 5. calculation of upwind fluxes. IF CFL > 1, the fluxes are the sum of
      !    integer-fluxes and a fractional flux. IF CFL <1 the fluxes are only
      !    comprised of the fractional flux. The fractional flux is calculated
      !    assuming a piecewise parabolic approx. for the subgrid distribution.
      !


      !
      ! 5a. First compute fluxes for the CFL<1 case for all grid points
      ! On the grid points where CFL>1, they will be overwritten afterwards 
      ! This part has been adopted from the restricted time step PPM-scheme.
      !
      DO jk = slevp1, elev

        ik   = MIN(jk,nlev)
        ! index of top half level
        ikm1 = jk -1

        DO jc = i_startidx, i_endidx

          ! if w < 0 , weta > 0 (physical downwelling)
          !
          ! note that the second coeff_grid factor in front of z_delta_p 
          ! is obsolete due to a compensating (-) sign emerging from the 
          ! computation of z_delta_p.
          z_delta_m = z_face_up(jc,ikm1) - z_face_low(jc,ikm1)
          z_a11     = p_cc(jc,ikm1,jb)                                  &
            &       - 0.5_wp * (z_face_low(jc,ikm1) + z_face_up(jc,ikm1))

          z_lext_1 = p_cc(jc,ikm1,jb)                                   &
            &  - (0.5_wp * z_delta_m * (1._wp - z_cflfrac_m(jc,jk,jb))) &
            &  - z_a11*(1._wp - 3._wp*z_cflfrac_m(jc,jk,jb)             &
            &  + 2._wp*z_cflfrac_m(jc,jk,jb)*z_cflfrac_m(jc,jk,jb))


          ! if w > 0 , weta < 0 (physical upwelling)
          !
          ! note that the second coeff_grid factor in front of z_delta_p 
          ! is obsolete due to a compensating (-) sign emerging from the 
          ! computation of z_delta_p.
          z_delta_p = z_face_up(jc,ik) - z_face_low(jc,ik)
          z_a12     = p_cc(jc,ik,jb)                                    &
            &       - 0.5_wp * (z_face_low(jc,ik) + z_face_up(jc,ik))

          z_lext_2 = p_cc(jc,ik,jb)                                     &
            &  + (0.5_wp * z_delta_p * (1._wp - z_cflfrac_p(jc,ik,jb))) &
            &  - z_a12*(1._wp - 3._wp*z_cflfrac_p(jc,ik,jb)             &
            &  + 2._wp*z_cflfrac_p(jc,ik,jb)*z_cflfrac_p(jc,ik,jb))

          !
          ! full flux
          !
          p_upflux(jc,jk,jb) =                                  &
            &  laxfr_upflux_v( p_mflx_contra_v(jc,jk,jb),       &
            &                z_lext_1, z_lext_2,                &
            &                coeff_grid )

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels


      !
      ! 5b. Now execute the special computations needed for CFL>1:
      !     Computation of integer fluxes and a fractional flux
      !
      IF (max_cfl_blk(jb) > 1._wp) THEN

        z_iflx_m(i_startidx:i_endidx,slev:nlevp1) = 0._wp
        z_iflx_p(i_startidx:i_endidx,slev:nlevp1) = 0._wp

        ! Loop over all lists (nlist will serve as shift index)
        DO nlist = 1, nlist_max

          IF (i_listdim_p(nlist,jb) == 0) CYCLE

          !
          ! loop over all cells in i_indlist_p
          !
          ! integer fluxes for w>0 (weta<0)
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ji_p=1,i_listdim_p(nlist,jb)

            ! get jc and jk index from precomputed list
            jc = i_indlist_p(ji_p,nlist,jb)
            jk = i_levlist_p(ji_p,nlist,jb)

            ! cycle if the model level is in a region where advection is 
            ! turned off for the present variable
            IF (jk < slevp1) CYCLE

            ! integer shift (depends on the applied list)
            jk_shift = jk + nlist - 1

            ! Integer flux (division by p_dtime is done at the end)
            z_iflx_p(jc,jk) = z_iflx_p(jc,jk) - coeff_grid * p_cc(jc,jk_shift,jb) &
              &              * p_cellmass_now(jc,jk_shift,jb)

          ENDDO  ! loop over cells in i_indlist_p

        ENDDO ! list loop


        ! Now use the first list (containing all points with CFL>1 and upward flux)
        ! to compute the corrected fluxes
        IF (i_listdim_p(1,jb) > 0) THEN
!CDIR NODEP,VOVERTAKE,VOB
          DO ji_p=1,i_listdim_p(1,jb)

            ! get jc and jk index from precomputed list
            jc = i_indlist_p(ji_p,1,jb)
            jk = i_levlist_p(ji_p,1,jb)

            ! cycle if the model level is in a region where advection is 
            ! turned off for the present variable
            IF (jk < slevp1) CYCLE

            ! fractional upward flux
            ! if w > 0 , weta < 0 (physical upwelling)
            ! note that the second coeff_grid factor in front of z_delta_p 
            ! is obsolete due to a compensating (-) sign emerging from the 
            ! computation of z_delta_p.
            z_delta_p = z_face_up(jc,jk_int_p(jc,jk,jb))                &
              &         - z_face_low(jc,jk_int_p(jc,jk,jb))
            z_a12     = p_cc(jc,jk_int_p(jc,jk,jb),jb)                 &
              &         - 0.5_wp * (z_face_low(jc,jk_int_p(jc,jk,jb))   &
              &         + z_face_up(jc,jk_int_p(jc,jk,jb)))


            ! fractional high order flux   
            z_flx_frac_high = ( - coeff_grid                                        &
              &         * p_cellmass_now(jc,jk_int_p(jc,jk,jb),jb)                  &
              &         * z_cflfrac_p(jc,jk,jb) *( p_cc(jc,jk_int_p(jc,jk,jb),jb)   &
              &         + (0.5_wp * z_delta_p * (1._wp - z_cflfrac_p(jc,jk,jb)))    &
              &         - z_a12*(1._wp - 3._wp*z_cflfrac_p(jc,jk,jb)                &
              &         + 2._wp*z_cflfrac_p(jc,jk,jb)**2) ) )                       &
              &         * rdtime

            ! full flux (integer- plus high order fractional flux)
            p_upflux(jc,jk,jb) = z_iflx_p(jc,jk)*rdtime + z_flx_frac_high

          ENDDO
        ENDIF


        DO nlist = 1, nlist_max

          IF (i_listdim_m(nlist,jb) == 0) CYCLE

          !
          ! loop over all cells in i_indlist_m
          !
          ! integer fluxes for w<0 (weta>0)
          !
!CDIR NODEP,VOVERTAKE,VOB
          DO ji_m=1,i_listdim_m(nlist,jb)

            ! get jc and jk index from precomputed list
            jc = i_indlist_m(ji_m,nlist,jb)
            jk = i_levlist_m(ji_m,nlist,jb)

            ! integer shift (depends on the applied list)
            jk_shift = jk - nlist

            ! cycle if the source model level is in a region where advection is 
            ! turned off for the present variable
            IF (jk_shift < slevp1) CYCLE

            ! Integer flux (division by p_dtime is done at the end)
            z_iflx_m(jc,jk) = z_iflx_m(jc,jk) + coeff_grid*p_cc(jc,jk_shift,jb) &
              &              * p_cellmass_now(jc,jk_shift,jb)

          ENDDO  ! loop over cells in i_indlist_m

        ENDDO  ! loop over index lists


        ! Now use the first list (containing all points with CFL>1 and downward flux)
        ! to compute the corrected fluxes
        IF (i_listdim_m(1,jb) > 0) THEN
!CDIR NODEP,VOVERTAKE,VOB
          DO ji_m=1,i_listdim_m(1,jb)

            ! get jc and jk index from precomputed list
            jc = i_indlist_m(ji_m,1,jb)
            jk = i_levlist_m(ji_m,1,jb)

            ! cycle if the model level is in a region where advection is 
            ! turned off for the present variable
            IF (jk < slevp1) CYCLE

            ! this is needed in addition in order to avoid accessing non-existing (uninitalized)
            ! source levels for tracers that are not advected on all model levels
            IF (jk_int_m(jc,jk,jb) < slev) CYCLE

            ! fractional downward flux
            ! if w < 0 , weta > 0 (physical downwelling)
            ! note that the second coeff_grid factor in front of z_delta_p 
            ! is obsolete due to a compensating (-) sign emerging from the 
            ! computation of z_delta_p.
            z_delta_m = z_face_up(jc,jk_int_m(jc,jk,jb))                &
              &         - z_face_low(jc,jk_int_m(jc,jk,jb))
            z_a11     = p_cc(jc,jk_int_m(jc,jk,jb),jb)                  &
              &         - 0.5_wp * (z_face_low(jc,jk_int_m(jc,jk,jb))   &
              &         + z_face_up(jc,jk_int_m(jc,jk,jb)))


            ! fractional high order flux           
            z_flx_frac_high = ( coeff_grid                                          &
              &         * p_cellmass_now(jc,jk_int_m(jc,jk,jb),jb)                  &
              &         * z_cflfrac_m(jc,jk,jb) * ( p_cc(jc,jk_int_m(jc,jk,jb),jb)  &
              &         - (0.5_wp * z_delta_m * (1._wp - z_cflfrac_m(jc,jk,jb)))    &
              &         - z_a11*(1._wp - 3._wp*z_cflfrac_m(jc,jk,jb)                &
              &         + 2._wp*z_cflfrac_m(jc,jk,jb)**2) ) )                       &
              &         * rdtime

            ! full flux (integer- plus fractional flux)
            p_upflux(jc,jk,jb) = z_iflx_m(jc,jk)*rdtime + z_flx_frac_high

          ENDDO
        ENDIF

      ENDIF  ! IF (max_cfl_blk(jb) > 1)


      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(nproma, p_upflux(:,slev+1,jb),    &! in
        &              p_mflx_contra_v(:,slev+1,jb),     &! in
        &              p_mflx_contra_v(:,slev  ,jb),     &! in
        &              p_iubc_adv, i_startidx, i_endidx, &! in
        &              zparent_topflx(:,jb),             &! in
        &              p_upflux(:,slev,jb),              &! out
        &              p_upflux(:,nlevp1,jb), llbc_adv)   ! out


      ! If desired, get edge value of advected quantity 
      IF ( l_out_edgeval ) THEN

        DO jk = slevp1, nlev
          DO jc = i_startidx, i_endidx
            p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb)/p_mflx_contra_v(jc,jk,jb)
          ENDDO
        ENDDO

      ENDIF

    ENDDO
#ifdef _OPENACC
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

    !
    ! 6. If desired, apply a flux limiter to limit computed fluxes.
    !    These flux limiters are based on work by Zalesak (1979)
    !
    IF ( iequations /= 3 ) THEN
      IF (p_itype_vlimit == ifluxl_vpd) THEN
        ! positive-definite (pd) limiter
        CALL vflx_limiter_pd_ha( p_patch, p_dtime, p_cc, p_upflux,  & !in,inout
          &                   opt_rlstart=i_rlstart,                & !in
          &                   opt_rlend=i_rlend, opt_slev=slev      ) !in
      ENDIF
    ELSE
      IF (p_itype_vlimit == ifluxl_vpd) THEN
        ! positive-definite (pd) limiter
        CALL vflx_limiter_pd( p_patch, p_dtime, p_cc, p_upflux,     & !in,inout
          &                   opt_rlstart=i_rlstart,                & !in
          &                   opt_rlend=i_rlend, opt_slev=slev      ) !in
      ENDIF
    ENDIF


#ifdef _OPENACC
! OpenACC is not available.  Skipping entire routine
    ENDIF
#else
    !
    ! If desired, print maximum vertical CFL number
    !
    IF ( ld_compute .AND. msg_level >= 10 .AND. lprint_cfl) THEN

      max_cfl_tot = MAXVAL(max_cfl_blk(i_startblk:i_endblk))

      ! Take maximum over all PEs
      IF (msg_level >= 13) THEN
        max_cfl_tot = global_max(max_cfl_tot)
      ELSE
        max_cfl_tot = global_max(max_cfl_tot, iroot=process_mpi_stdio_id)
      ENDIF
      IF (my_process_is_stdio() .OR. comm_lev>0 .AND. get_my_mpi_work_id() == get_glob_proc0() ) THEN
        ! otherwise it is possible that max_cfl_tot is undefined
        WRITE(message_text,'(a,e16.8)') 'maximum vertical CFL =',max_cfl_tot
        CALL message(TRIM(routine),message_text)
      ENDIF

      ! Add layer-wise diagnostic if the maximum CFL value is close to the stability limit
      IF (msg_level >= 13 .AND. max_cfl_tot > 4._wp) THEN
        DO jk = slevp1_ti, nlev
          max_cfl_lay_tot(jk) = MAXVAL(max_cfl_lay(jk,i_startblk:i_endblk))
        ENDDO

        max_cfl_lay_tot(slevp1_ti:nlev) = global_max(max_cfl_lay_tot(slevp1_ti:nlev), iroot=process_mpi_stdio_id)
        DO jk = slevp1_ti,nlev
          WRITE(message_text,'(a,i4,a,e16.8)') 'maximum vertical CFL in layer', jk,' =', max_cfl_lay_tot(jk)
          CALL message(TRIM(routine),message_text)
        ENDDO
      ENDIF

    END IF
#endif

    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays
      DEALLOCATE( i_indlist_p, i_indlist_m, i_levlist_p,           &
        &         i_levlist_m, i_listdim_p, i_listdim_m, jk_int_p, &
        &         jk_int_m, z_cflfrac_p, z_cflfrac_m, max_cfl_blk, & 
        &         STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine),                     &
          &  'deallocation for i_indlist_p, i_indlist_m, i_levlist_p, ' //  &
          &  'i_levlist_m, i_listdim_p, i_listdim_m, jk_int_p, '        //  &
          &  'jk_int_m, z_cflfrac_p, z_cflfrac_m, max_cfl_blk failed '      )
      ENDIF
    END IF


  END SUBROUTINE upwind_vflux_ppm_cfl



  !-------------------------------------------------------------------------
  !>
  !! Set upper and lower boundary condition for vertical transport
  !!
  !! Set upper and lower boundary condition for vertical transport.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-04-12)
  !!
  !
  SUBROUTINE set_bc_vadv(nproma, upflx_top_p1, mflx_top_p1, mflx_top, iubc_adv, &
    &                    i_start, i_end, parent_topflx, upflx_top,      &
    &                    upflx_bottom, llbc_adv )

!$ACC ROUTINE VECTOR

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_advection_vflux: set_ubc_adv'

    INTEGER, INTENT(IN)      :: nproma  ! To avoid the declaration of a module var.
    REAL(wp), INTENT(IN)     :: & !< computed tracer flux at second half level
      &  upflx_top_p1(nproma)
    REAL(wp), INTENT(IN)     :: & !< mass flux at second half level
      &  mflx_top_p1(nproma)
    REAL(wp), INTENT(IN)     :: & !< mass flux at upper boundary
      &  mflx_top(nproma)
    INTEGER, INTENT(IN)      :: & !< selects upper boundary condition
      &  iubc_adv
    INTEGER, INTENT(IN)      :: & !< start and end index
      & i_start, i_end
    REAL(wp), INTENT(IN)     :: & !< tracer flux at upper boundary, 
      &  parent_topflx(nproma)    !< interpolated from parent grid
    REAL(wp), INTENT(OUT)    :: & !< upper boundary condition
      &  upflx_top(nproma)
    REAL(wp), INTENT(OUT)    :: & !< lower boundary condition
      &  upflx_bottom(nproma)
    LOGICAL, INTENT(IN)      :: & !< apply lower boundary condition?
      &  llbc_adv

    !-------------------------------------------------------------------------

    ! 
    ! flux at top boundary
    ! 
    SELECT CASE (iubc_adv)
      CASE ( ino_flx )     ! no flux
        upflx_top(i_start:i_end) = 0._wp
 
      CASE ( izero_grad )  ! zero gradient
        upflx_top(i_start:i_end) = upflx_top_p1(i_start:i_end)      &
            &           * mflx_top(i_start:i_end)                   &
            &           / ( mflx_top_p1(i_start:i_end)              &
            &           + SIGN(dbl_eps, mflx_top_p1(i_start:i_end)))

      CASE ( iparent_flx ) ! interpolated flux from parent grid
        upflx_top(i_start:i_end) = parent_topflx(i_start:i_end)
    END SELECT

    !
    ! flux at bottom boundary
    !
    IF ( llbc_adv ) THEN
      upflx_bottom(i_start:i_end) = 0._wp
    END IF

  END SUBROUTINE set_bc_vadv



END MODULE mo_advection_vflux

