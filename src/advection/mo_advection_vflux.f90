!>
!! Computation of vertical tracer flux
!!
!! Vertical fluxes are calculated at triangle centers on half-levels.
!! Possible options for vertical flux calculation include
!! - first order Godunov method (UP1)
!! - third order PPM method without CFL restriction
!! - third order PSM method without CFL restriction
!!
!! Semi-monotone and monotone limiters are available for PPM
!!
!! These routines compute only the correct half level value of
!! 'c*w'. The vertical divergence is computed in step_advection.
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
!! Modification by Daniel Reinert, DWD (2015-11-26)
!! - included parabolic spline method (PSM)
!! Modification by Daniel Reinert, DWD (2016-03 ?)
!! - refactoring in upwind_vflux_ppm_cfl
!! Modification by Will Sawyer, CSCS (2016-07-15)
!! - added OpenACC support
!! Modification by Daniel Reinert, DWD (2019-02-18)
!! - remove restricted time step version of PPM and rename 
!!   upwind_vflux_ppm_cfl to upwind_vflux_ppm.
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
    &                               iup_v, ippm_v, ipsm_v, ippm4gpu_v,          &
    &                               islopel_vsm, islopel_vm, ifluxl_vpd,        &
    &                               ino_flx, izero_grad, iparent_flx
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_math_utilities,      ONLY: tdma_solver_vec
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level, lvert_nest, timers_level, iqtke
  USE mo_advection_config,    ONLY: advection_config, lcompute, lcleanup, t_trList 
  USE mo_advection_vlimit,    ONLY: v_limit_parabola_mo, v_limit_parabola_sm, &
   &                                vflx_limiter_pd,                          &
   &                                v_limit_slope_mo, v_limit_slope_sm,       &
   &                                v_limit_face_mo, v_limit_face_sm
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: global_max
  USE mo_mpi,                 ONLY: process_mpi_stdio_id, my_process_is_stdio, get_my_mpi_work_id, &
                                    get_glob_proc0, comm_lev
#ifdef _OPENACC
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, check_patch_array
  USE mo_mpi,                 ONLY: i_am_accel_node, my_process_is_work
#endif
  USE mo_timer,               ONLY: timer_adv_vert, timer_start, timer_stop


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: vert_upwind_flux
  PUBLIC :: upwind_vflux_up
  PUBLIC :: upwind_vflux_ppm
  PUBLIC :: upwind_vflux_ppm4gpu
  PUBLIC :: implicit_sedim_tracer

#if defined( _OPENACC )
#if defined(__ADVECTION_VFLUX_NOACC)
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
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !!
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !! using either
  !! - the first order Godunov method (UP1)
  !! - the third order PPM method
  !! - the third order PSM method
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
  ! see below
  !
  SUBROUTINE vert_upwind_flux( p_patch, p_cc, p_mflx_contra_v,                &
    &                      p_dtime, p_cellhgt_mc_now,                         &
    &                      p_cellmass_now, p_ivadv_tracer,                    &
    &                      p_itype_vlimit, p_ivlimit_selective,               &
    &                      p_iubc_adv, p_iadv_slev,                           &
    &                      lprint_cfl, p_upflux, opt_topflx_tra, opt_q_int,   &
    &                      opt_rlstart, opt_rlend  )

   CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux: vert_upwind_flux'

    TYPE(t_patch), INTENT(IN) ::  &  !< patch on which computation is 
      &  p_patch                             !< performed

    REAL(wp), INTENT(IN) ::  &      !< advected cell centered variable
      &  p_cc(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT) ::  &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)     !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    REAL(wp), INTENT(IN) ::  &      !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)    !< time step n (either \Delta p or \Delta z)
                                    !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN)::  &       !< NH: density times cell thickness at cell center
      &  p_cellmass_now(:,:,:)      !< at time step n [kg/m**2]
                                    !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN) ::   &      !< parameter to select numerical
      &  p_ivadv_tracer(:)          !< scheme for vertical transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< parameter to select the limiter
      &  p_itype_vlimit(:)          !< for vertical transport
                                    !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< avoids limiting of smooth extrema
      &  p_ivlimit_selective(:)     !< if activated

    INTEGER, INTENT(IN) ::   &      !< vertical start level for transport
      &  p_iadv_slev(:)             !< dim: (ntracer)

    INTEGER, INTENT(IN) ::   &      !< selects upper boundary condition
      &  p_iubc_adv

    LOGICAL, INTENT(IN) ::   &      !< determines if vertical CFL number shall be printed
      &  lprint_cfl                 !< in routine upwind_vflux_ppm

    REAL(wp), INTENT(INOUT) :: &    !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:,:)          !< dim: (nproma,nlevp1,nblks_c,ntracer)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:,:)          !< NH: [kg/m**2/s]
                                        !< dim: (nproma,nblks_c,ntracer)

    REAL(wp), INTENT(OUT), OPTIONAL :: & !< tracer value at upper boundary of child nest 
      &  opt_q_int(:,:,:)               !< NH: [kg/kg]
                                        !< dim: (nproma,nblks_c,ntracer)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    TYPE(t_trList), POINTER ::  &      !< pointer to tracer sublist
      &  trAdvect

    INTEGER :: jt, nt                  !< tracer index and loop index
    INTEGER :: jg                      !< patch ID
    INTEGER :: jc, jb                  !< cell and block loop index
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart_c, i_rlend_c
    INTEGER :: iadv_min_slev           !< scheme specific minimum slev

    REAL(wp) :: z_mflx_contra_v(nproma) !< auxiliary variable for computing vertical nest interface quantities

#ifdef _OPENACC
    LOGICAL  :: save_i_am_accel_node
#endif

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_mflx_contra_v
#endif
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


      CASE( ippm_v, ipsm_v )

        iadv_min_slev = advection_config(jg)%ppm_v%iadv_min_slev

#ifdef _OPENACC
! In GPU mode, copy data to HOST and perform upwind_vflux_ppm there, then update device
! NOTE: this is only for testing; use upwind_vflux_ppm4gpu for performance
        WRITE(message_text,'(a)') 'GPU mode: performing upwind_vflux_ppm on host; for performance use upwind_vflux_ppm4gpu'
        CALL message(TRIM(routine),message_text)
!$ACC UPDATE HOST( p_cc(:,:,:,jt), p_cellhgt_mc_now, p_cellmass_now  ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE HOST( p_upflux(:,:,:,jt), p_mflx_contra_v ), IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE HOST( opt_topflx_tra(:,:,jt) ), IF( i_am_accel_node .AND. acc_on .AND. PRESENT(opt_topflx_tra ) )
        save_i_am_accel_node = i_am_accel_node
        i_am_accel_node = .FALSE.                  ! deactivate GPUs throughout upwind_vflux_ppm
#endif
        ! CALL third order PPM/PSM (unrestricted timestep-version) (i.e. CFL>1)
        CALL upwind_vflux_ppm( p_patch, p_cc(:,:,:,jt), p_iubc_adv,        &! in
          &                  p_mflx_contra_v, p_dtime, lcompute%ppm_v(jt), &! in
          &                  lcleanup%ppm_v(jt), p_itype_vlimit(jt),       &! in
          &                  p_ivlimit_selective(jt),                      &! in
          &                  p_cellhgt_mc_now, p_cellmass_now,             &! in
          &                  lprint_cfl,                                   &! in
          &                  p_ivadv_tracer(jt),                           &! in
          &                  p_upflux(:,:,:,jt),                           &! out
          &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),        &! in
          &                  opt_slev=p_iadv_slev(jt),                     &! in
          &                  opt_ti_slev=iadv_min_slev,                    &! in
          &                  opt_rlstart=opt_rlstart,                      &! in
          &                  opt_rlend=i_rlend_c                           )! in
#ifdef _OPENACC
        i_am_accel_node =  save_i_am_accel_node    ! reactivate GPUs if appropriate
!$ACC UPDATE DEVICE( p_upflux(:,:,:,jt), p_mflx_contra_v ), IF( i_am_accel_node .AND. acc_on )
#endif


      CASE( ippm4gpu_v )

        iadv_min_slev = advection_config(jg)%ppm4gpu_v%iadv_min_slev

        ! CALL third order PPM (unrestricted timestep-version, optimized for GPU)
        CALL upwind_vflux_ppm4gpu( p_patch, p_cc(:,:,:,jt), p_iubc_adv,    &! in
            &                  p_mflx_contra_v, p_dtime, lcompute%ppm4gpu_v(jt), &! in
            &                  lcleanup%ppm4gpu_v(jt), p_itype_vlimit(jt),   &! in
            &                  p_ivlimit_selective(jt),                      &! in
            &                  p_cellhgt_mc_now, p_cellmass_now, lprint_cfl, &! in
            &                  p_upflux(:,:,:,jt),                           &! out
            &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),        &! in
            &                  opt_slev=p_iadv_slev(jt),                     &! in
            &                  opt_ti_slev=iadv_min_slev,                    &! in
            &                  opt_rlstart=opt_rlstart,                      &! in
            &                  opt_rlend=i_rlend_c                           )! in

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

      i_startblk = p_patch%cells%start_block(i_rlstart_c)
      i_endblk   = p_patch%cells%end_block(i_rlend_c)

!$OMP PARALLEL DO PRIVATE(jb,jt,jc,nt,i_startidx,i_endidx,z_mflx_contra_v) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

        ! Be sure to avoid division by zero
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          z_mflx_contra_v(jc) = SIGN( MAX(ABS(p_mflx_contra_v(jc,p_patch%nshift_child,jb)),dbl_eps), &
            &                                 p_mflx_contra_v(jc,p_patch%nshift_child,jb) )
        ENDDO
!$ACC END PARALLEL

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG
        DO nt = 1, trAdvect%len
          jt = trAdvect%list(nt)
          !$ACC LOOP VECTOR
          DO jc = i_startidx, i_endidx
            opt_q_int(jc,jb,jt) = p_upflux(jc,p_patch%nshift_child,jb,jt) / z_mflx_contra_v(jc)
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
!$OMP END PARALLEL DO

    ENDIF

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

    TYPE(t_patch), INTENT(IN) ::  & !< patch on which computation is performed
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
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    REAL(wp) ::  &                              !< necessary, to make this routine
      &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 
                                       
    INTEGER  :: slev                   !< vertical start level
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: jc, jk, jb             !< index of cell, vertical level and block
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend
    !-------------------------------------------------------------------------


!$ACC DATA CREATE( zparent_topflx ), PCOPYIN( p_cc, p_mflx_contra_v ), PCOPYOUT( p_upflux ), &
!$ACC IF( i_am_accel_node .AND. acc_on )

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: zparent_topflx
#endif

!$ACC UPDATE DEVICE( p_cc, p_mflx_contra_v ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev+1, nlev
        DO jc = i_startidx, i_endidx
          ! calculate vertical tracer flux   -- removed flaky laxfr macro
          p_upflux(jc,jk,jb) = p_mflx_contra_v(jc,jk,jb) *                    &
                               MERGE( p_cc(jc,jk,jb),p_cc(jc,jk-1,jb),        &
                                      p_mflx_contra_v(jc,jk,jb) .GE. 0.0_wp ) 
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
!$ACC END PARALLEL

      !
      ! set upper and lower boundary condition
      !
!
! With OpenACC this sometimes causes the compiler to crash, apparently due to the array syntax of one argument.  
!
      CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
        &              p_mflx_contra_v(:,slev+1,jb),     &! in
        &              p_mflx_contra_v(:,slev  ,jb),     &! in
        &              p_iubc_adv, i_startidx, i_endidx, &! in
        &              zparent_topflx(:,jb),             &! in
        &              p_upflux(:,slev,jb),              &! out
        &              p_upflux(:,nlevp1,jb), .TRUE.)     ! out

    ENDDO ! end loop over blocks

!$ACC UPDATE HOST( p_upflux ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
!$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE upwind_vflux_up





  !-------------------------------------------------------------------------
  !>
  !! The third order PPM/PSM scheme for large time steps (CFL>1)
  !!
  !! Calculation of time averaged vertical tracer fluxes or tracer edge 
  !! values using the third order PPM/PSM scheme. This scheme can handle 
  !! large time steps (i.e. CFL>1)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-01-14)
  !!
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201 (PPM)
  ! - Carpenter et al. (1989), MWR, 118, 586-612  (PPM)
  ! - Zerroukat et al. (2006), Int. J. Numer. Meth. Fluids, 51, 1297-1318 (PSM)
  ! - Lin et al (1994), MWR, 122, 1575-1593 (filtered reconstruction)
  ! - Lin and Rood (1996), MWR, 124, 2046-2070 (CFL-independent version)
  !
  SUBROUTINE upwind_vflux_ppm( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v,     &
    &                      p_dtime,  ld_compute, ld_cleanup, p_itype_vlimit,   &
    &                      p_ivlimit_selective,                                &
    &                      p_cellhgt_mc_now, p_cellmass_now,                   &
    &                      lprint_cfl, ivadv_tracer,                           &
    &                      p_upflux, opt_lout_edge, opt_topflx_tra, opt_slev,  &
    &                      opt_ti_slev, opt_rlstart, opt_rlend, opt_elev )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux:upwind_vflux_ppm'

    TYPE(t_patch), INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::  &    !< selects upper boundary condition
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

    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< NH: density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)    !< at time step n [kg/m**2]
                                  !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &    !< determines if vertical CFL number shall be written out
      &  lprint_cfl

    INTEGER, INTENT(IN) ::   &    !< type of vertical transport (PPM or PSM)
      &  ivadv_tracer

    REAL(wp), INTENT(INOUT) :: &  !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    LOGICAL, INTENT(IN), OPTIONAL ::  & !< optional: output edge value (.TRUE.),
      &  opt_lout_edge                  !< or the flux across the edge 
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
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    LOGICAL  :: l_out_edgeval     !< corresponding local variable; default 
                                  !< .FALSE. i.e. output flux across the edge

    REAL(wp) :: &                 !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< integer fluxes for w>0
      &  z_iflx_p(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< integer fluxes for w<0
      &  z_iflx_m(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< difference between upper and lower face value times 0.5
      &  z_delta_q(nproma,p_patch%nlev)
    REAL(wp) :: &                 !< 1/6 * a6,i (see Colella and Woodward (1984))
      &  z_a1(nproma,p_patch%nlev)

    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: ikm1, ikp1               !< vertical level minus and plus one
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
    INTEGER  :: i_rlstart, i_rlend

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

    INTEGER  :: jk_shift, jks
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

    REAL(wp) :: rdtime                  !< 1/dt


#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: z_face,z_face_up,z_face_low,z_iflx_p
!DIR$ ATTRIBUTES ALIGN : 64 :: z_iflx_m
!DIR$ ATTRIBUTES ALIGN : 64 :: z_cflfrac_m,max_cfl_blk,i_indlist_p
!DIR$ ATTRIBUTES ALIGN : 64 :: i_levlist_p,i_levlist_m,i_listdim_p
!DIR$ ATTRIBUTES ALIGN : 64 :: i_listdim_m,jk_int_p,jk_int_m
!DIR$ ATTRIBUTES ALIGN : 64 :: max_cfl,zparent_topflx,max_cfl_lay
!DIR$ ATTRIBUTES ALIGN : 64 :: z_aux_p,max_cfl_lay_tot,z_aux_m
#endif
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

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    !
    ! advection is done with an upwind scheme where a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    !
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

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jc,ik,ikm1,i_startidx,i_endidx,z_dummy,nlist_p,nlist_m,        &
!$OMP            counter_p,counter_m,counter_jip,counter_jim,max_cfl,                 &
!$OMP            z_aux_p,z_aux_m,ikp1,nlist,ji_p,                                     &
!$OMP            ji_m,jk_shift,jks,z_iflx_m,z_iflx_p,z_delta_q,z_a1,                  &
!$OMP            z_lext_1,z_lext_2,z_face,z_face_up,z_face_low,                       &
!$OMP            z_flx_frac_high) ICON_OMP_GUIDED_SCHEDULE
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
    !    for w<0 and w>0 and integer shift s
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
                &  .AND. p_mflx_contra_v(jc,jk,jb) > 0._wp ) THEN

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
                &  .AND. p_mflx_contra_v(jc,jk,jb) < 0._wp ) THEN

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
    ! 2. Edge value reconstruction
    !
    SELECT CASE(ivadv_tracer)
    CASE (IPPM_V)

      !
      ! PPM Reconstruction following Colella and Woodward (1984)
      !
      CALL compute_face_values_ppm( i_startidx       = i_startidx,               & !in
        &                           i_endidx         = i_endidx,                 & !in
        &                           slev             = slev,                     & !in
        &                           elev             = nlev,                     & !in
        &                           p_itype_vlimit   = p_itype_vlimit,           & !in
        &                           p_cc             = p_cc(:,:,jb),             & !in
        &                           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,jb), & !in
        &                           p_face           = z_face(:,:)               ) !inout


    CASE (IPSM_V)


      !
      ! PSM Reconstruction following Zerroukat et al. (2006)
      !
      CALL compute_face_values_psm( i_startidx       = i_startidx,               & !in
        &                           i_endidx         = i_endidx,                 & !in
        &                           slev             = slev,                     & !in
        &                           elev             = nlev,                     & !in
        &                           p_itype_vlimit   = p_itype_vlimit,           & !in
        &                           p_cc             = p_cc(:,:,jb),             & !in
        &                           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,jb), & !in
        &                           p_face           = z_face(:,:)               ) !inout


    END SELECT  ! ivadv_tracer


      !
      ! 4. Limitation/filtering of first guess parabola (which is based on z_face)
      ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
      ! the limitation procedure.
      ! Therefore 2 additional fields z_face_up and z_face_low are introduced.
      !
      SELECT CASE (p_itype_vlimit)
      CASE(ISLOPEL_VSM)

        ! semi-monotonic (sm) filter
        CALL v_limit_parabola_sm( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in

      CASE(ISLOPEL_VM)

        ! monotonic (mo) filter
        CALL v_limit_parabola_mo( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in


      CASE default

        ! simply copy face values to 'face_up' and 'face_low' arrays
        !
        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          z_face_up(i_startidx:i_endidx,jk)  = z_face(i_startidx:i_endidx,jk)
          z_face_low(i_startidx:i_endidx,jk) = z_face(i_startidx:i_endidx,ikp1)
        ENDDO

      END SELECT  ! p_itype_vlimit



      !
      ! 5. Computation of upwind fluxes. IF CFL > 1, the fluxes are the sum of
      !    integer-fluxes and a fractional flux. IF CFL <1 the fluxes are only
      !    comprised of the fractional flux. The fractional flux is calculated
      !    assuming a piecewise parabolic subgrid distribution.
      !

      ! 5a. Compute coefficients of reconstructed parabola as they are used at 
      !     various places below.
      !     Terminology follows Colella (1984)
      !     z_delta_q = 0.5*\Delta q
      !     z_a1 = 1/6*a_6
      !
      DO jk = slev, nlev
        DO jc = i_startidx, i_endidx
          z_delta_q(jc,jk) = 0.5_wp * (z_face_up(jc,jk) - z_face_low(jc,jk))
          z_a1(jc,jk)      = p_cc(jc,jk,jb) - 0.5_wp*(z_face_up(jc,jk) + z_face_low(jc,jk))
        ENDDO
      ENDDO

      !
      ! 5b. First compute fluxes for the CFL<1 case for all grid points
      ! On the grid points where CFL>1, they will be overwritten afterwards 
      ! This part has been adopted from the restricted time step PPM-scheme.
      !
      DO jk = slevp1, elev

        ik   = MIN(jk,nlev)
        ! index of top half level
        ikm1 = jk -1

#ifdef __INTEL_COMPILER
! HB: for some strange reason this loop introduces a decomposition dependency if
! vectorized... threfore inhibit vectorization here
!DIR$ NOVECTOR
#endif
        DO jc = i_startidx, i_endidx

          ! if w < 0 (physical downwelling)
          !
          z_lext_1 = p_cc(jc,ikm1,jb)                                   &
            &  - (z_delta_q(jc,ikm1) * (1._wp - z_cflfrac_m(jc,jk,jb))) &
            &  - z_a1(jc,ikm1)*(1._wp - 3._wp*z_cflfrac_m(jc,jk,jb)     &
            &  + 2._wp*z_cflfrac_m(jc,jk,jb)*z_cflfrac_m(jc,jk,jb))


          ! if w > 0 (physical upwelling)
          !
          z_lext_2 = p_cc(jc,ik,jb)                                     &
            &  + (z_delta_q(jc,ik) * (1._wp - z_cflfrac_p(jc,ik,jb)))   &
            &  - z_a1(jc,ik)*(1._wp - 3._wp*z_cflfrac_p(jc,ik,jb)       &
            &  + 2._wp*z_cflfrac_p(jc,ik,jb)*z_cflfrac_p(jc,ik,jb))

          !
          ! full flux  -- removed flaky laxfr macro
          !
          p_upflux(jc,jk,jb) = p_mflx_contra_v(jc,jk,jb)*               &
                               MERGE(z_lext_2,z_lext_1,p_mflx_contra_v(jc,jk,jb) >= 0.0_wp)

        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels


      !
      ! 5c. Now execute the special computations needed for CFL>1:
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
          ! integer fluxes for w>0
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
            z_iflx_p(jc,jk) = z_iflx_p(jc,jk) + p_cc(jc,jk_shift,jb) &
                             * p_cellmass_now(jc,jk_shift,jb)

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

            jks = jk_int_p(jc,jk,jb)

            ! fractional upward flux
            ! if w > 0 (physical upwelling)

            ! fractional high order flux   
            z_flx_frac_high = ( p_cellmass_now(jc,jks,jb)                           &
              &         * z_cflfrac_p(jc,jk,jb) *( p_cc(jc,jks,jb)                  &
              &         + (z_delta_q(jc,jks) * (1._wp - z_cflfrac_p(jc,jk,jb)))     &
              &         - z_a1(jc,jks)*(1._wp - 3._wp*z_cflfrac_p(jc,jk,jb)         &
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
          ! integer fluxes for w<0
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
            z_iflx_m(jc,jk) = z_iflx_m(jc,jk) - p_cc(jc,jk_shift,jb) &
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

            jks = jk_int_m(jc,jk,jb)

            ! this is needed in addition in order to avoid accessing non-existing (uninitalized)
            ! source levels for tracers that are not advected on all model levels
            IF (jks < slev) CYCLE

            ! fractional downward flux
            ! if w < 0 (physical downwelling)

            ! fractional high order flux           
            z_flx_frac_high = ( -1._wp * p_cellmass_now(jc,jks,jb)                  &
              &         * z_cflfrac_m(jc,jk,jb) * ( p_cc(jc,jks,jb)                 &
              &         - (z_delta_q(jc,jks) * (1._wp - z_cflfrac_m(jc,jk,jb)))     &
              &         - z_a1(jc,jks)*(1._wp - 3._wp*z_cflfrac_m(jc,jk,jb)         &
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
      CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
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

      !
      ! 6. If desired, apply positive-definite flux limiter to limit 
      !    computed fluxes (based on work by Zalesak (1979)).
      !
      IF (p_itype_vlimit == IFLUXL_VPD) THEN
        ! positive-definite (pd) flux limiter
        CALL vflx_limiter_pd( p_dtime         = p_dtime,                & !in
          &                   p_cc            = p_cc(:,:,jb),           & !in
          &                   p_rhodz_now     = p_cellmass_now(:,:,jb), & !in
          &                   p_mflx_tracer_v = p_upflux(:,:,jb),       & !inout
          &                   i_startidx      = i_startidx,             & !in
          &                   i_endidx        = i_endidx,               & !in
          &                   slev            = slev,                   & !in
          &                   elev            = elev_lim                ) !in
      ENDIF

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL



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
      IF (msg_level >= 13 .AND. max_cfl_tot > (nlist_max-1)) THEN
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

  END SUBROUTINE upwind_vflux_ppm





  !-------------------------------------------------------------------------
  !>
  !! The third order PPM scheme for large time steps (CFL>1)
  !! GPU-enabled version without index lists.
  !!
  !! Calculation of time averaged vertical tracer fluxes or tracer edge 
  !! values using the third order PPM scheme. Includes extension to large 
  !! time steps, i.e. this scheme can handle large time steps (CFL>1).
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-01-14)
  !! Modification by Daniel Reinert, DWD (2018-01-27)
  !! - optimized for GPU
  !!
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201
  ! - Carpenter et al. (1989), MWR, 118, 586-612
  ! - Lin et al (1994), MWR, 122, 1575-1593 (slope limiter)
  ! - Lin and Rood (1996), MWR, 124, 2046-2070 (CFL-independent version)
  !
  SUBROUTINE upwind_vflux_ppm4gpu( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v, &
    &                      p_dtime,  ld_compute, ld_cleanup, p_itype_vlimit,   &
    &                      p_ivlimit_selective,                                &
    &                      p_cellhgt_mc_now, p_cellmass_now, lprint_cfl,       &
    &                      p_upflux, opt_lout_edge, opt_topflx_tra, opt_slev,  &
    &                      opt_ti_slev, opt_rlstart, opt_rlend, opt_elev )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_advection_vflux:upwind_vflux_ppm4gpu'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(INOUT) ::  & !< contravariant vertical mass flux [kg/m**2/s]
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< time step [s]
      &  p_dtime

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. compute geometric terms
      &  ld_compute

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. clean up geometric terms
      &  ld_cleanup

    INTEGER, INTENT(IN)  ::  &    !< parameter to select the limiter for
      &  p_itype_vlimit           !< vertical transport

    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n [m]
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)    !< at time step n [kg/m**2]
                                  !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &    !< determines if vertical CFL number shall be written out
      &  lprint_cfl

    REAL(wp), INTENT(INOUT) :: &  !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    LOGICAL, INTENT(IN), OPTIONAL ::  & !< optional: output edge value (.TRUE.),
      &  opt_lout_edge                  !< or the flux across the edge 
                                        !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:)            !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level (for sedimentation)
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)



    ! local vars

    REAL(wp), ALLOCATABLE, SAVE  ::   &  !< sum of integer and fractional Courant number 
      &  z_cfl(:,:,:)                    !< Sign equals sign of w.

    REAL(wp) :: &                        !< absolute value of fractional Courant number
      &  z_cflfrac                       !< i.e. always >=0

    REAL(wp) :: &                        !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1)

    REAL(wp) :: &                        !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< monotonized slope
      &  z_slope(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< difference between upper and lower face value times 0.5
      &  z_delta_q(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< 1/6 * a6,i (see Colella and Woodward (1984))
      &  z_a1(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< mass crossing cell face during \Delta t [kg/m**2]
      &  z_mass                          !< can be positive, or negative, depending on the sign of w

    REAL(wp) :: &                        !< integrated value of q (from 0 to z_cflfrac)
      &  z_q_int

    REAL(wp) :: &                        !< integer flux for w>0 or w<0  [kg/m**2/s]
      &  z_iflx

    REAL(wp) :: p_cc_min, p_cc_max       !< 3-point max/min values

    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: ikm1, ikp1, ikp2         !< vertical level minus and plus one, plus two
    INTEGER  :: slev, slevp1             !< vertical start level and start level +1
    INTEGER  :: slev_ti, slevp1_ti       !< vertical start level (+1)  (tracer independent part)
    INTEGER  :: nlev, nlevp1             !< number of full and half levels

    INTEGER  :: jk_shift, jks            !< shifted vertical index

    INTEGER  :: js                       !< the shift itself (always positive), i.e. jks = jk \pm js

    ! JF: for treatment of sedimentation
    INTEGER  :: elev, elev_lim           !< vertical end level
    LOGICAL  :: llbc_adv                 !< apply lower boundary condition?

    LOGICAL  :: l_out_edgeval            !< corresponding local variable; default 
                                         !< .FALSE. i.e. output flux across the edge

    INTEGER  :: ist                      !< status variable
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend

    INTEGER  :: n                        !< loop index

    REAL(wp) :: wsign                    !< wind direction: introduced, in order to merge flux formula  
                                         !< for w>0 and w<0 into one.
                                         !< +1, if w >0
                                         !< -1, if w <0

    REAL(wp) ::   &                      !< absolute CFL at grid point
      &  abs_cfl(nproma)

    REAL(wp) ::   &                      !< maximum CFL for each layer (blockwise)
      &  max_cfl_lay(p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) ::   &                      !< maximum CFL for each layer  
      &  max_cfl_lay_tot(p_patch%nlevp1)

    REAL(wp) ::    &                     !< maximum vertical Courant number per block
      &  max_cfl_blk(p_patch%nblks_c)

    REAL(wp) ::   &                      !< domain-wide maximum CFL
      &  max_cfl_tot

    REAL(wp) ::  &                              !< necessary, to make this routine
      &  zparent_topflx(nproma,p_patch%nblks_c) !< compatible to the hydrost. core 

    REAL(wp) ::   &                      !< auxiliaries for optimization
      &  zfac, zfac_m1, zgeo1, zgeo2, zgeo3, zgeo4

    REAL(wp) :: rdtime                   !< 1/dt


    !-----------------------------------------------------------------------

    ! inverse of time step for computational efficiency
    rdtime = 1._wp/p_dtime

!$ACC UPDATE DEVICE( p_cc, p_cellhgt_mc_now, p_cellmass_now, p_mflx_contra_v, p_upflux ), &
!$ACC        IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

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

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    IF ( ld_compute ) THEN
      !
      ! allocate field for storing the density weighted Courant number 
      !
      ALLOCATE( z_cfl(nproma,nlevp1,p_patch%nblks_c), STAT=ist  )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'allocation for z_cfl failed')
      ENDIF
!$ACC ENTER DATA CREATE( z_cfl ),  IF( i_am_accel_node .AND. acc_on )
    END IF

!$ACC DATA CREATE( z_face, z_face_up, z_face_low, z_slope, z_delta_q, z_a1 ), &
!$ACC      PCOPYIN( p_cc, p_cellhgt_mc_now, p_cellmass_now ), PCOPY( p_mflx_contra_v, p_upflux ), &
!$ACC      IF( i_am_accel_node .AND. acc_on )


!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jc,ikm1,i_startidx,i_endidx,ikp1,ikp2,  &
!$OMP            jks,z_mass,p_cc_min,p_cc_max,jk_shift,js,n,   &
!$OMP            z_iflx,z_delta_q,z_a1,zfac,zfac_m1,           &
!$OMP            zgeo1,zgeo2,zgeo3,zgeo4,z_q_int,z_slope,      &
!$OMP            wsign,z_cflfrac,z_face,z_face_up,z_face_low  ) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      ! The contravariant mass flux should never exactly vanish
      !
      IF (l_out_edgeval) THEN
!$ACC PARALLEL  DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slevp1, elev
          DO jc = i_startidx, i_endidx
            p_mflx_contra_v(jc,jk,jb) =                            &
          &              p_mflx_contra_v(jc,jk,jb)               &
          &              + SIGN(dbl_eps,p_mflx_contra_v(jc,jk,jb))
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDIF

      !
      ! 1. Compute density weighted Courant number for w<0 and w>0. 
      !    It is the sum of the fractional Courant number and the integer shift s.
      !    Stored at cell faces
      !
      IF (ld_compute) THEN

        ! initialize Courant number
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
        z_cfl(i_startidx:i_endidx,slev_ti:nlevp1,jb) = 0._wp
!$ACC END KERNELS


        ! Split density-weighted Courant number into integer and fractional 
        ! part and store the sum in z_cfl (for w>0 and w<0)
        !
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( z_mass, jks ) COLLAPSE(2)
        DO jk = slevp1_ti, elev
          DO jc = i_startidx, i_endidx

            ! total mass crossing jk'th edge during \Delta t
            z_mass = p_dtime*p_mflx_contra_v(jc,jk,jb)

            !
            ! Case: w > 0
            !
            IF (z_mass > 0._wp) THEN

              jks = jk   ! initialize shifted index

              DO WHILE( (z_mass > p_cellmass_now(jc,jks,jb)) .AND. (jks <= nlev-1) )
                z_mass = z_mass - p_cellmass_now(jc,jks,jb)
                jks = jks+1
                ! update Courant number
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + 1._wp
              ENDDO

              ! now we add the fractional Courant number
              ! The MIN function is required here for the case that 
              ! we approach the lower boundary and exit the above loop 
              ! because of jks > nlev-1.
              z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + MIN(1._wp,z_mass/p_cellmass_now(jc,jks,jb))

            ELSE
            !
            ! Case w < 0
            !
              jks = jk-1   ! initialize shifted index

              DO WHILE( (ABS(z_mass) > p_cellmass_now(jc,jks,jb)) .AND. &
                &       (jks >= slevp1_ti) )
                z_mass = z_mass + p_cellmass_now(jc,jks,jb)
                jks = jks-1
                ! update Courant number
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) - 1._wp
              ENDDO

              ! now we add the fractional Courant number
              ! The MAX function is required here for the case that 
              ! we approach the upper boundary and exit the above loop 
              ! because of jks < slevp1_ti.
              z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + MAX(-1._wp,z_mass/p_cellmass_now(jc,jks,jb))

            ENDIF

          ENDDO  ! jc
        ENDDO  ! jk
!$ACC END PARALLEL

      END IF ! ld_compute



      !
      ! 2. Compute monotonized slope
      !

      ! Initialize z_slope for jk=slev
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
      z_slope(i_startidx:i_endidx,slev) = 0._wp
!$ACC END KERNELS

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikm1, ikp1, zfac_m1, zfac, p_cc_min, p_cc_max ) COLLAPSE(2)
      DO jk = slevp1, nlev
        DO jc = i_startidx, i_endidx
          ! index of top half level
          ikm1    = jk - 1
          ! index of bottom half level
          ikp1    = MIN( jk+1, nlev )

          zfac_m1 = 1._wp / (p_cellhgt_mc_now(jc,jk,jb) + p_cellhgt_mc_now(jc,ikm1,jb)) &
            &  * (p_cc(jc,jk,jb) - p_cc(jc,ikm1,jb))

          zfac = 1._wp / (p_cellhgt_mc_now(jc,ikp1,jb) + p_cellhgt_mc_now(jc,jk,jb)) &
            &  * (p_cc(jc,ikp1,jb) - p_cc(jc,jk,jb))

          z_slope(jc,jk) = ( p_cellhgt_mc_now(jc,jk,jb)                                          &
            &  / (p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)                      &
            &  + p_cellhgt_mc_now(jc,ikp1,jb)) )                                                 &
            &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1,jb) + p_cellhgt_mc_now(jc,jk,jb)) * zfac    &
            &  + (p_cellhgt_mc_now(jc,jk,jb) + 2._wp * p_cellhgt_mc_now(jc,ikp1,jb)) * zfac_m1 )

          ! equivalent formulation of Colella and Woodward (1984) slope limiter 
          ! following Lin et al (1994).
          p_cc_min = MIN(p_cc(jc,ikm1,jb),p_cc(jc,jk,jb),p_cc(jc,ikp1,jb))
          p_cc_max = MAX(p_cc(jc,ikm1,jb),p_cc(jc,jk,jb),p_cc(jc,ikp1,jb))
          z_slope(jc,jk) = SIGN(                                            &
            &  MIN( ABS(z_slope(jc,jk)), 2._wp*(p_cc(jc,jk,jb)-p_cc_min),   &
            &                            2._wp*(p_cc_max-p_cc(jc,jk,jb)) ), &
            &    z_slope(jc,jk))

        END DO  ! jc
      END DO  ! jk

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
!$ACC LOOP GANG VECTOR
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
!$ACC END PARALLEL

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikm1, ikp1, ikp2, zgeo1, zgeo2, zgeo3, zgeo4 ) COLLAPSE(2)
      DO jk = slevp1, nlev-2
        DO jc = i_startidx, i_endidx

          ! index of top half level
          ikm1 = jk - 1
          ! index of bottom half level
          ikp1 = jk + 1
          ikp2 = jk + 2

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
!$ACC END PARALLEL


      !
      ! 4. Limitation of first guess parabola (which is based on z_face)
      ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
      ! the limitation procedure.
      ! Therefore 2 additional fields z_face_up and z_face_low are
      ! introduced.
      !
      IF (p_itype_vlimit == ISLOPEL_VSM) THEN
        ! semi-monotonic (sm) filter
        CALL v_limit_parabola_sm( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in
      ELSE IF (p_itype_vlimit == ISLOPEL_VM) THEN
        ! monotonic (mo) filter
        CALL v_limit_parabola_mo( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in
      ENDIF


      IF (p_itype_vlimit /= islopel_vsm .AND. p_itype_vlimit /= islopel_vm) THEN
        ! simply copy face values to 'face_up' and 'face_low' arrays

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikp1 ) COLLAPSE(2)
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ! index of bottom half level
            ikp1 = jk + 1
            z_face_up(jc,jk)  = z_face(jc,jk)
            z_face_low(jc,jk) = z_face(jc,ikp1)
          ENDDO
        ENDDO
!$ACC END PARALLEL

      ENDIF


      !
      ! 5. calculation of upwind fluxes. For CFL>1, the total flux is the sum of
      !    integer-fluxes and a fractional flux. IF CFL<=1 the fluxes are only
      !    comprised of the fractional flux. The fractional flux is calculated
      !    by assuming a piecewise parabolic approx. for the subgrid distribution.
      !

      ! 5a. Compute coefficients of reconstructed parabola as they are used below.
      !     Terminology follows Colella (1984)
      !     z_delta_q = 0.5*\Delta q
      !     z_a1 = 1/6*a_6
      !
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, nlev
        DO jc = i_startidx, i_endidx
          z_delta_q(jc,jk) = 0.5_wp * (z_face_up(jc,jk) - z_face_low(jc,jk))
          z_a1(jc,jk)      = p_cc(jc,jk,jb) - 0.5_wp*(z_face_up(jc,jk) + z_face_low(jc,jk))
        ENDDO
      ENDDO
!$ACC END PARALLEL


      !
      ! 5b. First compute the fractional fluxes for all cell faces.
      !     For cell faces with CFL>1, integer fluxes will be added lateron.
      !
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikm1, js, z_cflfrac, jks, wsign, z_q_int ) COLLAPSE(2)
      DO jk = slevp1, elev

        DO jc = i_startidx, i_endidx
          ikm1 = jk-1
          ! get integer shift (always non-negative)
          js = FLOOR(ABS(z_cfl(jc,jk,jb)))

          ! get fractional part of Courant number (always non-negative)
          z_cflfrac = ABS(z_cfl(jc,jk,jb)) - REAL(js,wp)

          ! compute shifted cell index
          IF (z_cfl(jc,jk,jb) > 0._wp) THEN
            jks = MIN(jk,nlev)+js
            wsign = 1._wp
          ELSE
            jks = ikm1-js
            wsign = -1._wp
          ENDIF

          ! this is needed in addition in order to avoid accessing non-existing (uninitalized)
          ! source levels for tracers that are not advected on all model levels
          IF (jks < slev) CYCLE

          ! compute flux
          !
          ! flux formula differs between w>0 and w<0. 
          ! By using the coefficient 'wsign' we are able to merge 
          ! the two formula into one.
          z_q_int = p_cc(jc,jks,jb)                                 &
            &     + wsign*(z_delta_q(jc,jks) * (1._wp - z_cflfrac)) &
            &     - z_a1(jc,jks)*(1._wp - 3._wp*z_cflfrac + 2._wp*z_cflfrac*z_cflfrac)

          p_upflux(jc,jk,jb) = wsign * p_cellmass_now(jc,jks,jb)    &
            &                * z_cflfrac * z_q_int * rdtime
        ENDDO

      ENDDO ! end loop over vertical levels
!$ACC END PARALLEL


      !
      ! 5c. Now compute the integer fluxes and add them to the fractional flux
      !
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( js, z_iflx, jk_shift ) COLLAPSE(2)
      DO jk = slevp1, elev

        DO jc = i_startidx, i_endidx

          ! get integer shift (always non-negative)
          js = FLOOR(ABS(z_cfl(jc,jk,jb)))

          IF (js == 0) CYCLE   ! no work to do

          z_iflx = 0._wp

          ! case w > 0
          IF (z_cfl(jc,jk,jb) > 0._wp) THEN

            DO n = 1, js
              jk_shift = jk-1 + n
              ! Integer flux (division by p_dtime is done at the end)
              z_iflx = z_iflx + p_cc(jc,jk_shift,jb) * p_cellmass_now(jc,jk_shift,jb)
            ENDDO

          ! case w <= 0
          ELSE

            DO n = 1, js
              jk_shift = jk - n

              ! cycle if the source model level is in a region where advection is 
              ! turned off for the present variable
              IF (jk_shift < slevp1) CYCLE

              ! Integer flux (division by p_dtime is done at the end)
              z_iflx = z_iflx - p_cc(jc,jk_shift,jb) * p_cellmass_now(jc,jk_shift,jb)
            ENDDO

          ENDIF

          ! compute full (integer- plus high order fractional) flux
          p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb) + z_iflx*rdtime
        ENDDO  ! jc
      ENDDO ! jk
!$ACC END PARALLEL


      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(p_upflux(:,slev+1,jb),            &! in
        &              p_mflx_contra_v(:,slev+1,jb),     &! in
        &              p_mflx_contra_v(:,slev  ,jb),     &! in
        &              p_iubc_adv, i_startidx, i_endidx, &! in
        &              zparent_topflx(:,jb),             &! in
        &              p_upflux(:,slev,jb),              &! out
        &              p_upflux(:,nlevp1,jb), llbc_adv)   ! out


      ! If desired, get edge value of advected quantity 
      IF ( l_out_edgeval ) THEN

!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slevp1, nlev
          DO jc = i_startidx, i_endidx
            p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb)/p_mflx_contra_v(jc,jk,jb)
          ENDDO
        ENDDO
!$ACC END PARALLEL

      ENDIF


      !
      ! 6. If desired, apply positive-definite flux limiter to limit 
      !    computed fluxes (based on work by Zalesak (1979)).
      !
      IF (p_itype_vlimit == IFLUXL_VPD) THEN
        ! positive-definite (pd) flux limiter
        CALL vflx_limiter_pd( p_dtime         = p_dtime,                & !in
          &                   p_cc            = p_cc(:,:,jb),           & !in
          &                   p_rhodz_now     = p_cellmass_now(:,:,jb), & !in
          &                   p_mflx_tracer_v = p_upflux(:,:,jb),       & !inout
          &                   i_startidx      = i_startidx,             & !in
          &                   i_endidx        = i_endidx,               & !in
          &                   slev            = slev,                   & !in
          &                   elev            = elev_lim                ) !in
      ENDIF

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL




#ifndef _OPENACC
! These diagnostics are not viable on GPU
    !
    ! If desired, print maximum vertical CFL number
    !
    IF ( ld_compute .AND. msg_level >= 10 .AND. lprint_cfl ) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,abs_cfl)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        DO jk = slevp1_ti, nlev
          DO jc = i_startidx, i_endidx
            abs_cfl(jc) = ABS(z_cfl(jc,jk,jb))
          ENDDO  ! jc
          max_cfl_lay(jk,jb) = MAXVAL(abs_cfl(i_startidx:i_endidx))
        ENDDO  ! jk
        !
        max_cfl_blk(jb) = MAXVAL(max_cfl_lay(slevp1_ti:elev,jb))
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

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

!$ACC END DATA

    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays
!$ACC EXIT DATA DELETE( z_cfl ) IF( i_am_accel_node .AND. acc_on )
      DEALLOCATE( z_cfl, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish ( TRIM(routine), 'deallocation for z_cfl failed' )
      ENDIF
    END IF

!$ACC UPDATE HOST( p_mflx_contra_v, p_upflux ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

  END SUBROUTINE upwind_vflux_ppm4gpu



  !---------------------------------------------------------------
  !>
  !! Description:
  !!   solve the vertical flux advection equation for sedimentation
  !!   of scalar variables (a purely downward directed transport)
  !!   with the 2-point implicit scheme described in
  !!   COSMO Sci. doc. II, section 5.2.4.
  !!
  !! Method:
  !!   index convention for the sedimentation velocity :
  !!   v_new(i,j,k) = v(i,j,k-1/2)
  !!   sign: v_new, v_old > 0 ! (i.e. directed downward)
  !!
  !!   negative values in phi_new are clipped; this destroys
  !!   mass conservation.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2018-11-07)
  !
  SUBROUTINE implicit_sedim_tracer( tracer,                    &
    &                        rho, rho_inv,                     &
    &                        v_new, v_old,                     &
    &                        dt,                               &
    &                        p_patch, p_metrics,               &
    &                        i_rlstart, i_rlend,               &
    &                        rhoS )

    USE mo_parallel_config,         ONLY: nproma
    USE mo_nonhydro_types ,         ONLY: t_nh_metrics

    IMPLICIT NONE

    REAL (wp),     INTENT(INOUT) :: tracer(:,:,:) !advected cell centered variable
                                                       !< dim: (nproma, nlev, nblks_c)
    REAL (wp),     INTENT(IN)  :: rho    (:,:,:)  ! mass density  (in kg/m^3)
    REAL (wp),     INTENT(IN)  :: rho_inv(:,:,:)  ! 1/rho  (in m^3/kg)
    REAL (wp),     INTENT(IN)  :: v_new  (:,:,:)  ! sedimentation velocity at time level n+1 (in m/s)
    REAL (wp),     INTENT(IN)  :: v_old  (:,:,:)  ! sedimentation velocity at time level n   (in m/s)
    REAL (wp),     INTENT(IN)  :: dt     ! time step

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch      !< Patch on which computation is performed
    TYPE(t_nh_metrics), TARGET,  INTENT(IN) :: p_metrics    !< Metrical fields

    INTEGER, INTENT(IN) :: i_rlstart, i_rlend

    REAL (wp), INTENT(IN), OPTIONAL :: rhoS(:,:,:)

    INTEGER    :: jk, jb, jc
    REAL (wp)  :: h, dz, c, lambda_im

    REAL (wp) :: phi_old( nproma, p_patch%nlev, p_patch%nblks_c)
    REAL (wp) :: phi_new( nproma, p_patch%nlev, p_patch%nblks_c)

    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx

    LOGICAL :: is_rhoS_present

    is_rhoS_present = PRESENT( rhoS )   ! own variable may help in vectorisation for some compilers

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

#ifdef _OPENACC
    PRINT *, "Sorry: implicit_sedim_tracer not yet available for OpenACC"
    IF ( .FALSE. ) THEN
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,dz,c,lambda_im,h) ICON_OMP_GUIDED_SCHEDULE
#endif
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,        &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! calculate densities (for the following flux advection scheme)
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx
          phi_old(jc,jk,jb) = tracer(jc,jk,jb) * rho(jc,jk,jb)
        ENDDO ! jc
      ENDDO ! jk

      IF ( is_rhoS_present ) THEN

        ! top level
        jk = 1
        DO jc = i_startidx, i_endidx

          dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
          c = dt / ( 2.0_wp * dz );
          lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

          h = phi_old(jc,jk,jb) - c * ( v_old(jc,jk+1,jb) * phi_old(jc,jk,jb) )

          phi_new(jc,jk,jb) = MAX( lambda_im * ( h + rhoS(jc,jk,jb)*dt ), 0.0_wp)
        END DO ! jc

        DO jk=2, p_patch%nlev
          DO jc = i_startidx, i_endidx

            dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
            c = dt / ( 2.0_wp * dz );
            lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

            h = phi_old(jc,jk,jb) + c *                      &
              &  ( v_new(jc,jk  ,jb) * phi_new(jc,jk-1,jb)   &
              &  + v_old(jc,jk  ,jb) * phi_old(jc,jk-1,jb)   &
              &  - v_old(jc,jk+1,jb) * phi_old(jc,jk  ,jb)  )

            phi_new(jc,jk,jb) = MAX( lambda_im * ( h + rhoS(jc,jk,jb)*dt ), 0.0_wp)

          END DO ! jc
        END DO ! jk

      ELSE

        ! the same code as in the if-block before but without the source term rhoS:

        ! top level
        jk = 1
        DO jc = i_startidx, i_endidx

          dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
          c = dt / ( 2.0_wp * dz );
          lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

          h = phi_old(jc,jk,jb) - c * ( v_old(jc,jk+1,jb) * phi_old(jc,jk,jb) )

          phi_new(jc,jk,jb) = MAX( lambda_im * h, 0.0_wp)
        END DO ! jc

        DO jk=2, p_patch%nlev
          DO jc = i_startidx, i_endidx

            dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
            c = dt / ( 2.0_wp * dz );
            lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

            h = phi_old(jc,jk,jb) + c *                      &
              &  ( v_new(jc,jk  ,jb) * phi_new(jc,jk-1,jb)   &
              &  + v_old(jc,jk  ,jb) * phi_old(jc,jk-1,jb)   &
              &  - v_old(jc,jk+1,jb) * phi_old(jc,jk  ,jb)  )

            phi_new(jc,jk,jb) = MAX( lambda_im * h, 0.0_wp)
          END DO ! jc
        END DO ! jk

      END IF       ! IF ( is_rhoS_present )

      ! calculate back the specific mass:
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx
          tracer(jc,jk,jb) = phi_new(jc,jk,jb) * rho_inv(jc,jk,jb)
        ENDDO ! jc
      ENDDO ! jk

    END DO   ! jb
#ifdef _OPENACC
    END IF
#else
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

  END SUBROUTINE implicit_sedim_tracer



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
  SUBROUTINE set_bc_vadv(upflx_top_p1, mflx_top_p1, mflx_top, iubc_adv, &
    &                    i_start, i_end, parent_topflx, upflx_top,      &
    &                    upflx_bottom, llbc_adv )

!!$    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
!!$      &  routine = 'mo_advection_vflux: set_ubc_adv'

    REAL(wp), INTENT(IN)     :: & !< computed tracer flux at second half level
      &  upflx_top_p1(:)
    REAL(wp), INTENT(IN)     :: & !< mass flux at second half level
      &  mflx_top_p1(:)
    REAL(wp), INTENT(IN)     :: & !< mass flux at upper boundary
      &  mflx_top(:)
    INTEGER, INTENT(IN)      :: & !< selects upper boundary condition
      &  iubc_adv
    INTEGER, INTENT(IN)      :: & !< start and end index
      &  i_start, i_end
    REAL(wp), INTENT(IN)     :: & !< tracer flux at upper boundary, 
      &  parent_topflx(:)         !< interpolated from parent grid
    REAL(wp), INTENT(OUT)    :: & !< upper boundary condition
      &  upflx_top(:)
    REAL(wp), INTENT(OUT)    :: & !< lower boundary condition
      &  upflx_bottom(:)
    LOGICAL, INTENT(IN)      :: & !< apply lower boundary condition?
      &  llbc_adv

    !-------------------------------------------------------------------------

!$ACC DATA PCOPYIN( upflx_top_p1, mflx_top_p1, mflx_top, parent_topflx ), &
!$ACC      PCOPYOUT( upflx_top, upflx_bottom ), IF( i_am_accel_node .AND. acc_on )

    ! 
    ! flux at top boundary
    ! 
    SELECT CASE (iubc_adv)
      CASE ( ino_flx )     ! no flux
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
        upflx_top(i_start:i_end) = 0._wp
!$ACC END KERNELS
 
      CASE ( izero_grad )  ! zero gradient
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
        upflx_top(i_start:i_end) = upflx_top_p1(i_start:i_end)      &
            &           * mflx_top(i_start:i_end)                   &
            &           / ( mflx_top_p1(i_start:i_end)              &
            &           + SIGN(dbl_eps, mflx_top_p1(i_start:i_end)))
!$ACC END KERNELS

      CASE ( iparent_flx ) ! interpolated flux from parent grid
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
        upflx_top(i_start:i_end) = parent_topflx(i_start:i_end)
!$ACC END KERNELS
    END SELECT

    !
    ! flux at bottom boundary
    !
    IF ( llbc_adv ) THEN
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
      upflx_bottom(i_start:i_end) = 0._wp
!$ACC END KERNELS
    END IF

!$ACC END DATA

  END SUBROUTINE set_bc_vadv


  !-------------------------------------------------------------------------
  !>
  !! PSM Face value reconstruction after Zerroukat et al (2006)
  !!
  !! PSM Face value reconstruction after Zerroukat et al (2006)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-04-12)
  !!
  !
  SUBROUTINE compute_face_values_psm( i_startidx, i_endidx, slev, elev, &
    &                                 p_itype_vlimit, p_cc, p_cellhgt_mc_now, p_face )

    INTEGER, INTENT(IN) ::    &    !< horizontal start and end indices
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN) ::    &    !< vertical start and end levels
      &  slev, elev

    INTEGER, INTENT(IN) ::    &    !< selects the vertical limiter
      &  p_itype_vlimit

    REAL(wp), INTENT(IN) ::   &    !< cell centered variable (cell average)
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::   &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:)     !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT):: &    !< face values of transported field
      &  p_face(:,:)               !< dim: (nproma,nlev)


    ! local variables
    INTEGER :: jc, jk
    INTEGER :: ikp1                !< vertical level plus one
    INTEGER :: elevp1              !< end level + 1

    ! TDMA arrays
    REAL(wp) :: a(SIZE(p_face,1),SIZE(p_face,2))   !< sub-diagonal
    REAL(wp) :: b(SIZE(p_face,1),SIZE(p_face,2))   !< main diagonal
    REAL(wp) :: c(SIZE(p_face,1),SIZE(p_face,2))   !< super diagonal
    REAL(wp) :: rhs(SIZE(p_face,1),SIZE(p_face,2)) !< right hand side
    REAL(wp) :: dzfrac                             !< ratio of neighboring cell heights

    !-------------------------------------------------------------------------

    elevp1 = elev + 1

!$ACC DATA CREATE( a, b, c, rhs), PCOPYIN(p_cc, p_cellhgt_mc_now), &
!$ACC PCOPYOUT( p_face ), IF( i_am_accel_node .AND. acc_on )

    !
    ! 1. reconstruct face values at vertical half-levels using splines
    !
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
    ! top BC
    a  (i_startidx:i_endidx,slev) = 0._wp
    b  (i_startidx:i_endidx,slev) = 2._wp
    c  (i_startidx:i_endidx,slev) = 1._wp
    rhs(i_startidx:i_endidx,slev) = 3._wp*p_cc(i_startidx:i_endidx,slev)
    !
    ! bottom BC
    a  (i_startidx:i_endidx,elevp1) = 1._wp
    b  (i_startidx:i_endidx,elevp1) = 2._wp
    c  (i_startidx:i_endidx,elevp1) = 0._wp
    rhs(i_startidx:i_endidx,elevp1) = 3._wp*p_cc(i_startidx:i_endidx,elev)
!$ACC END KERNELS
    !
!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikp1, dzfrac ) COLLAPSE(2)
    DO jk=slev,elev-1
      DO jc = i_startidx, i_endidx
        ikp1  = jk+1 
        dzfrac = p_cellhgt_mc_now(jc,jk)/p_cellhgt_mc_now(jc,ikp1)
        a  (jc,ikp1) = 1._wp
        c  (jc,ikp1) = dzfrac
        b  (jc,ikp1) = 2._wp*(1._wp + dzfrac)
        rhs(jc,ikp1) = 3._wp*(dzfrac*p_cc(jc,ikp1) + p_cc(jc,jk))
      ENDDO
    ENDDO
!$ACC END PARALLEL

    ! solve tri-diagonal system
    CALL tdma_solver_vec(a       = a,          &
      &                  b       = b,          &
      &                  c       = c,          &
      &                  d       = rhs,        &
      &                  slev    = slev,       &
      &                  elev    = elevp1,     &
      &                  startidx= i_startidx, &
      &                  endidx  = i_endidx,   &
      &                  varout  = p_face      )  ! out



    ! 2. OPTIONAL: Limit face values
    !
    SELECT CASE (p_itype_vlimit)
    CASE(ISLOPEL_VSM)

      ! make sure that PSM face values lie within the range of values 
      ! in the neighbouring cells
      CALL v_limit_face_sm( p_cc       = p_cc(:,:),   & ! in
        &                   p_face     = p_face(:,:), & ! inout
        &                   i_startidx = i_startidx,  & ! in
        &                   i_endidx   = i_endidx,    & ! in
        &                   slev       = slev,        & ! in
        &                   elev       = elev-1       ) ! in

      ! top and bottom face
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
      p_face(i_startidx:i_endidx,slev)   = MAX(p_cc(i_startidx:i_endidx,slev),p_face(i_startidx:i_endidx,slev))
      p_face(i_startidx:i_endidx,elevp1) = MAX(p_cc(i_startidx:i_endidx,elev),p_face(i_startidx:i_endidx,elevp1))
!$ACC END KERNELS

    CASE(ISLOPEL_VM)

      ! make sure that PSM face values lie within the range of values 
      ! in the neighbouring cells
      CALL v_limit_face_mo( p_cc       = p_cc(:,:),   & ! in
        &                   p_face     = p_face(:,:), & ! inout
        &                   i_startidx = i_startidx,  & ! in
        &                   i_endidx   = i_endidx,    & ! in
        &                   slev       = slev,        & ! in
        &                   elev       = elev-1       ) ! in

      ! top and bottom face
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
      p_face(i_startidx:i_endidx,slev)   = p_cc(i_startidx:i_endidx,slev)
      p_face(i_startidx:i_endidx,elevp1) = p_cc(i_startidx:i_endidx,elev)
!$ACC END KERNELS

    CASE default
      ! do nothing
    END SELECT

!$ACC END DATA

  END SUBROUTINE compute_face_values_psm


  !-------------------------------------------------------------------------
  !>
  !! PPM Face value reconstruction after Colella and Woodward (1984)
  !!
  !! PPM Face value reconstruction after Colella and Woodward (1984)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-04-12)
  !!
  !
  SUBROUTINE compute_face_values_ppm( i_startidx, i_endidx, slev, elev, &
    &                                 p_itype_vlimit, p_cc, p_cellhgt_mc_now, p_face )

    INTEGER, INTENT(IN) ::    &    !< horizontal start and end indices
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN) ::    &    !< vertical start and end levels
      &  slev, elev

    INTEGER, INTENT(IN) ::    &    !< selects the vertical limiter
      &  p_itype_vlimit

    REAL(wp), INTENT(IN) ::   &    !< cell centered variable (cell average)
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::   &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:)     !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT):: &    !< face values of transported field
      &  p_face(:,:)               !< dim: (nproma,nlev)


    ! local variables
    INTEGER :: jc, jk
    INTEGER :: ikm1, ikp1, ikp2      !< vertical level minus and plus one, plus two
    INTEGER :: slevp1, elevp1        !< start/end level + 1

    REAL(wp) ::   &                  !< auxiliaries for optimization
      &  zfac, zfac_m1

#ifndef _OPENACC
    REAL(wp) ::   &                  !< auxiliary field for optimization
      &  zfac_n(nproma)
#endif

    REAL(wp) ::   &                  !< geometric factors
      &  zgeo1, zgeo2, zgeo3, zgeo4

    REAL(wp) :: &                     !< (monotonized) slope
      &  z_slope(SIZE(p_cc,1),SIZE(p_cc,2))

    !-------------------------------------------------------------------------

    slevp1 = slev + 1
    elevp1 = elev + 1

!$ACC DATA CREATE( z_slope), PCOPYIN(p_cc, p_cellhgt_mc_now), &
!$ACC PCOPYOUT( p_face ), IF( i_am_accel_node .AND. acc_on )

    !
    ! 1. Compute slope
    !
    ! Initialize z_slope for jk=slev
!$ACC KERNELS IF( i_am_accel_node .AND. acc_on )
    z_slope(i_startidx:i_endidx,slev) = 0._wp
!$ACC END KERNELS

#ifndef _OPENACC
    ! Initialize zfac_n for jk=slev
    zfac_n(i_startidx:i_endidx) = (p_cc(i_startidx:i_endidx,slevp1) - p_cc(i_startidx:i_endidx,slev)) &
      &                         /( p_cellhgt_mc_now(i_startidx:i_endidx,slevp1)                       &
      &                          + p_cellhgt_mc_now(i_startidx:i_endidx,slev) )
#endif

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikm1, ikp1, zfac, zfac_m1 ) COLLAPSE(2)
    DO jk = slevp1, elev
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, elev )

#ifndef _OPENACC
        zfac_m1 = zfac_n(jc)
#else
        zfac_m1 = (p_cc(jc,jk) - p_cc(jc,ikm1))  &
          &     / (p_cellhgt_mc_now(jc,jk) + p_cellhgt_mc_now(jc,ikm1))
#endif

        zfac = (p_cc(jc,ikp1) - p_cc(jc,jk)) &
          &  / (p_cellhgt_mc_now(jc,ikp1) + p_cellhgt_mc_now(jc,jk))

        z_slope(jc,jk) = ( p_cellhgt_mc_now(jc,jk)                                       &
          &  / (p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk)                      &
          &  + p_cellhgt_mc_now(jc,ikp1)) )                                              &
          &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk)) * zfac    &
          &  + (p_cellhgt_mc_now(jc,jk) + 2._wp * p_cellhgt_mc_now(jc,ikp1)) * zfac_m1)

#ifndef _OPENACC
        zfac_n(jc) = zfac
#endif
      END DO  ! jc

    END DO  ! jk
!$ACC END PARALLEL

    !
    ! 2. Optional: monotonize slope if necessary
    !     - only necessary, when using the monotonic or semi-monotonic 
    !       sub-grid scale filter (i.e. if the parabola is modified). 
    !
    IF (p_itype_vlimit == ISLOPEL_VSM) THEN
      CALL v_limit_slope_sm(p_cc(:,:), i_startidx, i_endidx, slevp1, elev, z_slope)
    ELSE IF (p_itype_vlimit == ISLOPEL_VM) THEN
      CALL v_limit_slope_mo(p_cc(:,:), i_startidx, i_endidx, slevp1, elev, z_slope)
    ENDIF



    !
    ! 3. face value reconstruction at vertical half-levels
    !

    ! Boundary values for two uppermost and lowermost half-levels
    !
    ! for faces k=slevp1 and k=nlevp1-1 face values are reconstructed by
    ! interpolating a quadratic (instead of quartic) polynomial through 3
    ! values of the indefinite integral A=\int_{z_{0}}^{z}q\,\mathrm{d}z
    !
    ! for faces k=slev and k=nlevp1 a zero gradient condition is assumed and the
    ! face values are set to the values of the corresponding cell centers
    !
!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx

      ! see transport documentation for derivation
      p_face(jc,slevp1) = p_cc(jc,slevp1)*(1._wp - (p_cellhgt_mc_now(jc,slevp1) &
        &       / p_cellhgt_mc_now(jc,slev))) + (p_cellhgt_mc_now(jc,slevp1)    &
        &       /(p_cellhgt_mc_now(jc,slev) + p_cellhgt_mc_now(jc,slevp1)))     &
        &       * ((p_cellhgt_mc_now(jc,slevp1) / p_cellhgt_mc_now(jc,slev))    &
        &       * p_cc(jc,slevp1) + p_cc(jc,slev))

      p_face(jc,elev) = p_cc(jc,elev)*( 1._wp                                &
        &       - (p_cellhgt_mc_now(jc,elev) / p_cellhgt_mc_now(jc,elev-1))) &
        &       + (p_cellhgt_mc_now(jc,elev)/(p_cellhgt_mc_now(jc,elev-1)    &
        &       + p_cellhgt_mc_now(jc,elev))) * ((p_cellhgt_mc_now(jc,elev)  &
        &       / p_cellhgt_mc_now(jc,elev-1)) * p_cc(jc,elev)               &
        &       + p_cc(jc,elev-1))

      p_face(jc,slev)   = p_cc(jc,slev)
      p_face(jc,elevp1) = p_cc(jc,elev)

    ENDDO
!$ACC END PARALLEL


!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR PRIVATE( ikm1, ikp1, ikp2, zgeo1, zgeo2, zgeo3, zgeo4 ) COLLAPSE(2)
    DO jk = slevp1, elev-2
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1 = jk - 1
        ! index of bottom half level
        ikp1 = jk + 1
        ikp2 = jk + 2

        zgeo1 = p_cellhgt_mc_now(jc,jk)                                      &
          &   / (p_cellhgt_mc_now(jc,jk) + p_cellhgt_mc_now(jc,ikp1))
        zgeo2 = 1._wp / (p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk) &
          &   + p_cellhgt_mc_now(jc,ikp1) + p_cellhgt_mc_now(jc,ikp2))
        zgeo3 = (p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk))        &
          &   / (2._wp*p_cellhgt_mc_now(jc,jk) + p_cellhgt_mc_now(jc,ikp1))
        zgeo4 = (p_cellhgt_mc_now(jc,ikp2) + p_cellhgt_mc_now(jc,ikp1))      &
          &   / (2._wp*p_cellhgt_mc_now(jc,ikp1) + p_cellhgt_mc_now(jc,jk))


        p_face(jc,ikp1) = p_cc(jc,jk)                                  &
          &  + zgeo1 * (p_cc(jc,ikp1) - p_cc(jc,jk))                   &
          &  + zgeo2 * ( (2._wp * p_cellhgt_mc_now(jc,ikp1) * zgeo1)   &
          &  * ( zgeo3 - zgeo4 ) * (p_cc(jc,ikp1) - p_cc(jc,jk))       &
          &  - zgeo3 * p_cellhgt_mc_now(jc,jk)   * z_slope(jc,ikp1)    &
          &  + zgeo4 * p_cellhgt_mc_now(jc,ikp1) * z_slope(jc,jk) )

      END DO  !jc

    END DO  !jk
!$ACC END PARALLEL

!$ACC END DATA

  END SUBROUTINE compute_face_values_ppm


END MODULE mo_advection_vflux
