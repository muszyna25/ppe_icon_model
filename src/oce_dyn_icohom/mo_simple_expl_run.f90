!>
!! Contains the stepping routine for the simplified explicit ocean core.
!!
!! @author Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-08
!!   - introduction of a simplified explicit core for test purposes
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_simple_expl_run
!
!
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
MODULE mo_simple_expl_run
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
USE mo_kind,                   ONLY: wp
USE mo_mpi,                    ONLY: my_process_is_stdio
USE mo_parallel_config,        ONLY: nproma
USE mo_io_units,               ONLY: filename_max
USE mo_impl_constants,         ONLY: sea_boundary, max_char_length,           &
  &                                  min_rlcell, min_rledge, toplev
USE mo_model_domain,           ONLY: t_patch
USE mo_model_domain_import,    ONLY: n_dom
USE mo_ocean_nml,              ONLY: n_zlev, iswm_oce
USE mo_dynamics_config,        ONLY: nold
USE mo_io_config,              ONLY: out_expname
USE mo_run_config,             ONLY: nsteps, dtime
USE mo_exception,              ONLY: message, message_text, finish, get_filename_noext
USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time
USE mo_loopindices,            ONLY: get_indices_c, get_indices_e
USE mo_oce_state,              ONLY: t_hydro_ocean_state
USE mo_oce_thermodyn,          ONLY: calc_density, calc_internal_press
USE mo_oce_math_operators,     ONLY: grad_fd_norm_oce_2d, grad_fd_norm_oce, div_oce
USE mo_advection_utils,        ONLY: laxfr_upflux, laxfr_upflux_v
USE mo_physical_constants,     ONLY: grav
USE mo_io_vlist,               ONLY: setup_vlist_oce, destruct_vlist_oce !, write_vlist_oce
USE mo_oce_index,              ONLY: c_b, c_i, c_k, ne_b, ne_i, nc_b, nc_i, form4ar, ldbg
USE mo_oce_physics,            ONLY: t_ho_params

IMPLICIT NONE
PRIVATE

INCLUDE 'cdi.inc'

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

!public interface
!
! subroutines
PUBLIC :: perform_expl_step_oce
!PUBLIC :: upwind_hflux_oce
!PUBLIC :: upwind_vflux_oce
!
!
!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Stepping routine to test basic characteristics of operators etc.
  !!
  !! Stepping routine to test basic characteristics of operators.
  !! A simplified model with explicit time stepping and upwind advection
  !! of tracers.
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-08)
  !
  !
  SUBROUTINE perform_expl_step_oce( ppatch, pstate_oce, datetime, n_io, n_file )

    TYPE(t_patch),             TARGET, INTENT(IN)     :: ppatch(n_dom)
    TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT)  :: pstate_oce(n_dom)
!   TYPE(t_ho_sfc_flx),                INTENT(INOUT)  :: p_sfc_flx
    TYPE (t_ho_params)                                :: p_phys_param
    TYPE(t_datetime), INTENT(INOUT)                   :: datetime
    INTEGER, INTENT(IN)                               :: n_io, n_file

    ! local variables

    REAL(wp) :: z_gradh_e  (nproma,           ppatch(1)%nblks_e)
    REAL(wp) :: z_vn_pred_e(nproma, n_zlev,   ppatch(1)%nblks_e)
    REAL(wp) :: z_w_pred_c (nproma, n_zlev+1, ppatch(1)%nblks_c)
    REAL(wp) :: z_div_c    (nproma, n_zlev,   ppatch(1)%nblks_c)
    REAL(wp) :: z_h_pred_c (nproma, 1,        ppatch(1)%nblks_c)  ! 3-dim for sbrt-calls
    REAL(wp) :: z_h_e      (nproma, 1,        ppatch(1)%nblks_e)  ! 3-dim for sbrt-calls
    REAL(wp) :: z_tracer_c (nproma, n_zlev,   ppatch(1)%nblks_c)
    REAL(wp) :: z_trtend_c (nproma, n_zlev,   ppatch(1)%nblks_c)  ! tracer tendency
    REAL(wp) :: z_upflux_e (nproma, n_zlev,   ppatch(1)%nblks_e)  ! horizontal mass flux
    REAL(wp) :: z_upflux_i (nproma, n_zlev+1, ppatch(1)%nblks_c)  ! vertical mass flux
    REAL(wp) :: z_divhor_c (nproma, n_zlev,   ppatch(1)%nblks_c)  ! horizontal tracer divergence
    REAL(wp) :: z_divver_c (nproma, n_zlev,   ppatch(1)%nblks_c)  ! vertical tracer divergence

    REAL(wp) :: delta_z

    INTEGER :: jb, jc, jk, jkp, jkp1, jt, jg
    INTEGER :: jstep, jfile, jlev
    INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: iswm_ocex

!    LOGICAL :: l_exist

    CHARACTER(LEN=filename_max)               :: outputfile
!    CHARACTER(LEN=filename_max)               :: gridfile
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &      routine = 'mo_simple_expl_oce:perform_expl_step_oce'

    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    ! file 1 is opened in control_model setup:
    jfile = 1

    ! #slo# due to nag -nan compiler-option set variables to zero to avoid undef values
    z_div_c(:,1,:)     = 0.0_wp
    z_gradh_e(:,:)     = 0.0_wp
    z_vn_pred_e(:,:,:) = 0.0_wp
    z_w_pred_c (:,:,:) = 0.0_wp
    z_h_e(:,1,:)       = 0.0_wp
    z_tracer_c (:,:,:) = 0.0_wp
    z_trtend_c (:,:,:) = 0.0_wp
    z_upflux_e (:,:,:) = 0.0_wp
    z_upflux_i (:,:,:) = 0.0_wp
    z_divhor_c (:,:,:) = 0.0_wp
    z_divver_c (:,:,:) = 0.0_wp

    ! set time level to old level (nold(1)=3)
    jt = nold(jg)

    !------------------------------------------------------------------
    ! call the simplified explicit core: start the time loop
    !------------------------------------------------------------------

    TIME_LOOP: DO jstep = 1, nsteps

      !WRITE(message_text,'(a,i10)') ' - begin of jstep =', jstep
      !CALL message (TRIM(routine), message_text)

      IF (my_process_is_stdio() .AND. ldbg) THEN

        jkp=c_k+1
        IF ( iswm_oce == 1 ) jkp=c_k
        WRITE(message_text,form4ar) &
          &   'Elevation h  =', pstate_oce(jg)%p_prog(jt)%h     (c_i,    c_b),   &
          &  '    Tracer 1  =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1), &
          &  '  Tracer(k+1) =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,jkp,c_b,1)
        CALL message (' ', message_text)
        WRITE(message_text,form4ar) &
          &  ' NormalVel.E1 =', pstate_oce(jg)%p_prog(jt)%vn(ne_i(1),c_k,ne_b(1)), &
          &  '  Ver.Vel(k)C =', pstate_oce(jg)%p_diag%w(c_i,c_k,c_b),            &
          &  '  Ver.Vel(k+1)=', pstate_oce(jg)%p_diag%w(c_i,jkp,c_b)
        CALL message (' ', message_text)

      END IF

      !*******************************************************************************
      ! Steps to be performed:
      !*******************************************************************************
      !
      ! Step 1: calculate density and pressure
      ! Step 2: calculate gradients of pressure and surface elevation
      ! Step 3: calculate horizontal velocity
      ! Step 4: diagnose vertical velocity
      ! Step 5: advect tracers
      ! Step 6: calculate final surface elevation
      ! Step 7: update prognositic variables for the next timestep
      ! Step 8: write output to files
      !*******************************************************************************


      !*******************************************************************************
      ! Step 1: calculate density and pressure
      !*******************************************************************************

      !CALL message (TRIM(routine), 'Integration Step 1: calculate density ')

      !------------------------------------------------------------------
      !  Not for shallow water mode
      IF ( iswm_oce /= 1 ) THEN

        !------------------------------------------------------------------
        ! calculate density using the temperature & salinity from timestep "n"(=old)
        !------------------------------------------------------------------
        !
        CALL calc_density ( ppatch(jg),                                  &
          &                 pstate_oce(jg)%p_prog(jt)%tracer(:,:,:,1:2), &
          &                 pstate_oce(jg)%p_diag%rho(:,:,:) )

        IF (my_process_is_stdio() .AND. ldbg) THEN
          WRITE(message_text,form4ar) &
            &   'Density    C =', pstate_oce(jg)%p_diag%rho(c_i,c_k,c_b),   &
            &  '  Dens(k+1) C =', pstate_oce(jg)%p_diag%rho(c_i,c_k+1,c_b), &
            &  '  Dens(k+2) C =', pstate_oce(jg)%p_diag%rho(c_i,c_k+2,c_b)
          CALL message (' ', message_text)
        END IF

        !------------------------------------------------------------------
        ! calculate hydrostatic pressure from density
        !------------------------------------------------------------------
        !
        CALL calc_internal_press ( ppatch(jg),p_phys_param,                &
          &                        pstate_oce(jg)%p_diag%rho(:,:,:),       &
          &                        pstate_oce(jg)%p_prog(jt)%h(:,:),       &
          &                        pstate_oce(jg)%p_diag%press_hyd(:,:,:) )
        IF (my_process_is_stdio() .AND. ldbg) THEN
          WRITE(message_text,form4ar) &
            &   'Hyd.Press. C =', pstate_oce(jg)%p_diag%press_hyd(c_i,c_k,c_b),   &
            &  '  Pres(k+1) C =', pstate_oce(jg)%p_diag%press_hyd(c_i,c_k+1,c_b), &
            &  '  Pres(k+2) C =', pstate_oce(jg)%p_diag%press_hyd(c_i,c_k+2,c_b)
          CALL message (' ', message_text)
          WRITE(message_text,form4ar) &
            &   'Hyd.Press.C1 =', pstate_oce(jg)%p_diag%press_hyd(nc_i(1),c_k,nc_b(1)), &
            &  '  Pressure C2 =', pstate_oce(jg)%p_diag%press_hyd(nc_i(2),c_k,nc_b(2)), &
            &  '  Pressure C3 =', pstate_oce(jg)%p_diag%press_hyd(nc_i(3),c_k,nc_b(3))
          CALL message (' ', message_text)
        END IF


      !*******************************************************************************
      ! Step 2: calculate gradients of pressure and surface elevation
      !*******************************************************************************

      !CALL message (TRIM(routine), 'Integration Step 2: calculate gradients ')

      !------------------------------------------------------------------
      ! calculate pressure gradients on edges in [m/s]
      !------------------------------------------------------------------

      ! calculate gradient of hydrostatic pressure in 3D
      CALL grad_fd_norm_oce( pstate_oce(jg)%p_diag%press_hyd(:,:,:),  &
        &                    ppatch(jg),                              &
        &                    pstate_oce(jg)%p_diag%press_grad(:,:,:))

      END IF

      !------------------------------------------------------------------
      ! calculate gradient of surface height
      !------------------------------------------------------------------
      !CALL message (TRIM(routine), 'call grad_fd_norm_oce_2d')

      ! #slo# use 3-d version!
      CALL grad_fd_norm_oce_2d( pstate_oce(jg)%p_prog(jt)%h, ppatch(jg), z_gradh_e)

      IF (my_process_is_stdio() .AND. ldbg) THEN
        WRITE(message_text,form4ar) &
          &   'Elev. at C1  =', pstate_oce(jg)%p_prog(jt)%h(nc_i(1),nc_b(1)), &
          &  '  Elevat. C2  =', pstate_oce(jg)%p_prog(jt)%h(nc_i(2),nc_b(2)), &
          &  '  Elevat. C3  =', pstate_oce(jg)%p_prog(jt)%h(nc_i(3),nc_b(3))
        CALL message (' ', message_text)
        WRITE(message_text,form4ar) &
          &   'Grad(Elv) E1 =', z_gradh_e(ne_i(1),ne_b(1)), &
          &  '   grad(h) E2 =', z_gradh_e(ne_i(2),ne_b(2)), &
          &  '   grad(h) E3 =', z_gradh_e(ne_i(3),ne_b(3))
        CALL message (' ', message_text)
      END IF


      !*******************************************************************************
      ! Step 3: calculate horizontal velocity
      !*******************************************************************************
      !
      !CALL message (TRIM(routine), 'Integration Step 3: calculate horizontal velocity')

      !------------------------------------------------------------------
      ! calculation of normal velocity vn on edges (c-grid)
      !  - pressure gradients multiplied by timestep determine predicted
      !    normal velocity z_vn_pred_e at edges
      !  - much simpler than using u and v on cells, no reconstructions needed here
      !  - in the shallow water case the contribution reduces implicitly
      !    to grav*grad(h)*dtime, which is calculated above
      !------------------------------------------------------------------

      ! values for the blocking
      rl_start = 1  !  without using neighbouring points
      rl_end = min_rledge
      i_startblk = ppatch(jg)%edges%start_blk(rl_start,1)
      i_endblk   = ppatch(jg)%edges%end_blk(rl_end,1)

      iswm_ocex=1
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ppatch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx,  &
          &                rl_start, rl_end)

#ifdef __SX__
!CDIR UNROLL=6
#endif
        DO jk = 1, n_zlev
          DO jc = i_startidx, i_endidx
!         IF ( iswm_oce == 1 ) THEN
          IF ( iswm_ocex == 1 ) THEN
              ! #slo# - z_gradh_e needed for swm only
            z_vn_pred_e(jc,jk,jb) = pstate_oce(jg)%p_prog(jt)%vn(jc,jk,jb) - &
              &                     grav * z_gradh_e(jc,jb) * dtime
          ELSE
              ! #slo# - 3-dim mode - grav*grad(h) is included in press_grad at jk=1
            z_vn_pred_e(jc,jk,jb) = pstate_oce(jg)%p_prog(jt)%vn(jc,jk,jb) - &
              &                     pstate_oce(jg)%p_diag%press_grad(jc,jk,jb) * dtime
          END IF
          END DO
        END DO
      END DO

      IF (my_process_is_stdio() .AND. ldbg) THEN
        WRITE(message_text,form4ar) &
          &   'Press.gradE1 =', pstate_oce(jg)%p_diag%press_grad(ne_i(1),c_k,ne_b(1)), &
          &  '   P_grad  E2 =', pstate_oce(jg)%p_diag%press_grad(ne_i(2),c_k,ne_b(2)), &
          &  '   P_grad  E3 =', pstate_oce(jg)%p_diag%press_grad(ne_i(3),c_k,ne_b(3))
        CALL message (' ', message_text)
        WRITE(message_text,form4ar) &
          &   'NormVeloc E1 =', pstate_oce(jg)%p_prog(jt)%vn(ne_i(1),c_k,ne_b(1)), &
          &  '  NorVeloc E2 =', pstate_oce(jg)%p_prog(jt)%vn(ne_i(2),c_k,ne_b(2)), &
          &  '  NorVeloc E3 =', pstate_oce(jg)%p_prog(jt)%vn(ne_i(3),c_k,ne_b(3))
        CALL message (' ', message_text)
        WRITE(message_text,form4ar) &
          &   'NorV_pred E1 =', z_vn_pred_e(ne_i(1),c_k,ne_b(1)), &
          &  '  NorV_pre E2 =', z_vn_pred_e(ne_i(2),c_k,ne_b(2)), &
          &  '  NorV_pre E3 =', z_vn_pred_e(ne_i(3),c_k,ne_b(3))
        CALL message (' ', message_text)
      END IF


      !*******************************************************************************
      ! Step 4: diagnose vertical velocity
      !*******************************************************************************
      !
      !CALL message (TRIM(routine), 'Integration Step 4: diagnose vertical velocity')

      !------------------------------------------------------------------
      ! store old vertical velocity in diagnostic component
      !------------------------------------------------------------------
      pstate_oce(jg)%p_diag%w_prev(:,:,:) = pstate_oce(jg)%p_diag%w(:,:,:)

      !------------------------------------------------------------------
      ! diagnose cell-based vertical velocity from the continuity equation and the
      ! incompressibility condition using predicted normal velocity on edges
      !------------------------------------------------------------------
      ! #slo# - change sequence in mo_math: call calc_vert_veloc(vn_e,ppatch,h_c,w_c)
      !PKCALL calc_vert_velocity(ppatch(jg), z_vn_pred_e, pstate_oce(jg)%p_prog(jt)%h, z_w_pred_c )

      IF (my_process_is_stdio() .AND. ldbg) THEN

        ! test call of divergence of hor.vel. as in calc_vert_velocity:
        CALL div_oce( z_vn_pred_e, ppatch(jg), z_divhor_c)
        jkp=c_k+1
        if ( iswm_oce == 1 ) jkp=1
        WRITE(message_text,form4ar) &
          &          'Div(NorV) C1 =', z_divhor_c(nc_i(1),c_k,nc_b(1)), &
          &         ' Div(NorV) C2 =', z_divhor_c(nc_i(2),c_k,nc_b(2)), &
          &         ' Div(NorV) C3 =', z_divhor_c(nc_i(3),c_k,nc_b(3))
        WRITE(message_text,form4ar) &
          &          'Div(NorV)  C =', z_divhor_c(c_i,c_k,c_b), &
          &         ' Div(NV)*dz C =', z_divhor_c(c_i,c_k,c_b)* &
          &   (pstate_oce(jg)%p_prog(jt)%h(c_i,c_b)+ppatch(jg)%patch_oce%del_zlev_m(c_k)), &
          &         ' Div(NV)lv2 C =', z_divhor_c(c_i,jkp,c_b)
        CALL message (' ', message_text)
      END IF


      !*******************************************************************************
      ! Step 5: calculate final surface elevation
      !*******************************************************************************
      !
      !CALL message (TRIM(routine), 'Integration Step 5: calculate surface elevation')

      !------------------------------------------------------------------
      ! calculate cell-based predicted surface elevation from vertical velocity at
      ! surface times timestep
      ! #slo# 2010-10-18: positive h for positive w
      !------------------------------------------------------------------

      z_h_pred_c(:,1,:) = pstate_oce(jg)%p_prog(jt)%h(:,:) + dtime * z_w_pred_c(:,1,:)

      IF (my_process_is_stdio() .AND. ldbg) THEN
        jkp=3
        if ( iswm_oce == 1 ) jkp=2
        WRITE(message_text,form4ar) &
          &     'w_pred(lv=1) =', z_w_pred_c(c_i,1,c_b), &
          &    '  w_pred(k=2) =', z_w_pred_c(c_i,2,c_b), &
          &    '  w_pred(k=3) =', z_w_pred_c(c_i,jkp,c_b)
        CALL message (' ', message_text)
        WRITE(message_text,form4ar) &
          &     'w_pred  C1   =', z_w_pred_c(ne_i(1),c_k,ne_b(1)), &
          &    ' w_pred  C2   =', z_w_pred_c(ne_i(2),c_k,ne_b(2)), &
          &    ' w_pred  C3   =', z_w_pred_c(ne_i(3),c_k,ne_b(3))
        CALL message (' ', message_text)
        WRITE(message_text,form4ar) &
          &     'h_pred  C    =', z_h_pred_c(c_i,    c_k,c_b),     &
          &    ' h_pred  C2   =', z_h_pred_c(nc_i(2),c_k,nc_b(2)), &
          &    ' h_pred  C3   =', z_h_pred_c(nc_i(3),c_k,nc_b(3))
        CALL message (' ', message_text)
      END IF


      !*******************************************************************************
      ! Step 6: advect tracers
      !*******************************************************************************
      !
      !------------------------------------------------------------------
      !  Not for shallow water mode
      !------------------------------------------------------------------
      IF ( iswm_oce /= 1 ) THEN

        !CALL message (TRIM(routine), 'Integration Step 6: advect tracers')

        !------------------------------------------------------------------
        !  map cell based surface elevation to edges
        !  #slo# - must be updated with dual_edge_length
        !------------------------------------------------------------------
        !PKcall map_cell2edge( ppatch(jg), z_h_pred_c, z_h_e, toplev, toplev)

        !------------------------------------------------------------------
        !  horizontal advection of tracers with upwind-method
        !  using predicted horizontal velocity and surface elevation at edges
        !------------------------------------------------------------------

        ! #slo# advect all tracers i=1,ntrac_oce

        z_tracer_c(:,:,:) = pstate_oce(jg)%p_prog(jt)%tracer(:,:,:,1)
        CALL upwind_hflux_oce( ppatch(jg), z_tracer_c, z_vn_pred_e, z_h_e(:,1,:), z_upflux_e )

        IF (my_process_is_stdio() .AND. ldbg) THEN
          WRITE(message_text,form4ar) &
            &     'Mapped H  E1 =', z_h_e(ne_i(1),c_k,ne_b(1)), &
            &    '  Mapped H E2 =', z_h_e(ne_i(2),c_k,ne_b(2)), &
            &    '  Mapped H E3 =', z_h_e(ne_i(3),c_k,ne_b(3))
          CALL message (' ', message_text)
          WRITE(message_text,form4ar) &
            &     'Tracer at C1 =', z_tracer_c(ne_i(1),c_k,ne_b(1)), &
            &    '  Tracer   C2 =', z_tracer_c(ne_i(2),c_k,ne_b(2)), &
            &    '  Tracer   C3 =', z_tracer_c(ne_i(3),c_k,ne_b(3))
          CALL message (' ', message_text)
        END IF

        !------------------------------------------------------------------
        !  divergence of horizontal upwind fluxes of tracers calculates tracer tendencies
        !------------------------------------------------------------------
        CALL div_oce( z_upflux_e, ppatch(jg), z_divhor_c)

        !------------------------------------------------------------------
        !  vertical advection of tracers with upwind-method
        !  using predicted vertical velocity
        !------------------------------------------------------------------
        CALL upwind_vflux_oce( ppatch(jg), z_tracer_c, z_w_pred_c, z_upflux_i )

        !------------------------------------------------------------------
        !  vertical divergence of upwind fluxes of tracers calculates tracer tendencies
        !------------------------------------------------------------------

        ! values for the blocking
        rl_start = 1  !  without using neighbouring points
        rl_end = min_rlcell
        i_startblk = ppatch(jg)%cells%start_blk(rl_start,1)
        i_endblk   = ppatch(jg)%cells%end_blk(rl_end,1)

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(ppatch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, &
            &                rl_start, rl_end)

          DO jk = 1, n_zlev

            ! index of intermediate level below main surface
            jkp1 = jk + 1

            ! depth at coordinate surfaces - elemental prism depth
            delta_z = ppatch(jg)%patch_oce%del_zlev_m(jk)

            DO jc = i_startidx, i_endidx

              ! at surface level add surface elevation z_h_pred_c
              IF (jk == toplev)  &
                &      delta_z = ppatch(jg)%patch_oce%del_zlev_m(jk) + z_h_pred_c(jc,1,jb)

              ! divide by delta_z at level k
              ! positive vertical divergence in direction of w (upward positive)
              z_divver_c(jc,jk,jb) = &
                &   ( z_upflux_i(jc,jk,jb)-z_upflux_i(jc,jkp1,jb) ) / delta_z

            END DO
          END DO
        END DO

        !------------------------------------------------------------------
        !  cumulate tracer tendencies from horizontal and vertical advection
        !------------------------------------------------------------------
        z_trtend_c(:,:,:) =                   - dtime * z_divhor_c(:,:,:)  &
          &                                   - dtime * z_divver_c(:,:,:)
        z_tracer_c(:,:,:) = z_tracer_c(:,:,:)*ppatch(jg)%patch_oce%wet_c(:,:,:) + &
          &                 z_trtend_c(:,:,:)

        IF (my_process_is_stdio() .AND. ldbg) THEN
          WRITE(message_text,form4ar) &
            &     'Horflx    E1 =', z_upflux_e(ne_i(1),c_k,ne_b(1)), &
            &    '  Horflx   E2 =', z_upflux_e(ne_i(2),c_k,ne_b(2)), &
           !&    '  Horflx   E3 =', z_upflux_e(ne_i(3),c_k,ne_b(3)), &
            &    '  Div(Hflx) C =', z_divhor_c(c_i,c_k,c_b)
          CALL message (' ', message_text)
          jkp=c_k+1
          WRITE(message_text,form4ar) &
            &     'Verflx(k)  C =', z_upflux_i(c_i,c_k,c_b), &
            &    '  Vflx(k+1) C =', z_upflux_i(c_i,jkp,c_b), &
            &    '  Div(Vflx) C =', z_divver_c(c_i,c_k,c_b)
          CALL message (' ', message_text)
          WRITE(message_text,form4ar) &
            &     'TrTend Hor C =',-z_divhor_c(c_i,c_k,c_b)*dtime, &
            &    ' TrTend Ver C =',-z_divver_c(c_i,c_k,c_b)*dtime, &
            &    ' TrTendency C =', z_trtend_c(c_i,c_k,c_b)
          CALL message (' ', message_text)
        END IF

      END IF


      !*******************************************************************************
      ! Step 7: update prognostic variables for the next timestep
      !*******************************************************************************
      !
      !CALL message (TRIM(routine), 'Integration Step 7: Update prognostic variables')
      pstate_oce(jg)%p_prog(jt)%h     (:,:)     = z_h_pred_c (:,1,:)
      pstate_oce(jg)%p_prog(jt)%vn    (:,:,:)   = z_vn_pred_e(:,:,:)
      pstate_oce(jg)%p_diag    %w     (:,:,:)   = z_w_pred_c (:,:,:)
      pstate_oce(jg)%p_prog(jt)%tracer(:,:,:,1) = z_tracer_c (:,:,:)


      !*******************************************************************************
      ! Step 8: write output to files
      !*******************************************************************************
      !
      !CALL message (TRIM(routine), 'Integration Step 8: write output ')

      !  - these procedures can be exported into subroutines for stepping and prepare_writing

      IF ( (jstep/=1 .AND. MOD(jstep-1,n_io)==0) .OR. jstep==nsteps ) THEN

        !CALL message (TRIM(routine),'Write output at:')
        !CALL print_datetime(datetime)

        !CALL write_vlist_oce( ppatch(jg), pstate_oce(jg), datetime )
        !CALL write_vlist_oce( ppatch(jg), pstate_oce(jg), p_sfc_flx, datetime )

      END IF

      ! close the current output file and trigger a new one
      ! #slo#: not synchronized with write_vlist - should be closed/renamed/opened
      !        before new writing timestep!
      IF (jstep/=1 .AND. (MOD(jstep-1,n_file)==0) .AND. jstep/=nsteps) THEN

        jlev = ppatch(jg)%level
        CALL destruct_vlist_oce( jg )

        ! contruct gridfile name once more as in control_model:
!         WRITE (gridfile,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',jlev,'-grid.nc'
!         INQUIRE (FILE=gridfile, EXIST=l_exist)
!         IF (.NOT. l_exist) CALL finish(TRIM(routine),' gridfile does not exist')

        ! contruct new outputfile name:      
        jfile = jfile +1
        WRITE (outputfile,'(a,a,a,a,i4.4,a)')  &
          &  TRIM(out_expname), '_', TRIM(get_filename_noext(ppatch(jg)%grid_filename)), &
          & '_', jfile, '.nc'
        WRITE(message_text,'(a,a)') 'New output file for setup_vlist_oce is ',TRIM(outputfile)
        CALL message(trim(routine),message_text)
        CALL setup_vlist_oce( ppatch(jg), TRIM(ppatch(jg)%grid_filename), TRIM(outputfile), jg )

      END IF

      ! One integration cycle finished on the lowest grid level (coarsest
      ! resolution). Set model time.

      CALL add_time(dtime,0,0,0,datetime)

      WRITE(message_text,'(a,i8,a)')  &
        &   'jStep =',jstep,' completed at:'
      CALL message (' ', message_text)
      !CALL message (TRIM(routine), message_text)
      CALL print_datetime(datetime)

    ENDDO TIME_LOOP


    IF (my_process_is_stdio() .AND. ldbg) THEN
      CALL message (TRIM(routine),' - after completion of time loop at cell C:')
      WRITE(message_text,form4ar) &
        &   'Elevation h  =', pstate_oce(jg)%p_prog(jt)%h (c_i,c_b), &
        &  '    Tracer 1  =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1), &
        &  '   Nor.Vel.E1 =', pstate_oce(jg)%p_prog(jt)%vn(ne_i(1),c_k,ne_b(1))
      CALL message (' ', message_text)
    END IF

  END SUBROUTINE perform_expl_step_oce

  !-----------------------------------------------------------------------
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
  SUBROUTINE upwind_hflux_oce( ppatch, pvar_c, pvn_e, ph_e, pupflux_e )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(IN)     :: pvn_e(:,:,:)       !< normal velocity on edges
    REAL(wp), INTENT(IN)     :: ph_e (:,:)         !< surface elevation on edges
    REAL(wp), INTENT(INOUT)  :: pupflux_e(:,:,:)   !< variable in which the upwind flux is stored

    ! local variables
    REAL(wp) :: zhorflx_e(nproma,n_zlev,ppatch%nblks_e) ! contravariant horizontal mass flux
    REAL(wp) :: delta_z

    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices

    INTEGER  :: nblks_c, npromz_c, nblks_e, npromz_e
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_rcstartlev, rl_start, rl_end
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block

    !-----------------------------------------------------------------------

    ! #slo# due to nag -nan compiler-option set variables to zero to avoid undef values
    !       and: values over land are set to zero
    zhorflx_e(:,:,:) = 0.0_wp

    ! values for the blocking
    ! #slo# - Attention - why not using ppatch%nblks_c etc. ?? - Daniel??
    nblks_c  = ppatch%nblks_int_c
    npromz_c = ppatch%npromz_int_c
    nblks_e  = ppatch%nblks_int_e
    npromz_e = ppatch%npromz_int_e

    !
    !  compute horizontal mass flux from horizontal velocity
    !
    !  loop over edges
    rl_start = 1
    rl_end = min_rledge
    i_startblk = ppatch%edges%start_blk(rl_start,1)
    i_endblk   = ppatch%edges%end_blk(rl_end,1)

    DO jb = i_startblk, i_endblk

    CALL get_indices_e(ppatch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, n_zlev
        DO je = i_startidx, i_endidx

          IF ( ppatch%patch_oce%lsm_oce_e(je,jk,jb) <= sea_boundary ) THEN

            ! depth at coordinate surfaces
            delta_z = ppatch%patch_oce%del_zlev_m(jk)

            ! at surface level add surface elevation h_e
            IF (jk == toplev) delta_z = ppatch%patch_oce%del_zlev_m(jk) + ph_e(je,jb)

            ! compute mass flux at edges
            zhorflx_e(je,jk,jb) = pvn_e(je,jk,jb) * delta_z

          END IF
        END DO
      END DO
    END DO

    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a picewise constant approx. of the cell centered values
    ! is used.
    !
    ! no ocean boundary treatment here
    !

    i_rcstartlev = 2
    i_startblk = ppatch%edges%start_blk(i_rcstartlev,1)

    ! line and block indices of two neighboring cells
    iilc => ppatch%edges%cell_idx
    iibc => ppatch%edges%cell_blk

    ! loop through all patch edges (and blocks)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx)
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ppatch, jb, i_startblk, nblks_e,   &
        &                i_startidx, i_endidx, i_rcstartlev)

#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = 1, n_zlev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that only the pvar_c*zhorflx_e value at cell edge is computed
          ! multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          !
          pupflux_e(je,jk,jb) =  &
            &  laxfr_upflux( zhorflx_e(je,jk,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                             pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )

        END DO  ! end loop over edges

      END DO  ! end loop over levels

    END DO  ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

    IF (my_process_is_stdio() .AND. ldbg) THEN
      WRITE(*,form4ar) &
        &  '  zhorflx  E1 =', zhorflx_e(ne_i(1),c_k,ne_b(1)), &
        & '         ph_e  =', ph_e(ne_i(1),ne_b(1))

      WRITE(*,form4ar) &
        &  '  zhorflx  E2 =', zhorflx_e(ne_i(2),c_k,ne_b(2)), &
        & '         ph_e  =', ph_e(ne_i(2),ne_b(2))

      WRITE(*,form4ar) &
        &  '  zhorflx  E3 =', zhorflx_e(ne_i(3),c_k,ne_b(3)), &
        & '         ph_e  =', ph_e(ne_i(3),ne_b(3))

      WRITE(*,form4ar) &
        &  '  pupflx_h E1 =', pupflux_e(ne_i(1),c_k,ne_b(1)), &
        &  '    pvar_c(1) =', pvar_c   (nc_i(1),c_k,nc_i(1)), &
        &  '    pvar_c(2) =', pvar_c   (nc_i(2),c_k,nc_i(2))
    END IF

  END SUBROUTINE upwind_hflux_oce

  !-------------------------------------------------------------------------
  !>
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
  SUBROUTINE upwind_vflux_oce( ppatch, pvar_c, pw_c, pupflux_i )

    TYPE(t_patch), TARGET, INTENT(IN) :: ppatch      !< patch on which computation is performed
    REAL(wp), INTENT(INOUT)  :: pvar_c(:,:,:)      !< advected cell centered variable
    REAL(wp), INTENT(IN)     :: pw_c(:,:,:)        !< vertical velocity on cells
    REAL(wp), INTENT(INOUT)  :: pupflux_i(:,:,:)   !< variable in which the upwind flux is stored
                                                   !< dim: (nproma,n_zlev+1,nblks_c)
    ! local variables
    !REAL(wp) :: zvertflx_i(nproma,n_zlev+1,ppatch%nblks_c) ! contravariant vertical mass flux
    REAL(wp) :: zcoeff_grid

   !INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER  :: nlen, nblks_c, npromz_c  !< values for the blocking
    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: jkm1                     !< jk - 1
    !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ppatch%nblks_int_c
    npromz_c = ppatch%npromz_int_c

    ! no fluxes at top boundary
    pupflux_i(:,1,:)        = 0.0_wp
    ! no fluxes at bottom boundary
    pupflux_i(:,n_zlev+1,:) = 0.0_wp

    ! height based but reversed (downward increasing depth) coordinate system,
    ! grid coefficient is negative (same as pressure based atmospheric coordinate system
    zcoeff_grid = -1.0_wp

    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    ! no ocean boundary treatment here
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,jkm1)
    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF

        DO jk = 2, n_zlev

          ! index of top half level
          jkm1 = jk - 1

          DO jc = 1, nlen

            ! calculate vertical tracer flux using upwind method
            pupflux_i(jc,jk,jb) =                 &
              &  laxfr_upflux_v( pw_c(jc,jk,jb),  &
              &                  pvar_c(jc,jkm1,jb), pvar_c(jc,jk,jb), zcoeff_grid )

            IF (my_process_is_stdio() .AND. ldbg) THEN
              IF ((jb==110).and.(jc==14).and.(jk==2)) THEN
                WRITE(*,form4ar) &
                  &  '  pupflux_v C =', pupflux_i(jc,jk,jb), &
                  &  '    pvar_c(1) =', pvar_c(jc,1,jb),     &
                  &  '    pvar_c(2) =', pvar_c(jc,2,jb)
              END IF
            END IF

          END DO ! end loop over cells
      END DO ! end loop over vertical levels
    END DO ! end loop over blocks
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE upwind_vflux_oce

  !-------------------------------------------------------------------------

END MODULE mo_simple_expl_run

