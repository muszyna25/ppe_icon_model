!>
!! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!!
!! @author Peter Korn, Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model V1.0.3 to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_hydro_ocean_run including updated reconstructions
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
MODULE mo_hydro_ocean_run
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
USE mo_impl_constants,         ONLY: max_char_length
USE mo_model_domain,           ONLY: t_patch
USE mo_model_domain_import,    ONLY: n_dom, nroot
USE mo_ocean_nml,              ONLY: n_zlev, iswm_oce, no_tracer
USE mo_dynamics_nml,           ONLY: nold,nnew
USE mo_io_nml,                 ONLY: out_expname
USE mo_run_nml,                ONLY: nsteps, dtime, ltimer
USE mo_exception,              ONLY: message, message_text, finish
USE mo_io_units,               ONLY: filename_max
USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time
USE mo_timer,                  ONLY: timer_total, timer_start, timer_stop
!USE mo_loopindices,            ONLY: get_indices_c, get_indices_e
USE mo_oce_ab_timestepping,    ONLY: solve_free_surface_eq_ab,            &
  &                                  calc_normal_velocity_ab  ,           &
  &                                  calc_vert_velocity
USE mo_oce_tracer_transport,   ONLY: advect_tracer_ab
USE mo_oce_state,              ONLY: t_hydro_ocean_state,                     &
  &                                  construct_hydro_ocean_state, destruct_hydro_ocean_state
USE mo_oce_physics,            ONLY: t_ho_params,t_ho_physics, construct_ho_physics,&
  &                                  init_ho_physics, destruct_ho_physics, construct_ho_params, &
  &                                  init_ho_params, destruct_ho_params
USE mo_io_vlist,               ONLY: setup_vlist_oce, write_vlist_oce, destruct_vlist_oce
USE mo_oce_index,              ONLY: c_b, c_i, c_k, ldbg, form4ar, init_index_test
USE mo_oce_forcing,            ONLY: t_ho_sfc_flx,construct_ho_sfcflx, update_ho_sfcflx
USE mo_interpolation,          ONLY: t_int_state
USE mo_oce_init,               ONLY: init_ho_testcases
!USE mo_oce_diagnostics,        ONLY: calculate_oce_diagnostics,&
!                                  & construct_oce_diagnostics,&
!                                  & destruct_oce_diagnostics, t_oce_timeseries
IMPLICIT NONE

PRIVATE

INCLUDE 'cdi.inc'

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

!public interface
!
! subroutines
PUBLIC :: perform_ho_stepping
PUBLIC :: prepare_ho_integration
PUBLIC :: finalise_ho_integration
!
!
!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Main stepping routine for call of hydrostatic ocean model
  !!
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE perform_ho_stepping( ppatch, pstate_oce,&
                                & datetime, n_io, n_file, p_int,&
                                & p_sfc_flx, p_phys_param)

  TYPE(t_patch),             TARGET, INTENT(IN)     :: ppatch(n_dom)
  TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT)  :: pstate_oce(n_dom)
  TYPE(t_datetime), INTENT(INOUT)                   :: datetime
  INTEGER, INTENT(IN)                               :: n_io, n_file
  TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL     :: p_int(n_dom)
  TYPE(t_ho_sfc_flx)                                :: p_sfc_flx
  TYPE (t_ho_params)                                :: p_phys_param 

  ! local variables
  INTEGER :: jt, jg, n_temp
  INTEGER :: jstep, jfile, jlev
  LOGICAL :: l_exist
  !TYPE(t_oce_timeseries), POINTER :: oce_ts

  CHARACTER(LEN=filename_max)  :: outputfile, gridfile
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_hydro_ocean_run:perform_ho_stepping'
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! no grid refinement allowed here so far
  !------------------------------------------------------------------

  IF (n_dom > 1 ) THEN
    CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
  END IF
  jg = n_dom

  ! set time level to old level (nold(1)=3)
  jt = nold(jg)

  ! file 1 is opened in control_model setup:
  jfile = 1

  CALL init_index_test( ppatch, pstate_oce )

  !CALL construct_oce_diagnostics( ppatch(jg), pstate_oce(jg), oce_ts)

  IF (ltimer) CALL timer_start(timer_total)

  !------------------------------------------------------------------
  ! call the dynamical core: start the time loop
  !------------------------------------------------------------------

  TIME_LOOP: DO jstep = 1, nsteps

    write(*,*)'-----------timestep:',jstep

    !In case of a time-varying forcing: 
    !CALL update_ho_sfcflx(ppatch(jg), p_sfc_flx)


    ! #slo# Put the steps from mo_ocean_semi_implicit_ab here!
    !  - avoid second definition of i_ost_idx/blk/ilv for tests

    ! solve for new free surface
    !CALL message (TRIM(routine), 'calculate free surface')

    CALL solve_free_surface_eq_ab (ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param, &
      &                            jstep, p_int(jg))

    ! Step 4: calculate final normal velocity from predicted horizontal velocity vn_pred
    !         and updated surface height
    CALL calc_normal_velocity_ab(ppatch(jg), pstate_oce(jg))

    IF ( iswm_oce /= 1 ) THEN

      ! Step 5: calculate vertical velocity from continuity equation under incompressiblity condition
      CALL calc_vert_velocity( ppatch(jg), pstate_oce(jg))

      ! Step 6 transport tracer and diffuse them
      IF(no_tracer>=1)THEN
        CALL advect_tracer_ab(ppatch(jg), pstate_oce(jg), p_phys_param,p_sfc_flx, jstep)
       ENDIF
    ENDIF

    ! Step 7: Swap time indices before output
    !         half time levels of semi-implicit Adams-Bashforth timestepping are
    !         stored in auxiliary arrays g_n and g_nimd of p_diag%aux
    n_temp    = nold(jg)
    nold(jg)  = nnew(jg)
    nnew(jg)  = n_temp

   !Actually diagnostics for 3D not implemented, PK March 2011 
   IF ( iswm_oce == 1 ) THEN
    !CALL calculate_oce_diagnostics( ppatch(jg),    &
    !                              & pstate_oce(jg),&
    !                              & p_sfc_flx,     &
    !                              & p_phys_param,  &
    !                              & jstep,         &
    !                              & oce_ts)
   ENDIF 
    ! Step 8: test output
!     IF (ldbg) WRITE(*,'(a,i5,2(a,g20.12))') '*** After jstep = ',jstep, &
!       &  '  Elevation h =', pstate_oce(jg)%p_prog(nold(jg))%h(c_i,c_b), &
!       &  '  Velocity  u =', pstate_oce(jg)%p_diag%u_pred(c_i,c_k,c_b)

     IF ( (jstep/=1 .AND. MOD(jstep-1,n_io)==0) .OR. jstep==nsteps ) THEN


      CALL message (TRIM(routine),'Write output at:')
      CALL print_datetime(datetime)

      CALL write_vlist_oce( ppatch(jg), pstate_oce(jg), p_sfc_flx, datetime )

    END IF

    ! close the current output file and trigger a new one
    ! #slo#: not synchronized with write_vlist - should be closed/renamed/opened
    !        before new writing timestep!
    IF (jstep/=1 .AND. (MOD(jstep-1,n_file)==0) .AND. jstep/=nsteps) THEN

      jlev = ppatch(jg)%level
      CALL destruct_vlist_oce( jg )

      ! contruct gridfile name once more as in control_model:
      WRITE (gridfile,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',jlev,'-grid.nc'
      INQUIRE (FILE=gridfile, EXIST=l_exist)
      IF (.NOT. l_exist) CALL finish(TRIM(routine),' gridfile does not exist')

      ! contruct new outputfile name:
      jfile = jfile +1
      WRITE (outputfile,'(a,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
           &  TRIM(out_expname), '_O.R', nroot, 'B', jlev, 'L', n_zlev, '_', jfile, '.nc'
      WRITE(message_text,'(a,a)') 'New output file for setup_vlist_oce is ',TRIM(outputfile)
      CALL message(trim(routine),message_text)
      CALL setup_vlist_oce( ppatch(jg), TRIM(gridfile), TRIM(outputfile), jg )

    END IF

    ! One integration cycle finished on the lowest grid level (coarsest
    ! resolution). Set model time.

    CALL add_time(dtime,0,0,0,datetime)

    CALL message (TRIM(routine),'Step completed at:')
    CALL print_datetime(datetime)

  ENDDO TIME_LOOP

  !CALL destruct_oce_diagnostics(oce_ts)

  IF (ldbg) THEN
    WRITE(*,*)  ' After run:'
    WRITE(*,form4ar) ' Elevation h at cell   =', pstate_oce(jg)%p_prog(jt)%h(c_i,c_b),  &
      &              '  Tracer 1 =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1)
  ENDIF

  IF (ltimer) CALL timer_stop(timer_total)

  END SUBROUTINE perform_ho_stepping

  !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
  !>
  !! Simple routine for preparing hydrostatic ocean model.
  !!
  !! Simple routine for preparing hydrostatic ocean model.
  !! Calls basic routines ...
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE prepare_ho_integration(ppatch, pstate_oce, p_sfc_flx, p_phys_param, p_physics_oce)

    TYPE(t_patch),             INTENT(INOUT)  :: ppatch(n_dom)
    TYPE(t_hydro_ocean_state), INTENT(INOUT)  :: pstate_oce(n_dom)
    TYPE(t_ho_sfc_flx)                        :: p_sfc_flx
    TYPE (t_ho_params),        INTENT(INOUT)  :: p_phys_param 
    TYPE(t_ho_physics),        INTENT(INOUT)  :: p_physics_oce

    ! local variables
    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &      routine = 'mo_test_hydro_ocean:prepare_ho_integration'

    !------------------------------------------------------------------
    ! no grid refinement allowed here so far
    !------------------------------------------------------------------

    IF (n_dom > 1 ) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF
    jg = n_dom

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! ppatch and pstate_oce have dimension n_dom
    CALL construct_hydro_ocean_state(ppatch, pstate_oce)

    CALL construct_ho_params(ppatch(jg), p_phys_param)
    CALL init_ho_params(p_phys_param)

    CALL construct_ho_physics(ppatch(jg), p_physics_oce)
    CALL init_ho_physics(p_physics_oce)



    !------------------------------------------------------------------
    ! construct ocean forcing and testcases
    !------------------------------------------------------------------
    CALL  construct_ho_sfcflx(ppatch(jg), p_sfc_flx)
    CALL init_ho_testcases(ppatch(jg), pstate_oce(jg), p_sfc_flx)

  END SUBROUTINE prepare_ho_integration

  !-------------------------------------------------------------------------
  !>
  !! Simple routine for finalising integration of hydrostatic ocean model.
  !!
  !! Simple routine for finalising integration of hydrostatic ocean model.
  !! Calls basic routines ...
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE finalise_ho_integration(p_os, p_phys_param, p_physics_oce)
    TYPE(t_hydro_ocean_state), INTENT(INOUT)  :: p_os(n_dom)
    TYPE (t_ho_params),        INTENT(INOUT)  :: p_phys_param 
    TYPE(t_ho_physics),        INTENT(INOUT)  :: p_physics_oce

    !------------------------------------------------------------------
    ! destruct ocean physics and forcing
    ! destruct ocean state is in control_model
    !------------------------------------------------------------------
     CALL destruct_hydro_ocean_state(p_os)
     CALL destruct_ho_params(p_phys_param)
    CALL destruct_ho_physics(p_physics_oce)

    !CALL destruct_ho_forcing

  END SUBROUTINE finalise_ho_integration


 

END MODULE mo_hydro_ocean_run

