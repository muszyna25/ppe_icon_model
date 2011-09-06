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
USE mo_model_domain_import,    ONLY: n_dom
USE mo_ocean_nml,              ONLY: iswm_oce, no_tracer,itestcase_oce, EOS_type
USE mo_dynamics_config,        ONLY: nold, nnew
USE mo_io_config,              ONLY: out_expname, istime4output, istime4newoutputfile
USE mo_run_config,             ONLY: nsteps, dtime, ltimer
USE mo_exception,              ONLY: message, message_text, finish
USE mo_ext_data,               ONLY: t_external_data
USE mo_io_units,               ONLY: filename_max
USE mo_datetime,               ONLY: t_datetime, print_datetime, add_time
USE mo_timer,                  ONLY: timer_total, timer_start, timer_stop
!USE mo_loopindices,            ONLY: get_indices_c, get_indices_e
USE mo_oce_ab_timestepping,    ONLY: solve_free_surface_eq_ab,            &
  &                                  calc_normal_velocity_ab  ,           &
  &                                  calc_vert_velocity
USE mo_oce_tracer_transport,   ONLY: advect_tracer_ab
USE mo_oce_state,              ONLY: t_hydro_ocean_state, t_hydro_ocean_base, &
  &                                  init_ho_base, v_base, &
  &                                  construct_hydro_ocean_base, destruct_hydro_ocean_base, &
  &                                  construct_hydro_ocean_state, destruct_hydro_ocean_state, &
  &                                  init_scalar_product_base, init_geo_factors_base
USE mo_oce_physics,            ONLY: t_ho_params, &
                                   & construct_ho_params, init_ho_params, &
                                   &destruct_ho_params, update_ho_params
USE mo_output,              ONLY: init_output_files, write_output
!USE mo_oce_index,              ONLY: c_b, c_i, c_k, ldbg, form4ar,
USE mo_oce_index,              ONLY: init_index_test

USE mo_interpolation,          ONLY: t_int_state
USE mo_oce_init,               ONLY: init_ho_testcases
USE mo_oce_diagnostics,        ONLY: calculate_oce_diagnostics,&
                                  & construct_oce_diagnostics,&
                                  & destruct_oce_diagnostics, t_oce_timeseries

USE mo_oce_forcing,            ONLY: construct_sfcflx , &
                                  & construct_atmos_for_ocean,&
                                  & destruct_atmos_for_ocean,&
                                  & construct_atmos_fluxes, destruct_atmos_fluxes,&
                                  & t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
USE mo_sea_ice,                ONLY:t_sea_ice, construct_sea_ice, destruct_sea_ice
USE mo_oce_bulk,               ONLY: update_sfcflx
USE mo_oce_thermodyn,          ONLY:calc_density_MPIOM_func, calc_density_lin_EOS_func,&
                                    &calc_density_JMDWFG06_EOS_func
!USE mo_physical_constants,     ONLY: kice

IMPLICIT NONE

PRIVATE
INTEGER, PARAMETER :: kice = 1

INCLUDE 'cdi.inc'

!VERSION CONTROL:
CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

! public interface
!
! public subroutines
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
  SUBROUTINE perform_ho_stepping( ppatch, pstate_oce, p_ext_data, &
                                & datetime, n_io, n_file, p_int,  &
                                & p_sfc_flx, p_phys_param, p_as, p_atm_f, p_ice)

  TYPE(t_patch),             TARGET, INTENT(IN)     :: ppatch(n_dom)
  TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT)  :: pstate_oce(n_dom)
  TYPE(t_external_data), TARGET, INTENT(IN)         :: p_ext_data(n_dom)
  TYPE(t_datetime), INTENT(INOUT)                   :: datetime
  INTEGER, INTENT(IN)                               :: n_io, n_file
  TYPE(t_int_state),TARGET,INTENT(IN), OPTIONAL     :: p_int(n_dom)
  TYPE(t_sfc_flx)                                   :: p_sfc_flx
  TYPE (t_ho_params)                                :: p_phys_param 
  TYPE(t_atmos_for_ocean),  INTENT(INOUT)           :: p_as
  TYPE(t_atmos_fluxes ),       INTENT(INOUT)        :: p_atm_f
  TYPE (t_sea_ice),             INTENT(INOUT)       :: p_ice



  ! local variables
  INTEGER :: jstep, jt, jg, n_temp
  INTEGER :: jfile
  TYPE(t_oce_timeseries), POINTER :: oce_ts

  !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
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

  CALL init_index_test( ppatch, pstate_oce, p_ext_data )
  IF ( iswm_oce == 1 ) THEN
    CALL construct_oce_diagnostics( ppatch(jg), pstate_oce(jg), p_ext_data(jg), oce_ts)
  ENDIF

  IF (ltimer) CALL timer_start(timer_total)

  !------------------------------------------------------------------
  ! call the dynamical core: start the time loop
  !------------------------------------------------------------------
  TIME_LOOP: DO jstep = 1, nsteps

    write(*,*)'-----------timestep:',jstep

    IF(itestcase_oce==28)THEN
      CALL advect_tracer_ab(ppatch(jg), pstate_oce(jg), p_phys_param,p_sfc_flx, jstep)
    ELSE

    !In case of a time-varying forcing: 
    CALL update_sfcflx(ppatch(jg), pstate_oce(jg), p_as, p_ice, p_atm_f, p_sfc_flx, jstep)

    SELECT CASE (EOS_TYPE)
      CASE(1)

        CALL update_ho_params(ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param,&
                            & calc_density_lin_EOS_func)

     CASE(2)
       CALL update_ho_params(ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param,&
                             & calc_density_MPIOM_func)

    CASE(3)
       CALL update_ho_params(ppatch(jg), pstate_oce(jg), p_sfc_flx, p_phys_param,&
                             & calc_density_JMDWFG06_EOS_func)

    CASE DEFAULT

   END SELECT


    ! solve for new free surface
    CALL solve_free_surface_eq_ab (ppatch(jg), pstate_oce(jg), p_ext_data(jg), &
      &                            p_sfc_flx, p_phys_param, jstep, p_int(jg))

    ! Step 4: calculate final normal velocity from predicted horizontal velocity vn_pred
    !         and updated surface height
    CALL calc_normal_velocity_ab(ppatch(jg), pstate_oce(jg), p_ext_data(jg), p_phys_param)

    IF ( iswm_oce /= 1 ) THEN
      ! Step 5: calculate vertical velocity from continuity equation under incompressiblity condition
      CALL calc_vert_velocity( ppatch(jg), pstate_oce(jg))
    ENDIF

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
    CALL calculate_oce_diagnostics( ppatch(jg),    &
                                  & pstate_oce(jg),&
                                  & p_sfc_flx,     &
                                  & p_phys_param,  &
                                  & jstep,         &
                                  & oce_ts)
   ENDIF 
    ! Step 8: test output
!     IF (ldbg) WRITE(*,'(a,i5,2(a,g20.12))') '*** After jstep = ',jstep, &
!       &  '  Elevation h =', pstate_oce(jg)%p_prog(nold(jg))%h(c_i,c_b), &
!       &  '  Velocity  u =', pstate_oce(jg)%p_diag%u_pred(c_i,c_k,c_b)
     IF ( istime4output(jstep) .OR. jstep==nsteps ) THEN
      CALL message (TRIM(routine),'Write output at:')
      CALL print_datetime(datetime)

      CALL write_output( datetime )
    END IF

    ! close the current output file and trigger a new one
    ! #slo#: not synchronized with write_vlist - should be closed/renamed/opened
    !        before new writing timestep!
    IF (istime4newoutputfile(jstep) .AND. jstep/=nsteps) THEN

      jfile = jfile +1
      CALL init_output_files(jfile,lclose=.TRUE.)

    END IF

    ! One integration cycle finished on the lowest grid level (coarsest
    ! resolution). Set model time.
    CALL add_time(dtime,0,0,0,datetime)

!     CALL message (TRIM(routine),'Step completed at:')
!     CALL print_datetime(datetime)

  ENDDO TIME_LOOP
  IF ( iswm_oce == 1 ) THEN
    CALL destruct_oce_diagnostics(oce_ts)
   ENDIF

!   IF (ldbg) THEN
!     WRITE(*,*)  ' After run:'
!     WRITE(*,form4ar) ' Elevation h at cell   =', pstate_oce(jg)%p_prog(jt)%h(c_i,c_b),  &
!       &              '  Tracer 1 =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1)
!   ENDIF

  IF (ltimer) CALL timer_stop(timer_total)

  END SUBROUTINE perform_ho_stepping
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
  SUBROUTINE prepare_ho_integration(ppatch, pstate_oce, p_ext_data, p_sfc_flx, &
                                  & p_phys_param, p_as,&
                                  & p_atm_f, p_ice)

    TYPE(t_patch),                INTENT(INOUT)  :: ppatch(n_dom)
    TYPE(t_hydro_ocean_state),    INTENT(INOUT)  :: pstate_oce(n_dom)
    TYPE(t_external_data),        INTENT(INOUT)  :: p_ext_data(n_dom)
    TYPE(t_sfc_flx),              INTENT(INOUT)  :: p_sfc_flx
    TYPE (t_ho_params),           INTENT(INOUT)  :: p_phys_param 
    !TYPE(t_ho_physics),           INTENT(INOUT)  :: p_physics_oce 
    TYPE(t_atmos_for_ocean ),     INTENT(INOUT)  :: p_as
    TYPE(t_atmos_fluxes ),        INTENT(INOUT)  :: p_atm_f
    TYPE (t_sea_ice),             INTENT(INOUT)  :: p_ice

    ! local variables
    !TYPE(t_hydro_ocean_base)                  :: p_base
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

    ! hydro_ocean_base contains the 3-dimensional structures for the ocean state
    CALL construct_hydro_ocean_base(ppatch(jg), v_base)
    CALL init_ho_base(ppatch(jg), p_ext_data(jg), v_base)
    CALL init_scalar_product_base( ppatch(jg), v_base )
    CALL init_geo_factors_base( ppatch(jg), v_base )

    !------------------------------------------------------------------
    ! construct ocean state and physics
    !------------------------------------------------------------------

    ! ppatch and pstate_oce have dimension n_dom
    CALL construct_hydro_ocean_state(ppatch, pstate_oce)

    CALL construct_ho_params(ppatch(jg), p_phys_param)
    CALL init_ho_params(p_phys_param)

    !CALL construct_ho_physics(ppatch(jg), p_physics_oce)
    !CALL init_ho_physics(p_physics_oce)

    !------------------------------------------------------------------
    ! construct ocean forcing and testcases
    !------------------------------------------------------------------

    CALL  construct_sfcflx(ppatch(jg), p_sfc_flx)

    CALL construct_sea_ice(ppatch(jg), p_ice, kice)
    CALL construct_atmos_for_ocean(ppatch(jg), p_as)
    CALL construct_atmos_fluxes(ppatch(jg), p_atm_f, kice)

    CALL init_ho_testcases(ppatch(jg), pstate_oce(jg), p_ext_data(jg), p_sfc_flx)
  ! CALL init_ho_testcases(ppatch(jg), pstate_oce(jg), p_sfc_flx)

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
  SUBROUTINE finalise_ho_integration(p_os, p_phys_param, p_as, p_atm_f, p_ice)
    TYPE(t_hydro_ocean_state), INTENT(INOUT) :: p_os(n_dom)
    TYPE (t_ho_params),        INTENT(INOUT) :: p_phys_param 
    TYPE(t_atmos_for_ocean),   INTENT(INOUT) :: p_as
    TYPE(t_atmos_fluxes ),     INTENT(INOUT) :: p_atm_f
    TYPE (t_sea_ice),          INTENT(INOUT) :: p_ice


    !TYPE(t_ho_physics),        INTENT(INOUT)  :: p_physics_oce

    !------------------------------------------------------------------
    ! destruct ocean physics and forcing
    ! destruct ocean state is in control_model
    !------------------------------------------------------------------
     CALL destruct_hydro_ocean_state(p_os)
     CALL destruct_ho_params(p_phys_param)

     CALL destruct_sea_ice(p_ice)
     CALL destruct_atmos_for_ocean(p_as)
     CALL destruct_atmos_fluxes(p_atm_f)


  END SUBROUTINE finalise_ho_integration


 

END MODULE mo_hydro_ocean_run

