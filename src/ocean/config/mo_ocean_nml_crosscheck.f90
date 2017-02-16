!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ocean_nml_crosscheck

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish, warning
  USE mo_grid_config,       ONLY: init_grid_configuration
  USE mo_parallel_config,   ONLY: check_parallel_configuration, p_test_run, l_fast_sum
  USE mo_run_config,        ONLY: nsteps, dtime, nlev
  USE mo_time_config,       ONLY: time_config, dt_restart
  USE mo_io_config,         ONLY: dt_checkpoint, write_initial_state
  USE mo_grid_config,       ONLY: grid_rescale_factor, use_duplicated_connectivity
  USE mo_ocean_nml
  USE mo_master_config,     ONLY: isRestart
  USE mo_time_management,   ONLY: compute_timestep_settings,                        &
    &                             compute_restart_settings,                         &
    &                             compute_date_settings

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ocean_crosscheck

CONTAINS
  

!<Optimize:inUse>
  SUBROUTINE check_thicknesses

    ! ensure, that all used thicknesses are non-zero in case of non-shallow-water.
    !For shallow-water this test makes no sense since dzlev==0 
    IF(iswm_oce==0)THEN
      IF (MINVAL(dzlev_m(1:n_zlev)) <= 0.0_wp) THEN
        CALL finish("check_thicknesses","Found zero or negative thicknesses")
      END IF
    ENDIF  
  END SUBROUTINE check_thicknesses

!<Optimize:inUse>
  SUBROUTINE ocean_crosscheck()

    CHARACTER(len=*), PARAMETER :: method_name =  'mo_ocean_nml_crosscheck:ocean_crosscheck'

    CALL check_parallel_configuration()

    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("oce", dt_restart, nsteps)

    CALL init_grid_configuration

    ! set the patch-related nlev variable to the ocean setup n_ zlev
    nlev = n_zlev

    IF (p_test_run .AND. l_fast_sum ) THEN                      
       CALL warning(method_name, "p_test_run sets l_fast_sum=.f alse.")
       l_fast_sum = .false.                                     
    ENDIF                                                       
    
    SELECT CASE (select_solver)
      CASE (select_gmres)

      CASE (select_restart_gmres, select_restart_mixedPrecision_gmres)

        IF (p_test_run .OR. .NOT. l_fast_sum ) THEN
           CALL warning(method_name, "p_test_run .OR. .NOT. l_fast_sum cannot be used by the restart gmres solver")
           CALL message(method_name, "Using the standard gmres solver")
           select_solver = select_gmres
!        ELSE
!           use_absolute_solver_tolerance = .true.
        ENDIF

      CASE default
        CALL finish(method_name, "Unknown solver")

    END SELECT
    
    
    IF (no_tracer < 1) THEN
      CALL warning("ocean_crosscheck", "no_tracer < 1, use_constant_mixing")
      PPscheme_type = PPscheme_Constant_type
    ENDIF

    CALL check_thicknesses

    IF  (RichardsonDiffusion_threshold < convection_InstabilityThreshold) &
      CALL finish (method_name, "RichardsonDiffusion_threshold < convection_InstabilityThreshold")

     
    IF (l_rigid_lid .AND. iswm_oce /= 1) THEN
      CALL finish(method_name, "l_rigid_lid .AND. iswm_oce /= 1")
    ENDIF

    IF (use_duplicated_connectivity) THEN
      use_duplicated_connectivity = .FALSE.
      CALL message(method_name, "Set use_duplicated_connectivity to FALSE")
    ENDIF
    
    IF (isRestart() .AND. write_initial_state) THEN
      CALL warning(method_name, "write_initial_state is disbaled for restarts")
      write_initial_state = .false.
    ENDIF

    IF ((VelocityDiffusion_order == 21 .or. VelocityDiffusion_order == 213) .and. .not. laplacian_form == 1) &
      CALL finish(method_name,"harmonic+biharmonic velocity diffusion requires curl-curl form")

  END SUBROUTINE ocean_crosscheck


END MODULE mo_ocean_nml_crosscheck
