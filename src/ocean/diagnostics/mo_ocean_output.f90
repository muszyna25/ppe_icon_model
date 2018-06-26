!>
!! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!!
!! @author Peter Korn, Stephan Lorenz, MPI
!!
!! @par Revision History
!!  Initial version by Stephan Lorenz (MPI-M), (2010-04).
!!   - renaming and adjustment of hydrostatic ocean model V1.0.3 to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2010-10
!!   - new module mo_ocean_output including updated reconstructions
!
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_output
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_grid_subset,            ONLY: get_index_range
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, gibraltar, lhamocc
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: timeSteps_per_outputStep
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  !USE mo_io_units,               ONLY: filename_max
  USE mo_ocean_types,              ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog
  USE mo_ocean_state,              ONLY: ocean_restart_list
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface
  !USE mo_ocean_physics,            ONLY: t_ho_params
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_ocean_diagnostics,        ONLY: calc_slow_oce_diagnostics, calc_fast_oce_diagnostics, &
    & destruct_oce_diagnostics, calc_moc, calc_psi
  USE mo_linked_list,            ONLY: t_list_element
  USE mo_var_list,               ONLY: print_var_list, find_list_element
  USE mo_mpi,                    ONLY: my_process_is_stdio
  USE mo_statistics
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mo_ocean_statistics
  USE mtime,                     ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToPosixString
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: output_ocean

  INTEGER :: nsteps_since_last_output = 0

  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! @par Revision History
  !! Developed by Peter Korn, MPI-M  (2008-2010).
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
!<Optimize:inUse>
  SUBROUTINE output_ocean( &
    & patch_3d,        &
    & ocean_state,     &
    & this_datetime,   &
    & surface_fluxes,  &
    & sea_ice,         &
    & jstep, jstep0,   &
    & force_output)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: surface_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: sea_ice
    INTEGER,   INTENT(in)                            :: jstep, jstep0
    LOGICAL, OPTIONAL                                :: force_output
   
    ! local variables
    LOGICAL :: use_force_output
    INTEGER :: jg, jtrc, out_step
    INTEGER :: ocean_statistics
    !LOGICAL                         :: l_outputtime
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    INTEGER, POINTER :: dolic(:,:)
    REAL(wp), POINTER :: prism_thickness(:,:,:)
    
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datestring
    CHARACTER(len=32) :: fmtstr
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_output:output_ocean'
    !------------------------------------------------------------------

    patch_2D      => patch_3d%p_patch_2d(1)
    jg = 1
    nsteps_since_last_output = nsteps_since_last_output + 1
    !------------------------------------------------------------------
    use_force_output = .false.
    IF (PRESENT(force_output)) &
      use_force_output = force_output

    out_step = jstep
    IF (use_force_output) THEN
      out_step = jstep - nsteps_since_last_output + timeSteps_per_outputStep
    ENDIF

   !write(0,*) "out_step=", jstep, nsteps_since_last_output, timeSteps_per_outputStep, out_step
   IF (.not. istime4name_list_output(jstep) )  RETURN
   !write(0,*) "write ....."

    !------------------------------------------------------------------
   !CALL calc_slow_oce_diagnostics( patch_3d       , &
   !  &                             ocean_state(jg), &
   !  &                             surface_fluxes      , &
   !  &                             sea_ice          , &
   !  &                             jstep-jstep0   , &
   !  &                             this_datetime) ! , &
          ! &                             oce_ts)
  
!   IF (diagnostics_level > 0 ) THEN
!     IF (no_tracer>=2) THEN
!       CALL calc_moc (patch_2d,patch_3d, ocean_state(jg)%p_diag%w(:,:,:), this_datetime)
!       CALL calc_moc (patch_2d,patch_3d, ocean_state(jg)%p_acc%w(:,:,:), this_datetime)
!     ENDIF
!   ENDIF

    IF (output_mode%l_nml) CALL write_name_list_output(out_step)

    fmtstr = '%Y-%m-%d %H:%M:%S'
    call datetimeToPosixString(this_datetime, datestring, fmtstr)
    WRITE(message_text,'(a,a)') 'Write output at:', TRIM(datestring)
    CALL message (TRIM(routine),message_text)
  END SUBROUTINE output_ocean
  !-------------------------------------------------------------------------

END MODULE mo_ocean_output
