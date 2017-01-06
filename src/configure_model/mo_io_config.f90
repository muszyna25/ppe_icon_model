!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
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
MODULE mo_io_config

  USE mo_kind,           ONLY: wp
  USE mo_run_config,     ONLY: dtime
  USE mo_io_units,       ONLY: filename_max
  USE mo_cdi,            ONLY: FILETYPE_NC2
  USE mo_impl_constants, ONLY: max_dom

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Derived type
  !--------------------------------------------------------------------------

  ! from namelist

  LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
  REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
  REAL(wp):: gust_interval(max_dom)     ! time interval over which maximum wind gusts are taken
  REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

  INTEGER :: inextra_2d                 ! number of extra output fields for debugging
  INTEGER :: inextra_3d                 ! number of extra output fields for debugging

  LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated
                                        ! from the beginning of the run
                                        ! if .TRUE. the output fluxex are average values
                                        ! from the beginning of the run, except of
                                        ! TOT_PREC that would be accumulated

  INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure
  INTEGER :: itype_rh                   ! Specifies method for computation of relative humidity

  CHARACTER(LEN=filename_max) :: &
    &        output_nml_dict,    &      !< maps variable names onto the internal ICON names.
    &        netcdf_dict                !< maps internal variable names onto names in output file (NetCDF only).

  LOGICAL :: lnetcdf_flt64_output       !< if .TRUE. floating point valued NetCDF output
                                        !  is written in 64-bit instead of 32-bit accuracy

  ! derived variables
  !
  INTEGER, PARAMETER :: read_netcdf_broadcast_method  = 1
  INTEGER, PARAMETER :: read_netcdf_distribute_method = 2
  INTEGER :: default_read_method = 2

  INTEGER :: restart_file_type = FILETYPE_NC2

  LOGICAL :: write_initial_state = .true.
  LOGICAL :: write_last_restart  = .false.
  INTEGER :: timeSteps_per_outputStep    = 0

  INTEGER :: n_chkpt           ! number of timesteps between successive checkpoint events
  INTEGER :: n_diag            ! number of timesteps between successive tot_int diag events


  ! currently used by hydrostatic model only
  LOGICAL :: l_outputtime      ! if .true., output is written at the end of the time step.
  LOGICAL :: l_diagtime        ! if .true., diagnostic output is computed and written at the end of the time step.

  LOGICAL :: lmask_boundary    ! flag: true, if interpolation zone should be masked *in output*

CONTAINS

  !>
  !! Set up derived components of the I/O config state
  !!
  !! Set up derived components of the I/O config state. This routine is
  !! called, after all namelists have been read and a synoptic consistency
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-11-28)
  !!
  SUBROUTINE configure_io()

    !-----------------------------------------------------------------------

    ! number of timesteps between successive checkpoint events
    n_chkpt = n_checkpoints()

    ! number of timesteps between successive tot_int diag events
    n_diag  = n_diags()

  END SUBROUTINE configure_io



  !----------------------------------------------------------------------------------
   FUNCTION n_checkpoints()

     INTEGER :: n_checkpoints

     n_checkpoints = NINT(dt_checkpoint/dtime)  ! write restart files
   END FUNCTION n_checkpoints
  !----------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------
   FUNCTION n_diags()

     INTEGER :: n_diags

     n_diags  = MAX(1,NINT(dt_diag/dtime)) ! number of: diagnose of total integrals
   END FUNCTION n_diags
  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------

   FUNCTION is_checkpoint_time(current_step, n_checkpoints, n_steps) RESULT(l_checkpoint)
     INTEGER, INTENT(IN)            :: current_step, n_checkpoints
     INTEGER, INTENT(IN), OPTIONAL  :: n_steps

     LOGICAL              :: l_checkpoint

     IF (n_checkpoints == 0) THEN
       l_checkpoint = .FALSE.
     ELSE
       IF (PRESENT(n_steps)) THEN
         IF ( MOD(current_step,n_checkpoints)==0 .AND. current_step/=n_steps ) THEN
           l_checkpoint = .TRUE.
         ELSE
           l_checkpoint = .FALSE.
         END IF
       ELSE
         IF ( MOD(current_step,n_checkpoints)==0 ) THEN
           l_checkpoint = .TRUE.
         ELSE
           l_checkpoint = .FALSE.
         END IF
       END IF
     END IF
   END FUNCTION is_checkpoint_time
  !----------------------------------------------------------------------------------

  !>
  !! Decides about diagnostic computation of total integrals
  !!
  !! Decides about diagnostic computation of total integrals, which
  !! is performed in "supervise_total_integrals_nh"
  !! Total integrals are computed
  !! - at the first time step (or the first time step after restart)
  !! - if (MOD(current_step,n_diag) == 0)
  !! - at the very last time step
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-12-01)
  !!
  FUNCTION is_totint_time(current_step, restart_step, n_diag, n_steps)

    INTEGER, INTENT(IN) :: current_step !< current time step number
    INTEGER, INTENT(IN) :: restart_step !< time step for which the restart file was produced
                                        !< rfile_step+1: first step after restart
    INTEGER, INTENT(IN) :: n_steps      !< total number of time steps
    INTEGER, INTENT(IN) :: n_diag       !< number of timesteps between successive calls

    LOGICAL :: is_totint_time           ! Result

    ! local variables
    INTEGER :: kstep                    ! time step number relative to restart step

    kstep = current_step - restart_step

    is_totint_time = ((kstep == 1)                        .OR. &
      &              (MOD(current_step,n_diag) == 0)      .OR. &
      &              (kstep==n_steps))                    .AND.&
      &              (current_step > 0 )
  END FUNCTION is_totint_time

END MODULE mo_io_config
