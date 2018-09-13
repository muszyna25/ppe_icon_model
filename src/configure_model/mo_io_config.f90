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
  USE mo_cdi,                     ONLY: FILETYPE_NC2
  USE mo_exception,               ONLY: finish, message
  USE mo_impl_constants,          ONLY: max_dom
  USE mo_io_units,                ONLY: filename_max
  USE mo_kind,                    ONLY: wp
  USE mo_parallel_config,         ONLY: num_restart_procs
  USE mo_run_config,              ONLY: dtime
  USE mo_util_string,             ONLY: int2string

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

  LOGICAL :: linvert_dict               !< inverts columns in output_nml_dict (allows using the same dictionary file as for input)

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

  CHARACTER(LEN = 256) :: restart_write_mode

  ! When using the restart write mode "dedicated proc mode", it is
  ! possible to split the restart output into several files, as if
  ! "nrestart_streams" * "num_io_procs" restart processes were
  ! involved. This speeds up the read-in process, since all the files
  ! may then be read in parallel.
  INTEGER :: nrestart_streams

  ! constants to communicate which restart writing MODULE to USE
  ENUM, BIND(C)
    ENUMERATOR :: kSyncRestartModule = 1, kAsyncRestartModule, kMultifileRestartModule
  END ENUM

  CHARACTER(*), PARAMETER :: modname = "mo_io_config"

  INTEGER, PARAMETER :: ALL_WORKERS_INVOLVED = -1


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
     IF (n_checkpoints == 0) n_checkpoints = HUGE(1)
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

  ! Inquire the parameters for restart writing.
  !
  ! opt_dedicatedProcCount: The number of processes that are split off
  !                         the MPI_Communicator to serve as dedicated
  !                         restart processes.
  !
  ! opt_restartProcCount:   The number of processes that actually
  !                         perform restart writing. This IS never zero,
  !                         especially NOT when opt_dedicatedProcCount
  !                         IS zero.
  !
  ! opt_restartModule:      The id of the MODULE that handles the restart
  !                         writing. Possible values are
  !                         kSyncRestartModule, kAsyncRestartModule, AND
  !                         kMultifileRestartModule.
  !
  ! opt_lDedicatedProcMode: Whether the restart processes are split
  !                         off the MPI_Communicator OR are a subset
  !                         of mpi_comm_work.
  !
  SUBROUTINE restartWritingParameters(opt_dedicatedProcCount, &
    &                                 opt_restartProcCount,   &
    &                                 opt_restartModule,      &
    &                                 opt_lDedicatedProcMode, &
    &                                 opt_nrestart_streams)

    INTEGER, INTENT(OUT), OPTIONAL :: opt_dedicatedProcCount, &
      &                               opt_restartProcCount,   &
      &                               opt_restartModule,      &
      &                               opt_nrestart_streams
    LOGICAL, INTENT(OUT), OPTIONAL :: opt_lDedicatedProcMode

    LOGICAL, SAVE :: cacheValid = .FALSE., lDedicatedProcMode = .FALSE.
    INTEGER, SAVE :: dedicatedProcCount = -1, &
      &              restartProcCount   = -1, &
      &              restartModule      = -1, &
      &              nrestartStreams    = -1
    CHARACTER(:), ALLOCATABLE :: errorMessage
    CHARACTER(*), PARAMETER :: routine = modname//":restartWritingParameters"

    ! If this IS the first CALL of this FUNCTION, analyze the namelist
    ! parameters AND cache the RESULT, sanity checking the settings.
    IF(.NOT.cacheValid) THEN
      IF(num_restart_procs < 0) THEN
        ! No matter what, negative process counts are illegal.
        errorMessage = "illegal value of namelist parameter num_restart_procs: value must not be negative"
        
      ELSE IF(restart_write_mode == "") THEN
        ! No restart_write_mode given, so we fall back to the old
        ! behavior of switching between sync/async restart mode based
        ! on num_restart_procs for backwards compatibility.
        dedicatedProcCount = num_restart_procs
        restartProcCount   = MAX(1, num_restart_procs)
        lDedicatedProcMode = num_restart_procs > 0
        nrestartStreams    = 1
        IF(lDedicatedProcMode) THEN
          restartModule = kAsyncRestartModule
        ELSE
          restartModule = kSyncRestartModule
        END IF
        
      ELSE IF(restart_write_mode == "sync") THEN
        IF(num_restart_procs /= 0) THEN
          errorMessage = "inconsistent namelist parameters: num_restart_procs must be zero OR unset &
            &for the given restart_write_mode"
        END IF
        dedicatedProcCount = 0
        restartProcCount   = 1
        restartModule      = kSyncRestartModule
        lDedicatedProcMode = .FALSE.
        nrestartStreams    = 1
        
      ELSE IF(restart_write_mode == "async") THEN
        IF(num_restart_procs == 0) THEN
          errorMessage = "inconsistent namelist parameters: num_restart_procs must be non-zero for &
            &the given restart_write_mode"
        END IF
        dedicatedProcCount = num_restart_procs
        restartProcCount   = num_restart_procs
        restartModule      = kAsyncRestartModule
        lDedicatedProcMode = .TRUE.
        nrestartStreams    = 1

      ELSE IF(restart_write_mode == "joint procs multifile") THEN
        dedicatedProcCount = 0
        
        ! if not set otherwise (num_restart_procs=0), set the
        ! number of restart PEs to the number of worker
        ! PEs... this is done later, since the number of workers
        ! is not yet available.
        IF (num_restart_procs == 0) THEN
          restartProcCount = ALL_WORKERS_INVOLVED
        ELSE
          restartProcCount = num_restart_procs
        END IF
        restartModule      = kMultifileRestartModule
        lDedicatedProcMode = .FALSE.
        nrestartStreams    = 1
        
      ELSE IF(restart_write_mode == "dedicated procs multifile") THEN
        IF(num_restart_procs == 0) THEN
          errorMessage = "inconsistent namelist parameters: num_restart_procs must be non-zero for &
            &the given restart_write_mode"
        END IF
        dedicatedProcCount = num_restart_procs
        restartProcCount   = num_restart_procs
        restartModule      = kMultifileRestartModule
        lDedicatedProcMode = .TRUE.
        nrestartStreams    = nrestart_streams
        
      ELSE
        errorMessage = "illegal value of namelist parameter restart_write_mode: expected one of 'sync', 'async', &
          &'joint procs multifile', or 'dedicated procs multifile'"
        
      END IF
      IF(ALLOCATED(errorMessage)) THEN
        CALL finish(routine, errorMessage//" (got restart_write_mode = '"//TRIM(restart_write_mode)// &
          &"', num_restart_procs = "//TRIM(int2string(num_restart_procs))//")")
      END IF
      cacheValid = .TRUE.
    END IF
    
    ! Ok, the cache IS up to date. What info did the caller want again?
    IF(PRESENT(opt_dedicatedProcCount))  opt_dedicatedProcCount = dedicatedProcCount;
    IF(PRESENT(opt_restartProcCount))    opt_restartProcCount   = restartProcCount;
    IF(PRESENT(opt_restartModule))       opt_restartModule      = restartModule;
    IF(PRESENT(opt_lDedicatedProcMode))  opt_lDedicatedProcMode = lDedicatedProcMode;
    IF(PRESENT(opt_nrestart_streams))    opt_nrestart_streams   = nrestartStreams

  END SUBROUTINE restartWritingParameters

END MODULE mo_io_config
