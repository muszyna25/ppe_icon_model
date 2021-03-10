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
  USE mo_impl_constants,          ONLY: max_dom, max_echotop
  USE mo_io_units,                ONLY: filename_max
  USE mo_kind,                    ONLY: wp
  USE mo_parallel_config,         ONLY: num_restart_procs
  USE mo_run_config,              ONLY: dtime
  USE mo_util_string,             ONLY: int2string
  USE mtime,                      ONLY: max_timedelta_str_len
  USE mo_name_list_output_config, ONLY: first_output_name_list, &
                                        is_variable_in_output, is_variable_in_output_dom
  USE mo_lnd_nwp_config,          ONLY: groups_smi

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Derived type
  !--------------------------------------------------------------------------

  ! from namelist

  LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
  REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
  REAL(wp):: gust_interval(max_dom)     ! time interval [seconds] over which maximum wind gusts are taken
  REAL(wp):: celltracks_interval(max_dom)  ! time interval [seconds] over which extrema of cell track vars are taken
                                           !  (LPI_MAX, UH_MAX, VORW_CTMAX, W_CTMAX, DBZ_CTMAX)
  CHARACTER(len=max_timedelta_str_len) :: precip_interval(max_dom)   ! time interval over which precipitation variables are accumulated
  CHARACTER(len=max_timedelta_str_len) :: runoff_interval(max_dom)   ! time interval over which runoff variables are accumulated
  CHARACTER(len=max_timedelta_str_len) :: maxt_interval(max_dom)     ! time interval for tmax_2m, tmin_2m
  REAL(wp):: dt_lpi                     ! calling frequency [seconds] of lpi diagnosis for hourly maximum calculation
  REAL(wp):: dt_celltracks              ! calling frequency [seconds] of celltrack diagnosis for hourly maximum calculation
                                        ! this pertains to the following variables: tcond_max/tcond10_max, uh_max, vorw_ctmax, w_ctmax
  REAL(wp):: dt_radar_dbz               ! calling frequency [seconds] of radar reflectivity diagnosis for hourly maximum calculation

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

  TYPE t_echotop_meta
    REAL(wp) :: time_interval           !< time interval [seconds] over which echotops are maximized/minimized
                                        !   (ECHOTOP, ECHOTOPinM) for a domain
    INTEGER  :: nechotop                !< number of actual echotop levels for domain
    REAL(wp) :: dbzthresh(max_echotop)  !< actual echotop levels from namelist
  END TYPE t_echotop_meta

  TYPE(t_echotop_meta) :: echotop_meta(max_dom)

  INTEGER :: bvf2_mode                  !< computation mode for square of Brunt-Vaisala frequency
  INTEGER :: parcelfreq2_mode           !< computation mode for square of general air parcel oscillation frequency

  ! Derived type to collect logical variables indicating if optional diagnostics are requested for output
  TYPE t_var_in_output
    LOGICAL :: pres_msl    = .FALSE. !< Flag. TRUE if computation of mean sea level pressure desired
    LOGICAL :: omega       = .FALSE. !< Flag. TRUE if computation of vertical velocity desired
    LOGICAL :: rh          = .FALSE. !< Flag. TRUE if computation of relative humidity desired
    LOGICAL :: pv          = .FALSE. !< Flag. TRUE if computation of potential vorticity desired
    LOGICAL :: sdi2        = .FALSE. !< Flag. TRUE if computation of supercell detection index desired
    LOGICAL :: lpi         = .FALSE. !< Flag. TRUE if computation of lightning potential index desired
    LOGICAL :: lpi_max     = .FALSE. !< Flag. TRUE if computation of max. of lightning potential index desired
    LOGICAL :: ceiling     = .FALSE. !< Flag. TRUE if computation of ceiling height desired
    LOGICAL :: hbas_sc     = .FALSE. !< Flag. TRUE if computation of height of base from shallow convection desired
    LOGICAL :: htop_sc     = .FALSE. !< Flag. TRUE if computation of height of top  from shallow convection desired
    LOGICAL :: twater      = .FALSE. !< Flag. TRUE if computation of total column integrated water desired
    LOGICAL :: q_sedim     = .FALSE. !< Flag. TRUE if computation of specific content of precipitation particles desired
    LOGICAL :: dbz         = .FALSE. !< Flag. TRUE if computation of radar reflectivity is desired
    LOGICAL :: dbz850      = .FALSE. !< Flag. TRUE if computation of radar reflectivity in approx. 850 hPa is desired
    LOGICAL :: dbzcmax     = .FALSE. !< Flag. TRUE if computation of radar reflectivity column maximum is desired
    LOGICAL :: dbzctmax    = .FALSE. !< Flag. TRUE if computation of radar reflectivity column and time maximum is desired
    LOGICAL :: echotop     = .FALSE. !< Flag. TRUE if computation of echo tops in hPa of radar reflectivity is desired
    LOGICAL :: echotopinm  = .FALSE. !< Flag. TRUE if computation of echo tops in m MSL of radar reflectivity is desired
    LOGICAL :: smi         = .FALSE. !< Flag. TRUE if computation of soil moisture index desired
    LOGICAL :: tcond_max   = .FALSE. !< Flag. TRUE if computation of total column-integrated condensate desired
    LOGICAL :: tcond10_max = .FALSE. !< Flag. TRUE if computation of total column-integrated condensate above z(T=-10 degC) desired
    LOGICAL :: uh_max      = .FALSE. !< Flag. TRUE if computation of updraft helicity desired
    LOGICAL :: vorw_ctmax  = .FALSE. !< Flag. TRUE if computation of maximum rotation amplitude desired
    LOGICAL :: w_ctmax     = .FALSE. !< Flag. TRUE if computation of maximum updraft track desired
    LOGICAL :: vor_u       = .FALSE. !< Flag. TRUE if computation of zonal component of relative vorticity desired
    LOGICAL :: vor_v       = .FALSE. !< Flag. TRUE if computation of meridional component of relative vorticity desired
    LOGICAL :: bvf2        = .FALSE. !< Flag. TRUE if computation of square of Brunt-Vaisala frequency desired
    LOGICAL :: parcelfreq2 = .FALSE. !< Flag. TRUE if computation of square of general parcel oscillation frequency desired
  END TYPE t_var_in_output

  TYPE(t_var_in_output), ALLOCATABLE :: var_in_output(:)
  
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
  !! Precomputation of derived type collecting logical variables indicating whether
  !! optional diagnostics are requested in the output namelists
  !!
  !! Replaces repeated calculations of the same that used to be scattered around various places in the model code
  !!
  !! @par Revision History
  !! Initial revision by Guenther Zaengl, DWD (2020-02-13)
  !!
  SUBROUTINE init_var_in_output(n_dom, lnwp)

    INTEGER, INTENT(in)  :: n_dom  ! number of model domains
    LOGICAL, INTENT(in)  :: lnwp   ! true if ICON runs in NWP mode, implying that the full set of variables 
                                   ! needs to be computed

    INTEGER :: jg, jgr

    ALLOCATE(var_in_output(n_dom))

    DO jg=1,n_dom
      var_in_output(jg)%pres_msl = is_variable_in_output(first_output_name_list, var_name="pres_msl") .OR. &
        &                          is_variable_in_output(first_output_name_list, var_name="psl_m")
      var_in_output(jg)%omega    = is_variable_in_output(first_output_name_list, var_name="omega")    .OR. &
        &                          is_variable_in_output(first_output_name_list, var_name="wap_m")
      var_in_output(jg)%vor_u    = is_variable_in_output_dom(first_output_name_list, var_name="vor_u", jg=jg)
      var_in_output(jg)%vor_v    = is_variable_in_output_dom(first_output_name_list, var_name="vor_v", jg=jg)
      var_in_output(jg)%bvf2     = is_variable_in_output_dom(first_output_name_list, var_name="bvf2", jg=jg)
      var_in_output(jg)%parcelfreq2 = is_variable_in_output_dom(first_output_name_list, var_name="parcelfreq2", jg=jg)
    END DO


    IF (lnwp) THEN
      DO jg=1,n_dom
        var_in_output(jg)%rh          = is_variable_in_output_dom(first_output_name_list, var_name="rh", jg=jg)
        var_in_output(jg)%pv          = is_variable_in_output_dom(first_output_name_list, var_name="pv", jg=jg)
        var_in_output(jg)%sdi2        = is_variable_in_output_dom(first_output_name_list, var_name="sdi2", jg=jg)
        var_in_output(jg)%lpi         = is_variable_in_output_dom(first_output_name_list, var_name="lpi", jg=jg)
        var_in_output(jg)%lpi_max     = is_variable_in_output_dom(first_output_name_list, var_name="lpi_max", jg=jg)
        var_in_output(jg)%ceiling     = is_variable_in_output_dom(first_output_name_list, var_name="ceiling", jg=jg)
        var_in_output(jg)%hbas_sc     = is_variable_in_output_dom(first_output_name_list, var_name="hbas_sc", jg=jg)
        var_in_output(jg)%htop_sc     = is_variable_in_output_dom(first_output_name_list, var_name="htop_sc", jg=jg)
        var_in_output(jg)%twater      = is_variable_in_output_dom(first_output_name_list, var_name="twater", jg=jg)
        var_in_output(jg)%q_sedim     = is_variable_in_output_dom(first_output_name_list, var_name="q_sedim", jg=jg)
        var_in_output(jg)%tcond_max   = is_variable_in_output_dom(first_output_name_list, var_name="tcond_max", jg=jg)
        var_in_output(jg)%tcond10_max = is_variable_in_output_dom(first_output_name_list, var_name="tcond10_max", jg=jg)
        var_in_output(jg)%uh_max      = is_variable_in_output_dom(first_output_name_list, var_name="uh_max", jg=jg)
        var_in_output(jg)%vorw_ctmax  = is_variable_in_output_dom(first_output_name_list, var_name="vorw_ctmax", jg=jg)
        var_in_output(jg)%w_ctmax     = is_variable_in_output_dom(first_output_name_list, var_name="w_ctmax", jg=jg)
        var_in_output(jg)%dbz         = is_variable_in_output_dom(first_output_name_list, var_name="dbz", jg=jg)
        var_in_output(jg)%dbz850      = is_variable_in_output_dom(first_output_name_list, var_name="dbz_850", jg=jg)
        var_in_output(jg)%dbzcmax     = is_variable_in_output_dom(first_output_name_list, var_name="dbz_cmax", jg=jg)
        var_in_output(jg)%dbzctmax    = is_variable_in_output_dom(first_output_name_list, var_name="dbz_ctmax", jg=jg)
        var_in_output(jg)%echotop     = is_variable_in_output_dom(first_output_name_list, var_name="echotop", jg=jg)
        var_in_output(jg)%echotopinm  = is_variable_in_output_dom(first_output_name_list, var_name="echotopinm", jg=jg)
        var_in_output(jg)%smi         = is_variable_in_output_dom(first_output_name_list, var_name="smi", jg=jg)

        ! Check for special case: SMI is not in one of the output lists but it is part of a output group.
        ! In this case, the group can not be checked, as the connection between SMI and the group will be
        ! established during the add_var call. However, add_var for SMI will only be called if l_smi =.true.
        ! As a crutch, a character array containing the output groups of SMI from mo_lnd_nwp_config is used
        ! here and also at the add_var call.
        ! The check loops through the output groups. It has to be checked if l_smi is already .true., to not
        ! overwrite an existing .true. with a false. 

        IF (.NOT. var_in_output(jg)%smi) THEN 
          ! Check for output groups containing SMI
          DO jgr = 1,SIZE(groups_smi)
            IF (.NOT. var_in_output(jg)%smi) THEN
              var_in_output(jg)%smi = is_variable_in_output_dom(first_output_name_list, &
                                      var_name='group:'//TRIM(groups_smi(jgr)), jg=jg)
            END IF
          END DO
        END IF
      END DO
    END IF

  END SUBROUTINE init_var_in_output


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
