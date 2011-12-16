!>
!! @brief workflow for the ICON atmospheric hydrostatic model
!!
!! @author Hui Wan (MPI-M)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_atmo_hydrostatic

  USE mo_exception,         ONLY: message
  USE mo_impl_constants,    ONLY: iecham, ildf_echam
  USE mo_timer,             ONLY: print_timer

  USE mo_master_control,  ONLY: is_restart_run
  USE mo_time_config,       ONLY: time_config
  USE mo_run_config,        ONLY: dtime, nsteps, ltestcase, ltimer,iforcing, nlev, &
    & msg_level
  USE mo_ha_testcases,      ONLY: ctest_name
  USE mo_io_config,         ONLY: n_diags, n_checkpoints,n_files,n_ios,lwrite_initial

  USE mo_model_domain,        ONLY: p_patch
  USE mo_intp_data_strc,      ONLY: p_int_state, p_int_state_lonlat
  USE mo_grf_intp_data_strc,  ONLY: p_grf_state

  USE mo_vertical_coord_table,ONLY: vct_a, vct_b, ceta
  USE mo_icoham_dyn_memory,   ONLY: p_hydro_state, destruct_icoham_dyn_state
  USE mo_ha_stepping,         ONLY: prepare_ha_dyn, initcond_ha_dyn, &
                                    perform_ha_stepping

  USE mo_echam_phy_init,      ONLY: prepare_echam_phy, initcond_echam_phy, &
                                    additional_restart_init
  USE mo_echam_phy_cleanup,   ONLY: cleanup_echam_phy

  USE mo_io_restart,           ONLY: read_restart_files
  USE mo_io_restart_attributes,ONLY: get_restart_attribute
  USE mo_output,               ONLY: init_output_files, close_output_files,&
                                     write_output
  USE mo_parallel_config,      ONLY: use_icon_comm
  USE mo_icon_comm_interface,  ONLY: construct_icoham_communication, &
    & destruct_icoham_communication

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: atmo_hydrostatic

CONTAINS
  !>
  !!
  !!
  SUBROUTINE atmo_hydrostatic

    LOGICAL :: l_have_output
    INTEGER :: n_io, n_file, n_diag, n_chkpt
    INTEGER :: jfile
#ifndef NOMPI
    INTEGER :: jg
#endif

    CHARACTER(*), PARAMETER :: routine = "atmo_hydrostatic"

    !------------------------------------------------------------------
    ! Initialize parameters and solvers;
    ! Allocate memory for model state vectors.
    !------------------------------------------------------------------

    CALL prepare_ha_dyn( p_patch(1:) )

    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
      CALL prepare_echam_phy( p_patch(1:), ltestcase, ctest_name, &
                            & nlev, vct_a, vct_b, ceta )
    END IF

    !-------------------------------------------------------------------
    ! Initialize icon_comm_lib
    !-------------------------------------------------------------------
    IF (use_icon_comm) THEN
      CALL construct_icoham_communication()
    ENDIF
    
    !------------------------------------------------------------------
    ! Set initial conditions for time integration.
    !------------------------------------------------------------------
    ! Here constant values are also initialized,
    ! such as anayltical topography.
    ! It should be caleed even in reastart mode
    CALL initcond_ha_dyn( p_patch(1:), p_int_state(1:),  &
                        & p_grf_state(1:), p_hydro_state )

    IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) &
      CALL initcond_echam_phy( p_patch(1:),p_hydro_state, ltestcase, ctest_name )

    IF (is_restart_run()) THEN
    ! This is an resumed integration. Read model state from restart file(s).

#ifdef NOMPI
      CALL read_restart_files
#else
      jg = 1
     !DO jg = n_dom_start,n_dom
        CALL read_restart_files( p_patch(jg) )
     !END DO
#endif
      CALL message(TRIM(routine),'normal exit from read_restart_files')

      ! Re-initialize SST, sea ice and glacier for certain experiments; 
      ! Initialize logical variables in echam physics state.
      ! The latter is necessary for now because logical arrays can not yet be
      ! written into restart files.

      IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) THEN
        CALL additional_restart_init( p_patch(1:), ltestcase, ctest_name )
      END IF
!     ELSE
!     ! This is an initial run (cold start). Compute initial condition for
!     ! test cases, or read externally given initial conditions.
! 
!       CALL initcond_ha_dyn( p_patch(1:), p_int_state(1:),  &
!                           & p_grf_state(1:), p_hydro_state )
! 
!       IF (iforcing==IECHAM.OR.iforcing==ILDF_ECHAM) &
!       CALL initcond_echam_phy( p_patch(1:),p_hydro_state, ltestcase, ctest_name )

    END IF ! is_restart_run()

    !------------------------------------------------------------------
    ! The most primitive event handling algorithm:
    ! compute time step interval for taking a certain action
    !------------------------------------------------------------------

    n_io    = n_ios()        ! number of: write output
    n_file  = n_files()        ! number of: trigger new output file
    n_chkpt = n_checkpoints()       ! number of: write restart files
    n_diag  = n_diags() ! number of: diagnose of total integrals

    !------------------------------------------------------------------
    ! Initialize output file if necessary;
    ! Write out initial conditions.
    !------------------------------------------------------------------

    IF (.NOT.is_restart_run()) THEN
    ! Initialize the first output file which will contain also the
    ! initial conditions.

      jfile = 1
      CALL init_output_files(jfile, lclose=.FALSE.)
      IF (lwrite_initial) CALL write_output( time_config%cur_datetime )
      l_have_output = .TRUE.

    ELSE
    ! No need to write out the initial condition, thus no output
    ! during the first integration step. This run will produce
    ! output if n_io <= integration_length.

      CALL get_restart_attribute('next_output_file',jfile)

      IF (n_io.le.(nsteps-1)) THEN
         CALL init_output_files(jfile, lclose=.FALSE.)
         l_have_output = .TRUE.
      ELSE
         l_have_output = .FALSE.
      END IF

    END IF ! (not) is_restart_run()

    !------------------------------------------------------------------
    ! Time integraion
    !------------------------------------------------------------------

    CALL perform_ha_stepping( p_patch(1:), p_int_state(1:), &
                            & p_grf_state(1:),                               &
                            & p_hydro_state, time_config%cur_datetime,       &
                            & n_io, n_file, n_chkpt, n_diag, jfile,          &
                            & l_have_output                                  )

    IF (ltimer) CALL print_timer

    !---------------------------------------------------------------------
    ! Integration finished. Start to clean up.
    !---------------------------------------------------------------------
    CALL message(TRIM(routine),'start to clean up')

    CALL destruct_icoham_dyn_state
    
    IF (use_icon_comm) THEN
      CALL destruct_icoham_communication()
    ENDIF

    IF (iforcing==IECHAM .OR. iforcing==ILDF_ECHAM) THEN
      CALL cleanup_echam_phy
    ENDIF

    IF (msg_level > 5) CALL message(TRIM(routine),'echam_phy clean up is done')
    
    IF (l_have_output) CALL close_output_files
    
    IF (msg_level > 5) CALL message(TRIM(routine),'close_output_files is done')

  END SUBROUTINE atmo_hydrostatic
  !-------------

END MODULE mo_atmo_hydrostatic

