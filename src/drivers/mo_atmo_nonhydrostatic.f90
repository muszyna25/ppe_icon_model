!>
!! @brief branch for the non-hydrostatic ICON workflow
!!
!! @author Kristina Froehlich, MPI-M (2011-07-19)
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
MODULE mo_atmo_nonhydrostatic

USE mo_kind,                 ONLY: wp
USE mo_exception,            ONLY: message, finish
USE mo_impl_constants,       ONLY: SUCCESS, max_dom
USE mo_timer,                ONLY: print_timer, timers_level, timer_start, &
  &                                timer_stop, timer_model_init
USE mo_master_control,       ONLY: is_restart_run
USE mo_output,               ONLY: init_output_files, close_output_files,&
  &                                write_output
USE mo_intp_rbf,             ONLY: rbf_vec_interpol_cell
USE mo_intp,                 ONLY: edges2cells_scalar
USE mo_time_config,          ONLY: time_config      ! variable
USE mo_io_restart,           ONLY: read_restart_files
USE mo_io_restart_attributes,ONLY: get_restart_attribute
USE mo_io_config,            ONLY: dt_data,dt_file,dt_diag,dt_checkpoint, &
  &                                lwrite_pzlev
USE mo_parallel_config,      ONLY: nproma
USE mo_nh_pzlev_config,      ONLY: configure_nh_pzlev
USE mo_run_config,           ONLY: &
  &                               dtime, dtime_adv,     & !    namelist parameter
  &                               ltestcase,            &
  &                               nsteps,               & !    :
  &                               ltimer,               & !    :
  &                               iforcing                !    namelist parameter
USE mo_dynamics_config,      ONLY: nnow, nnow_rcf
USE mo_impl_constants,       ONLY: inwp
USE mo_lnd_nwp_config,       ONLY: configure_lnd_nwp, nsfc_subs, p_tiles, &
  &                                destruct_tiles_arrays
! Horizontal grid
USE mo_model_domain,         ONLY: p_patch
USE mo_grid_config,          ONLY: n_dom
! to break circular dependency KF???
USE mo_intp_data_strc,       ONLY: p_int_state, p_int_state_lonlat
USE mo_grf_intp_data_strc,   ONLY: p_grf_state
! NH-namelist state
USE mo_atm_phy_nwp_config,   ONLY: configure_atm_phy_nwp, atm_phy_nwp_config
! NH-Model states
USE mo_nonhydro_state,       ONLY: p_nh_state, construct_nh_state, destruct_nh_state
USE mo_opt_diagnostics,      ONLY: p_nh_opt_diag, construct_opt_diag,          &
  &                                destruct_opt_diag
USE mo_nwp_phy_state,        ONLY: construct_nwp_phy_state,                    &
  &                                destruct_nwp_phy_state, prm_diag
USE mo_nwp_lnd_state,        ONLY: p_lnd_state, construct_nwp_lnd_state,       &
  &                                destruct_nwp_lnd_state
USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
! Time integration
USE mo_nh_stepping,          ONLY: prepare_nh_integration, perform_nh_stepping
! Initialization with real data
USE mo_prepicon_utils,      ONLY: init_prepicon, prepicon, copy_prepicon2prog, &
  &                               compute_coord_fields,  deallocate_prepicon
USE mo_prepicon_config,     ONLY: i_oper_mode, l_sfc_in
USE mo_nh_vert_interp,      ONLY: vertical_interpolation
USE mo_ext_data_state,      ONLY: ext_data
! meteogram output
USE mo_meteogram_output,    ONLY: meteogram_init, meteogram_finalize
USE mo_meteogram_config,    ONLY: meteogram_output_config
USE mo_name_list_output_config, ONLY: first_output_name_list
USE mo_name_list_output,    ONLY: init_name_list_output,  &
  &                               write_name_list_output, &
  &                               close_name_list_output
USE mo_nh_pzlev_config,     ONLY: nh_pzlev_config
USE mo_var_list,            ONLY: nvar_lists, var_lists
USE mo_pp_scheduler,        ONLY: t_simulation_status, new_simulation_status, &
  &                               pp_scheduler_init, pp_scheduler_process, pp_scheduler_finalize

!-------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_nonhydrostatic

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_nonhydrostatic"


    INTEGER :: jg, jfile, n_file, ist, n_diag, n_chkpt, ntl, ntlr
    LOGICAL :: l_have_output, l_realcase
    INTEGER :: pat_level(n_dom)
    TYPE(t_simulation_status) :: simulation_status

    IF (timers_level > 3) CALL timer_start(timer_model_init)

    DO jg=1,n_dom
       pat_level(jg)= p_patch(jg)%level
    ENDDO

    IF(iforcing == inwp) THEN
     CALL configure_atm_phy_nwp(n_dom, pat_level(:), ltestcase, dtime_adv )
    ENDIF

    IF (.NOT. ltestcase .AND. iforcing == inwp) THEN
      l_realcase = .TRUE.
    ELSE
      l_realcase = .FALSE.
    ENDIF
 
    !---------------------------------------------------------------------
    ! 4.c Non-Hydrostatic / NWP
    !---------------------------------------------------------------------

    ALLOCATE (p_nh_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_nh_state failed')
    ENDIF

    ! Note(GZ): Land state now needs to be allocated even if physics is turned
    ! off because ground temperature is included in feedback since r8133
    ! However, setting inwp_surface = 0 effects that only a few 2D fields are allocated
    ALLOCATE (p_lnd_state(n_dom), stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'allocation for p_lnd_state failed')
    ENDIF

    CALL configure_lnd_nwp(p_patch(1:), n_dom, nproma)
    IF(iforcing /= inwp) atm_phy_nwp_config(:)%inwp_surface = 0

      !------------------------------------------------------------------
      ! Prepare initial conditions for time integration.
      !------------------------------------------------------------------

    ! Initialize model with real atmospheric data if appropriate switches are set
    ! This has been shifted before the allocation of the NH and physics state
    ! in order to reduce the maximum amount of memory consumed
    IF (l_realcase .AND. .NOT. is_restart_run()) THEN

      CALL message(TRIM(routine),'Real-data mode: perform initialization with IFS2ICON data')

      ! Allocate prepicon data type
      ALLOCATE (prepicon(n_dom), stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'allocation for prepicon failed')
      ENDIF

      i_oper_mode = 3 ! For the time being, the only option that works when called
                      ! from the main ICON program

      IF (atm_phy_nwp_config(1)%inwp_surface > 0 .AND. .NOT. l_sfc_in) THEN
        CALL finish(TRIM(routine),'A real-data run with surface scheme requires &
                                   &surface input data')
      ENDIF
 
      ! allocate memory for topography and coordinate fields,
      ! read topo data from netCDF file, 
      ! optionally smooth topography data,
      ! and, in case of nesting, perform topography blending and feedback
      CALL init_prepicon (p_int_state(1:), p_grf_state(1:), prepicon, ext_data)

    ENDIF

    ! Now allocate memory for the states
    CALL construct_nh_state(p_patch(1:), p_nh_state, n_timelevels=2)

    ! Add optional diagnostic variable lists (might remain empty)
    CALL construct_opt_diag(p_patch(1:))

    IF(iforcing == inwp) THEN
      CALL construct_nwp_phy_state( p_patch(1:) )
    ENDIF

    CALL construct_nwp_lnd_state( p_patch(1:),p_lnd_state,n_timelevels=2 )

    !---------------------------------------------------------------------
    ! 5. Perform time stepping
    !---------------------------------------------------------------------
      !------------------------------------------------------------------
      ! Prepare for time integration
      !------------------------------------------------------------------


    CALL prepare_nh_integration(p_patch(1:), p_nh_state, p_int_state(1:), p_grf_state(1:))


    !
    ! Read restart files (if necessary)
    !
    IF (is_restart_run()) THEN
      ! This is a resumed integration. Read model state from restart file(s).

#ifdef NOMPI
      CALL read_restart_files
#else
      jg = 1
      !DO jg = n_dom_start,n_dom
      CALL read_restart_files( p_patch(jg) )
      !END DO
#endif      
      CALL message(TRIM(routine),'normal exit from read_restart_files')

    ENDIF


    !---------------------------------------------------------------------
    !     Setup of meteogram output
    !---------------------------------------------------------------------

    IF (.NOT. ltestcase) THEN
      DO jg =1,n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_init(meteogram_output_config(jg), jg, p_patch(jg), &
            &                ext_data(jg), p_nh_state(jg), prm_diag(jg),    &
            &                p_lnd_state(jg), iforcing)
        END IF
      END DO
    END IF

    ! Continue operations for real-data initialization
    IF (l_realcase .AND. .NOT. is_restart_run()) THEN

      ! Compute the 3D coordinate fields
      CALL compute_coord_fields(p_int_state(1:), prepicon)

      ! Perform vertical interpolation from intermediate IFS2ICON grid to ICON grid
      ! and convert variables to the NH set of prognostic variables
      CALL vertical_interpolation(p_patch(1:), p_int_state(1:), p_grf_state(1:), prepicon)

      ! Finally copy the results to the prognostic model variables
      CALL copy_prepicon2prog(prepicon, p_nh_state, p_lnd_state, ext_data(:))

      ! Deallocate prepicon data type
      CALL deallocate_prepicon(prepicon)
      DEALLOCATE (prepicon, stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'deallocation for prepicon failed')
      ENDIF

    ENDIF

    !---------------------------------------------------------
    ! The most primitive event handling algorithm: 
    ! compute time step interval for taking a certain action
    !--------------------------------------------------------- 
 
    ! writing output is now controlled via 'istime4output'
    n_file  = NINT(dt_file/dtime)        ! trigger new output file
    n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
    n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals


    !------------------------------------------------------------------
    ! Prepare output file
    !------------------------------------------------------------------
    !
    ! if output on z and/or p-levels is required do some config
    IF (lwrite_pzlev) THEN
      DO jg = 1, n_dom
        CALL configure_nh_pzlev(jg, nproma, p_patch(jg)%npromz_c,  &
          &                     p_patch(jg)%nblks_c, lwrite_pzlev  )
      ENDDO
    ENDIF

    IF (.NOT.is_restart_run()) THEN
      ! Initialize the first output file which will contain also the 
      ! initial conditions.

      jfile = 1
      CALL init_output_files(jfile, lclose=.FALSE.)

    ELSE
    ! No need to write out the initial condition, thus no output
    ! during the first integration step. This run will produce
    ! output if n_io <= integration_length. 

      CALL get_restart_attribute('next_output_file',jfile)

!DR may need to be fine-tuned later on
!DR      IF (n_io.le.(nsteps-1)) THEN
         CALL init_output_files(jfile, lclose=.FALSE.)
         l_have_output = .TRUE.  
!DR      ELSE
!DR         l_have_output = .FALSE.
!DR      END IF

    END IF

    ! setup of post-processing job queue, e.g. setup of optional
    ! diagnostic quantities like pz-level interpolation
    IF (.NOT. ltestcase) THEN
      CALL pp_scheduler_init(p_patch(1:), p_nh_state, prm_diag, p_nh_opt_diag, &
        &                    nh_pzlev_config(:), first_output_name_list,   & 
        &                    var_lists, nvar_lists )     
    ENDIF
    ! If async IO is in effect, init_name_list_output is a collective call
    ! with the IO procs and effectively starts async IO
    CALL init_name_list_output

    !------------------------------------------------------------------
    !  get and write out some of the inital values
    !------------------------------------------------------------------
    IF (.NOT.is_restart_run()) THEN

      ! diagnose u and v to have meaningful initial output
      ! For real-case runs, also diagnose pressure and temperature
      ! because these variables are needed for initializing the physics parameterizations

      DO jg = 1, n_dom

        ! time levels
        ntl  = nnow(jg)
        ntlr = nnow_rcf(jg)

        SELECT CASE (p_patch(jg)%cell_type)
          CASE (3)
          CALL rbf_vec_interpol_cell(p_nh_state(jg)%prog(ntl)%vn,p_patch(jg),&
            & p_int_state(jg),p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)
          CASE (6)
          CALL edges2cells_scalar(p_nh_state(jg)%prog(ntl)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_east,p_nh_state(jg)%diag%u)
          CALL edges2cells_scalar(p_nh_state(jg)%prog(ntl)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_north,p_nh_state(jg)%diag%v)
        END SELECT

        IF (l_realcase) THEN ! for test cases, diagnose_pres_temp is currently 
                             ! called in nh_testcases
          CALL diagnose_pres_temp(p_nh_state(jg)%metrics, p_nh_state(jg)%prog(ntl),     &
            &                  p_nh_state(jg)%prog(ntlr), p_nh_state(jg)%diag,          &
            &                  p_patch(jg), opt_calc_temp=.TRUE., opt_calc_pres=.TRUE., &
            &                  lnd_prog=p_lnd_state(jg)%prog_lnd(ntlr) )
        ENDIF

      ENDDO

      !--------------------------------------------------------------------------
      ! loop over the list of internal post-processing tasks, e.g.
      ! interpolate selected fields to p- and/or z-levels
      simulation_status = new_simulation_status(l_first_step =.TRUE., l_output_step=.TRUE.)
      IF (.NOT. ltestcase) &
        CALL pp_scheduler_process(simulation_status)

      ! Note: here the derived output variables are not yet available
      ! (divergence, vorticity)
      !
      CALL write_output( time_config%cur_datetime )
      l_have_output = .TRUE.
      IF (.NOT. ltestcase) &
        CALL write_name_list_output( time_config%cur_datetime, 0._wp, .FALSE. )

    END IF ! not is_restart_run()

    IF (timers_level > 3) CALL timer_stop(timer_model_init)

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

    CALL perform_nh_stepping( p_patch, p_int_state, p_grf_state, p_nh_state, &
      &                       time_config%cur_datetime,       &
      &                       n_file, jfile, n_chkpt, n_diag, l_have_output  )
 
    IF (ltimer) CALL print_timer

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

    CALL message(TRIM(routine),'start to clean up')

    ! Destruction of post-processing job queue
    IF (.NOT. ltestcase) &
      CALL pp_scheduler_finalize()

    ! Delete optional diagnostics
    CALL destruct_opt_diag()
   
    ! Delete state variables

    CALL destruct_nh_state( p_nh_state )
    DEALLOCATE (p_nh_state, STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
    ENDIF

    IF (iforcing == inwp) THEN
      CALL destruct_nwp_phy_state
      CALL destruct_tiles_arrays(p_tiles)
    ENDIF
    CALL destruct_nwp_lnd_state( p_lnd_state )

    ! Delete output variable lists
    IF (l_have_output) CALL close_output_files
    CALL close_name_list_output

    ! finalize meteogram output
    IF (.NOT. ltestcase) THEN
      DO jg = 1, n_dom
        IF (meteogram_output_config(jg)%lenabled) THEN
          CALL meteogram_finalize(jg)
        END IF
      END DO
      DO jg = 1, max_dom
        DEALLOCATE(meteogram_output_config(jg)%station_list)
      END DO
    END IF

    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE atmo_nonhydrostatic
  
END MODULE mo_atmo_nonhydrostatic

