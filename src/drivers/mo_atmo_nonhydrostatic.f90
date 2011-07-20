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

USE mo_exception,            ONLY: message, finish
USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH
USE mo_timer,                ONLY: print_timer
USE mo_master_nml,           ONLY: lrestart
USE mo_output,               ONLY: init_output_files, close_output_files,&
  &                                write_output
USE mo_interpolation,        ONLY: rbf_vec_interpol_cell,       &
  &                                edges2cells_scalar
USE mo_time_config,          ONLY: time_config      ! variable
USE mo_io_restart,           ONLY: read_restart_files
USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute
USE mo_io_config,            ONLY: dt_data,dt_file,dt_diag,dt_checkpoint
USE mo_run_config,           ONLY: &
  &                               dtime,                & !    namelist parameter
  &                               ltestcase,            &
  &                               nsteps,               & !    :
  &                               ltimer,               & !    :
  &                               iforcing                !    namelist parameter
USE mo_impl_constants,       ONLY: inoforcing,           & !    :
  &                                inwp
! Horizontal grid
USE mo_atmo_control,         ONLY: p_patch_subdiv, p_patch
USE mo_grid_config,          ONLY: n_dom
! to break circular dependency KF???
USE mo_intp_data_strc,       ONLY: p_int_state_global, p_int_state_subdiv, p_int_state
USE mo_grf_intp_data_strc,   ONLY: p_grf_state_global, p_grf_state_subdiv, p_grf_state
! NH-namelist state
USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config,&  !> namelist state
  &                                configure_atm_phy_nwp !> subroutine
! NH-Model states
USE mo_nonhydro_state,       ONLY: p_nh_state,   &
  &                                destruct_nh_state
USE mo_nwp_phy_state,        ONLY: construct_nwp_phy_state, &
  &                                destruct_nwp_phy_state
USE mo_lnd_nwp_nml,          ONLY: setup_nwp_lnd
USE mo_nwp_lnd_state,        ONLY: construct_nwp_lnd_state,   &
  &                                destruct_nwp_lnd_state, p_lnd_state
! Time integration
USE mo_nh_stepping,          ONLY: prepare_nh_integration, perform_nh_stepping


!-------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE

PUBLIC :: atmo_nonhydrostatic

CONTAINS

  !---------------------------------------------------------------------
  SUBROUTINE atmo_nonhydrostatic
    
    CHARACTER(*), PARAMETER :: routine = "mo_atmo_nonhydrostatic"


    INTEGER :: n_io, jg, jfile, n_file, ist, n_diag, n_chkpt
    LOGICAL :: l_have_output
   

    IF(iforcing == inwp) THEN
     CALL configure_atm_phy_nwp(n_dom, ltestcase)
    ENDIF
 
    !---------------------------------------------------------------------
    ! 4.c Non-Hydrostatic / NWP
    !---------------------------------------------------------------------

     ALLOCATE (p_nh_state(n_dom), stat=ist)
     IF (ist /= success) THEN
       CALL finish(TRIM(routine),'allocation for p_nh_state failed')
     ENDIF

     IF(iforcing == inwp) THEN
       CALL construct_nwp_phy_state( p_patch(1:) )
       
       IF (atm_phy_nwp_config(1)%inwp_surface > 0 ) THEN
         ALLOCATE (p_lnd_state(n_dom), stat=ist)
         IF (ist /= success) THEN
           CALL finish(TRIM(routine),'allocation for p_lnd_state failed')
         ENDIF
         CALL construct_nwp_lnd_state( p_patch(1:),p_lnd_state,n_timelevels=2 )
       ENDIF

     ENDIF

    !---------------------------------------------------------------------
    ! 5. Perform time stepping
    !---------------------------------------------------------------------
      !------------------------------------------------------------------
      ! Prepare for time integration
      !------------------------------------------------------------------

      !------------------------------------------------------------------
      ! Set initial conditions for time integration.
      !------------------------------------------------------------------
     
    IF (lrestart) THEN
      ! This is an resumed integration. Read model state from restart file(s).

      CALL read_restart_files
      CALL message(TRIM(routine),'normal exit from read_restart_files')

    ENDIF
                                                                                       
   !--------------------

     CALL prepare_nh_integration(p_patch(1:), p_nh_state, p_int_state(1:), p_grf_state(1:))

     !---------------------------------------------------------
     ! The most primitive event handling algorithm: 
     ! compute time step interval for taking a certain action
     !--------------------------------------------------------- 
 
     n_io    = NINT(dt_data/dtime)        ! write output
     n_file  = NINT(dt_file/dtime)        ! trigger new output file
     n_chkpt = NINT(dt_checkpoint/dtime)  ! write restart files
     n_diag  = MAX(1,NINT(dt_diag/dtime)) ! diagnose of total integrals

     !------------------------------------------------------------------
     ! Prepare output file
     !------------------------------------------------------------------
     IF (.NOT.lrestart) THEN
       ! Initialize the first output file which will contain also the 
       ! initial conditions.

       jfile = 1
       CALL init_output_files(jfile, lclose=.FALSE.)

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

    END IF

    !------------------------------------------------------------------
    !  get and write out some of the inital values
    !------------------------------------------------------------------
    IF (.NOT.lrestart) THEN

    ! diagnose u and v to have meaningful initial output
    
    DO jg = 1, n_dom

        SELECT CASE (p_patch(jg)%cell_type)
        CASE (3)
          CALL rbf_vec_interpol_cell(p_nh_state(jg)%prog(1)%vn,p_patch(jg),&
            & p_int_state(jg),p_nh_state(jg)%diag%u,p_nh_state(jg)%diag%v)
        CASE (6)
          CALL edges2cells_scalar(p_nh_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_east,p_nh_state(jg)%diag%u)
          CALL edges2cells_scalar(p_nh_state(jg)%prog(1)%vn,p_patch(jg), &
            & p_int_state(jg)%hex_north,p_nh_state(jg)%diag%v)
        END SELECT

    ENDDO
    
    ! Note: here the derived output variables are not yet available
    ! (divergence, vorticity)
    CALL write_output( time_config%cur_datetime )
    l_have_output = .TRUE.

    END IF ! not lrestart

    !------------------------------------------------------------------
    ! Now start the time stepping:
    ! The special initial time step for the three time level schemes
    ! is executed within process_grid_level
    !------------------------------------------------------------------

      CALL perform_nh_stepping( p_patch, p_int_state, p_grf_state, p_nh_state,   &
                              & time_config%cur_datetime,                        &
                              & n_io, n_file, n_chkpt, n_diag, l_have_output     )
 
    IF (ltimer) CALL print_timer

    !---------------------------------------------------------------------
    ! 6. Integration finished. Clean up.
    !---------------------------------------------------------------------

   CALL message(TRIM(routine),'start to clean up')
   
   ! Delete state variables

   CALL destruct_nh_state( p_nh_state )
   DEALLOCATE (p_nh_state, STAT=ist)
   IF (ist /= SUCCESS) THEN
     CALL finish(TRIM(routine),'deallocation for p_nh_state failed')
   ENDIF

   IF (iforcing == inwp) THEN
     CALL destruct_nwp_phy_state
     CALL destruct_nwp_lnd_state(p_lnd_state)
   ENDIF

   ! Delete output variable lists
   IF (l_have_output) CALL close_output_files

    CALL message(TRIM(routine),'clean-up finished')
    
  END SUBROUTINE atmo_nonhydrostatic
  
END MODULE mo_atmo_nonhydrostatic

