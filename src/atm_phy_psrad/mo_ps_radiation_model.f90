!>
!! @page pagecontrolmodelf901 Main program for the ICON ps_radiation_model 
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! @author L. Linardakis, MPI-M, Hamburg
!!
!! @remarks
!!  

MODULE mo_ps_radiation_model

  USE mo_exception,               ONLY: message, message_text, finish
  USE mo_psrad_interface_namelist,ONLY: configure_ps_radiation, number_of_levels
  USE mo_time_config,             ONLY: time_config
!  USE mo_grid_config,             ONLY: n_dom

  USE mo_model_domain,            ONLY: t_patch, p_patch
  USE mo_run_config,              ONLY: ltimer, nshift
  USE mo_build_decomposition,     ONLY: build_decomposition
!  USE mo_icon_comm_interface,     ONLY: construct_icon_communication
  USE mo_icon_comm_interface,     ONLY: destruct_icon_communication

  
  USE mo_psrad_interface_memory,  ONLY: construct_psrad_interface_memory, destruct_psrad_interface_memory
 
  USE mo_echam_phy_config,        ONLY: echam_phy_tc
  USE mo_echam_rad_config       , ONLY: echam_rad_config
  USE mo_bc_aeropt_kinne         ,ONLY: read_bc_aeropt_kinne
  USE mo_bc_aeropt_stenchikov    ,ONLY: read_bc_aeropt_stenchikov
  USE mo_bc_aeropt_splumes,       ONLY: setup_bc_aeropt_splumes

  USE mtime,                      ONLY: datetime, timedelta, datetimeToString,      &
    &   OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*), OPERATOR(<),            &
    &   ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=)
  USE mo_mpi,                     ONLY: set_mpi_work_communicators  
  USE mo_parallel_config,         ONLY: p_test_run, l_test_openmp, num_io_procs,              &
    &                                   num_restart_procs
  USE mo_timer,                   ONLY: init_timer, timer_start, timer_stop,                  &
    &                                   print_timer
!  USE mo_timer,                   ONLY: timer_model_init

  USE mo_atmo_psrad_interface,    ONLY: psrad_concurrent_interface, finalize_psrad_concurrent

  USE mo_psrad_communication,     ONLY: setup_psrad_2_atmo_communication
  USE mo_master_control,          ONLY: ps_radiation_process
  USE mo_timer,                   ONLY: ltimer, timer_start, timer_stop, timer_total

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ps_radiation_model

  CONTAINS

    !--------------------------------------------------------------------------
    !>
    SUBROUTINE ps_radiation_model(ps_rad_namelist_filename,shr_namelist_filename)

      CHARACTER(LEN=*), INTENT(in) :: ps_rad_namelist_filename,shr_namelist_filename

      !-------------------------------------------------------------------
      CALL construct_ps_radiation_model(ps_rad_namelist_filename,shr_namelist_filename)
 


      !------------------------------------------------------------------
      ! Initialize output file if necessary;
      ! Write out initial conditions.
      !------------------------------------------------------------------
!     IF (output_mode%l_nml) THEN
!       CALL parse_variable_groups()
!       ! compute sim_start, sim_end
!       CALL datetimeToString(time_config%tc_exp_startdate, sim_step_info%sim_start)
!       CALL datetimeToString(time_config%tc_stopdate, sim_step_info%sim_end)
!       CALL datetimeToString(time_config%tc_startdate, sim_step_info%run_start)
!       CALL datetimeToString(time_config%tc_stopdate, sim_step_info%restart_time)
! 
!       sim_step_info%dtime      = dtime
!       jstep0 = 0
! 
!       restartAttributes => getAttributesForRestarting()
!       IF (ASSOCIATED(restartAttributes)) THEN
! 
!         ! get start counter for time loop from restart file:
!         jstep0 = restartAttributes%getInteger("jstep")
!       END IF
!       sim_step_info%jstep0    = jstep0
!       CALL init_statistics_stream
!       CALL init_name_list_output(sim_step_info, opt_lprintlist=.TRUE.,opt_l_is_ps_radiation_model=.TRUE.)
!       CALL create_mipz_level_selections(output_file)
!     ENDIF

    !------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_total)

    CALL run_ps_radiation_model()

    IF (ltimer) CALL timer_stop(timer_total)
    !------------------------------------------------------------------
    !  cleaning up process
    CALL destruct_ps_radiation_model()


    CALL print_timer()

  END SUBROUTINE ps_radiation_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE run_ps_radiation_model()

    TYPE(datetime),  POINTER        :: mtime_current     ! current datetime (mtime)
    TYPE(t_patch),   POINTER        :: patch    ! Patch
    TYPE(timedelta), POINTER        :: radiation_time_step => NULL()
    INTEGER                         :: timestep
    CHARACTER(LEN=32)               :: datestring, datestring2
    CHARACTER(LEN=*), PARAMETER     :: method_name="ps_radiation_model"

    CALL message ("", "-----------------------------------------------------------")
    CALL message (method_name, " starts...")

    mtime_current => time_config%tc_current_date
!     radiation_time_step => time_config%tc_dt_model  ! NOTE: here the model timestep in the namelists is radiation timestep
    radiation_time_step => echam_phy_tc(1)%dt_rad

    CALL datetimeToString(mtime_current,           datestring)
    CALL datetimeToString(time_config%tc_stopdate, datestring2)
    WRITE(message_text,'(4a)') ' start datetime: ', datestring, '; end datetime: ', datestring2
    CALL message (method_name, message_text)

    patch => p_patch(1)
    timestep = 0
    ! loop 
    DO WHILE (mtime_current < time_config%tc_stopdate) 

      CALL datetimeToString(mtime_current, datestring)
      WRITE(message_text,'(a,i10,2a)') '  Begin of timestep = ',timestep,'  datetime: ', datestring
      CALL message (method_name, message_text)
      
      CALL ps_rad_run_bc(mtime_current, patch)

      CALL psrad_concurrent_interface(mtime_current)
 
      mtime_current = mtime_current + radiation_time_step
      timestep = timestep + 1
   
    ENDDO

    CALL finalize_psrad_concurrent

    CALL message (method_name, " ended")
    CALL message ("", "-----------------------------------------------------------")

  END SUBROUTINE run_ps_radiation_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE ps_rad_run_bc(mtime_current, patch)
    TYPE(datetime),  POINTER, INTENT(in)  :: mtime_current     ! current datetime (mtime)
    TYPE(t_patch),   POINTER, INTENT(in)  :: patch    ! Patch

    CHARACTER(LEN=*), PARAMETER     :: method_name="ps_rad_run_bc"

    IF (.not. associated(mtime_current)) CALL finish(method_name,".not. associated(mtime_current)")
!     write(0,*) method_name, "..."
!   CALL to get pp_sfc, tk_hl
!   SUBROUTINE calculate_temperatur_pressure (                &
!       & jg            ,&!< in  domain index
!       & jb            ,&!< in  block index
!       & kproma        ,&!< in  end index for loop over block
!       & kbdim         ,&!< in  dimension of block over cells
!       & klev          ,&!< in  number of full levels = number of layers
!       & klevp1        ,&!< in  number of half levels = number of layer interfaces
!       & pp_hl         ,&! in
!       & pp_fl         ,&
!       & tk_fl         ,&
!       & tk_sfc        ,&
!       & pp_sfc        ,& ! out
!       & tk_hl   )

      ! tropospheric aerosol optical properties
      IF (echam_rad_config(1)%irad_aero == 13) THEN
        CALL read_bc_aeropt_kinne(mtime_current, patch) 
      END IF
      !
      ! stratospheric aerosol optical properties
      IF (echam_rad_config(1)%irad_aero == 14) THEN
        CALL read_bc_aeropt_stenchikov(mtime_current, patch)
      END IF
      !
      ! tropospheric and stratospheric aerosol optical properties
      IF (echam_rad_config(1)%irad_aero == 15) THEN
        CALL read_bc_aeropt_kinne     (mtime_current, patch)
        CALL read_bc_aeropt_stenchikov(mtime_current, patch)
      END IF
      ! tropospheric background aerosols (Kinne) and stratospheric
      ! aerosols (Stenchikov) + simple plumes (analytical, nothing to be read
      ! here, initialization see init_echam_phy (mo_echam_phy_init)) 
      IF (echam_rad_config(1)%irad_aero == 18) THEN
        CALL read_bc_aeropt_kinne     (mtime_current, patch)
        CALL read_bc_aeropt_stenchikov(mtime_current, patch)
      END IF

!     write(0,*) method_name, " done."

  END SUBROUTINE ps_rad_run_bc
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  !>
  SUBROUTINE destruct_ps_radiation_model()

    CHARACTER(*), PARAMETER :: method_name = "mo_ps_radiation_model:destruct_ps_radiation_model"

    !------------------------------------------------------------------
    !  cleaning up process
    !------------------------------------------------------------------
    CALL message(TRIM(method_name),'start to clean up')

    CALL destruct_psrad_interface_memory()
  
    CALL destruct_icon_communication()

    CALL message(TRIM(method_name),'clean-up finished')

  END SUBROUTINE destruct_ps_radiation_model
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  SUBROUTINE construct_ps_radiation_model(ps_rad_namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: ps_rad_namelist_filename,shr_namelist_filename

    CHARACTER(*), PARAMETER :: method_name = "mo_ps_radiation_model:construct_ps_radiation_model"
!    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datetime_string

    !-------------------------------------------------------------------
    CALL message(TRIM(method_name),'start construction')

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    CALL configure_ps_radiation(ps_rad_namelist_filename,shr_namelist_filename)

    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------
!     CALL ps_radiation_model_crosscheck()

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs, &
                                    num_restart_procs, &
                                    opt_comp_id=ps_radiation_process)

    !-------------------------------------------------------------------
    ! 3.2 Initialize various timers
    !-------------------------------------------------------------------
    CALL init_timer

!     IF (ltimer) CALL timer_start(timer_model_init)

    !-------------------------------------------------------------------
    ! 4. Setup IO procs
    !-------------------------------------------------------------------
    ! If we belong to the I/O PEs just call xxx_io_main_proc before
    ! reading patches.  This routine will never return
!     CALL init_io_processes()

    ! 4. Import patches
    !-------------------------------------------------------------------
!     CALL build_decomposition(number_of_levels,nshift, is_ps_radiation_model_decomposition =.TRUE., &
!       & patch_3d=ps_radiation_model_patch_3d)
    CALL build_decomposition(number_of_levels, nshift, is_ocean_decomposition = .false.)
!     CALL construct_icon_communication(p_patch, n_dom)

    !-------------------------------------------------------------------
    ! 5. set up the communication between atmo and psrad processes
    !-------------------------------------------------------------------
    CALL setup_psrad_2_atmo_communication()

    !--------------------------------------------
    ! Setup the information for the physical patches
!     CALL setup_phys_patches

    !------------------------------------------------------------------
    ! step 5b: allocate state variables
    !------------------------------------------------------------------
    CALL construct_psrad_interface_memory( p_patch )    

    !--------------------------------------------
    IF (echam_rad_config(1)%irad_aero == 18) THEN
      CALL setup_bc_aeropt_splumes
    END IF
    !--------------------------------------------
 
    CALL message(TRIM(method_name),'end construction')

  END SUBROUTINE construct_ps_radiation_model
  !--------------------------------------------------------------------------


END MODULE mo_ps_radiation_model

