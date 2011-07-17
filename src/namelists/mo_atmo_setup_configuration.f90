!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
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
MODULE mo_atmo_setup_configuration

  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_master_nml,          ONLY: lrestart
  USE mo_namelist,            ONLY: open_nml_output, close_nml_output

  USE mo_time_nml,            ONLY: read_time_namelist
  USE mo_parallel_nml,        ONLY: read_parallel_namelist

  USE mo_run_nml,             ONLY: read_run_namelist
  USE mo_ha_testcases,        ONLY: read_ha_testcase_namelist 
  USE mo_nh_testcases,        ONLY: read_nh_testcase_namelist

  USE mo_dynamics_nml,        ONLY: read_dynamics_namelist
  USE mo_nonhydrostatic_nml,  ONLY: read_nonhydrostatic_namelist
  USE mo_ha_dyn_nml,          ONLY: read_ha_dyn_namelist 
  USE mo_diffusion_nml,       ONLY: read_diffusion_namelist 
  USE mo_io_nml,              ONLY: read_io_namelist
  USE mo_extpar_nml,          ONLY: read_extpar_namelist
  USE mo_advection_nml,       ONLY: read_transport_namelist
  USE mo_gridref_nml,         ONLY: read_gridref_namelist
  USE mo_echam_phy_nml,       ONLY: read_echam_phy_namelist
  USE mo_vdiff_nml,           ONLY: read_vdiff_namelist
  USE mo_echam_conv_nml,      ONLY: read_echam_conv_namelist
  USE mo_atm_phy_nwp_nml,     ONLY: read_nwp_phy_namelist
  USE mo_radiation_nml,       ONLY: read_radiation_namelist
  USE mo_gw_hines_nml,        ONLY: read_gw_hines_namelist
  USE mo_lnd_nwp_nml,         ONLY: read_nwp_lnd_namelist
  USE mo_sleve_nml,           ONLY: read_sleve_namelist
  USE mo_grid_nml,            ONLY: read_grid_namelist
  USE mo_interpol_nml,        ONLY: read_interpol_namelist

  USE mo_io_restart_namelist,  ONLY: read_restart_namelists
  USE mo_io_restart_attributes,ONLY: read_restart_attributes, get_restart_attribute
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_atmo_namelists !, setup_atmo_configuration
  
CONTAINS
  !>
  !! Read namelists;
  !! Create a new file in which all the namelist variables and their
  !! actual values used in the model run will be stored.
  !!
  SUBROUTINE read_atmo_namelists(atm_namelist_filename,shr_namelist_filename)
    
    CHARACTER(LEN=*), INTENT(in) :: atm_namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    IF(my_process_is_stdio()) CALL open_nml_output('NAMELIST_ICON_output_atm')

    ! Shared with the ocean model when performing a coupled simulation

    CALL read_time_namelist       (TRIM(shr_namelist_filename))
    CALL read_parallel_namelist   (TRIM(shr_namelist_filename))

    ! General

    CALL read_run_namelist        (TRIM(atm_namelist_filename))
    CALL read_ha_testcase_namelist(TRIM(atm_namelist_filename))
    CALL read_nh_testcase_namelist(TRIM(atm_namelist_filename))
    CALL read_io_namelist         (TRIM(atm_namelist_filename))

    ! Grid, dynamics, and transport

    CALL read_grid_namelist       (TRIM(atm_namelist_filename))
    CALL read_gridref_namelist    (TRIM(atm_namelist_filename))
    CALL read_interpol_namelist   (TRIM(atm_namelist_filename))

    CALL read_dynamics_namelist   (TRIM(atm_namelist_filename))
    CALL read_ha_dyn_namelist     (TRIM(atm_namelist_filename))
    CALL read_nonhydrostatic_namelist(TRIM(atm_namelist_filename))
    CALL read_sleve_namelist      (TRIM(atm_namelist_filename))

    CALL read_diffusion_namelist  (TRIM(atm_namelist_filename))

    CALL read_transport_namelist  (TRIM(atm_namelist_filename))

    ! Physics

    CALL read_extpar_namelist     (TRIM(atm_namelist_filename))
    CALL read_radiation_namelist  (TRIM(atm_namelist_filename))
    CALL read_vdiff_namelist      (TRIM(atm_namelist_filename))
    CALL read_echam_conv_namelist (TRIM(atm_namelist_filename))
    CALL read_gw_hines_namelist   (TRIM(atm_namelist_filename))
    CALL read_nwp_lnd_namelist    (TRIM(atm_namelist_filename))
    CALL read_nwp_phy_namelist    (TRIM(atm_namelist_filename))
    CALL read_echam_phy_namelist  (TRIM(atm_namelist_filename))
      
    IF (my_process_is_stdio()) THEN
      CALL close_nml_output
    END IF
        
  END SUBROUTINE read_atmo_namelists
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE old_read_atmo_namelists !(namelist_filename)
!     
!     CHARACTER(LEN=*), INTENT(in) :: namelist_filename
! 
!     CHARACTER(*), PARAMETER :: method_name = "old_read_atmo_namelists"
! 
!     CHARACTER(LEN=MAX_CHAR_LENGTH) :: grid_file_name 
!     INTEGER :: n_io, jg, jfile, n_file, ist, n_diag, n_chkpt
!     LOGICAL :: lsuccess, l_have_output
!    
! 
!     !-------------------------------------------------------------------
!     ! 1. Open the atmosphere-specific namelist file and create a 
!     !    new file in which all the namelist variables and their
!     !    actual values used in the model run will be stored.
!     !-------------------------------------------------------------------
!     CALL open_nml(TRIM(namelist_filename))
!     IF(p_pe == p_io) CALL open_nml_output('NAMELIST_ICON_output_atm')
! 
!     ! The namelists ('run_nml' and 'testcase_ctl') are read in seperate
!     ! subroutines. The validity of the user-specified configurations is
!     ! checked right after each namelist is read.
!     ! The two 'setup' subroutines above must be called before grid and patch
!     ! import because during the import procedure the topography and land-sea
!     ! mask will be initialized. They are related to the atmos/ocean switch
!     ! as well as the selected test case.
!DONE     
!DONE     CALL run_nml_setup
!DONE     
!DONE     !-------------------------------------------------------------------
!DONE     ! parallel_nml_setup must be called after setup_run since it needs
!DONE     ! some variables read in setup_run
!DONE     
!DONE     CALL parallel_nml_setup
!     
!     !-------------------------------------------------------------------
!     ! Initialize test case setup 
!     !-------------------------------------------------------------------
!     IF (ltestcase) THEN
!       
!       SELECT CASE (iequations)
!       CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
!         CALL setup_testcase
!         
!       CASE (inh_atmosphere)
!         CALL setup_nh_testcase
!         
!       CASE DEFAULT
!         CALL finish( TRIM(method_name),' invalid value for iequations!' )
!       END SELECT
!       
!     ENDIF
!     
!DONE     !-------------------------------------------------------------------
!DONE     ! step 2: import computational domain of the model
!DONE     !-------------------------------------------------------------------
!DONE     ! 2a) read namelist 'grid_ctl' from the namelist file (already
!DONE     !     opened above) and set up the grid/patch configuration.
!DONE     !-------------------------------------------------------------------
!DONE     
!DONE     CALL grid_nml_setup
!DONE         
!DONE     !------------------------------------------------------------------
!DONE     ! step 3a: Read namelist 'dynamics_ctl'. This has to be done
!DONE     !          before the interpolation state is computed.
!DONE     !------------------------------------------------------------------
!DONE     CALL dynamics_nml_setup(n_dom)
!DONE     
!DONE     !------------------------------------------------------------------
!DONE     ! Read specific dynamics namelists for the nonhydrost. dynamical core
!DONE     !------------------------------------------------------------------
!DONE     CALL nonhydrostatic_nml_setup
!     
!     !------------------------------------------------------------------
!     ! step 3a-a: ! read nwp physics namelist, ...
!     !------------------------------------------------------------------
!       IF ( iforcing == inwp) THEN
!         CALL setup_nwp_phy( p_patch_global(1:) )  ! read Namelist, ...
!         IF (inwp_surface > 0) CALL setup_nwp_lnd
!       ENDIF
!     
!     !------------------------------------------------------------------
!     ! step 3b: Read namelist 'transport_ctl',
!     !------------------------------------------------------------------
!     CALL transport_nml_setup
!     
!     !------------------------------------------------------------------
!     ! step 4: construct interpolation and the grid refinement state
!     !------------------------------------------------------------------
!     ! - Allocate memory for the interpolation state;
!     ! - Calculate interpolation coefficients.
!     !------------------------------------------------------------------
!     
!     CALL gridref_nml_setup
!     
!     ! interpolation state not used for ocean model
!     ! #slo# - temporarily switched on for comparison with rbf-reconstruction
!     
! !     CALL interpol_nml_setup(p_patch_global)
!     
!     !-------------------------------------------------------------------------
!     ! set up horizontal diffusion after the vertical coordinate is configured
!     !-------------------------------------------------------------------------
! !     CALL diffusion_nml_setup(n_dom,parent_id)
!     
!    
!     !------------------------------------------------------------------
!     ! step 6: initialize output
!     !------------------------------------------------------------------    
!     CALL io_nml_setup
!      
!     !------------------------------------------------------------------
!     ! Create and optionally read external data fields
!     !------------------------------------------------------------------
!         
!     CALL close_nml
!     IF (p_pe == p_io) THEN
!       CALL close_nml_output
!     END IF
        
  END SUBROUTINE old_read_atmo_namelists
  !-------------------------------------------------------------------------
  
  
  
END MODULE mo_atmo_setup_configuration

