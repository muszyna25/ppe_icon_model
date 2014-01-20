!>
!!        Contains the variables to set up the ocean model.
!!
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3814)
!!   Modification by Constantin Junk (2010-03-18)
!!     - separated namelist mpiom_phy_nml, ocean_nml und octst_nml
!!       from mo_global_variables
!!     - therefore, added mo_ocean_nml module
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_ocean_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_ocean_config,       ONLY: config_ignore_land_points => ignore_land_points
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
#ifndef __NO_ICON_ATMO__
  USE mo_coupling_config,    ONLY: is_coupled_run
#endif
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC


  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters for mpiom_phy_nml
  !     mpiom forcing (right hand side)
  ! ------------------------------------------------------------------------

  ! switches for parameterizations (examples)
  LOGICAL  :: lmpiom_radiation
  LOGICAL  :: lmpiom_convection
  LOGICAL  :: lmpiom_gentmcwill

  NAMELIST/mpiom_phy_nml/ lmpiom_radiation, lmpiom_convection, lmpiom_gentmcwill


  ! ------------------------------------------------------------------------
  ! 2.0 Namelist variables and auxiliary parameters for ocean_nml
  !      - contains all default values to minimize ocean namelist (SLO, 2012/03)
  ! ------------------------------------------------------------------------

  INTEGER  :: n_zlev        ! number of ocean levels
  REAL(wp) :: dzlev_m(100)  ! namelist input of layer thickness

  INTEGER, PARAMETER :: toplev    = 1   ! surface ocean level

  ! parameterized forcing for ocean model:
  INTEGER            :: iforc_oce                 =  0  ! index of parameterized forcing
  INTEGER, PARAMETER :: NO_FORCING                = 10
  INTEGER, PARAMETER :: ANALYT_FORC               = 11
  INTEGER, PARAMETER :: FORCING_FROM_FILE_FLUX    = 12  ! OMIP or NCEP type forcing
  INTEGER, PARAMETER :: FORCING_FROM_FILE_FIELD   = 13  ! not yet
  INTEGER, PARAMETER :: FORCING_FROM_COUPLED_FLUX = 14  ! parameter for a coupled atmosphere-ocean run
  INTEGER, PARAMETER :: FORCING_FROM_COUPLED_FIELD= 15  ! not yet

  ! read time varying OMIP or NCEP flux forcing from file:
                      ! 1: read wind stress (records 1, 2) and temperature (record 3)
                      ! 2: read full OMIP dataset for bulk formula in mo_oce_bulk (12 records)
                      ! 3: as 1; read surface heat (record 4) and freshwater flux (record 5) add.
                      ! 4: as 1; read 4 parts of heat flux, precip/evap flux additionally
                      ! 5: read full NCEP datasets; read monthly mean data of consecutive years
  INTEGER            :: iforc_type = 10

  ! length of time varying flux forcing: 12: read 12 months, other: read daily values
  INTEGER            :: forcing_timescale = 1

  ! switch for stationary forcing for special testcases of ocean model:
  INTEGER            :: iforc_stat_oce = 3
!TODO

  LOGICAL  :: forcing_enable_freshwater            = .FALSE.    ! .TRUE.: apply freshwater forcing boundary condition
  LOGICAL  :: forcing_set_runoff_to_zero           = .FALSE.    ! .TRUE.: set river runoff to zero for comparion to MPIOM
  LOGICAL  :: use_new_forcing                      = .FALSE.
  INTEGER  :: forcing_fluxes_type                  = 0
  INTEGER  :: forcing_windstress_u_type            = 0
  INTEGER  :: forcing_windstress_v_type            = 0
  REAL(wp) :: forcing_windstress_zonal_waveno      = 3.0_wp  ! For the periodic analytic forcing (wind)
  REAL(wp) :: forcing_windstress_meridional_waveno = 3.0_wp
  REAL(wp) :: forcing_wind_u_amplitude             = 1.0_wp
  REAL(wp) :: forcing_wind_v_amplitude             = 1.0_wp
  REAL(wp) :: analytic_wind_amplitude              = 1.0_wp


  INTEGER            :: relax_temp_type           = 0 ! will cover parts of init_oce_relax
  INTEGER            :: relax_temp_param          = 0 ! replacement for relaxation_param
  INTEGER            :: relax_sal_type            = 0 ! will cover parts of init_oce_relax
  INTEGER            :: relax_sal_param           = 0 ! replacement for relax_2d_mon_S

  ! switch for reading relaxation data: 1: read from file
  INTEGER            :: init_oce_relax = 0


  ! switch for ocean diagnostics - 0: no diagnostics; 1: write to stderr
  INTEGER            :: diagnostics_level      = 0

  ! switch for ocean stream function (not yet activated):
  !                   ! 0: no output
                      ! 1: write barotropic velocity in output file
                      ! 2: write barotropic stream function on regular grid
  INTEGER            :: idiag_psi      = 0

  ! parameterized velocity boundary conditions
                      ! Velocity boundary condition: Currently only no-slip is supported !!
                      ! i_bc_veloc_lateral = 0: boundary condition for velocity is no-slip: normal
                      !                         and tangential velocity components at lateral 
                      !                         boundaries are set to zero
                      ! i_bc_veloc_lateral = 1: boundary condition for velocity is free-slip: 
                      !                         normal velocity components at lateral boundariea is
                      !                         set to zero, tangential not.
  INTEGER            :: i_bc_veloc_lateral = 0   

  INTEGER            :: i_bc_veloc_top = 1  !Top boundary condition for velocity: 
                                            ! i_bc_veloc_top =0 : zero value at top boundary,no wind
                                            ! i_bc_veloc_top =1 : forced by wind field that is
                                            !                     stored in p_os%p_aux%bc_top_veloc
                                            ! i_bc_veloc_top =2 : forced by difference between wind
                                            !                     field in p_os%p_aux%bc_top_veloc 
                                            !                     and ocean velocity at top layer
  INTEGER            :: i_bc_veloc_bot = 0  !Bottom boundary condition for velocity: 
                                            ! i_bc_veloc_bot =0 : zero value at bottom boundary 
                                            ! i_bc_veloc_bot =1 : bottom boundary friction
                                            ! i_bc_veloc_bot =2 : bottom friction plus topographic
                                            !                     slope (not implemented yet)
  ! parameterized choice of tracer transport scheme
  INTEGER, PARAMETER :: UPWIND = 1
  INTEGER, PARAMETER :: CENTRAL= 2
  INTEGER, PARAMETER :: MIMETIC= 3
  INTEGER, PARAMETER :: MIMETIC_MIURA= 4
  !INTEGER, PARAMETER :: ADPO   = 5
  INTEGER, PARAMETER :: FLUX_CORR_TRANSP_horz = 5
  INTEGER, PARAMETER :: FLUX_CORR_TRANSP_vert_adpo = 6
  !INTEGER            :: FLUX_CALCULATION_HORZ = MIMETIC  ! consistent with l_edge_based=.FALSE.
  !INTEGER            :: FLUX_CALCULATION_VERT = UPWIND   ! consistent with l_edge_based=.FALSE.
  INTEGER            :: FLUX_CALCULATION_HORZ = MIMETIC_MIURA  ! consistent with l_edge_based=.TRUE.
  INTEGER            :: FLUX_CALCULATION_VERT = MIMETIC_MIURA  ! consistent with l_edge_based=.TRUE.
  LOGICAL            :: l_adpo_flowstrength   = .TRUE.     ! .TRUE.: activate second condition for adpo weight


  !this distinction is no longer used: INTEGER  :: i_sfc_forcing_form        = 0
  !=0: surface forcing applied as top boundary condition to vertical diffusion
  !=1: surface forcing applied as volume forcing at rhs, i.e.part of explicit term in momentum and tracer eqs. 
  !in this case, top boundary ondition of vertical Laplacians are homogeneous
  !INTEGER            :: vbc_zero_cond   =   0   ! no or zero boundary condition

  ! parameterized shallow water mode in the ocean model
  INTEGER            :: iswm_oce        =   0  ! switch for shallow water mode (1 = on, 0 = 3dim)
  INTEGER            :: discretization_scheme    =   1  ! discretization scheme: 1 for mimetic, 
                                               ! 2 for RBF-type of discretization

  ! parameters for Adams-Bashforth semi-implicit time stepping scheme
  ! are set according to Marshall et al paper
  REAL(wp) :: ab_const              = 0.1_wp     ! Adams-Bashforth constant
  REAL(wp) :: ab_beta               = 0.6_wp     ! Parameter in semi-implicit timestepping
  REAL(wp) :: ab_gam                = 0.6_wp     ! Parameter in semi-implicit timestepping


  ! parameters for gmres solver
  REAL(wp) :: solver_tolerance                 = 1.e-11_wp   ! Maximum value allowed for solver absolute tolerance
!  REAL(wp) :: solver_start_tolerance           = -1.0_wp
  INTEGER  :: solver_max_restart_iterations    = 100       ! For restarting gmres
  INTEGER  :: solver_max_iter_per_restart      = 200       ! For inner loop after restart
!  REAL(wp) :: solver_tolerance_decrease_ratio  = 0.1_wp    ! For restarting gmres, must be < 1
  LOGICAL  :: use_absolute_solver_tolerance    = .true.   ! Maximum value allowed for solver tolerance
  INTEGER, PARAMETER :: select_gmres = 1
  INTEGER, PARAMETER :: select_restart_gmres = 2
  INTEGER :: select_solver = select_restart_gmres
  LOGICAL :: use_continuity_correction = .false.
  LOGICAL :: use_edges2edges_viacell_fast = .false.


  ! physical parameters for  aborting the ocean model
  REAL(wp) :: dhdtw_abort           =  3.17e-11_wp  ! abort criterion for gmres solution (~1mm/year)
  REAL(wp) :: threshold_vn          = 10.0_wp    ! abort criterion for absolute velocity maximum
  REAL(wp) :: threshold_min_T       = -4.0_wp    ! abort criterion for temperature minimum
  REAL(wp) :: threshold_max_T       = 100._wp    ! abort criterion for temperature minimum
  REAL(wp) :: threshold_min_S       =  0.0_wp    ! abort criterion for salinity minimum
  REAL(wp) :: threshold_max_S       = 60.0_wp    ! abort criterion for salinity minimum

  INTEGER  :: EOS_TYPE              = 2          ! 1=linear EOS,2=(nonlinear, from MPIOM)
                                                 ! 3=nonlinear Jacket-McDoudgall-formulation (not yet recommended)
  INTEGER  :: density_computation   = 1          ! 1 = calc_density_MPIOM_func,2 = calc_density_MPIOM_elemental,
                                                 ! 3 = calc_density_MPIOM_elemental_wrap
  INTEGER  :: no_tracer             = 2          ! number of tracers 

  ! more ocean parameters, not yet well placed
  INTEGER  :: expl_vertical_velocity_diff = 1    ! 0=explicit, 1 = implicit  
  INTEGER  :: expl_vertical_tracer_diff   = 1    ! 0=explicit, 1 = implicit
  INTEGER  :: HORZ_VELOC_DIFF_TYPE  = 1          ! 0=no hor.diff; 1=constant Laplacian coefficients
                                                 ! 2=constant coefficients satisfying Munk criterion
                                                 ! 3=variable coefficients satisfying Munk criterion
  INTEGER  :: N_POINTS_IN_MUNK_LAYER = 1
  INTEGER  :: veloc_diffusion_order = 1          !order of friction/diffusion in velocity eq.: 1=laplacian, 2=biharmonic
  INTEGER  :: veloc_diffusion_form  = 1          !form of friction/diffusion operator
                                                 !1: Laplace=curlcurl-graddiv
                                                 !2: Laplace=div k  grad
                                                 !For the corresponding biharmonic choice the laplacian in their form 1 or 2 are iterated
  REAL(wp) :: k_veloc_h             = 1.0E+5_wp  ! horizontal diffusion coefficient
  REAL(wp) :: k_veloc_v             = 1.0E-3_wp  ! vertical diffusion coefficient
  REAL(wp) :: k_pot_temp_h          = 1.0E+3_wp  ! horizontal mixing coefficient for pot. temperature
  REAL(wp) :: k_pot_temp_v          = 1.0E-4_wp  ! vertical mixing coefficient for pot. temperature
  REAL(wp) :: k_sal_h               = 1.0E+3_wp  ! horizontal diffusion coefficient for salinity
  REAL(wp) :: k_sal_v               = 1.0E-4_wp  ! vertical diffusion coefficient for salinity
  REAL(wp) :: MAX_VERT_DIFF_VELOC   = 0.0_wp     ! maximal diffusion coefficient for velocity
  REAL(wp) :: MAX_VERT_DIFF_TRAC    = 0.0_wp     ! maximal diffusion coefficient for tracer
  REAL(wp) :: biharmonic_diffusion_factor = 5.0E12_wp! factor for adjusting the biharmonic diffusion coefficient
                                      !has to be adjusted for each resolution, the bigger this number 
                                      !the smaller becomes the effect of biharmonic diffusion.The appropriate
                                      !size of this number depends also on the position of the biharmonic diffusion coefficient
                                      !within the biharmonic operator. Currently the coefficient is placed in front of the operator.
  LOGICAL  :: l_smooth_veloc_diffusion = .TRUE.

  REAL(wp) :: richardson_veloc      = 0.5E-2_wp  ! Factor with which the richarseon related part of the vertical 
                                                 ! diffusion is multiplied before it is added to the background 
                                                 ! vertical diffusion ! coeffcient for the velocity. See usage in
                                                 ! mo_oce_physics.f90, update_ho_params, variable z_av0
  REAL(wp) :: richardson_tracer     = 0.5E-2_wp  ! see above, valid for tracer instead velocity, see variable z_dv0 in update_ho_params
  LOGICAL  :: l_constant_mixing     = .FALSE.    ! .TRUE.: the vertical mixing coefficients for velocity and tracer
                                                 ! are kept constant over time and are set to the background values; no convection
!  LOGICAL  :: l_convection          = .TRUE.     ! .FALSE.: the vertical mixing coefficients for velocity and tracer
!                                                 ! are unchanged in case of instable stratification
!  LOGICAL  :: l_pp_scheme           = .TRUE.     ! .FALSE.: the vertical mixing coefficients for velocity and tracer
                                                 ! are set to the background values in case of stable stratification
  LOGICAL  :: l_wind_mixing         = .FALSE.    ! .TRUE.: activate wind mixing part of Marsland et al. (2003)
  REAL(wp) :: bottom_drag_coeff     = 2.5E-3_wp  ! chezy coefficient for bottom friction
  REAL(wp) :: wstress_coeff         = 0.3_wp     ! windstress coefficient for analytical wind forcing
                                                 ! 2-dimensional surface relaxation of temperature and salinity
  INTEGER  :: temperature_relaxation= 0          ! 0=no relax.; 1=on for some testcases; 2=use OMIP-file
                                                 ! 3=use initialized values for temperature relaxation
  REAL(wp) :: relaxation_param      = 1.0_wp     ! strength of 2-dim temperatuere relaxation in months
  INTEGER  :: irelax_2d_S           = 0          ! 0=no relax.; 3=use initialized values for relaxation
  REAL(wp) :: relax_2d_mon_S        = 1.0_wp     ! strength of 2-dim salinity relaxation in months
                                                 ! 3-dimensional relaxation of temperature and salinity
  INTEGER  :: irelax_3d_T           = 0          ! 0: no 3-dim relax.,  3: use initial T read with use_file_initialConditions=1
  REAL(wp) :: relax_3d_mon_T        = 1.0_wp     ! strength of 3-dim relaxation for temperature in months
  INTEGER  :: irelax_3d_S           = 0          ! 0: no 3-dim relax.,  3: use initial S read with use_file_initialConditions=1
  REAL(wp) :: relax_3d_mon_S        = 1.0_wp     ! strength of 3-dim relaxation for salinity in months
  LOGICAL  :: limit_elevation       = .FALSE.    ! .TRUE.: balance sea level elevation
  REAL(wp) :: seaice_limit          = 0.5_wp     ! limit sea ice to fraction of surface layer thickness (1.0: no limit)

  INTEGER  :: coriolis_type         = 1          ! 0=zero Coriolis, the non-rotating case
                                                 ! 1=full varying Coriolis
                                                 ! 2=beta-plane (linear) approximation to Coriolis
                                                 ! 3=f-plane (constant) approximation to Coriolis
  ! The variables below are used to set up in basin configuration the Coriolis (f/beta-plane) and
  !   to adjust the analytic wind forcing, units are degrees
  REAL(wp) :: basin_center_lat      = 30.0_wp    ! lat coordinate of basin center
  REAL(wp) :: basin_center_lon      =  0.0_wp    ! lon coordinate of basin center
  REAL(wp) :: basin_width_deg       = 60.0_wp    ! basin extension in zonal direction
  REAL(wp) :: basin_height_deg      = 60.0_wp    ! basin extension in meridional direction
  REAL(wp) :: CWA                   = 5.0E-4_wp  ! Tuning parameters for vertical mixing
  REAL(wp) :: CWT                   = 5.0E-4_wp  !   of tracer and velocity

  LOGICAL  :: lviscous              = .TRUE.    
  LOGICAL  :: l_rigid_lid           = .FALSE.    ! include friction or not
  LOGICAL  :: l_inverse_flip_flop   = .FALSE.    ! true=complete discrete scalarproduct (slow)
                                                 ! false=use a shortcut (faster)
  LOGICAL  :: l_max_bottom          = .TRUE.     ! wet cell: true=if bathy deeper than top
                                                 !           false=bathy deeper mid of cell
  LOGICAL  :: l_edge_based          = .TRUE.     ! mimetic discretization based on edges (true) or cells (false)
  LOGICAL  :: l_partial_cells       = .FALSE.    ! partial bottom cells=true: local varying bottom depth
  LOGICAL  :: l_staggered_timestep  = .FALSE.    ! TRUE=staggering between thermodynamic and dynamic part,
                                                 !   offset of half timestep between dynamic and thermodynamic variables;
                                                 !   thermodynamic and dynamic variables are colocated in time
  INTEGER  :: i_apply_bulk          = 0          ! 0=no bulk formula; 1=apply bulk formula without sea ice
  INTEGER  :: i_sea_ice             = 1          ! 0 = no sea ice; 1=apply sea ice model using sea_ice_nml
  LOGICAL  :: l_relaxsal_ice        = .TRUE.     ! TRUE: relax salinity below sea ice
                                                 ! false = salinity is relaxed under sea ice completely

  LOGICAL  :: l_skip_tracer         = .FALSE.    ! TRUE: no advection and diffusion (incl. convection) of tracer
  LOGICAL  :: use_tracer_x_height   = .FALSE.    ! use the tracer_x_height to calculate advection, in order to minimize round-off errors
  LOGICAL  :: l_with_horz_tracer_diffusion = .TRUE.  ! FALSE: no horizontal tracer diffusion
  LOGICAL  :: l_with_vert_tracer_diffusion = .TRUE.  ! FALSE: no vertical tracer diffusion
  LOGICAL  :: l_with_horz_tracer_advection = .TRUE.  ! FALSE: no horizontal tracer advection
  LOGICAL  :: l_with_vert_tracer_advection = .TRUE.  ! FALSE: no vertical tracer advection
  LOGICAL  :: l_horz_limiter_advection     = .TRUE.  ! FALSE: no horizontal limiter for tracer advection
  LOGICAL  :: l_vert_limiter_advection     = .TRUE.  ! FALSE: no vertical limiter for tracer advection
                                                     ! Note that only in vertical ppm-scheme a limiter is used


  ! special diagnostics configuration
  !
  ! vertex list of throughflows
  INTEGER :: denmark_strait(100)         = -1
  INTEGER :: gibraltar(100)              = -1
  INTEGER :: drake_passage(100)          = -1
  INTEGER :: indonesian_throughflow(100) = -1
  INTEGER :: scotland_iceland(100)       = -1


  NAMELIST/ocean_dynamics_nml/ n_zlev, dzlev_m, discretization_scheme,              &
    &                 iswm_oce, l_staggered_timestep,                      &
    &                 i_bc_veloc_lateral,i_bc_veloc_top,i_bc_veloc_bot,    &
    &                 ab_const, ab_beta, ab_gam, solver_tolerance,         &
    &                 l_RIGID_LID, lviscous, l_inverse_flip_flop,          &
    &                 l_edge_based, i_apply_bulk, l_max_bottom,            &
    &                 l_partial_cells, l_skip_tracer,                      &
    &                 coriolis_type, basin_center_lat, basin_center_lon,   &
    &                 basin_width_deg,basin_height_deg,                    &
    &                 expl_vertical_velocity_diff,                         &
    &                 expl_vertical_tracer_diff,                           &
    &                 veloc_diffusion_order,veloc_diffusion_form,          &
    &                 FLUX_CALCULATION_HORZ, FLUX_CALCULATION_VERT,        &
    &                 dhdtw_abort, l_adpo_flowstrength,                    &
    &                 threshold_min_T, threshold_max_T, threshold_vn,      &
    &                 threshold_min_S, threshold_max_S,                    &
    &                 solver_max_restart_iterations,                       &
    &                 solver_max_iter_per_restart,                         &
    &                 select_solver, use_continuity_correction,            &
    &                 use_edges2edges_viacell_fast


  NAMELIST/ocean_physics_nml/EOS_TYPE, density_computation,                &
    &                 no_tracer, HORZ_VELOC_DIFF_TYPE,                     &
    &                 N_POINTS_IN_MUNK_LAYER,                              &
    &                 k_veloc_h, k_veloc_v,  k_pot_temp_h, k_pot_temp_v,   &
    &                 k_sal_h, k_sal_v,                                    &
    &                 MAX_VERT_DIFF_VELOC, MAX_VERT_DIFF_TRAC,             &
    &                 CWA, CWT,  bottom_drag_coeff, wstress_coeff,         &
    &                 i_sea_ice,                                           &
    &                 biharmonic_diffusion_factor,                         &
    &                 l_smooth_veloc_diffusion,                            &
    &                 richardson_veloc, richardson_tracer,                 &
    &                 l_constant_mixing, l_wind_mixing,                    &
    &                 l_with_vert_tracer_diffusion,                        &
    &                 l_with_horz_tracer_diffusion,                        &
    &                 l_with_vert_tracer_advection,                        &
    &                 l_with_horz_tracer_advection,                        &
    &                 l_horz_limiter_advection, l_vert_limiter_advection,  &
    &                 use_tracer_x_height


  NAMELIST/ocean_forcing_nml/iforc_oce, iforc_type, forcing_timescale,    &
    &                 iforc_stat_oce, init_oce_relax,       &
    &                 l_relaxsal_ice,            &
    &                 temperature_relaxation, relaxation_param,            &
    &                 irelax_2d_S, relax_2d_mon_S,&!relax_2d_T, relax_2d_mon_T, &
    &                 irelax_3d_S, relax_3d_mon_S, irelax_3d_T, relax_3d_mon_T, &
    &                 forcing_enable_freshwater, limit_elevation, seaice_limit,        &
    &                 forcing_set_runoff_to_zero,         &
    &                 use_new_forcing,                    &
    &                 forcing_fluxes_type,                &
    &                 forcing_windstress_u_type,          &
    &                 forcing_windstress_v_type,          &
    &                 forcing_windstress_zonal_waveno,    &
    &                 forcing_windstress_meridional_waveno,&
    &                 forcing_wind_u_amplitude ,          &
    &                 forcing_wind_v_amplitude ,          &
    &                 analytic_wind_amplitude

  !----------------------------------------------------------------------------
  ! initial conditions
  LOGICAL            :: use_file_initialConditions  = .false.
  REAL(wp) :: initial_temperature_top           = 16.0_wp    ! reference temperature used for initialization in testcase 46
  REAL(wp) :: initial_temperature_bottom        = 16.0_wp    ! reference temperature used for initialization in testcase 46
  REAL(wp) :: initial_salinity_top              = 35.0_wp    ! reference salinity used for initialization in testcase 46
  REAL(wp) :: initial_salinity_bottom           = 35.0_wp    ! reference salinity used for initialization in testcase 46
!  INTEGER  :: scatter_levels(10)    = 0          ! levels for possible scattering of the constant tracer fields
!  REAL(wp) :: scatter_t             = 20.0_wp    ! temperature value for scattering
!  REAL(wp) :: scatter_s             = 10.0_wp    ! salinity value for scattering
  INTEGER  :: topography_type  = 0               ! >= 200 analytic bathymetry
  REAL(wp) :: topography_height_reference = 0.0_wp         ! used if topography_type >= 200
  INTEGER  :: sea_surface_height_type  = 0       ! >= 200 sea_surface_height
  INTEGER  :: initial_salinity_type    = 0
  INTEGER  :: initial_temperature_type = 0
  CHARACTER(LEN=max_char_length) :: initial_sst_type = 'sst1'

  ! test cases for ocean model; for the index see run scripts
  INTEGER            :: itestcase_oce  = 0
  NAMELIST/ocean_initialConditions_nml/ use_file_initialConditions, itestcase_oce,      &
    &  initial_temperature_bottom, initial_temperature_top, &
    & initial_salinity_top, initial_salinity_bottom, &
    &  topography_type, topography_height_reference, sea_surface_height_type, &
    &  initial_salinity_type, initial_temperature_type, initial_sst_type
  !----------------------------------------------------------------------------

  NAMELIST/ocean_diagnostics_nml/ diagnostics_level, denmark_strait,drake_passage,gibraltar,  &
    &                 indonesian_throughflow, scotland_iceland

  ! ------------------------------------------------------------------------
  ! 3.0 Namelist variables and auxiliary parameters for octst_nml
  !     This namelists mainly exists during the development of the ocean model
  ! ------------------------------------------------------------------------

  ! location of single cell for input of test values
  INTEGER  :: i_ocv_blk = 1       ! input test block
  INTEGER  :: i_ocv_idx = 1       ! input test index
  INTEGER  :: i_ocv_ilv = 1       ! input test level
  REAL(wp) :: h_val     = 0.0_wp  ! input test value for elevation
  REAL(wp) :: t_val     = 0.0_wp  ! input test value for temperature
  !REAL(wp) :: s_val     = 0.0_wp  ! input  test value for salinity

  ! latitude/longitude location of single cell output for debugging
  REAL(wp) :: rlat_in   = 0.0_wp  ! latitude of cell for debug output
  REAL(wp) :: rlon_in   = 0.0_wp  ! longitude of cell for debug output

  CHARACTER(len=3) :: str_proc_tst(10)   ! namelist string of source processes to print

  NAMELIST/octst_nml/  h_val, t_val, rlat_in, rlon_in

  CONTAINS

 !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
 !>
 !! Initialization of variables that set up the configuration of the ocean model
 !!
 !!               Initialization of variables that set up the configuration
 !!               of the ocean using values read from
 !!               namelist 'ocean_nml' and 'octst_nml'.
 !!
 !! @par Revision History
 !!   Modification by Constantin Junk, MPI-M (2010-02-22)
 !!    - separated subroutine ocean_nml_setup from the original
 !!      setup_run subroutine (which is moved to mo_run_nml)
 !!
 SUBROUTINE setup_ocean_nml( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename

    LOGICAL  :: ignore_land_points = .false.

    NAMELIST/ocean_run_nml/ ignore_land_points

    INTEGER :: i_status, istat
    INTEGER :: iunit

    CHARACTER(len=max_char_length), PARAMETER :: &
            routine = 'mo_ocean_nml/setup_ocean_nml:'

    CALL message(TRIM(routine),'running the hydrostatic ocean model')

     !------------------------------------------------------------
     ! 4.0 set up the default values for ocean_nml
     !------------------------------------------------------------

     ! default values when namelist is not present and no default on definition

     n_zlev            = 5
     dzlev_m(:)        = -99.99_wp!

     dzlev_m(1:n_zlev) =  (/ 50.0_wp, 150.0_wp, 500.0_wp, 1300.0_wp, 2500.0_wp  /)
     !  lower level of layers:  50       200       700       2000       4500
     !  surface coord. levels:  25       125       450       1350       3250

     ! maximal diffusion coefficient for tracer used in implicit vertical tracer diffusion,
     !   if stability criterion is met
     MAX_VERT_DIFF_TRAC  = 100.0_wp * k_veloc_v
     MAX_VERT_DIFF_VELOC = 100.0_wp * k_pot_temp_v

     !------------------------------------------------------------
     ! 5.0 Read ocean_nml namelist
     !------------------------------------------------------------
     ! (done so far by all MPI processes)

     CALL open_nml(TRIM(filename))

     ! setup for the ocean_run_nml
     ignore_land_points = config_ignore_land_points
     CALL position_nml ('ocean_run_nml', status=i_status)
     IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, ocean_run_nml)   ! write defaults to temporary text file
     END IF
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_run_nml, iostat=istat)                           ! overwrite default settings
       IF (my_process_is_stdio()) THEN
         iunit = temp_settings()
         WRITE(iunit, ocean_run_nml)   ! write settings to temporary text file
       END IF
     END SELECT
     config_ignore_land_points = ignore_land_points

     CALL position_nml ('ocean_dynamics_nml', status=i_status)
     IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, ocean_dynamics_nml) ! write defaults to temporary text file
     END IF
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_dynamics_nml, iostat=istat)                         ! overwrite default settings
       IF (my_process_is_stdio()) THEN
         iunit = temp_settings()
         WRITE(iunit, ocean_dynamics_nml) ! write settings to temporary text file
       END IF
     END SELECT

     CALL position_nml ('ocean_physics_nml', status=i_status)
     IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, ocean_physics_nml)    ! write defaults to temporary text file
     END IF
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_physics_nml, iostat=istat)                            ! overwrite default settings
       IF (my_process_is_stdio()) THEN
         iunit = temp_settings()
         WRITE(iunit, ocean_physics_nml)    ! write settings to temporary text file
       END IF
     END SELECT

     CALL position_nml ('ocean_forcing_nml', status=i_status)
     IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, ocean_forcing_nml)  ! write defaults to temporary text file
     END IF
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_forcing_nml, iostat=istat)                          ! overwrite default settings
       IF (my_process_is_stdio()) THEN
         iunit = temp_settings()
         WRITE(iunit, ocean_forcing_nml)  ! write settings to temporary text file
       END IF
     END SELECT

     CALL position_nml ('ocean_initialConditions_nml', status=i_status)
     IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, ocean_initialConditions_nml)  ! write defaults to temporary text file
     END IF
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_initialConditions_nml, iostat=istat)                          ! overwrite default settings
       IF (my_process_is_stdio()) THEN
         iunit = temp_settings()
         WRITE(iunit, ocean_initialConditions_nml)  ! write settings to temporary text file
       END IF
     END SELECT

     CALL position_nml ('ocean_diagnostics_nml', status=i_status)
     IF (my_process_is_stdio()) THEN
       iunit = temp_defaults()
       WRITE(iunit, ocean_diagnostics_nml)   ! write defaults to temporary text file
     END IF
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_diagnostics_nml, iostat=istat)                           ! overwrite default settings
       IF (my_process_is_stdio()) THEN
         iunit = temp_settings()
         WRITE(iunit, ocean_diagnostics_nml)   ! write settings to temporary text file
       END IF
     END SELECT

     !------------------------------------------------------------
     ! 6.0 check the consistency of the parameters
     !------------------------------------------------------------

     IF( iswm_oce == 1 .AND. n_zlev > 1 ) THEN
       CALL message(TRIM(routine),'WARNING, shallow water model (ocean): n_zlev set to 1')
       n_zlev = 1
     ENDIF

     IF(discretization_scheme == 1)THEN
       CALL message(TRIM(routine),'You have choosen the mimetic dicretization')
     ELSEIF(discretization_scheme == 2)THEN
       CALL message(TRIM(routine),'You have choosen the RBF dicretization')
     ELSE
       CALL finish(TRIM(routine), 'wrong parameter for discretization scheme')
     ENDIF

     IF( FLUX_CALCULATION_HORZ > 6 .OR. FLUX_CALCULATION_HORZ <1 ) THEN
       CALL finish(TRIM(routine), 'wrong parameter for advection scheme; use 1-6')
     ENDIF

     IF( FLUX_CALCULATION_VERT > 6 .OR. FLUX_CALCULATION_VERT <1 ) THEN
       CALL finish(TRIM(routine), 'wrong parameter for advection scheme; use 1-6')
     ENDIF

     IF( FLUX_CALCULATION_HORZ == UPWIND .AND. l_horz_limiter_advection ) THEN
       CALL message(TRIM(routine),'WARNING, limiter for horizontal upwind advection set to false')
       l_horz_limiter_advection = .FALSE.
     ENDIF

     IF( FLUX_CALCULATION_HORZ == FLUX_CORR_TRANSP_horz .AND. l_horz_limiter_advection ) THEN
       CALL message(TRIM(routine),'WARNING, limiter for horizontal flux corrected transport set to false')
       l_horz_limiter_advection = .FALSE.
     ENDIF
     
     ! not necessary, limiter within ppm advection (Miura) only
     IF( FLUX_CALCULATION_VERT == UPWIND .AND. l_vert_limiter_advection ) THEN
       CALL message(TRIM(routine),'WARNING, UPWIND advection, limiter for vertical upwind advection set to false')
       l_vert_limiter_advection = .FALSE.
     ENDIF

     ! not necessary, limiter within ppm advection (Miura) only
     IF( FLUX_CALCULATION_VERT == FLUX_CORR_TRANSP_vert_adpo .AND. l_vert_limiter_advection ) THEN
       CALL message(TRIM(routine),'WARNING, ADPO advection, limiter for vertical upwind advection set to false')
       l_vert_limiter_advection = .FALSE.
     ENDIF

     IF(i_bc_veloc_lateral/= 0) THEN
       CALL finish(TRIM(routine), &
         &  'free-slip boundary condition for velocity currently not supported')
     ENDIF
     IF(i_bc_veloc_top < 0.OR.i_bc_veloc_top > 2) THEN
       CALL finish(TRIM(routine), &
         &  'top boundary condition for velocity currently not supported: choose = 0 or =1 or =2')
     ENDIF
     IF(i_bc_veloc_bot < 0 .OR. i_bc_veloc_bot>1) THEN
       CALL finish(TRIM(routine), &
         &  'bottom boundary condition for velocity currently not supported: choose = 0 or =1')
     ENDIF

     IF(no_tracer == 1 .OR. no_tracer < 0 .OR. no_tracer > 2) THEN
       IF(no_tracer == 1) THEN
         CALL message(TRIM(routine), 'WARNING - You have chosen tracer temperature only')
         CALL message(TRIM(routine), ' - this generates error in mo_varlist/mo_oce_state')
         CALL finish(TRIM(routine),  'no_tracer=1 not supported - choose =0 or =2')
       ENDIF
       CALL finish(TRIM(routine),  'no_tracer not supported - choose =0 or =2')
     ENDIF


!     IF (solver_start_tolerance <= 0.0_wp) THEN
!       solver_start_tolerance = solver_tolerance
!       solver_tolerance_decrease_ratio  = 0.1_wp ! must be < 1
!     ENDIF

     IF (forcing_enable_freshwater) THEN
       limit_elevation = .TRUE.
       CALL message(TRIM(routine),'WARNING, limit_elevation set to .TRUE. with forcing_enable_freshwater=.TRUE.')
     END IF

     IF (forcing_set_runoff_to_zero) THEN
       CALL message(TRIM(routine),'WARNING, forcing_set_runoff_to_zero is .TRUE. - forcing with river runoff is set to zero')
     END IF

#ifndef __NO_ICON_ATMO__
     IF ( is_coupled_run() ) THEN
       iforc_oce = FORCING_FROM_COUPLED_FLUX
       CALL message(TRIM(routine),'WARNING, iforc_oce set to 14 for coupled experiment')
       limit_elevation = .FALSE.
       CALL message(TRIM(routine),'WARNING, limit_elevation set to .FALSE. for coupled experiment')
       seaice_limit = 1.0_wp
       CALL message(TRIM(routine),'WARNING, seaice_limit set to 1.0 - no limit for coupled experiment')
     END IF
#endif

     ! write the contents of the namelist to an ASCII file
     IF(my_process_is_stdio()) THEN
       WRITE(nnml_output,nml=ocean_dynamics_nml)
       WRITE(nnml_output,nml=ocean_physics_nml)
       WRITE(nnml_output,nml=ocean_forcing_nml)
       WRITE(nnml_output,nml=ocean_initialConditions_nml)
       WRITE(nnml_output,nml=ocean_diagnostics_nml)
     ENDIF
     !------------------------------------------------------------
     ! 6.0 Read octst_nml namelist
     !------------------------------------------------------------
     ! (done so far by all MPI processes)

     ! 3-char string with marked processes to be printed out for debug purposes
     str_proc_tst =  (/  & 
       &  'all', &  ! initiate print messages in all routines
       &  'abm', &  ! main timestepping routines       in mo_oce_ab_timestepping (mimetic/rbf)
       &  'vel', &  ! velocity advection and diffusion in mo_oce_veloc_advection
       &  'dif', &  ! diffusion                        in mo_oce_diffusion
       &  'trc', &  ! tracer advection and diffusion   in mo_oce_tracer_transport
       &  '   ', &  ! ...
       &  '   ', &
       &  '   ', &
       &  '   ', &
       &  '   '  /)

     CALL position_nml ('octst_nml', status=i_status)
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, octst_nml)
     END SELECT
     CALL close_nml

     ! write the contents of the namelist to an ASCII file
     IF(my_process_is_stdio()) WRITE(nnml_output,nml=octst_nml)

END SUBROUTINE setup_ocean_nml

END MODULE mo_ocean_nml
