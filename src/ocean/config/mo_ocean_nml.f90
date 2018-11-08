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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_ocean_nml
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp, sp
  USE mo_exception,          ONLY: message, warning, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_units,           ONLY: filename_max
  USE mo_physical_constants, ONLY: a_T, b_S,rho_ref, grav, sitodbar
  USE mo_param1_bgc,         ONLY: n_bgctra, ntraad

#ifndef __NO_ICON_ATMO__
  USE mo_coupling_config,    ONLY: is_coupled_run
#endif
  IMPLICIT NONE

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
  INTEGER, PARAMETER :: max_allocated_levels = 1024
  REAL(wp) :: dzlev_m(max_allocated_levels)  ! namelist input of layer thickness

  INTEGER, PARAMETER :: toplev    = 1   ! surface ocean level

  INTEGER :: surface_module = 3         ! 1: old mo_ocean_bulk      - obsolete and no longer technically functional
                                        ! 2: old mo_ocean_surface   - 1st restructuring of sea ice thermodynamics (Feature #5913, Stephan, 03/2015)
                                        ! 3: new mo_ocean_surface_refactor - 2nd restructuring of sea ice thermodynamics (Feature #5936, Vladimir, 03/2017)

  ! switch for reading relaxation data: 1: read from file
  INTEGER :: init_oce_relax = 0
  INTEGER :: relax_analytical_type     = 0 ! special setup for analytic testases, replacement for itestcase_oce in the

!  LOGICAL :: l_time_marching    = .TRUE.  !=.TRUE. is default, the time loop is entered
!                                          !=.FALSE. the time loop is NOT entered and tests with stationary fields can
!                                          !be run. Example: tracer tests with prescribed time-invariant velocity and height.
  ! ----------------------------------------------------------------------------
  ! DIAGNOSTICS
  ! switch for ocean diagnostics - 0: no diagnostics; 1: write to stderr
  INTEGER            :: diagnostics_level      = 1

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
  INTEGER,PARAMETER  :: i_bc_veloc_lateral_noslip   = 0                         
  INTEGER,PARAMETER  :: i_bc_veloc_lateral_freeslip = 1                           
  INTEGER            :: i_bc_veloc_lateral = i_bc_veloc_lateral_noslip   

  INTEGER            :: i_bc_veloc_top = 1  !Top boundary condition for velocity: 
                                            ! i_bc_veloc_top =0 : zero value at top boundary,no wind stress
                                            ! i_bc_veloc_top =1 : forced by wind stress 
                                            !                     stored in p_os%p_aux%bc_top_veloc
                                            ! i_bc_veloc_top =2 : forced by difference between wind
                                            !                     field in p_os%p_aux%bc_top_veloc 
                                            !                     and ocean velocity at top layer
  INTEGER            :: i_bc_veloc_bot = 1  !Bottom boundary condition for velocity:
                                            ! i_bc_veloc_bot =0 : zero value at bottom boundary 
                                            ! i_bc_veloc_bot =1 : bottom boundary friction
                                            ! i_bc_veloc_bot =2 : bottom friction plus topographic
                                            !                     slope (not implemented yet)
  ! Parameters for tracer transport scheme
  !
  !identifiers for different modes of updating the tracer state (substep
  INTEGER, PARAMETER :: use_none    = 0
  INTEGER, PARAMETER :: use_all     = 32767    ! =7FFF

  INTEGER, PARAMETER :: cell_based    = 1
  INTEGER, PARAMETER :: edge_based    = 2

  INTEGER :: tracer_update_mode = use_all
  INTEGER :: tracer_HorizontalAdvection_type = cell_based

  !Options for non-linear corilois term in vector invariant velocity equations
  INTEGER, PARAMETER :: no_coriolis                   = 0
  INTEGER, PARAMETER :: nonlinear_coriolis_dual_grid  = 200 !Default
  INTEGER, PARAMETER :: nonlinear_coriolis_primal_grid= 201
  INTEGER            :: NonlinearCoriolis_type            = nonlinear_coriolis_dual_grid
  INTEGER, PARAMETER :: NoKineticEnergy               = 0
  INTEGER, PARAMETER :: KineticEnergy_onDualGrid      = 200
  INTEGER, PARAMETER :: KineticEnergy_onPrimalGrid    = 201
  INTEGER            :: KineticEnergy_type            = KineticEnergy_onDualGrid
  
  LOGICAL :: l_ANTICIPATED_VORTICITY  = .FALSE.

  !Identifiers for advection schemes
  INTEGER, PARAMETER :: upwind                     = 1
  INTEGER, PARAMETER :: central                    = 2  
  INTEGER, PARAMETER :: lax_friedrichs             = 3
  INTEGER, PARAMETER :: miura_order1               = 4
  INTEGER, PARAMETER :: horz_flux_twisted_vec_recon                   = 5
  INTEGER, PARAMETER :: fct_vert_adpo              = 6
  INTEGER, PARAMETER :: fct_vert_ppm               = 7  
  INTEGER, PARAMETER :: fct_vert_minmod            = 8      
  INTEGER, PARAMETER :: fct_vert_zalesak           = 9      
  !Additional parameters for FCT: High and low order flux calculations can be
  !chosen from list above. Below is the default option
  INTEGER            :: fct_high_order_flux= central
  INTEGER            :: fct_low_order_flux = upwind 
  !Limiters for Flux-Corrected Transport
  INTEGER, PARAMETER :: fct_limiter_horz_zalesak   =100
  INTEGER, PARAMETER :: fct_limiter_horz_minmod    =101
  INTEGER, PARAMETER :: fct_limiter_horz_posdef    =102  
  
  !The default setting concerning tracer advection
  !horizontal
  INTEGER            :: flux_calculation_horz      = horz_flux_twisted_vec_recon      
  INTEGER            :: fct_limiter_horz           = fct_limiter_horz_zalesak!! 
  !vertical
  INTEGER            :: flux_calculation_vert      = fct_vert_ppm !fct_vert_ppm

  
  LOGICAL            :: l_adpo_flowstrength        = .FALSE.   ! .TRUE.: activate second condition for adpo weight

  LOGICAL            :: l_LAX_FRIEDRICHS          =.FALSE.  !Additional LAX-Friedich for horizontal tracer advection of full mimetic scheme horz_flux_twisted_vec_recon
  LOGICAL            :: l_GRADIENT_RECONSTRUCTION = .FALSE. !Additional Gradient-Reconstruction for horizontal tracer advection of full mimetic scheme horz_flux_twisted_vec_recon
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
  REAL(wp) :: solver_tolerance                   = 1.e-14_wp   ! Maximum value allowed for solver absolute tolerance
  REAL(wp) :: MassMatrix_solver_tolerance        = 1.e-11_wp   ! Maximum value allowed for solver absolute tolerance
  !  REAL(wp) :: solver_start_tolerance          = -1.0_wp
  INTEGER  :: solver_max_restart_iterations      = 100       ! For restarting gmres
  INTEGER  :: solver_max_iter_per_restart        = 200       ! For inner loop after restart
  INTEGER  :: solver_max_iter_per_restart_sp     = 200       ! For inner loop after restart
  REAL(sp) :: solver_tolerance_sp                = 1.e-11_sp   ! Maximum value allowed for solver absolute tolerance
  LOGICAL  :: use_absolute_solver_tolerance      = .true.   ! Maximum value allowed for solver tolerance
  INTEGER, PARAMETER :: select_gmres             = 1
  INTEGER, PARAMETER :: select_restart_gmres     = 2
  INTEGER, PARAMETER :: select_restart_mixedPrecision_gmres     = 3
  INTEGER :: select_solver                       = select_restart_gmres
  INTEGER :: solver_FirstGuess                   = 0

  INTEGER, PARAMETER :: select_lhs_operators = 1
  INTEGER, PARAMETER :: select_lhs_matrix    = 2
  INTEGER :: select_lhs = select_lhs_operators




  LOGICAL :: use_continuity_correction           = .true.  
  INTEGER :: fast_performance_level              = 50 ! 5  ! 0= most safe, bit identical results, should be fast_sum = .false.
                                                      ! 1 = no optimized calls
                                                      ! 5 = standard (use of gmres restart)
                                                      ! > 10 = latest performnce optimizations
!   LOGICAL :: use_edges2edges_viacell_fast        = .false.

  INTEGER           :: MASS_MATRIX_INVERSION_TYPE =0  
  INTEGER,PARAMETER :: NO_INVERSION               =0
  INTEGER,PARAMETER :: MASS_MATRIX_INVERSION_ADVECTION =1
  INTEGER,PARAMETER :: MASS_MATRIX_INVERSION_ALLTERMS  =2


  ! physical parameters for  aborting the ocean model
  REAL(wp) :: dhdtw_abort           =  3.17e-11_wp  ! abort criterion for gmres solution (~1mm/year)
  REAL(wp) :: threshold_vn          = 10.0_wp    ! abort criterion for absolute velocity maximum
  REAL(wp) :: threshold_min_T       = -4.0_wp    ! abort criterion for temperature minimum
  REAL(wp) :: threshold_max_T       = 100._wp    ! abort criterion for temperature minimum
  REAL(wp) :: threshold_min_S       =  0.0_wp    ! abort criterion for salinity minimum
  REAL(wp) :: threshold_max_S       = 60.0_wp    ! abort criterion for salinity minimum

  REAL(wp) :: tracer_threshold_min(8), tracer_threshold_max(8) ! as above but with indexes
  CHARACTER(LEN=8) :: namelist_tracer_name(12) 

  INTEGER  :: no_tracer             = 2          ! number of tracers

  ! more ocean parameters, not yet well placed
!   INTEGER  :: expl_vertical_velocity_diff = 1    ! 0=explicit, 1 = implicit
!   INTEGER, PARAMETER :: explicit_diffusion = 2
!   INTEGER, PARAMETER :: implicit_diffusion = 1
!   INTEGER  :: vertical_tracer_diffusion_type   = 1   ! not used !
  INTEGER  :: VelocityDiffusion_order = 1         !1=laplacian, 2=biharmonic, 21=laplacian+biharmonic
  INTEGER  :: laplacian_form  = 1          !form of friction/diffusion operator
                                                 !1: Laplace=curlcurl-graddiv
                                                 !2: Laplace=div k  grad
                                                 !For the corresponding biharmonic choice the laplacian in their form 1 or 2 are iterated

  REAL(wp) :: bottom_drag_coeff     = 2.5E-3_wp  ! chezy coefficient for bottom friction
                                                 ! 2-dimensional surface relaxation of temperature and salinity:
  INTEGER  :: type_surfRelax_Temp  = 0           ! 0=no relax.; 1=on for some testcases; 2=use OMIP-file
                                                 ! 3=use initialized values for temperature relaxation
  REAL(wp) :: para_surfRelax_Temp   = 1.0_wp     ! strength of 2-dim temperatuere relaxation in months
  INTEGER  :: type_surfRelax_Salt   = 0          ! 0=no relax.; 3=use initialized values for relaxation
  REAL(wp) :: para_surfRelax_Salt   = 1.0_wp     ! strength of 2-dim salinity relaxation in months
                                                 ! 3-dimensional relaxation of temperature and salinity
  INTEGER  :: type_3dimRelax_Temp   = 0          ! 0: no 3-dim relax.,  3: use initial T read with use_file_initialConditions=1
  REAL(wp) :: para_3dimRelax_Temp   = 1.0_wp     ! strength of 3-dim relaxation for temperature in months
  INTEGER  :: type_3dimRelax_Salt   = 0          ! 0: no 3-dim relax.,  3: use initial S read with use_file_initialConditions=1
  REAL(wp) :: para_3dimRelax_Salt   = 1.0_wp     ! strength of 3-dim relaxation for salinity in months
  LOGICAL  :: limit_elevation       = .FALSE.    ! .TRUE.: balance sea level elevation
  LOGICAL  :: limit_seaice          = .TRUE.     ! .TRUE.: set a cutoff limit to sea ice thickness
  REAL(wp) :: seaice_limit          = 0.4_wp     ! limit sea ice thickness to fraction of surface layer thickness

  INTEGER  :: coriolis_type         = 1          ! 0=zero Coriolis, the non-rotating case
                                                 ! 1=full varying Coriolis
                                                 ! 2=beta-plane (linear) approximation to Coriolis
                                                 ! 3=f-plane (constant) approximation to Coriolis
  REAL(wp) :: coriolis_fplane_latitude = 0.0_wp
  ! The variables below are used to set up in basin configuration the Coriolis (f/beta-plane) and
  !   to adjust the analytic wind forcing, units are degrees
  REAL(wp) :: basin_center_lat      = 30.0_wp    ! lat coordinate of basin center
  REAL(wp) :: basin_center_lon      =  0.0_wp    ! lon coordinate of basin center
  REAL(wp) :: basin_width_deg       = 60.0_wp    ! basin extension in zonal direction
  REAL(wp) :: basin_height_deg      = 60.0_wp    ! basin extension in meridional direction

  LOGICAL  :: lviscous              = .TRUE.
  LOGICAL  :: l_rigid_lid           = .FALSE.    ! include friction or not
  LOGICAL  :: l_inverse_flip_flop   = .FALSE.    ! true=complete discrete scalarproduct (slow)
                                                 ! false=use a shortcut (faster)
  LOGICAL  :: l_max_bottom          = .TRUE.     ! wet cell: true=if bathy deeper than top
                                                 !           false=bathy deeper mid of cell
  LOGICAL  :: l_edge_based          = .TRUE.     ! mimetic discretization based on edges (true) or cells (false)
  LOGICAL  :: l_partial_cells       = .FALSE.    ! partial bottom cells=true: local varying bottom depth

  INTEGER  :: i_sea_ice             = 1          ! 0 = no sea ice; 1=apply sea ice model using sea_ice_nml
  LOGICAL  :: l_relaxsal_ice        = .TRUE.     ! TRUE: relax salinity below sea ice
                                                 ! false = salinity is relaxed under sea ice completely

  LOGICAL  :: use_tracer_x_height          = .FALSE. ! use the tracer_x_height to calculate advection, in order to minimize round-off errors
  LOGICAL  :: l_with_horz_tracer_diffusion = .TRUE.  ! FALSE: no horizontal tracer diffusion
  LOGICAL  :: l_with_vert_tracer_diffusion = .TRUE.  ! FALSE: no vertical tracer diffusion
  LOGICAL  :: l_with_horz_tracer_advection = .TRUE.  ! FALSE: no horizontal tracer advection
  LOGICAL  :: l_with_vert_tracer_advection = .TRUE.  ! FALSE: no vertical tracer advection

  ! cfl related
  LOGICAL  :: cfl_check             = .TRUE.
  LOGICAL  :: cfl_write             = .FALSE.
  LOGICAL  :: cfl_stop_on_violation = .FALSE.
  REAL(wp) :: cfl_threshold         = 1.0_wp

  INTEGER, PARAMETER :: VerticalAdvection_None           = 0
  INTEGER, PARAMETER :: VerticalAdvection_MimeticRotationalForm = 1
  INTEGER, PARAMETER :: VerticalAdvection_DivergenceForm = 2
  INTEGER, PARAMETER :: VerticalAdvection_RotationalForm = 3
  INTEGER :: HorizonatlVelocity_VerticalAdvection_form = VerticalAdvection_MimeticRotationalForm

  LOGICAL  :: use_smooth_ocean_boundary  = .TRUE.
  LOGICAL  :: createSolverMatrix  = .FALSE.

  ! adjoint run setup
  ! number of checkpoints
  INTEGER :: ncheckpoints           = -1
  INTEGER, PARAMETER :: RUN_FORWARD            = 0
  INTEGER, PARAMETER :: RUN_ADJOINT            = 1
  INTEGER :: run_mode               = 0

  NAMELIST/ocean_dynamics_nml/&
    &                 ab_beta                      , &
    &                 ab_const                     , &
    &                 ab_gam                       , &
    &                 basin_center_lat             , &
    &                 basin_center_lon             , &
    &                 basin_height_deg             , &
    &                 basin_width_deg              , &
    &                 cfl_check                    , &
    &                 cfl_write                    , &
    &                 cfl_stop_on_violation        , &
    &                 cfl_threshold                , &
    &                 coriolis_type                , &
    &                 coriolis_fplane_latitude       , &
    &                 dhdtw_abort                  , &
    &                 discretization_scheme        , &
    &                 dzlev_m                      , &
    &                 i_bc_veloc_bot               , &
    &                 i_bc_veloc_lateral           , &
    &                 i_bc_veloc_top               , &
    &                 iswm_oce                     , &
    &                 l_RIGID_LID                  , &
    &                 l_edge_based                 , &
    &                 l_inverse_flip_flop          , &
    &                 l_max_bottom                 , &
    &                 l_partial_cells              , &
    &                 lviscous                     , &
    &                 n_zlev                       , &
    &                 select_solver                , &
    &                 use_absolute_solver_tolerance, &
    &                 solver_max_iter_per_restart  , &
    &                 solver_max_restart_iterations, &
    &                 solver_tolerance             , &
    &                 solver_max_iter_per_restart_sp, &
    &                 solver_tolerance_sp          , &
    &                 select_lhs                   , &
    &                 MassMatrix_solver_tolerance  , &
    &                 threshold_vn                 , &
    &                 surface_module               , &
    &                 use_continuity_correction    , &
    &                 fast_performance_level       , &
    &                 MASS_MATRIX_INVERSION_TYPE   , &
    &                 NonlinearCoriolis_type       , &
    &                 KineticEnergy_type           , &
    &                 HorizonatlVelocity_VerticalAdvection_form, &
    &                 solver_FirstGuess            , &
    &                 use_smooth_ocean_boundary    , &
    &                 run_mode                     , &
    &                 ncheckpoints                 , &
    &                 createSolverMatrix


  NAMELIST/ocean_tracer_transport_nml/&
    &                 no_tracer                    , &  
    &                 flux_calculation_horz        , &
    &                 flux_calculation_vert        , & 
    &                 fct_high_order_flux          , &
    &                 fct_low_order_flux           , &
    &                 fct_limiter_horz             , & 
    &                 l_adpo_flowstrength          , &     
    &                 l_with_horz_tracer_advection , &
    &                 l_with_horz_tracer_diffusion , &
    &                 l_with_vert_tracer_advection , &
    &                 l_with_vert_tracer_diffusion , &
    &                 use_tracer_x_height          , &           
    &                 threshold_max_S              , &
    &                 threshold_max_T              , &
    &                 threshold_min_S              , &
    &                 threshold_min_T              , &
    &                 tracer_update_mode           , &
    &                 tracer_HorizontalAdvection_type, &
    &                 l_LAX_FRIEDRICHS             , &
    &                 l_GRADIENT_RECONSTRUCTION


  ! tracer horizontal diffusion
  REAL(wp) :: Temperature_HorizontalDiffusion_Background = 0.0_wp
  REAL(wp) :: Temperature_HorizontalDiffusion_Reference  = 1.0E+3_wp
  REAL(wp) :: Salinity_HorizontalDiffusion_Background    = 0.0_wp
  REAL(wp) :: Salinity_HorizontalDiffusion_Reference     = 1.0E+3_wp
  INTEGER  :: TracerHorizontalDiffusion_scaling          = 1 ! 1= constant, 5=scale with edge (dual) **3
  REAL(wp) :: TracerHorizontalDiffusion_ScaleWeight      = 1.0_wp
  REAL(wp) :: Tracer_HorizontalDiffusion_PTP_coeff       = 1.0E+3_wp  ! horizontal mixing coefficient for ptp 

  REAL(wp) :: TracerDiffusion_LeithWeight = 0.0_wp ! if Leith is active then the Leith coeff*this id added to the tracer diffusion coeff
  REAL(wp) :: max_turbulenece_TracerDiffusion = 4.0_wp ! max tracer diffusion amplification from turbulenece on top of the standard one

  ! tracer vertical diffusion
  INTEGER, PARAMETER  :: PPscheme_Constant_type   = 0  ! are kept constant over time and are set to the background values; no convection
  INTEGER, PARAMETER  :: PPscheme_ICON_type       = 1
  INTEGER, PARAMETER  :: PPscheme_MPIOM_type      = 2
  INTEGER, PARAMETER  :: PPscheme_ICON_Edge_type  = 3
  INTEGER, PARAMETER  :: PPscheme_ICON_Edge_vnPredict_type = 4
  INTEGER  :: PPscheme_type = PPscheme_MPIOM_type

  REAL(wp) :: tracer_convection_MixingCoefficient        = 0.1_wp     ! convection diffusion coefficient for tracer, used in PP scheme
  REAL(wp) :: convection_InstabilityThreshold            = -5.0E-8_wp ! used in PP scheme
  REAL(wp) :: RichardsonDiffusion_threshold              =  5.0E-8_wp ! used in PP scheme
  REAL(wp) :: tracer_RichardsonCoeff                     =  1.5E-3  ! factor for vertical diffusion coefficient in PP schemes
  REAL(wp) :: Temperature_VerticalDiffusion_background   = 1.5E-5   ! vertical mixing coefficient for pot. temperature
  REAL(wp) :: Salinity_VerticalDiffusion_background      = 1.5E-5   ! vertical diffusion coefficient for salinity
  REAL(wp) :: Salinity_ConvectionRestrict = 0.0_wp ! do not change !
  REAL(wp) :: richardson_tracer     = 0.5E-2_wp  ! see above, valid for tracer instead velocity, see variable z_dv0 in update_ho_params
  REAL(wp) :: lambda_wind           = 0.05_wp     ! 0.03_wp for 20km omip   !  wind mixing stability parameter, eq. (16) of Marsland et al. (2003)
  REAL(wp) :: wma_diff              = 5.0e-4_wp  !  wind mixing amplitude for diffusivity
  REAL(wp) :: wma_visc              = 5.0e-4_wp  !  wind mixing amplitude for viscosity
  LOGICAL  :: use_wind_mixing = .FALSE.          ! .TRUE.: wind mixing parametrization switched on
  LOGICAL  :: use_reduced_mixing_under_ice = .TRUE. ! .TRUE.: reduced wind mixing under sea ice in pp-scheme
  REAL(wp) :: tracer_TopWindMixing   = 2.5E-4_wp
  REAL(wp) :: velocity_TopWindMixing = 2.5E-4_wp
  REAL(wp) :: WindMixingDecayDepth  = 40.0

 
  ! velocity diffusion
  REAL(wp) :: VerticalViscosity_TimeWeight = 0.0_wp
  REAL(wp) :: velocity_VerticalDiffusion_background      = 1.0E-3_wp  ! vertical diffusion coefficient
  REAL(wp) :: velocity_RichardsonCoeff      = 0.5E-2_wp  ! Factor with which the richarseon related part of the vertical
                                                 ! diffusion is multiplied before it is added to the background
                                                 ! vertical diffusion coeffcient for the velocity. See usage in
                                                 ! mo_pp_scheme.f90

  REAL(wp) :: BiharmonicViscosity_background = 0.0_wp
  REAL(wp) :: BiharmonicViscosity_reference  = 0.0_wp ! 1.0E-2_wp
  INTEGER  :: BiharmonicViscosity_scaling = 1         ! 6
  REAL(wp) :: BiharmonicVort_weight = 1.0_wp
  REAL(wp) :: BiharmonicDiv_weight = 1.0_wp

  REAL(wp) :: HarmonicViscosity_background = 0.0_wp
  REAL(wp) :: HarmonicViscosity_reference  = 0.0_wp ! 4.0E-2_wp
  INTEGER  :: HarmonicViscosity_scaling = 1         ! 2
  REAL(wp) :: HarmonicVort_weight = 1.0_wp
  REAL(wp) :: HarmonicDiv_weight = 1.0_wp

  INTEGER  :: N_POINTS_IN_MUNK_LAYER = 1

  REAL(wp) :: LeithHarmonicViscosity_background = 0.0_wp
  REAL(wp) :: LeithHarmonicViscosity_reference = 0.0_wp   !3.82E-12_wp
  INTEGER  :: LeithHarmonicViscosity_scaling = 1          ! 6
  REAL(wp) :: LeithBiharmonicViscosity_background = 0.0_wp
  REAL(wp) :: LeithBiharmonicViscosity_reference  = 0.0_wp   ! 3.82E-12_wp
  INTEGER  :: LeithBiharmonicViscosity_scaling = 1           ! 6
  INTEGER  :: LeithClosure_order = 0 ! 1=laplacian, 2= biharmonc, 21 =laplacian+biharmonc
  INTEGER  :: LeithClosure_form = 0  ! 1=vort, 2=vort+div
!   REAL(wp) :: LeithClosure_gamma = 0.25_wp !dimensionless constant for Leith closure, not used
  INTEGER  :: HorizontalViscosity_SmoothIterations = 0
  REAL(wp) :: HorizontalViscosity_SpatialSmoothFactor = 0.5_wp
  REAL(wp) :: HorizontalViscosity_ScaleWeight = 0.5_wp

 NAMELIST/ocean_horizontal_diffusion_nml/&
    & &! define harmonic and biharmonic parameters !
    &  laplacian_form,                  & ! 1=curlcurl-graddiv, 2=div k  grad
    &  VelocityDiffusion_order        , & ! 1=harmonic, 2=biharmonic, 21=harmonic+biharmonc
    &  N_POINTS_IN_MUNK_LAYER         , &
    &  BiharmonicViscosity_scaling,     & ! the scaling type for the biharmonic viscosity
    &  HarmonicViscosity_scaling,       & ! the scaling type for the harmonic viscosity
    &  HarmonicViscosity_background,    & ! the harmonic viscosity background value (not scaled)
    &  HarmonicViscosity_reference,     & ! the harmonic viscosity parameter for scaling
    &  HarmonicVort_weight,             &
    &  HarmonicDiv_weight,              &
    &  BiharmonicViscosity_background,  & ! the biharmonic viscosity background value (not scaled)
    &  BiharmonicViscosity_reference,   & ! the biharmonic viscosity parameter for scaling
    &  BiharmonicVort_weight,           &
    &  BiharmonicDiv_weight,            &
    &  HorizontalViscosity_SmoothIterations,      & ! smoothing iterations for the scaled viscosity (both harmonic and biharmonic)
    &  HorizontalViscosity_SpatialSmoothFactor,   & ! the weight of the neigbors during the smoothing
    & &
    & &  ! define tracer horizontal diffusion parameters !
    &  TracerHorizontalDiffusion_scaling,         & ! the scaling type for the  tracer diffusion
    &  Temperature_HorizontalDiffusion_Background,&
    &  Temperature_HorizontalDiffusion_Reference, &
    &  Salinity_HorizontalDiffusion_Background,   &
    &  Salinity_HorizontalDiffusion_Reference,    &
    & &
    & & ! define Leith parameters
    &  LeithClosure_order,              & ! 1=harmonic, 2=biharmonic, 21=biharmonic+harmonic
    &  LeithClosure_form,               & ! 1=rotation only, 2=rot+div, 4=rot+div using a div grad harmonic form
    &  LeithHarmonicViscosity_background,       &
    &  LeithHarmonicViscosity_reference,        &
    &  LeithHarmonicViscosity_scaling,          &
    &  LeithBiharmonicViscosity_background,     &
    &  LeithBiharmonicViscosity_reference,      &
    &  LeithBiharmonicViscosity_scaling,        &
    &  TracerDiffusion_LeithWeight,             &     ! experimental, do not use!
    &  max_turbulenece_TracerDiffusion, &  ! experimental, do not use!
    & &
    & & ! other
    & Tracer_HorizontalDiffusion_PTP_coeff

  NAMELIST/ocean_vertical_diffusion_nml/&
    &  PPscheme_type               ,&         !2=as in MPIOM, 4=used for higher resolutions
    &  VerticalViscosity_TimeWeight,&         ! timeweight of the vertical viscosity calculated from the previous velocity (valid only with PPscheme_type=4)
    &  Temperature_VerticalDiffusion_background, &
    &  Salinity_VerticalDiffusion_background,    &
    &  velocity_VerticalDiffusion_background,    &
    &  bottom_drag_coeff           ,&
    &  velocity_RichardsonCoeff    ,&
    &  tracer_RichardsonCoeff,      &
    &  lambda_wind                 ,&
    &  wma_visc                    ,&
    &  use_reduced_mixing_under_ice,&
    &  use_wind_mixing,             &
    &  tracer_TopWindMixing,        &
    &  WindMixingDecayDepth,        &
    &  velocity_TopWindMixing,      &
    &  tracer_convection_MixingCoefficient ,    &
    &  convection_InstabilityThreshold, &
    &  RichardsonDiffusion_threshold,   &
    &  Salinity_ConvectionRestrict

  !Parameters for GM-Redi configuration
  REAL(wp) :: k_tracer_dianeutral_parameter   = 1.0E-4_wp  !dianeutral tracer diffusivity for GentMcWilliams-Redi parametrization
  REAL(wp) :: k_tracer_isoneutral_parameter   = 600.0_wp  !isoneutral tracer diffusivity for GentMcWilliams-Redi parametrization
  REAL(wp) :: k_tracer_GM_kappa_parameter     = 600.0_wp  !kappa parameter in GentMcWilliams parametrization
  INTEGER            :: GMRedi_configuration=0
  INTEGER, PARAMETER :: Cartesian_Mixing=0
  INTEGER, PARAMETER :: GMRedi_combined =1   !both parametrizations active
  INTEGER, PARAMETER :: GM_only         =2
  INTEGER, PARAMETER :: Redi_only       =3
  !Parameters for tapering configuration
  INTEGER, PARAMETER :: tapering_DanaMcWilliams=1
  INTEGER, PARAMETER :: tapering_Large=2
  INTEGER, PARAMETER :: tapering_Griffies=3
  INTEGER            :: tapering_scheme=tapering_DanaMcWilliams
  LOGICAL            :: switch_off_diagonal_vert_expl=.TRUE.
  LOGICAL            :: GMREDI_COMBINED_DIAGNOSTIC=.TRUE.
  LOGICAL            :: GM_INDIVIDUAL_DIAGNOSTIC=.TRUE.  
  LOGICAL            :: REDI_INDIVIDUAL_DIAGNOSTIC=.TRUE.    
  LOGICAL            :: TEST_MODE_GM_ONLY=.FALSE.
  LOGICAL            :: TEST_MODE_REDI_ONLY=.FALSE.
  LOGICAL            :: SWITCH_OFF_TAPERING=.FALSE.
  LOGICAL            :: SWITCH_ON_REDI_BALANCE_DIAGONSTIC=.FALSE.
  LOGICAL            :: SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION=.FALSE.
  LOGICAL            :: SLOPE_CALC_VIA_TEMPERTURE_SALINITY=.FALSE.
  LOGICAL            :: BOLUS_VELOCITY_DIAGNOSTIC=.FALSE.
  LOGICAL            :: REVERT_VERTICAL_RECON_AND_TRANSPOSED=.FALSE.  
  LOGICAL            :: INCLUDE_SLOPE_SQUARED_IMPLICIT=.TRUE.
  !Parameters for tapering schemes
  LOGICAL  :: GMRedi_usesRelativeMaxSlopes = .true. ! the slopes are defined relatively the the grid slopes: dz/dx
  REAL(wp) :: S_max      = 1.0_wp                   !maximally allowed slope
  REAL(wp) :: S_critical = 0.25_wp                   !critical value at which tapering reduces slope by 50%
  REAL(wp) :: S_d        = 0.01_wp                   !width of transition zone from untapered to tapered
  REAL(wp) :: c_speed    = 2.0_wp                   !aproximation to first baroclinic wave speed. Used in tapering scheme "Large"
  REAL(wp) :: RossbyRadius_min = 15000.0_wp
  REAL(wp) :: RossbyRadius_max =100000.0_wp


 NAMELIST/ocean_GentMcWilliamsRedi_nml/&
    &  GMRedi_configuration           ,&
    &  tapering_scheme                ,&
    &  GMRedi_usesRelativeMaxSlopes   ,&
    &  S_max                          ,&
    &  S_d                            ,&
    &  S_critical                     ,&
    &  c_speed                        ,&
    &  k_tracer_dianeutral_parameter  ,&
    &  k_tracer_isoneutral_parameter  ,&
    &  k_tracer_GM_kappa_parameter    ,&
    &  RossbyRadius_min               ,&
    &  RossbyRadius_max               ,&
    &  switch_off_diagonal_vert_expl, &
    &  GMREDI_COMBINED_DIAGNOSTIC,    &
    & TEST_MODE_GM_ONLY,              &
    & TEST_MODE_REDI_ONLY,            &
    & SWITCH_OFF_TAPERING,            &
    & SWITCH_ON_REDI_BALANCE_DIAGONSTIC,&
    & SWITCH_ON_TAPERING_HORIZONTAL_DIFFUSION,&
    & SLOPE_CALC_VIA_TEMPERTURE_SALINITY,&
    & BOLUS_VELOCITY_DIAGNOSTIC,         &
    & REVERT_VERTICAL_RECON_AND_TRANSPOSED,&
    & INCLUDE_SLOPE_SQUARED_IMPLICIT 
  
  INTEGER  :: EOS_TYPE              = 2          ! 1=linear EOS,2=(nonlinear, from MPIOM)
                                                 ! 3=nonlinear Jacket-McDoudgall-formulation (not yet recommended)
                                                 ! 10 = EOS10, note that the GMRedo needs to be updated to use the EOS10
  REAL(wp) :: LinearThermoExpansionCoefficient = a_T
  REAL(wp) :: LinearHalineContractionCoefficient=b_S
  REAL(wp) :: OceanReferenceDensity = rho_ref
  REAL(wp) :: OceanReferenceDensity_inv
  REAL(wp) :: ReferencePressureIndbars
  
  ! ist : todo move into different nml
  LOGICAL  :: lhamocc=.FALSE.
  LOGICAL  :: lbgcadv=.FALSE.
  INTEGER  :: nbgctra, nbgcadv 
                                 
  
  NAMELIST/ocean_physics_nml/&
    &  EOS_TYPE                    , &
    &  i_sea_ice                   , &
    &  LinearThermoExpansionCoefficient,  &
    &  LinearHalineContractionCoefficient,&
    &  OceanReferenceDensity,       &
    &  lhamocc, lbgcadv

  ! ------------------------------------------------------------------------
  ! FORCING {
  ! iforc_oce: parameterized forcing for ocean model:
  INTEGER, PARAMETER :: No_Forcing                = 10
  INTEGER, PARAMETER :: Analytical_Forcing        = 11
  INTEGER, PARAMETER :: OMIP_FluxFromFile         = 12  ! OMIP or NCEP type forcing
  INTEGER, PARAMETER :: Atmo_FluxFromFile         = 13  ! not yet
  INTEGER, PARAMETER :: Coupled_FluxFromAtmo      = 14  ! parameter for a coupled atmosphere-ocean run
  INTEGER, PARAMETER :: Coupled_FluxFromFile      = 15  ! not yet
  INTEGER            :: iforc_oce                 =  0  ! index of parameterized forcing

  ! new/renamed switches
  ! length of time varying flux forcing: 12: read 12 months, other: read daily values
  INTEGER  :: forcing_timescale                    = 1
  LOGICAL  :: forcing_enable_freshwater            = .TRUE.    ! .TRUE.: apply freshwater forcing boundary condition
  LOGICAL  :: forcing_set_runoff_to_zero           = .FALSE.   ! .TRUE.: set river runoff to zero for comparion to MPIOM
  LOGICAL  :: zero_freshwater_flux                 = .FALSE.   ! .TRUE.: zero freshwater fluxes but salt-change possible
  LOGICAL  :: use_new_forcing                      = .FALSE.
  ! _type variables range
  !    0    : not used
  !   1:100 : file based input
  ! 101:200 : analytic setup
  ! forcing_windstress_(u|v|fluxes)_type values
  ! 1   : omip input
  ! 5   : ncep input
  INTEGER  :: forcing_fluxes_type                  = 0
  INTEGER  :: forcing_windstress_u_type            = 0
  INTEGER  :: forcing_windstress_v_type            = 0

  REAL(wp) :: forcing_windstress_zonal_waveno      = 3.0_wp  ! For the periodic analytic forcing (wind)
  REAL(wp) :: forcing_windstress_zonalWavePhase    = 0.0_wp
!DR  REAL(wp) :: forcing_windstress_meridional_waveno = 3.0_wp
  REAL(wp) :: forcing_windstress_merid_waveno      = 3.0_wp
  REAL(wp) :: forcing_windStress_u_amplitude       = 1.0_wp
  REAL(wp) :: forcing_windStress_v_amplitude       = 1.0_wp
  REAL(wp) :: forcing_center                       = 0.0_wp
  INTEGER  :: forcing_smooth_steps                 = 1
  REAL(wp) :: forcing_windStress_weight            = 1.0_wp
  INTEGER  :: forcing_windspeed_type               = 0
  REAL(WP) :: forcing_windspeed_amplitude          = 1.0_wp
  REAL(wp) :: relax_temperature_min                = 10.0_wp  ! in cases of analytic relaxation
  REAL(wp) :: relax_temperature_max                = 10.0_wp  ! in cases of analytic relaxation
  REAL(wp) :: relax_width           = 1.5_wp     ! the spacial width in degrees where relaxation is applied
  REAL(wp) :: forcing_HeatFlux_amplitude = 10._wp
  REAL(wp) :: forcing_HeatFlux_base      = 0.0_wp
  REAL(wp) :: forcing_temperature_poleLat          = 90.0_wp  ! place the pole at this latitude
                                                              ! for temperature forcing (degrees)
  INTEGER  :: atmos_flux_analytical_type           = 0        ! type of atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_SWnet_const                    = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_LWnet_const                    = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_lat_const                      = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_sens_const                     = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_SWnetw_const                   = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_LWnetw_const                   = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_latw_const                     = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_sensw_const                    = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_precip_const                   = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
  REAL(wp) :: atmos_evap_const                     = 0.0_wp   ! constant atmospheric fluxes for analytical forcing
!   INTEGER  :: windstress_smoothIterations          = 0
!   REAL(wp) :: windstress_smoothWeight              = 0.0_wp
                                                              

  NAMELIST/ocean_forcing_nml/&
    &                 forcing_center                      , &
    &                 forcing_enable_freshwater           , &
    &                 zero_freshwater_flux                , &
    &                 forcing_fluxes_type                 , &
    &                 forcing_set_runoff_to_zero          , &
    &                 forcing_timescale                   , &
    &                 forcing_windStress_u_amplitude      , &
    &                 forcing_windStress_v_amplitude      , &
!DR    &                 forcing_windstress_meridional_waveno, &
    &                 forcing_windstress_merid_waveno     , &
    &                 forcing_windstress_u_type           , &
    &                 forcing_windstress_v_type           , &
    &                 forcing_windstress_zonal_waveno     , &
    &                 forcing_windstress_zonalWavePhase   , &
    &                 forcing_windspeed_type              , &
    &                 forcing_windspeed_amplitude         , &
!     &                 windstress_smoothIterations         , &
!     &                 windstress_smoothWeight            , &
    &                 forcing_HeatFlux_amplitude   , &
    &                 forcing_HeatFlux_base        , &
    &                 iforc_oce                           , &
    &                 init_oce_relax                      , &
    &                 type_surfRelax_Salt                 , &
    &                 type_3dimRelax_Temp                 , &
    &                 type_3dimRelax_Salt                 , &
    &                 l_relaxsal_ice                      , &
    &                 limit_elevation                     , &
    &                 para_surfRelax_Temp                 , &
    &                 para_surfRelax_Salt                 , &
    &                 para_3dimRelax_Temp                 , &
    &                 para_3dimRelax_Salt                 , &
    &                 relax_analytical_type               , &
    &                 limit_seaice                        , &
    &                 seaice_limit                        , &
    &                 type_surfRelax_Temp                 , &
    &                 relax_temperature_min               , &
    &                 relax_temperature_max               , &
    &                 relax_width,                          &
    &                 atmos_flux_analytical_type          , &
    &                 atmos_SWnet_const                   , &
    &                 atmos_LWnet_const                   , &
    &                 atmos_lat_const                     , &
    &                 atmos_sens_const                    , &
    &                 atmos_SWnetw_const                  , &
    &                 atmos_LWnetw_const                  , &
    &                 atmos_latw_const                    , &
    &                 atmos_sensw_const                   , &
    &                 atmos_precip_const                  , &
    &                 atmos_evap_const                    , &
    &                 forcing_temperature_poleLat         , &
    &                 forcing_smooth_steps                , &
    &                 forcing_windStress_weight           , &
    &                 use_new_forcing                    
  ! } END FORCING

  !----------------------------------------------------------------------------
  ! initial conditions
  LOGICAL  :: use_file_initialConditions  = .false.
  REAL(wp) :: initial_temperature_top     = 16.0_wp    ! reference temperature used for initialization in testcase 46
  REAL(wp) :: initial_temperature_bottom  = 16.0_wp    ! reference temperature used for initialization in testcase 46
  REAL(wp) :: initial_temperature_shift   =  0.0_wp    ! temperature used for adding meridional temperature gradient
  REAL(wp) :: initial_temperature_north   =  5.0_wp    ! reference temperature used for initialization Atlantic channel
  REAL(wp) :: initial_temperature_south   =  0.0_wp    ! reference temperature used for initialization Atlantic channel
  REAL(wp) :: initial_temperature_scale_depth = 1000.0_wp ! reference scale depth used for initialization Atlantic channel
  REAL(wp) :: initial_temperature_VerticalGradient = -8.2E-3_wp ! as in Danilov 2012, Ocean Modelling 47
  REAL(wp) :: initial_salinity_top        = 35.0_wp    ! reference salinity used for initialization in testcase 46
  REAL(wp) :: initial_salinity_bottom     = 35.0_wp    ! reference salinity used for initialization in testcase 46
  !  INTEGER  :: scatter_levels(10)       = 0          ! levels for possible scattering of the constant tracer fields
  !  REAL(wp) :: scatter_t                = 20.0_wp    ! temperature value for scattering
  !  REAL(wp) :: scatter_s                = 10.0_wp    ! salinity value for scattering
  INTEGER  :: topography_type             = 0          ! >= 200 analytic bathymetry
  REAL(wp) :: topography_height_reference = 0.0_wp     ! used if topography_type >= 200
  INTEGER  :: sea_surface_height_type     = 0          ! >= 200 sea_surface_height
  INTEGER  :: initial_salinity_type       = 0
  INTEGER  :: initial_temperature_type    = 0
  CHARACTER(LEN= max_char_length) :: initial_sst_type = 'sst1'
  INTEGER  :: initial_velocity_type       = 0
  REAL(wp) :: initial_velocity_amplitude  = 0.0_wp
  CHARACTER(filename_max) :: InitialState_InputFileName = "initial_state.nc"   !< file name for reading in
  REAL(wp) :: smooth_initial_height_weights(2)  = 0.0_wp   ! if > 0, initial height is smoothed by these weights, 1st=this, 2nd=neigbors
  REAL(wp) :: smooth_initial_salinity_weights(2)  = 0.0_wp   ! if > 0, initial height is smoothed by these weights, 1st=this, 2nd=neigbors
  REAL(wp) :: smooth_initial_temperature_weights(2)  = 0.0_wp   ! if > 0, initial height is smoothed by these weights, 1st=this, 2nd=neigbors
  INTEGER  :: smooth_initial_salinity_iterations = 0
  INTEGER  :: smooth_initial_temperature_iterations = 0
  INTEGER  :: smooth_initial_height_iterations = 0
  REAL(wp) :: initial_perturbation_waveNumber = 2.0_wp 
  REAL(wp) :: initial_perturbation_max_ratio  = 0.05_wp 
  LOGICAL  :: initialize_fromRestart = .false.

  ! test cases for ocean model; for the index see run scripts
  INTEGER            :: itestcase_oce  = 0
  NAMELIST/ocean_initialConditions_nml/ &
    & use_file_initialConditions , &
    & initial_temperature_bottom , &
    & initial_temperature_top    , &
    & initial_temperature_shift  , &
    & initial_temperature_north  , &
    & initial_temperature_south  , &
    & initial_temperature_VerticalGradient, &
    & initial_salinity_top       , &
    & initial_salinity_bottom    , &
    & initial_salinity_type      , &
    & initial_temperature_type   , &
    & initial_sst_type           , &
    & initial_velocity_type      , &
    & initial_velocity_amplitude , &
    & topography_type            , &
    & topography_height_reference, &
    & sea_surface_height_type    , &
    & InitialState_InputFileName , &
    & smooth_initial_height_weights, &
    & smooth_initial_height_iterations, &
    & smooth_initial_salinity_weights, &
    & smooth_initial_salinity_iterations, &
    & smooth_initial_temperature_weights, &
    & smooth_initial_temperature_iterations, &
    & initial_temperature_scale_depth, &
    & initial_perturbation_waveNumber, & 
    & initial_perturbation_max_ratio,  &
    & initialize_fromRestart
  !----------------------------------------------------------------------------
  ! vertex list of throughflows
  INTEGER :: denmark_strait(100)         = -1
  INTEGER :: gibraltar(100)              = -1
  INTEGER :: drake_passage(100)          = -1
  INTEGER :: florida_strait(100)          = -1
  INTEGER :: indonesian_throughflow(100) = -1
  INTEGER :: scotland_iceland(100)       = -1
  INTEGER :: mozambique(100)             = -1
  INTEGER :: framStrait(100)             = -1
  INTEGER :: beringStrait(100)           = -1
  INTEGER :: barentsOpening(100)         = -1
  INTEGER :: agulhas(100)                = -1
  INTEGER :: agulhas_long(100)           = -1
  INTEGER :: agulhas_longer(100)         = -1
  LOGICAL :: diagnose_for_horizontalVelocity = .false.
  
  NAMELIST/ocean_diagnostics_nml/ diagnostics_level, &
    & florida_strait, &
    & denmark_strait, &
    & drake_passage, &
    & gibraltar,  &
    & indonesian_throughflow, &
    & scotland_iceland, &
    & mozambique, &
    & framStrait, &
    & beringStrait, &
    & barentsOpening, &
    & agulhas, &
    & agulhas_long, &
    & agulhas_longer, &
    & diagnose_for_horizontalVelocity
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

  ! namelist diagnostics
  LOGICAL :: use_omip_windstress, use_omip_fluxes, use_omip_forcing

  CONTAINS

 !-------------------------------------------------------------------------

 !-------------------------------------------------------------------------
 !>
 !! Initialization of variables that set up the configuration of the ocean model
 !!
 !!               Initialization of variables that set up the configuration
 !!               of the ocean using values read from
 !!               namelist 'ocean_nml' and 'octst_nml'.
!<Optimize:inUse>
 SUBROUTINE read_ocean_namelist( filename )

  CHARACTER(LEN=*), INTENT(IN) :: filename

  ! LOGICAL  :: ignore_land_points = .false.

  ! NAMELIST/ocean_run_nml/ ignore_land_points

  INTEGER :: i_status, istat
  INTEGER :: iunit

  CHARACTER(*), PARAMETER :: &
          method_name = 'mo_ocean_nml/read_ocean_namelist:'

  CALL message(method_name,'running the hydrostatic ocean model')

    !------------------------------------------------------------
    ! 4.0 set up the default values for ocean_nml
    !------------------------------------------------------------

    ! default values when namelist is not present and no default on definition

    n_zlev            = -1 ! 5
    dzlev_m(:)        = -1.0_wp

    !dzlev_m(1:n_zlev) =  (/ 50.0_wp, 150.0_wp, 500.0_wp, 1300.0_wp, 2500.0_wp  /)
    !  lower level of layers:  50       200       700       2000       4500
    !  surface coord. levels:  25       125       450       1350       3250

    ! maximal diffusion coefficient for tracer used in implicit vertical tracer diffusion,
    !   if stability criterion is met
    tracer_convection_MixingCoefficient  = 100.0_wp * Temperature_VerticalDiffusion_background

    !------------------------------------------------------------
    ! 5.0 Read ocean_nml namelist
    !------------------------------------------------------------
    ! (done so far by all MPI processes)

    CALL open_nml(TRIM(filename))
    !==================================================================
    ! NOTE: DO NOT USE STATUS FLAG in READ(nnml) WITHOUT CHECKING IT  !
    ! This will result undetected unread namelists                    !
    !==================================================================

    CALL position_nml ('ocean_dynamics_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_dynamics_nml) ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_dynamics_nml)                         ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_dynamics_nml) ! write settings to temporary text file
      END IF
    END SELECT

    ! physics
    ! this is a stupid way to see if we read the OceanReferenceDensity from the namelist
    OceanReferenceDensity = -1.0_wp
    CALL position_nml ('ocean_physics_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_physics_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_physics_nml)                            ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_physics_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    IF (OceanReferenceDensity == -1.0_wp) THEN ! replace by default value
      SELECT CASE(eos_type)
      CASE(3)
        OceanReferenceDensity = 1035.0_wp
      CASE(10)
        OceanReferenceDensity = 1035.0_wp
        CALL warning("ocean_physics_nml","EOS10 requires absolute salinity and conservetive temperature")
      CASE default
        OceanReferenceDensity = rho_ref
      END SELECT
    ENDIF
    OceanReferenceDensity_inv = 1.0_wp/OceanReferenceDensity
    ReferencePressureIndbars = OceanReferenceDensity * grav * sitodbar

    CALL position_nml ('ocean_horizontal_diffusion_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_horizontal_diffusion_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_horizontal_diffusion_nml)                            ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_horizontal_diffusion_nml)    ! write settings to temporary text file
      END IF
    END SELECT

    CALL position_nml ('ocean_vertical_diffusion_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_vertical_diffusion_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_vertical_diffusion_nml)                            ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_vertical_diffusion_nml)    ! write settings to temporary text file
      END IF
    END SELECT

    CALL position_nml ('ocean_GentMcWilliamsRedi_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_GentMcWilliamsRedi_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_GentMcWilliamsRedi_nml)                            ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_GentMcWilliamsRedi_nml)    ! write settings to temporary text file
      END IF
!     CASE default
!       call finish("","ocean_GentMcWilliamsRedi_nml not positioned")
    END SELECT

    CALL position_nml ('ocean_tracer_transport_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_tracer_transport_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_tracer_transport_nml)                            ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_tracer_transport_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    
    CALL position_nml ('ocean_forcing_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ocean_forcing_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (i_status)
    CASE (positioned)
      READ (nnml, ocean_forcing_nml)                          ! overwrite default settings
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
      READ (nnml, ocean_initialConditions_nml)                          ! overwrite default settings
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
      READ (nnml, ocean_diagnostics_nml)                           ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ocean_diagnostics_nml)   ! write settings to temporary text file
      END IF
    END SELECT

    !------------------------------------------------------------
    ! 6.0 check the consistency of the parameters
    !------------------------------------------------------------
    If (laplacian_form /= 1 .and. laplacian_form /= 2) THEN
      CALL finish(method_name, 'wrong laplacian_form parameter')
    ENDIF
      

    IF( n_zlev < 1 ) &
      & CALL finish(method_name,  'n_zlev < 1')
    IF( iswm_oce == 1 .AND. n_zlev > 1 ) THEN
      CALL message(method_name,'WARNING, shallow water model (ocean): n_zlev set to 1')
      n_zlev = 1
    ENDIF

    IF(discretization_scheme == 1)THEN
      CALL message(method_name,'You have choosen the mimetic dicretization')
    !ELSEIF(discretization_scheme == 2)THEN
    !  CALL message(method_name,'You have choosen the RBF dicretization')
    ELSE
      CALL finish(method_name, 'wrong parameter for discretization scheme')
    ENDIF

    !consistency check for horizontal advection in edge_based configuration
    IF(l_edge_based)THEN
      CALL message(TRIM(method_name),'You are using the EDGE_BASED discretization')
      IF( flux_calculation_horz > horz_flux_twisted_vec_recon .OR. flux_calculation_horz <upwind ) THEN
        CALL finish(TRIM(method_name), 'wrong parameter for horizontal advection scheme; use 1-5')
      ENDIF
      !the fct case requires suitable choices of high- and low order fluxes and of limiter
      IF( flux_calculation_horz == horz_flux_twisted_vec_recon) THEN
        !high and low order flux check
        IF(fct_low_order_flux/=upwind)THEN
          CALL finish(method_name, 'wrong parameter for low order advection scheme in horizontal fct')
        ENDIF
        IF(fct_high_order_flux/= central.AND.fct_high_order_flux/=lax_friedrichs.AND.fct_high_order_flux/=miura_order1)THEN
          CALL finish(method_name, 'wrong parameter for high order advection scheme in horizontal fct')
        ENDIF
        !limiter check
        IF(      fct_limiter_horz/=fct_limiter_horz_zalesak&
          &.AND.fct_limiter_horz/=fct_limiter_horz_minmod &
          &.AND.fct_limiter_horz/=fct_limiter_horz_posdef)THEN
          CALL finish(method_name, 'wrong parameter for limiter in horizontal fct')
        ENDIF
    
      ENDIF
    !consistency check for horizontal advection in cell_based configuration       
    ELSEIF(.NOT.l_edge_based)THEN
      CALL message(TRIM(method_name),'You are using the CELL_BASED discretization')
      IF( flux_calculation_horz > horz_flux_twisted_vec_recon&
       & .OR. flux_calculation_horz <upwind.OR.flux_calculation_horz==lax_friedrichs ) THEN
        CALL finish(TRIM(method_name), 'wrong parameter for horizontal advection scheme; use 1-5 without 3')
      ENDIF     
      IF( flux_calculation_horz == horz_flux_twisted_vec_recon) THEN
        !high and low order flux check
        IF(fct_low_order_flux/=upwind .AND. fct_low_order_flux/=miura_order1)THEN
          CALL finish(TRIM(method_name), 'wrong parameter for low order advection scheme in horizontal fct')
        ENDIF
        !there is no option for high- or low order fluxes in cell_based config, this is all prescribed.
        !a wrong option has no effect.
        !limiter check
        IF(     fct_limiter_horz/=fct_limiter_horz_zalesak &
          &.AND.fct_limiter_horz/=fct_limiter_horz_minmod &
          &.AND.fct_limiter_horz/=fct_limiter_horz_posdef)THEN
          CALL finish(method_name, 'wrong parameter for limiter in horizontal fct')
        ENDIF     
      ENDIF     
    ENDIF
    
    !check for vertical advection
    IF(      flux_calculation_vert/=upwind        &
      &.AND.flux_calculation_vert/=fct_vert_ppm  &
      &.AND.flux_calculation_vert/=fct_vert_adpo &
      &.AND.flux_calculation_vert/=fct_vert_zalesak&
      &.AND.flux_calculation_vert/=fct_vert_minmod)THEN
      CALL finish(method_name, 'wrong parameter for vertical advection')
    ENDIF

    IF(i_bc_veloc_lateral/= 0) THEN
      CALL finish(method_name, &
        &  'free-slip boundary condition for velocity currently not supported')
    ENDIF
    IF(i_bc_veloc_top < 0 .OR. (i_bc_veloc_top > 1 .and. i_bc_veloc_top /= 4)) THEN
    !  option >1 disabled due to unphysical difference of stress minus velocity
    !  see method_name top_bound_cond_horz_veloc (#slo#, 2014-04)
      CALL finish(method_name, &
        &  'top boundary condition for velocity currently not supported: choose = 0,1')
    ENDIF
    IF(i_bc_veloc_bot < 0 .OR. i_bc_veloc_bot>2) THEN
      CALL finish(method_name, &
        &  'bottom boundary condition for velocity currently not supported: choose = 0, 1, 2')
    ENDIF

!      IF(no_tracer == 1 .OR. no_tracer < 0 .OR. no_tracer > 2) THEN
!        IF(no_tracer == 1) THEN
!          CALL message(method_name, 'WARNING - You have chosen tracer temperature only')
!          CALL message(method_name, ' - this generates error in mo_varlist/mo_ocean_state')
!          CALL finish(method_name,  'no_tracer=1 not supported - choose =0 or =2')
!        ENDIF
!        CALL finish(method_name,  'no_tracer not supported - choose =0 or =2')
!      ENDIF


!     IF (solver_start_tolerance <= 0.0_wp) THEN
!       solver_start_tolerance = solver_tolerance
!       solver_tolerance_decrease_ratio  = 0.1_wp ! must be < 1
!     ENDIF

    !IF (forcing_enable_freshwater) THEN
    !  !limit_elevation = .TRUE.
    !  CALL message(method_name,'WARNING, limit_elevation set to .TRUE. with forcing_enable_freshwater=.TRUE.')
    !END IF

    IF (forcing_set_runoff_to_zero) THEN
      CALL message(method_name,'WARNING, forcing_set_runoff_to_zero is .TRUE. - forcing with river runoff is set to zero')
    END IF

#ifndef __NO_ICON_ATMO__
    IF ( is_coupled_run() ) THEN
      iforc_oce = Coupled_FluxFromAtmo
      CALL message(method_name,'WARNING, iforc_oce set to 14 for coupled experiment')
 !!!  limiters can now be set by namelist
 !!!  limit_elevation = .FALSE.
 !!!  CALL message(method_name,'WARNING, limit_elevation set to .FALSE. for coupled experiment')
 !!!  limit_seaice = .FALSE.
 !!!  CALL message(method_name,'WARNING, limit_seaice set to .FALSE. - no limit for coupled experiment')
    END IF
#endif

    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) THEN
      WRITE(nnml_output,nml=ocean_dynamics_nml)
      WRITE(nnml_output,nml=ocean_physics_nml) 
      WRITE(nnml_output,nml=ocean_horizontal_diffusion_nml)
      WRITE(nnml_output,nml=ocean_vertical_diffusion_nml)
      WRITE(nnml_output,nml=ocean_tracer_transport_nml)
      WRITE(nnml_output,nml=ocean_GentMcWilliamsRedi_nml)
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
      &  'all', &  ! initiate print messages in all method_names
      &  'abm', &  ! main timestepping method_names       in mo_ocean_ab_timestepping (mimetic/rbf)
      &  'vel', &  ! velocity advection and diffusion in mo_ocean_velocity_advection
      &  'dif', &  ! diffusion                        in mo_ocean_diffusion
      &  'trc', &  ! tracer advection and diffusion   in mo_ocean_tracer_transport
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

    tracer_threshold_min(1) = threshold_min_T
    tracer_threshold_max(1) = threshold_max_T
    tracer_threshold_min(2) = threshold_min_S
    tracer_threshold_max(2) = threshold_max_S
    namelist_tracer_name(1) = "Temperature"
    namelist_tracer_name(2) = "Salinity"
     
    use_omip_windstress = ( forcing_windstress_u_type == 1 ) .AND. (forcing_windstress_v_type == 1)
    use_omip_fluxes     = ( forcing_fluxes_type == 1 )
    use_omip_forcing    = use_omip_windstress .OR. use_omip_fluxes

    nbgctra = 0
    nbgcadv = 0
    if(lhamocc) then 
       nbgctra = n_bgctra
       if(lbgcadv) nbgcadv =  ntraad
    endif

END SUBROUTINE read_ocean_namelist

END MODULE mo_ocean_nml
