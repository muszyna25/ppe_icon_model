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
  ! ------------------------------------------------------------------------

  INTEGER  :: n_zlev        ! number of ocean levels
  REAL(wp) :: dzlev_m(100)  ! namelist input of layer thickness

  INTEGER, PARAMETER :: toplev    = 1   ! surface ocean level

  ! parameterized forcing for ocean model:
  INTEGER            :: iforc_oce                 =  0   ! index of parameterized forcing
  INTEGER, PARAMETER :: NO_FORCING                = 10
  INTEGER, PARAMETER :: ANALYT_FORC               = 11
  INTEGER,PARAMETER  :: FORCING_FROM_FILE_FLUX    = 12
  INTEGER,PARAMETER  :: FORCING_FROM_FILE_FIELD   = 13
  INTEGER, PARAMETER :: FORCING_FROM_COUPLED_FLUX = 14
  INTEGER, PARAMETER :: FORCING_FROM_COUPLED_FIELD= 15

  ! type of time varying OMIP forcing: 1: u/v-stress/SST; 2: full OMIP for sea-ice
  INTEGER            :: iforc_omip     = 0

  ! length of time varying flux forcing: 1: read 12 months, 2: read 365 steps
  INTEGER            :: iforc_len      = 1

  ! stationary forcing for ocean model:
  INTEGER            :: iforc_stat_oce = 0

  ! switch for reading prognostic variables:
  INTEGER            :: init_oce_prog  = 0

  ! test cases for ocean model; for the index see run scripts
  INTEGER            :: itestcase_oce  = 0

  ! switch for ocean diagnostics - 0: no diagnostics; 1: write to stderr
  INTEGER            :: idiag_oce      = 0

  ! parameterized velocity boundary conditions
                      ! Velocity boundary condition: Currently only no-slip is supported !!
                      ! i_bc_veloc_lateral = 0: boundary condition for velocity is no-slip: normal
                      !                         and tangential velocity components at lateral 
                      !                         boundaries are set to zero
                      ! i_bc_veloc_lateral = 1: boundary condition for velocity is free-slip: 
                      !                         normal velocity components at lateral boundariea is
                      !                         set to zero, tangential not.
  INTEGER            :: i_bc_veloc_lateral = 0   

  INTEGER            :: i_bc_veloc_top = 0  !Top boundary condition for velocity: 
                                            ! i_bc_veloc_top =0 :zero value at top boundary,no wind
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

  !INTEGER            :: vbc_zero_cond   =   0   ! no or zero boundary condition

  ! parameterized time stepping scheme
  !  #slo# - might later be splitted into stepping scheme and call of core
  ! in mo_test_hydro_ocean:
  !   i_oce_stepping = 1 => call test_expl_core_oce
  !   i_oce_stepping = 2 => call perform_ho_stepping
  !   i_oce_stepping = 2 => call perform_core201006
  INTEGER            :: i_oce_stepping  =   1  ! switch for time stepping scheme


  ! parameterized shallow water mode in the ocean model
  INTEGER            :: iswm_oce        =   0  ! switch for shallow water mode (1 = on, 0 = 3dim)
  INTEGER            :: idisc_scheme    =   0  ! discretization scheme: 1 for mimetic, 
                                               ! 2 for RBF-type of discretization
 
  ! parameters for Adams-Bashforth semi-implicit time stepping scheme
  ! are set according to Marshall et al paper
  REAL(wp) :: ab_const              =  0.1_wp  ! Adams-Bashforth constant
  REAL(wp) :: ab_beta               =  0.0_wp  ! Parameter in semi-implicit timestepping
  REAL(wp) :: ab_gam                =  0.0_wp  ! Parameter in semi-implicit timestepping
  REAL(wp) :: solver_tolerance      =  0.0_wp

  INTEGER :: EOS_TYPE               = 0        ! 1=linear EOS,2=(nonlinear) 
                                               ! Jacket-McDoudgall-formulation
  INTEGER :: no_tracer              = 0        ! number of tracers 

  ! more ocean parameters, not yet well placed
  INTEGER  :: expl_vertical_velocity_diff=0    ! 0=explicit, 1 = implicit  
  INTEGER  :: expl_vertical_tracer_diff  = 0   ! 0=explicit, 1 = implicit
  INTEGER  :: HORZ_VELOC_DIFF_TYPE  = 0
  REAL(wp) :: k_veloc_h             = 0.0_wp   ! horizontal diffusion coefficient
  REAL(wp) :: k_veloc_v             = 0.0_wp   ! vertical diffusion coefficient
  REAL(wp) :: k_pot_temp_h          = 0.0_wp   ! horizontal mixing coefficient for pot. temperature
  REAL(wp) :: k_pot_temp_v          = 0.0_wp   ! vertical mixing coefficient for pot. temperature
  REAL(wp) :: k_sal_h               = 0.0_wp   ! horizontal diffusion coefficient for salinity
  REAL(wp) :: k_sal_v               = 0.0_wp   ! vertical diffusion coefficient for salinity

  REAL(wp) :: t_ref                 = 0.0_wp   ! reference temperature for initialization
  REAL(wp) :: s_ref                 = 0.0_wp   ! reference salinity for initialization
  REAL(wp) :: bottom_drag_coeff     = 0.002_wp ! chezy coefficient for bottom friction
  REAL(wp) :: wstress_coeff         = 1.e-4_wp ! windstress coefficient
  INTEGER  :: temperature_relaxation= 1
  REAL(wp) :: relaxation_param      = 0.0_wp   ! parameter for relaxation of tracer surface fluxes (same value for all tracers)
                                               !this value is divided by number of seconds per month (=30*24*3600)

  INTEGER  :: coriolis_type         = 1        ! 0=zero Coriolis, the non-rotating case
                                               ! 1=full varying Coriolis
                                               ! 2=beta-plane (linear) approximation to Coriolis
                                               ! 3=f-plane (constant) approximation to Coriolis
                                               ! The variables below are used to set up in basin configuration
                                               ! the Coriolis (f/beta-plane) and to adjust the analytic wind forcing
  REAL(wp) :: basin_center_lat     = 0.0_wp    ! lat-lon coordinate of basin center, used 
  REAL(wp) :: basin_center_lon     = 0.0_wp    ! in (non-global) basin configuration such as the Stommel-type tests
  REAL(wp) :: basin_width_deg      = 0.0_wp    ! basin extension in x-direction, units are degrees
  REAL(wp) :: basin_height_deg     = 0.0_wp    ! basin extension in y-direction, units are degrees

  REAL(wp) :: MAX_VERT_DIFF_VELOC = 0.5_wp     !maximal diffusion coeff, used in implicit vertical velocity diffusion, 
                                               !if stability criterion is met
  REAL(wp) :: MAX_VERT_DIFF_TRAC = 0.5_wp      !maximal diffusion coeff, used in implicit vertical tracer diffusion, 
                                               ! if stability criterion is met
  REAL(wp) :: CWA = 0.5E-3_wp                  !Tuning parameters for vertical mixing of tracer and velocity
  REAL(wp) :: CWT = 0.5E-3_wp

                                               
  LOGICAL  :: lviscous              = .TRUE.
  LOGICAL  :: l_RIGID_LID           = .FALSE.  ! include friction or not
  LOGICAL  :: l_inverse_flip_flop   = .FALSE.  ! true=complete discrete scalarproduct (slow)
                                               ! false=use a shortcut (faster)

  NAMELIST/ocean_dynamics_nml/ n_zlev, dzlev_m, idisc_scheme,              &
    &                 iswm_oce, i_oce_stepping,                            &
    &                 i_bc_veloc_lateral,i_bc_veloc_top,i_bc_veloc_bot,    &
    &                 ab_const, ab_beta, ab_gam, solver_tolerance,         &
    &                 l_RIGID_LID, lviscous, l_inverse_flip_flop,          &
    &                 coriolis_type, basin_center_lat, basin_center_lon,   &
    &                 basin_width_deg,basin_height_deg,                    &
    &                 expl_vertical_velocity_diff,                         &
    &                 expl_vertical_tracer_diff 
 


  NAMELIST/ocean_physics_nml/EOS_TYPE, no_tracer, HORZ_VELOC_DIFF_TYPE,    &
    &                 k_veloc_h, k_veloc_v,  k_pot_temp_h, k_pot_temp_v,   &
    &                 k_sal_h, k_sal_v,                                    &
    &                 MAX_VERT_DIFF_VELOC, MAX_VERT_DIFF_TRAC,             &
    &                 CWA, CWT,  bottom_drag_coeff, wstress_coeff       


  NAMELIST/ocean_forcing_and_init_nml/iforc_oce, iforc_omip, iforc_len,    &
    &                 iforc_stat_oce, init_oce_prog, itestcase_oce,        &
    &                 idiag_oce, temperature_relaxation, relaxation_param


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

  ! switches for level of debugging the ocean core
  INTEGER  :: i_dbg_oce = 0       ! different levels of debug output (1-5, 0: no output)
  INTEGER  :: i_dbg_inx = 0       ! different levels of debug output of values at indices given below

  ! longitude/latitude location of single cell output for debugging
  REAL(wp) :: rlat_in   = 0.0_wp  ! latitude of cell for debug output
  REAL(wp) :: rlon_in   = 0.0_wp  ! longitude of cell for debug output

  ! block/index location of cell output for debugging
  INTEGER  :: i_oct_blk = 0       ! output test block
  INTEGER  :: i_oct_idx = 0       ! output test index
  INTEGER  :: i_oct_ilv = 1       ! output test level

  CHARACTER(len=3) :: str_proc_tst(10)   ! namelist string of source processes to print

  NAMELIST/octst_nml/   i_oct_blk, i_oct_idx, i_oct_ilv,     &
    &                   i_ocv_blk, i_ocv_idx, i_ocv_ilv,     &
    &                   i_dbg_oce, i_dbg_inx, str_proc_tst,  &
    &                   h_val, t_val, rlat_in, rlon_in

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

     INTEGER :: i_status

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


     !------------------------------------------------------------
     ! 5.0 Read ocean_nml namelist
     !------------------------------------------------------------
     ! (done so far by all MPI processes)

     CALL open_nml(TRIM(filename))
     CALL position_nml ('ocean_dynamics_nml', status=i_status)
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_dynamics_nml)
     END SELECT

     CALL position_nml ('ocean_physics_nml', status=i_status)
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_physics_nml)
     END SELECT

     CALL position_nml ('ocean_forcing_and_init_nml', status=i_status)
     SELECT CASE (i_status)
     CASE (positioned)
       READ (nnml, ocean_forcing_and_init_nml)
     END SELECT


     !------------------------------------------------------------
     ! 6.0 check the consistency of the parameters
     !------------------------------------------------------------

     IF( iswm_oce == 1 .AND. n_zlev > 1 ) THEN
       CALL message(TRIM(routine),'WARNING, shallow water model (ocean): n_zlev set to 1')
       n_zlev = 1
     ENDIF

     IF(idisc_scheme == 1)THEN
       CALL message(TRIM(routine),'You have choosen the mimetic dicretization')
     ELSEIF(idisc_scheme == 2)THEN
       CALL message(TRIM(routine),'You have choosen the RBF dicretization')
     ELSE
       CALL finish(TRIM(routine), 'wrong parameter for discretization scheme')
     ENDIF
    
    
     IF(i_bc_veloc_lateral/= 0) THEN
       CALL finish(TRIM(routine), &
         &  'free-slip boundary condition for velocity currently not supported')
     ENDIF
     IF(i_bc_veloc_top < 0.AND.i_bc_veloc_top > 2) THEN
       CALL finish(TRIM(routine), &
         &  'top boundary condition for velocity currently not supported: choose = 0 or =1 or =2')
     ENDIF
     IF(i_bc_veloc_bot < 0 .AND. i_bc_veloc_bot>1) THEN
       CALL finish(TRIM(routine), &
         &  'bottom boundary condition for velocity currently not supported: choose = 0 or =1')
     ENDIF
 
     ! write the contents of the namelist to an ASCII file
     IF(my_process_is_stdio()) WRITE(nnml_output,nml=ocean_dynamics_nml)
     IF(my_process_is_stdio()) WRITE(nnml_output,nml=ocean_physics_nml)
     IF(my_process_is_stdio()) WRITE(nnml_output,nml=ocean_forcing_and_init_nml)
     !------------------------------------------------------------
     ! 6.0 Read octst_nml namelist
     !------------------------------------------------------------
     ! (done so far by all MPI processes)

     ! 3-char string with marked processes to be printed out for debug purposes
     str_proc_tst =  (/  & 
       &  'abm', &  ! main timestepping routines       in mo_oce_ab_timestepping (mimetic/rbf)
       &  'vel', &  ! velocity advection and diffusion in mo_oce_veloc_advection
       &  'dif', &  ! diffusion                        in mo_oce_diffusion
       &  'trc', &  ! tracer advection and diffusion   in mo_oce_tracer_transport
       &  '   ', &  ! ...
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
