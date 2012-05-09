!>
!!  Contains the data structures for the hydrostatic ocean model.
!!
!!  Contains the data structures to store the hydrostatic & boussinesq ocean model state.
!!  Implementation is based on ICON-Shallow-Water model
!!  to store the shallow water model state and other auxiliary variables.
!!  Constructors and destructors for these data structures are also defined here.
!!
!! @par Revision History
!!  Initial version by Peter Korn (MPI-M), (2006).
!!  Big recoding by P. Korn (MPI-M), (2009/2010)
!!  Modification by Stephan Lorenz, MPI-M, (2010-03-19):
!!   - renaming and adjustment to ocean domain and patch_oce
!!  Modification by Stephan Lorenz, MPI-M, 2011-07
!!   - 3-dim ocean structures moved from patch_oce to hydro_ocean_base
!
!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_oce_state
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,                ONLY: wp
! USE mo_mpi,                 ONLY: p_pe_work
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    &                               success, max_char_length, min_rledge, min_rlcell,  &
    &                               min_rlvert, min_rlcell_int,                        &
    &                               full_coriolis, beta_plane_coriolis,                &
    &                               f_plane_coriolis, zero_coriolis
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer,                        &
    &                               CORIOLIS_TYPE, basin_center_lat, basin_height_deg
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_config,         ONLY: n_dom
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_math_utilities,      ONLY: gc2cc,t_cartesian_coordinates,      &
    &                               t_geographical_coordinates, &!vector_product, &
    &                               arc_length
  USE mo_math_constants,      ONLY: deg2rad,rad2deg
  USE mo_physical_constants,  ONLY: re, omega, rho_ref
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, sync_patch_array, global_sum_array
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: add_var,                  &
    &                               new_var_list,             &
    &                               delete_var_list,          &
    &                               default_var_list_settings,&
    &                               add_ref
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range


  IMPLICIT NONE
  PRIVATE

! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !public interface
  !
  ! subroutines
  PUBLIC :: construct_hydro_ocean_base
  PUBLIC :: destruct_hydro_ocean_base
  PUBLIC :: construct_hydro_ocean_state
  PUBLIC :: destruct_hydro_ocean_state
  PUBLIC :: set_lateral_boundary_values
  PUBLIC :: init_ho_base
  PUBLIC :: init_ho_basins
  PUBLIC :: init_coriolis_oce
  PUBLIC :: set_del_zlev, set_zlev
  PUBLIC :: is_initial_timestep
  PUBLIC :: init_oce_config

  !
  ! types
  PUBLIC :: t_hydro_ocean_base
  PUBLIC :: t_hydro_ocean_state
  PUBLIC :: t_hydro_ocean_prog
  PUBLIC :: t_hydro_ocean_aux
  PUBLIC :: t_hydro_ocean_diag
  PUBLIC :: t_ptr3d
  PUBLIC :: t_oce_config

  !
  !constructors
  PRIVATE :: construct_hydro_ocean_diag
  PRIVATE :: construct_hydro_ocean_prog
  PRIVATE :: construct_hydro_ocean_aux
  !destructors
  PRIVATE :: destruct_hydro_ocean_diag
  PRIVATE :: destruct_hydro_ocean_aux
  !INTEGER, PRIVATE :: i_cell_type=3

!
!! basis types for constructing 3-dim ocean state
!
  TYPE t_hydro_ocean_base

    !! The ocean uses z-coordinates in meters in the vertical.
    !! The following data are required:
    !!
    !! n_zlev: number of z-coordinate surfaces
    !! n_zlvp: number of intermediate levels (+1)
    !! n_zlvm: number of z-coordinate distances (-1)
    INTEGER :: n_zlev, n_zlvp, n_zlvm

    !! del_zlev_m: thickness (height) of elemental prism, defined as the 
    !!             distance between top and bottom of elemental prism, 
    !!             i.e. the distance between two intermediate z-coordinate 
    !!             surfaces. These data are provided by the user, all other 
    !!             vertical information is calculated from this array of 
    !!             thicknesses.
    !!             Dimension: n_zlev
    REAL(wp), ALLOCATABLE :: del_zlev_m(:)

    !! zlev_m    : position of the vertical cell centers, i.e. below zero surface;
    !!             Numbering starts from surface and increases downwards to bottom.
    !!             Dimension: n_zlev
    !!             At these surfaces the horizontal velocities, vorticity, divergence
    !!             and scalar variables are evaluated.
    REAL(wp), ALLOCATABLE :: zlev_m(:)


    !! zlev_i    : vertical position of the UPPER BORDER of the vertical cell
    !!             i.e. the position of the top of elemental prisms.
    !!             Position of first surface is 0.
    !!             Dimension: n_zlvp = n_zlev + 1
    !!             The vertical velocities are evaluated at such surfaces.
    REAL(wp), ALLOCATABLE :: zlev_i(:)

    !! del_zlev_i: distance between two z-coordinate surfaces. The first is 
    !!             the distance from the ocean surface = zlev_m(1)
    !!             Dimension: n_zlev
    REAL(wp), ALLOCATABLE :: del_zlev_i(:)


    ! land-sea-mask for ocean has 3 dimensions (the 2nd is the number of 
    ! vertical levels)
    ! sea=-2, sea_boundary=-1, boundary (edges only)=0, land_boundary=1, land=2
    !
    ! land-sea-mask for cell centers
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_c
    INTEGER, ALLOCATABLE :: lsm_oce_c(:,:,:)
    ! land-sea-mask for cell edges
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e
    INTEGER, ALLOCATABLE :: lsm_oce_e(:,:,:)
    ! land-sea-mask for cell vertices
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_v
    ! INTEGER, ALLOCATABLE :: lsm_oce_v(:,:,:)


    ! To simplify the acess to the required information within these loops 
    ! we store an cell and edge based version of the deepest ocean layer 
    ! in column. dolic_e(edge1) and dolic_c(cell1) are identical if 'edge1' 
    ! is one of the edges of 'cell1'.
    ! If the ocean bottom is flat dolic_c and dolic_e are identical and equal 
    ! to the number of z-coodinate surfaces.

    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: dolic_c(:,:)
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: dolic_e(:,:)

    ! For diagnosis like stream functions and area calculations we add surface arrays
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER,  ALLOCATABLE :: basin_c(:,:)  ! basin information Atlantic/Indian/Pacific
    INTEGER,  ALLOCATABLE :: regio_c(:,:)  ! area information like tropical Atlantic etc.
    REAL(wp), ALLOCATABLE :: rbasin_c(:,:) ! real for output
    REAL(wp), ALLOCATABLE :: rregio_c(:,:) ! real for output

    ! To simply set land points to zero we store additional 3-dim wet points
    ! dimensions as in lsm_oce:
    REAL(wp), ALLOCATABLE :: wet_c(:,:,:)  ! cell centers
    REAL(wp), ALLOCATABLE :: wet_e(:,:,:)  ! cell edges
    !REAL(wp), ALLOCATABLE :: wet_i(:,:,:)  ! vertical velocity points 
    !                                       ! on intermediate levels


!!$    ! Arrays that describe vertical connectivity of the triangles.
!!$    ! The indices of triangles are stored in a whole vertical column.
!!$    ! index1=1,nproma, index2=1,nblks_c, index3=1,n_zlev
!!$    INTEGER, ALLOCATABLE :: neighbor_c(:,:,:)
!!$    ! index1=1,nproma, index2=1,nblks_e, index3=1,n_zlev
!!$    INTEGER, ALLOCATABLE :: neighbor_e(:,:,:)
  END TYPE t_hydro_ocean_base

  TYPE t_ptr3d
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr3d
!
!! prognostic variables
!
  TYPE t_hydro_ocean_prog

    REAL(wp), POINTER ::    &
      &  h(:,:)                ,& ! height of the free surface. Unit: [m]
                                  ! dimension:(nproma, nblks_c)
      &  vn(:,:,:)             ,& ! velocity component normal to cell edge. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  tracer(:,:,:,:)          ! tracer concentration.
                                  ! dimension: (nproma, n_zlev, nblks_c, no_tracer)
                                  ! Ordering of tracers:
                                  !   1) pot_temp:= potential temperature, Unit: [deg C]
                                  !   2) salinity:= salinity, Unit [psu]

    TYPE(t_ptr3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
  END TYPE t_hydro_ocean_prog

!
!! diagnostic variables
!
  TYPE t_hydro_ocean_diag

    REAL(wp), POINTER ::        &
      &  vt(:,:,:)             ,& ! tangential velocity component at edges. Unit [m/s].
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  rho(:,:,:)            ,& ! density. Unit: [kg/m^3]
                                  ! dimension: (nproma,n_zlev, nblks_c)
      &  h_e(:,:)              ,& ! surface height at cell edges. Unit [m].
                                  ! dimension: (nproma, nblks_e)
      &  thick_c(:,:)          ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma, nblks_c)
      &  thick_e(:,:)          ,& ! individual fluid column thickness at edges. Unit [m].
                                  ! dimension: (nproma, nblks_e)
      &  mass_flx_e(:,:,:)     ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  div_mass_flx_c(:,:,:) ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma,n_zlev, nblks_c)
      &  cons_thick_c(:,:,:)   ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma,n_zlev, nblks_c)
      &  depth_c(:,:,:)        ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma,n_zlev, nblks_c)
      &  w(:,:,:)              ,& ! vertical velocity. Unit [m/s].
                                  ! dimension: (nproma, n_zlev+1, nblks_c)
      &  w_old(:,:,:)          ,& ! vertical velocity from previous timestep. Unit [m/s].
                                  ! dimension: (nproma, n_zlev+1, nblks_c)
      &  w_e(:,:,:)            ,& ! vertical velocity at edges. Unit [m/s]
                                  ! dimension: (nproma, n_zlev+1, nblks_e)
      &  wtemp(:,:,:)              ,& ! vertical velocity. Unit [m/s].
      &  w_prev(:,:,:)         ,& ! vertical velocity at cells, from previous timestep. Unit [m/s]
                                  ! dimension: (nproma, n_zlev+1, nblks_c)
      &  u(:,:,:)              ,& ! reconstructed zonal velocity component. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  v(:,:,:)              ,& ! reconstructed meridional velocity component. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  psi(:,:)              ,& ! barotropic stream function. Unit [kg/s]
                                  ! dimension: (nproma, nblks_c)
      &  ptp_vn(:,:,:)         ,& ! normal velocity after mapping P^T P
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_pred(:,:,:)        ,& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_impl_vert_diff(:,:,:),& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_time_weighted(:,:,:),&  ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  w_time_weighted(:,:,:),& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vort(:,:,:)           ,& ! vorticity at triangle vertices. Unit [1/s]
                                  ! dimension: (nproma, n_zlev, nblks_v)
      &  vort_e(:,:,:)         ,& ! vorticity interpolated to triangle edges. Unit [1/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  kin(:,:,:)            ,& ! kinetic energy. Unit [m/s].
                                  ! (nproma, n_zlev, nblks_c)
      &  veloc_adv_horz(:,:,:) ,& ! horizontal velocity advection
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  veloc_adv_vert(:,:,:) ,& ! vertical velocity advection
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  laplacian_horz(:,:,:) ,& ! horizontal diffusion of horizontal velocity
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  laplacian_vert(:,:,:) ,& ! vertical diffusion of horizontal velocity
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  grad(:,:,:)           ,& ! gradient of kinetic energy. Unit [m/s]
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  div(:,:,:)            ,& ! divergence. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  press_hyd(:,:,:)      ,& ! hydrostatic pressure. Unit [m]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  press_grad(:,:,:)     ,& ! hydrostatic pressure gradient term. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      & temp_insitu(:,:,:)
    TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_vn(:,:,:)              ! reconstructed velocity at cell center in cartesian coordinates
                                  ! dimension: (nproma, n_zlev, nblks_c)
    TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_vn_dual(:,:,:)         ! reconstructed velocity at vertex in cartesian coordinates
                                  ! dimension: (nproma, n_zlev, nblks_v)
    TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_vn_mean(:,:,:)         ! reconstructed velocity at vertex in cartesian coordinates
                                  ! dimension: (nproma, n_zlev, nblks_v). For mimetic miura scheme

  END TYPE t_hydro_ocean_diag

!
!! auxiliary data
!
  TYPE t_hydro_ocean_aux

    REAL(wp), POINTER ::       &
      &  g_n(:,:,:)           ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                 ! at timelevel n
                                 ! dimension: (nproma, n_zlev, nblks_e)
      &  g_nm1(:,:,:)         ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                 ! at timelevel n-1
                                 ! dimension: (nproma, n_zlev, nblks_e)
      &  g_nimd(:,:,:)           ! explicit velocity term in Adams-Bashford time marching routines,
                                 ! located at intermediate timelevel

      ! Variables for each tracer
    REAL(wp), POINTER ::       &
      &  g_n_c_h(:,:,:,:)        ! explicit tracer term in Adams-Bashford time marching routines,
                                 ! at timelevel n for each tracer, horizontal
                                 ! dimension: (nproma, n_zlev, nblks_c, no_tracer )
    TYPE(t_ptr3d),ALLOCATABLE :: g_n_c_h_tracer_ptr(:)   !< pointer array: one pointer for each tracer

    REAL(wp), POINTER ::       &
      &  g_nm1_c_h(:,:,:,:)      ! explicit tracer term in Adams-Bashford time marching routines,
                                 ! at timelevel n-1 for each tracer, horizontal
                                 ! dimension: (nproma, n_zlev, nblks_c, no_tracer)
    TYPE(t_ptr3d),ALLOCATABLE :: g_nm1_c_h_tracer_ptr(:)   !< pointer array: one pointer for each tracer

    REAL(wp), POINTER ::       &
      &  g_nimd_c_h(:,:,:,:)     ! explicit tracer term in Adams-Bashford time marching routines,
                                 ! located at intermediate timelevel for each tracer, horizontal
                                 ! dimension: (nproma, n_zlev, nblks_c,no_tracer )
    TYPE(t_ptr3d),ALLOCATABLE :: g_nimd_c_h_tracer_ptr(:)   !< pointer array: one pointer for each tracer

    REAL(wp), POINTER ::       &
      &  g_n_c_v(:,:,:,:)        ! explicit tracer term in Adams-Bashford time marching routines,
                                 ! at timelevel n for each tracer, vertical
                                 ! dimension: (nproma, n_zlev, nblks_c, no_tracer )
    TYPE(t_ptr3d),ALLOCATABLE :: g_n_c_v_tracer_ptr(:)   !< pointer array: one pointer for each tracer

    REAL(wp), POINTER ::       &
      &  g_nm1_c_v(:,:,:,:)      ! explicit tracer term in Adams-Bashford time marching routines,
                                 ! at timelevel n-1 for each tracer, vertical
                                 ! dimension: (nproma, n_zlev, nblks_c, no_tracer)
    TYPE(t_ptr3d),ALLOCATABLE :: g_nm1_c_v_tracer_ptr(:)   !< pointer array: one pointer for each tracer

    REAL(wp), POINTER ::       &
      &  g_nimd_c_v(:,:,:,:)     ! explicit tracer term in Adams-Bashford time marching routines,
                                 ! located at intermediate timelevel for each tracer, vertical
                                 ! dimension: (nproma, n_zlev, nblks_c,no_tracer )
    TYPE(t_ptr3d),ALLOCATABLE :: g_nimd_c_v_tracer_ptr(:)  !< pointer array: one pointer for each tracer

    REAL(wp), POINTER ::       &
      &  bc_top_vn(:,:)       ,& ! normal velocity boundary condition at surface
                                 ! dimension: (nproma,nblks_e)
      &  bc_bot_vn(:,:)       ,& ! normal velocity boundary condition at bottom
                                 ! dimension: (nproma,nblks_c)
      &  bc_top_u(:,:)        ,& ! zonal velocity boundary condition at surface
                                 ! dimension: (nproma,nblks_c)
      &  bc_top_v(:,:)        ,& ! meridional velocity boundary condition at surface
                                 ! dimension: (nproma,nblks_c)
      &  bc_bot_u(:,:)        ,& ! zonal velocity boundary condition at bottom
                                 ! dimension: (nproma,nblks_c)
      &  bc_bot_v(:,:)        ,& ! meridional velocity boundary condition at bottom
                                 ! dimension: (nproma,nblks_c)
      &  bc_top_w(:,:)        ,& ! vertical velocity boundary condition at surface
                                 ! dimension: (nproma,nblks_c)
      &  bc_bot_w(:,:)        ,& ! vertical velocity boundary condition at bottom
      &  bc_top_tracer(:,:,:) ,& ! vertical velocity boundary condition at surface
                                 ! dimension: (nproma,nblks_c)
      &  bc_bot_tracer(:,:,:) ,& ! vertical velocity boundary condition at bottom
      &  p_rhs_sfc_eq(:,:)!,   & ! right hand side of surface equation
                                         ! dimension: (nproma,nblks_c)
     TYPE(t_cartesian_coordinates), POINTER :: bc_top_veloc_cc(:,:), &
                                  &                bc_bot_veloc_cc(:,:)
     TYPE(t_ptr3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

    ! Variables for 3-dim tracer relaxation:
    REAL(wp), POINTER ::         &
      &  relax_3d_data_T(:,:,:), & ! 3-dim temperature relaxation data (T*)
                                   ! dimension: (nproma,n_zlev,nblks_c)
      &  relax_3d_forc_T(:,:,:), & ! 3-dim temperature relaxation forcing (1/tau*(T-T*))
                                   ! dimension: (nproma,n_zlev,nblks_c)
      &  relax_3d_data_S(:,:,:), & ! 3-dim salinity relaxation data (T*)
                                   ! dimension: (nproma,n_zlev,nblks_c)
      &  relax_3d_forc_S(:,:,:)    ! 3-dim salinity relaxation forcing (1/tau*(T-T*))
                                   ! dimension: (nproma,n_zlev,nblks_c)

  END TYPE t_hydro_ocean_aux

!
!! array of states
!
  TYPE t_hydro_ocean_state

    TYPE(t_hydro_ocean_prog), POINTER :: p_prog(:)    ! time array of prognostic states at
                                                        ! different time levels
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_aux)  :: p_aux

  END TYPE t_hydro_ocean_state

  INTEGER, PARAMETER               :: max_tracers = 2
  TYPE t_oce_config
    CHARACTER(len=max_char_length) :: tracer_names(max_tracers)
    CHARACTER(len=max_char_length) :: tracer_longnames(max_tracers)
    CHARACTER(len=max_char_length) :: tracer_units(max_tracers)
    CHARACTER(len=max_char_length) :: tracer_tags(max_tracers)
    INTEGER                        :: tracer_codes(max_tracers)
  END TYPE t_oce_config

  ! variables
  TYPE(t_var_list)         , PUBLIC                      :: ocean_var_list
  TYPE(t_hydro_ocean_state), PUBLIC, TARGET, ALLOCATABLE :: v_ocean_state(:)
  TYPE(t_hydro_ocean_base) , PUBLIC, TARGET              :: v_base
  TYPE(t_oce_config)       , PUBLIC                      :: oce_config

!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
!

!>
!! Constructor for hydrostatic ocean state + diagnostic and auxiliary  states.
!!
!! Constructor for hydrostatic ocean state
!! It calls  constructors to single time level
!! auxiliary and diagnostic states. Then it constructs state array,
!! whose components (representing multiple time levels).
!! Initialization of all components with zero.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2007).
!!  Modification by Stephan Lorenz, MPI-M, (2010-06-01) - no temporary memory array
!
!
  SUBROUTINE construct_hydro_ocean_state( p_patch, p_os )

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)
    TYPE(t_hydro_ocean_state), TARGET :: p_os(n_dom)

    !INTEGER, OPTIONAL, INTENT(IN)             :: prog_length
    !INTEGER, OPTIONAL, INTENT(IN)             :: k_no_temp_mem

    ! local variables

    INTEGER           :: jg

    INTEGER           :: i_status, jp, prlength ! local prognostic array length
    !INTEGER           :: no_temp_memory         ! no of temporary memory elements
    CHARACTER(len=max_char_length) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_state'

    CALL message(TRIM(routine), 'start to construct hydro_ocean state' )

    !IF (PRESENT(prog_length)) THEN
    !  prlength = prog_length
    !ELSE
    !  prlength =2
    !END IF

    !
    ! Using Adams-Bashforth semi-implicit timestepping with 3 prognostic time levels:
    prlength = 3

    ! #slo# preliminary without dimensioning of prog/diag/aux

    !create state array for each domain
    DO jg = 1, n_dom

      ALLOCATE(p_os(jg)%p_prog(1:prlength), STAT=i_status)
      IF (i_status/=SUCCESS) THEN
         CALL finish(TRIM(routine), 'allocation of progn. state array failed')
      END IF

      ! construction loop: create components of state array
      ! !TODO organize var_lists for the multiple timesteps of prog. state
      WRITE(listname,'(a)')  'ocean_restart_var_list'
      CALL new_var_list(ocean_var_list, listname, patch_id=p_patch(jg)%id)
      CALL default_var_list_settings( ocean_var_list,            &
                                    & lrestart=.TRUE.,           &
                                    & restart_type=FILETYPE_NC2, &
                                    & model_type='oce' )
      DO jp = 1, prlength
         CALL construct_hydro_ocean_prog(p_patch(jg), p_os(jg)%p_prog(jp),jp)
      END DO

      CALL construct_hydro_ocean_diag(p_patch(jg), p_os(jg)%p_diag)

      CALL construct_hydro_ocean_aux(p_patch(jg), p_os(jg)%p_aux)

      CALL message(TRIM(routine),'construction of hydrostatic ocean state finished')

    END DO

  END SUBROUTINE construct_hydro_ocean_state

!-------------------------------------------------------------------------
!>
!!               Destructor for hydrostatic ocean state.
!
!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2006).
!!
  SUBROUTINE destruct_hydro_ocean_state(p_os)
    TYPE(t_hydro_ocean_state), TARGET,INTENT(inout)   :: p_os(n_dom)

    ! local variables

    INTEGER                                   :: jg, prlength, ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_state'


!-------------------------------------------------------------------------

    CALL message(TRIM(routine), 'start to destruct hydro ocean state ')

    prlength = SIZE(p_os(1)%p_prog)

    IF (prlength==0) THEN
      CALL finish(TRIM(routine),'prog array has length zero')
    END IF

    CALL delete_var_list(ocean_var_list)

    DO jg = 1, n_dom


      CALL destruct_hydro_ocean_diag(p_os(jg)%p_diag)
      CALL destruct_hydro_ocean_aux (p_os(jg)%p_aux)

      ! destruct state array
      ist = 1
      DEALLOCATE(p_os(jg)%p_prog, STAT=ist)
      IF (ist/=SUCCESS) THEN
         CALL finish(TRIM(routine),'deallocation of state array failed')
      END IF

    END DO

    CALL message(TRIM(routine),'destruction of hydrostatic ocean state finished')

  END SUBROUTINE destruct_hydro_ocean_state

!-------------------------------------------------------------------------
!>
!! Allocation of basic 3-dimensional structure components of hydrostatic ocean state.
!! Initialization of components with zero.
!
!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2011/06).
!!

  SUBROUTINE construct_hydro_ocean_base(p_patch, v_base)

    TYPE(t_patch), TARGET, INTENT(IN)          :: p_patch
    TYPE(t_hydro_ocean_base), INTENT(INOUT)    :: v_base

    ! local variables

    INTEGER :: ist
    INTEGER :: nblks_c, nblks_e, nblks_v, n_zlvp, n_zlvm!, ie
  ! INTEGER ::  jc,jb,jk, rl_start, rl_end
  ! INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_base'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine), 'start to construct basic hydro ocean state')

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v
    n_zlvp = n_zlev + 1
    n_zlvm = n_zlev - 1

    ! allocate and set vertical level thickness from the namelist
    ALLOCATE(v_base%del_zlev_m(n_zlev),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating del_zlev_m failed')
    ENDIF
    ALLOCATE(v_base%zlev_m(n_zlev),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating zlev_m failed')
    ENDIF
    ALLOCATE(v_base%zlev_i(n_zlvp),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating zlev_i failed')
    ENDIF
    ALLOCATE(v_base%del_zlev_i(n_zlev),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating del_zlev_i failed')
    ENDIF

    !
    !! 3-dim land-sea-mask at cells, edges and vertices
    !
    ! cells
    ALLOCATE(v_base%lsm_oce_c(nproma,n_zlev,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating lsm_oce_c failed')
    ENDIF
    ! edges
    ALLOCATE(v_base%lsm_oce_e(nproma,n_zlev,nblks_e),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating lsm_oce_e failed')
    ENDIF
    ! deepest ocean layer in column
    ALLOCATE(v_base%dolic_c(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating dolic_c failed')
    ENDIF
    ALLOCATE(v_base%dolic_e(nproma,nblks_e),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating dolic_e failed')
    ENDIF
    ! 2-dim basins and areas
    ALLOCATE(v_base%basin_c(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating basin_c failed')
    ENDIF
    ALLOCATE(v_base%regio_c(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating regio_c failed')
    ENDIF
    ALLOCATE(v_base%rbasin_c(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating basin_c failed')
    ENDIF
    ALLOCATE(v_base%rregio_c(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating regio_c failed')
    ENDIF
    ! 3-dim real land-sea-mask
    ! cells
    ALLOCATE(v_base%wet_c(nproma,n_zlev,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating wet_c failed')
    ENDIF
    ! edges
    ALLOCATE(v_base%wet_e(nproma,n_zlev,nblks_e),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating wet_e failed')
    ENDIF

    v_base%del_zlev_m = 0._wp
    v_base%del_zlev_i = 0._wp
    v_base%zlev_m     = 0._wp
    v_base%zlev_i     = 0._wp

    v_base%lsm_oce_c = 0
    v_base%lsm_oce_e = 0
    v_base%dolic_c = 0
    v_base%dolic_e = 0
    v_base%basin_c = 0
    v_base%regio_c = 0

    v_base%rbasin_c = 0.0_wp
    v_base%rregio_c = 0.0_wp

    v_base%wet_c = 0.0_wp
    v_base%wet_e = 0.0_wp

  END SUBROUTINE construct_hydro_ocean_base

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Stephan Lorenz, MPI-M (2011/06).
  !!
  SUBROUTINE destruct_hydro_ocean_base(v_base)

    TYPE(t_hydro_ocean_base), INTENT(INOUT) :: v_base

    ! local variables

    INTEGER :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_base'

    !CALL message(TRIM(routine),' start to destruct hydrostatic ocean basic state')

    DEALLOCATE(v_base%zlev_m,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating zlev_m failed')
    ENDIF

  END SUBROUTINE destruct_hydro_ocean_base

!-------------------------------------------------------------------------
!>
!!               Allocation of components of hydrostatic ocean prognostic state.
!!               Initialization of components with zero.
!
!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2006).
!!
  SUBROUTINE construct_hydro_ocean_prog(p_patch, p_os_prog, timelevel)

    TYPE(t_patch), INTENT(in), TARGET         :: p_patch
    TYPE(t_hydro_ocean_prog), INTENT(inout)   :: p_os_prog
    INTEGER, INTENT(IN)                       :: timelevel

    INTEGER  :: nblks_c, nblks_e !, nblks_v
    INTEGER  :: jtrc
    INTEGER, PARAMETER             :: max_oce_tracer = 2
    CHARACTER(len=max_char_length) :: oce_tracer_names(max_oce_tracer),&
    &                                 oce_tracer_units(max_oce_tracer),&
    &                                 oce_tracer_longnames(max_oce_tracer)
    INTEGER                        :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(len=max_char_length) :: var_suffix


    !-------------------------------------------------------------------------
    WRITE(var_suffix,'(a,i2.2)') '_TL',timelevel

    !-------------------------------------------------------------------------
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! height
    CALL add_var(ocean_var_list, 'h'//TRIM(var_suffix), p_os_prog%h , &
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, &
      &          t_cf_var('h', 'm', 'surface elevation at cell center'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,nblks_c/))

    !! normal velocity component
    CALL add_var(ocean_var_list, 'vn'//TRIM(var_suffix), p_os_prog%vn , GRID_UNSTRUCTURED_EDGE, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vn', 'm/s', 'normale velocity on edge,m'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL set_oce_tracer_info(max_oce_tracer      , &
        &                      oce_tracer_names    , &
        &                      oce_tracer_longnames, &
        &                      oce_tracer_codes    , &
        &                      oce_tracer_units,     &
        &                      var_suffix)
      CALL add_var(ocean_var_list, 'tracers'//TRIM(var_suffix), p_os_prog%tracer , &
      &            GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
      &            t_cf_var('tracers'//TRIM(var_suffix), '', '1:temperature 2:salinity'),&
      &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &            ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
      &            lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! Reference to individual tracer, for I/O

      ALLOCATE(p_os_prog%tracer_ptr(no_tracer))
      DO jtrc = 1,no_tracer
!      write(0,*)'jtrc:',jtrc

        CALL add_ref( ocean_var_list, 'tracers'//TRIM(var_suffix),              &
                    & oce_tracer_names(jtrc),                 &
                    & p_os_prog%tracer_ptr(jtrc)%p,                             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA,            &
                    & t_cf_var(oce_tracer_names(jtrc), &
                    &          oce_tracer_units(jtrc), &
                    &          oce_tracer_longnames(jtrc)), &
                    & t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
                    & ldims=(/nproma,n_zlev,nblks_c/))

      END DO
    ENDIF ! no_tracer > 0
  END SUBROUTINE construct_hydro_ocean_prog

  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean diagnostic state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!

  SUBROUTINE construct_hydro_ocean_diag(p_patch,p_os_diag)

    TYPE(t_patch), TARGET, INTENT(IN)          :: p_patch
    TYPE(t_hydro_ocean_diag), INTENT(INOUT)    :: p_os_diag

    ! local variables

    INTEGER :: ist
    INTEGER :: nblks_c, nblks_e, nblks_v
    INTEGER ::  jc,jb,jk, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_diag'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine), 'start to construct diagnostic hydro ocean state')

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    CALL add_var(ocean_var_list, 'rho', p_os_diag%rho , GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('rho', 'kg/m^3', 'density'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))

    CALL add_var(ocean_var_list, 'vt', p_os_diag%vt, GRID_UNSTRUCTURED_EDGE, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vt','m/s','tangential velocity at edges'),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE,GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    CALL add_var(ocean_var_list, 'h_e', p_os_diag%h_e, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_SURFACE, &
    &            t_cf_var('h_e','m','surface height ar edges'),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE,GRID_EDGE),&
    &            ldims=(/nproma,nblks_e/))

    ! thicknesses
    CALL add_var(ocean_var_list, 'thick_c', p_os_diag%thick_c,  &
    &            GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, &
    &            t_cf_var('thick_c','m','fluid column thickness at cells'),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE,GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list, 'thick_e', p_os_diag%thick_e, &
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, &
    &            t_cf_var('thick_e','m','fluid column thickness at edges'),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE,GRID_EDGE),&
    &            ldims=(/nproma,nblks_e/))

    CALL add_var(ocean_var_list, 'div_mass_flx_c', p_os_diag%div_mass_flx_c,&
    &            GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('div mass flux','','divergence mass flux at cells'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))
    CALL add_var(ocean_var_list, 'cons_thick_c', p_os_diag%cons_thick_c, &
    &            GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('cons thick','','consistent thickness at cells'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))
    CALL add_var(ocean_var_list, 'depth_c', p_os_diag%depth_c, &
    &            GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('depth_c','','time dependent depth at cells'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))
    CALL add_var(ocean_var_list, 'mass flux', p_os_diag%mass_flx_e, &
    &            GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('mass flux','',' mass flux'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! velocities
    CALL add_var(ocean_var_list, 'w', p_os_diag%w, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('w','m/s','vertical velocity at cells'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev+1,nblks_c/))
    CALL add_var(ocean_var_list, 'wtemp', p_os_diag%wtemp, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('wtemp','m/s','vertical velocity at cells'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev+1,nblks_c/))
    CALL add_var(ocean_var_list, 'w_old', p_os_diag%w_old, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA,&
    &            t_cf_var('w_old','m/s','vertical velocity at cells'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev+1,nblks_c/))
    CALL add_var(ocean_var_list, 'w_e', p_os_diag%w_e, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('w_e','m/s','vertical velocity at edges'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev+1,nblks_e/))
    CALL add_var(ocean_var_list, 'w_prev', p_os_diag%w_prev, &
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('w_prev','m/s','vertical velocity at edges'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev+1,nblks_c/),lrestart=.FALSE.)
    ! reconstructed u velocity component
    CALL add_var(ocean_var_list, 'u', p_os_diag%u, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('u','m/s','u velocity component'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))
    ! reconstructed v velocity component
    CALL add_var(ocean_var_list, 'v', p_os_diag%v, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('v','m/s','v velocity component'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))
    ! reconstrcuted velocity in cartesian coordinates
!   CALL add_var(ocean_var_list, 'p_vn', p_os_diag%p_vn, GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
!   &            t_cf_var('p_vn','m/s','normal velocity in cartesian coordinates'),&
!   &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
!   &            ldims=(/nproma,n_zlev,nblks_c/))
    ! integrated barotropic stream function
    CALL add_var(ocean_var_list, 'psi', p_os_diag%psi, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_SURFACE, &
    &            t_cf_var('psi','kg/s','barotropic stream function'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list, 'ptp_vn', p_os_diag%ptp_vn, &
    &            GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('ptp_vn','m/s','normal velocity in cartesian coordinates'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))
    ! predicted vn normal velocity component
    CALL add_var(ocean_var_list, 'vn_pred', p_os_diag%vn_pred, &
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vn_pred','m/s','predicted vn normal velocity component'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))
    ! predicted vn normal velocity component
    CALL add_var(ocean_var_list, 'vn_impl_vert_diff', p_os_diag%vn_impl_vert_diff,&
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vn_impl_vert_diff','m/s','predicted vn normal velocity component'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    CALL add_var(ocean_var_list, 'vn_time_weighted', p_os_diag%vn_time_weighted, &
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vn_pred','m/s','average vn normal velocity component'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    CALL add_var(ocean_var_list, 'w_time_weighted', p_os_diag%w_time_weighted, &
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('w_prev','m/s','vertical velocity at edges'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev+1,nblks_c/),lrestart=.FALSE.)

    ! vorticity
    CALL add_var(ocean_var_list, 'vort', p_os_diag%vort, &
    &            GRID_UNSTRUCTURED_VERT, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vort','1/s','vorticity'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_VERTEX),&
    &            ldims=(/nproma,n_zlev,nblks_v/))
    CALL add_var(ocean_var_list, 'vort_e', p_os_diag%vort_e, &
    &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('vort_e','1/s','vorticity at edges'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! kinetic energy component
    CALL add_var(ocean_var_list, 'kin', p_os_diag%kin, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('kin','J','kinetic energy'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))

    ! gradient term
    CALL add_var(ocean_var_list, 'grad', p_os_diag%grad, GRID_UNSTRUCTURED_EDGE, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('grad','','gradient'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! divergence component
    CALL add_var(ocean_var_list, 'div', p_os_diag%div, GRID_UNSTRUCTURED_CELL, &
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('div','','divergence'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_c/))

    ! pressures
    CALL add_var(ocean_var_list, 'press_hyd', p_os_diag%press_hyd, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('press_hyd','','hydrostatic pressure'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c/))
    CALL add_var(ocean_var_list, 'press_grad', p_os_diag%press_grad, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('press_grad','',' pressure gradient'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! mass flux
    !CALL add_var(ocean_var_list, 'flux_mass', p_os_diag%flux_mass, GRID_UNSTRUCTURED_EDGE,&
    !&            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('flux_mass','','mass flux at edges'),&
    !&            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    !&            ldims=(/nproma,n_zlev,nblks_e/))
    !CALL add_var(ocean_var_list, 'flux_tracer', p_os_diag%flux_tracer, GRID_UNSTRUCTURED_EDGE,&
    !&            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('flux_tracer','','tracers flux at edges'),&
    !&            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    !&            ldims=(/nproma,n_zlev,nblks_e/))

    ! horizontal velocity advection
    CALL add_var(ocean_var_list, 'veloc_adv_horz', p_os_diag%veloc_adv_horz, &
    &            GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('veloc_adv_horz','','horizontal velocity advection'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! vertical velocity advection
    CALL add_var(ocean_var_list, 'veloc_adv_vert', p_os_diag%veloc_adv_vert, &
    &            GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('veloc_adv_vert','','vertical velocity advection'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! horizontal diffusion
    CALL add_var(ocean_var_list, 'laplacian_horz', p_os_diag%laplacian_horz, &
    &            GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('laplacian_horz','','horizontal diffusion'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))
    ! vertical diffusion
    CALL add_var(ocean_var_list, 'laplacian_vert', p_os_diag%laplacian_vert, &
    &            GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('laplacian_vert','','vertical diffusion'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    ! initialize density with reference value rather than with zero
    !  - mainly for plotting purpose
    p_os_diag%rho(:,:,:) = rho_ref

    !reconstrcuted velocity in cartesian coordinates
    ALLOCATE(p_os_diag%p_vn(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at cells failed')
    END IF
    rl_start = 1
    rl_end = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
        &                  i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_os_diag%p_vn(jc,jk,jb)%x(:)=0.0_wp
        END DO
      END DO
    END DO

    ALLOCATE(p_os_diag%p_vn_dual(nproma,n_zlev,nblks_v), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at verts failed')
    END IF
    rl_start = 1
    rl_end = min_rlvert

    i_startblk = p_patch%verts%start_blk(rl_start,1)
    i_endblk   = p_patch%verts%end_blk(rl_end,1)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk
        CALL get_indices_v(p_patch, jb, i_startblk, i_endblk,&
        &                  i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_os_diag%p_vn_dual(jc,jk,jb)%x(:)=0.0_wp
        END DO
      END DO
    END DO

   !reconstrcuted velocity at edges in cartesian coordinates
    ALLOCATE(p_os_diag%p_vn_mean(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn_mean at edges failed')
    END IF
    rl_start = 1
    rl_end = min_rledge

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,1)
    DO jk=1,n_zlev
      DO jb = i_startblk, i_endblk
        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,&
        &                  i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_os_diag%p_vn_mean(jc,jk,jb)%x(:)=0.0_wp
        END DO
      END DO
    END DO

    !remapped velocity at cell edges
    ALLOCATE(p_os_diag%ptp_vn(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for ptp_vn at edges failed')
    END IF
    ! initialize all components with zero (this is preliminary)
    p_os_diag%ptp_vn    = 0.0_wp

    CALL add_var(ocean_var_list, 'temp_insitu', p_os_diag%temp_insitu, GRID_UNSTRUCTURED_CELL, &
      &          ZAXIS_DEPTH_BELOW_SEA, &
      &          t_cf_var('temp_insitu', 'K', 'in situ temperature'),&
      &          t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c/))

    !CALL message(TRIM(routine), 'construction of hydrostatic ocean diagnostic state finished')

  END SUBROUTINE construct_hydro_ocean_diag


  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
  SUBROUTINE destruct_hydro_ocean_diag(p_os_diag)

    TYPE(t_hydro_ocean_diag), INTENT(INOUT) :: p_os_diag

    ! local variables

    INTEGER :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_diag'

    DEALLOCATE(p_os_diag%p_vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn failed')
    END IF
    DEALLOCATE(p_os_diag%p_vn_dual, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn_dual failed')
    END IF
    DEALLOCATE(p_os_diag%ptp_vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for ptp_vn failed')
    END IF

  END SUBROUTINE destruct_hydro_ocean_diag
  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean auxiliary state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
  SUBROUTINE construct_hydro_ocean_aux(p_patch, p_os_aux)

    TYPE(t_patch),TARGET, INTENT(IN)                :: p_patch
    TYPE(t_hydro_ocean_aux), TARGET,INTENT(INOUT)   :: p_os_aux

    ! local variables

    INTEGER ::  ist, jc,jb, rl_start, rl_end, jtrc
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER ::  nblks_c, nblks_e, nblks_v

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_aux'
    !-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start to construct hydro ocean auxiliary state')

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! allocation for Adam-Bashford time stepping
    CALL add_var(ocean_var_list,'g_n',p_os_aux%g_n, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('g_n','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/),loutput=.TRUE.)
    CALL add_var(ocean_var_list,'g_nm1',p_os_aux%g_nm1, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('g_nm1','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/),loutput=.TRUE.)
    CALL add_var(ocean_var_list,'g_nimd',p_os_aux%g_nimd, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, t_cf_var('g_nimd','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/),loutput=.TRUE.)

    !-------------------------------------------------------------------------
    ! time stepping tracers go into the restart file. 4d has to be handles as 3D-references
    CALL add_var(ocean_var_list, 'g_n_c_h', p_os_aux%g_n_c_h , &
    &            GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('g_n_c_h', '', ''),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
    &            lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(p_os_aux%g_n_c_h_tracer_ptr(no_tracer))
    DO jtrc = 1,no_tracer
      CALL add_ref(ocean_var_list,'g_n_c_h',&
        &          'g_n_c_h_'//TRIM(oce_config%tracer_names(jtrc)),&
        &           p_os_aux%g_n_c_h_tracer_ptr(jtrc)%p, &
        &           GRID_UNSTRUCTURED_CELL,&
        &           ZAXIS_DEPTH_BELOW_SEA, &
        &           t_cf_var('g_n_c_h'//TRIM(oce_config%tracer_names(jtrc)),'',''),&
        &           t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &           ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END DO
    !-------------------------------------------------------------------------
    CALL add_var(ocean_var_list,'g_nm1_c_h',p_os_aux%g_nm1_c_h,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
      &          t_cf_var('g_nm1_c_h','',''),&
      &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(p_os_aux%g_nm1_c_h_tracer_ptr(no_tracer))
    DO jtrc = 1,no_tracer
      CALL add_ref(ocean_var_list,'g_nm1_c_h',&
        &          'g_nm1_c_h_'//TRIM(oce_config%tracer_names(jtrc)),&
        &          p_os_aux%g_nm1_c_h_tracer_ptr(jtrc)%p, &
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, &
        &          t_cf_var('g_nm1_c_h'//TRIM(oce_config%tracer_names(jtrc)),'',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END DO
    !-------------------------------------------------------------------------
    CALL add_var(ocean_var_list,'g_nimd_c_h',p_os_aux%g_nimd_c_h, &
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
      &          t_cf_var('g_nimd_c_h','',''),&
      &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(p_os_aux%g_nimd_c_h_tracer_ptr(no_tracer))
    DO jtrc = 1,no_tracer
      CALL add_ref(ocean_var_list,'g_nimd_c_h',&
        &          'g_nimd_c_h_'//TRIM(oce_config%tracer_names(jtrc)),&
        &          p_os_aux%g_nimd_c_h_tracer_ptr(jtrc)%p, &
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, &
        &          t_cf_var('g_nimd_c_h'//TRIM(oce_config%tracer_names(jtrc)),'',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END DO
    !-------------------------------------------------------------------------
    CALL add_var(ocean_var_list,'g_n_c_v',p_os_aux%g_n_c_v,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
      &          t_cf_var('g_n_c_v','',''),&
      &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(p_os_aux%g_n_c_v_tracer_ptr(no_tracer))
    DO jtrc = 1,no_tracer
      CALL add_ref(ocean_var_list,'g_n_c_v',&
        &          'g_n_c_v_'//TRIM(oce_config%tracer_names(jtrc)),&
        &          p_os_aux%g_n_c_v_tracer_ptr(jtrc)%p, &
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, &
        &          t_cf_var('g_n_c_v'//TRIM(oce_config%tracer_names(jtrc)),'',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END DO
    !-------------------------------------------------------------------------
    CALL add_var(ocean_var_list,'g_nm1_c_v', p_os_aux%g_nm1_c_v,&
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA,&
      &          t_cf_var('g_nm1_c_v','',''),&
      &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(p_os_aux%g_nm1_c_v_tracer_ptr(no_tracer))
    DO jtrc = 1,no_tracer
      CALL add_ref(ocean_var_list,'g_nm1_c_v',&
        &          'g_nm1_c_v_'//TRIM(oce_config%tracer_names(jtrc)),&
        &          p_os_aux%g_nm1_c_v_tracer_ptr(jtrc)%p, &
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, &
        &          t_cf_var('g_nm1_c_h'//TRIM(oce_config%tracer_names(jtrc)),'',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END DO
    !-------------------------------------------------------------------------
    CALL add_var(ocean_var_list,'g_nimd_c_v',&
      &          p_os_aux%g_nimd_c_v, GRID_UNSTRUCTURED_CELL,&
      &          ZAXIS_DEPTH_BELOW_SEA, t_cf_var('g_nimd_c_v','',''),&
      &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
      &          ldims=(/nproma,n_zlev,nblks_c,no_tracer/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
    ALLOCATE(p_os_aux%g_nimd_c_v_tracer_ptr(no_tracer))
    DO jtrc = 1,no_tracer
      CALL add_ref(ocean_var_list,'g_nimd_c_v',&
        &          'g_nimd_c_v_'//TRIM(oce_config%tracer_names(jtrc)),&
        &          p_os_aux%g_nimd_c_v_tracer_ptr(jtrc)%p, &
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, &
        &          t_cf_var('g_nimd_c_v'//TRIM(oce_config%tracer_names(jtrc)),'',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END DO
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------

    CALL add_var(ocean_var_list,'p_rhs_sfc_eq',p_os_aux%p_rhs_sfc_eq, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('p_rhs_sfc_eq','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))

    ! allocation for boundary conditions
    CALL add_var(ocean_var_list,'bc_top_u',p_os_aux%bc_top_u, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_top_u','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list,'bc_top_v',p_os_aux%bc_top_v, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_top_v','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    !CALL add_var(ocean_var_list,'bc_top_veloc_cc',p_os_aux%bc_top_veloc_cc, GRID_UNSTRUCTURED_CELL,&
    !&            ZAXIS_SURFACE, t_cf_var('bc_top_veloc_cc','',''),&
    !&            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL,ZAXIS_SURFACE),&
    !&            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list,'bc_top_vn',p_os_aux%bc_top_vn, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_SURFACE, t_cf_var('bc_top_vn','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,nblks_e/))

    CALL add_var(ocean_var_list,'bc_bot_u',p_os_aux%bc_bot_u, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_bot_u','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list,'bc_bot_v',p_os_aux%bc_bot_v, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_bot_v','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    !CALL add_var(ocean_var_list,'bc_bot_veloc_cc',p_os_aux%bc_bot_veloc_cc, GRID_UNSTRUCTURED_CELL,&
    !&            ZAXIS_SURFACE, t_cf_var('bc_bot_veloc_cc','',''),&
    !&            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL,ZAXIS_SURFACE),&
    !&            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list,'bc_bot_vn',p_os_aux%bc_bot_vn, GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_SURFACE, t_cf_var('bc_bot_vn','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,nblks_e/))

    CALL add_var(ocean_var_list,'bc_bot_w',p_os_aux%bc_bot_w, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_bot_w','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list,'bc_top_w',p_os_aux%bc_top_w, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_top_w','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c/))
    CALL add_var(ocean_var_list,'bc_bot_tracer',p_os_aux%bc_bot_tracer, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_bot_tracer','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c,no_tracer/),&
    &            lrestart=.FALSE.)
    CALL add_var(ocean_var_list,'bc_top_tracer',p_os_aux%bc_top_tracer, GRID_UNSTRUCTURED_CELL,&
    &            ZAXIS_SURFACE, t_cf_var('bc_top_tracer','',''),&
    &            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    &            ldims=(/nproma,nblks_c,no_tracer/),&
    &            lrestart=.FALSE.)

    ! allocation for divergence of fluxes
    !CALL add_var(ocean_var_list,'p_div_flux_horiz_act',p_os_aux%p_div_flux_horiz_act, GRID_UNSTRUCTURED_CELL,&
    !&            ZAXIS_DEPTH_BELOW_SEA, t_cf_vat('p_div_flux_horiz_act','',''),&
    !&            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    !&            ldims=(/nproma,n_zlev,nblks_c/))
    !CALL add_var(ocean_var_list,'p_div_flux_vert_act',p_os_aux%p_div_flux_vert_act, GRID_UNSTRUCTURED_CELL,&
    !&            ZAXIS_DEPTH_BELOW_SEA, t_cf_vat('p_div_flux_vert_act','',''),&
    !&            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    !&            ldims=(/nproma,n_zlev,nblks_c/))
    !CALL add_var(ocean_var_list,'p_div_flux_horiz_prev',p_os_aux%p_div_flux_horiz_prev, GRID_UNSTRUCTURED_CELL,&
    !&            ZAXIS_DEPTH_BELOW_SEA, t_cf_vat('p_div_flux_horiz_prev','',''),&
    !&            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    !&            ldims=(/nproma,n_zlev,nblks_c/))
    !CALL add_var(ocean_var_list,'p_div_flux_vert_prev',p_os_aux%p_div_flux_vert_prev, GRID_UNSTRUCTURED_CELL,&
    !&            ZAXIS_DEPTH_BELOW_SEA, t_cf_vat('p_div_flux_vert_prev','',''),&
    !&            t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
    !&            ldims=(/nproma,n_zlev,nblks_c/))

    ALLOCATE(p_os_aux%bc_top_veloc_cc(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    ALLOCATE(p_os_aux%bc_bot_veloc_cc(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    rl_start = 1
    rl_end = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        p_os_aux%bc_top_veloc_cc(jc,jb)%x = 0.0_wp
        p_os_aux%bc_bot_veloc_cc(jc,jb)%x = 0.0_wp
      END DO
   END DO

    ! allocation of 3-dim tracer relaxation:
    IF (no_tracer >= 1) THEN
      CALL add_var(ocean_var_list,'relax_3d_data_T',p_os_aux%relax_3d_data_T,&
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, t_cf_var('relax_3d_data_T','',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.FALSE.)
      CALL add_var(ocean_var_list,'relax_3d_forc_T',p_os_aux%relax_3d_forc_T,&
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, t_cf_var('relax_3d_forc_T','',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END IF
    IF (no_tracer == 2) THEN
      CALL add_var(ocean_var_list,'relax_3d_data_S',p_os_aux%relax_3d_data_S,&
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, t_cf_var('relax_3d_data_S','',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.FALSE.)
      CALL add_var(ocean_var_list,'relax_3d_forc_S',p_os_aux%relax_3d_forc_S,&
        &          GRID_UNSTRUCTURED_CELL,&
        &          ZAXIS_DEPTH_BELOW_SEA, t_cf_var('relax_3d_forc_S','',''),&
        &          t_grib2_var(255,255,255,16,GRID_REFERENCE, GRID_CELL),&
        &          ldims=(/nproma,n_zlev,nblks_c/),loutput=.TRUE.)
    END IF

  END SUBROUTINE construct_hydro_ocean_aux

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of auxilliary hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
  SUBROUTINE destruct_hydro_ocean_aux(p_os_aux)

    TYPE(t_hydro_ocean_aux), INTENT(INOUT)      :: p_os_aux

    ! local variables

    INTEGER                                   :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_aux'

!-------------------------------------------------------------------------


    DEALLOCATE(p_os_aux%bc_top_veloc_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond cc failed')
    END IF

   DEALLOCATE(p_os_aux%bc_bot_veloc_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of bot boundary cond cc failed')
    END IF

  END SUBROUTINE destruct_hydro_ocean_aux
!-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Sbr set boundary values for velocity field.
  !!
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2011).
  !
  !!  mpi parallelized LL (no sync required)
  SUBROUTINE  set_lateral_boundary_values( p_patch, vn )

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch
    REAL(wp)                  :: vn(:,:,:)

    ! local variables
    INTEGER :: jb, je, jk !, il_v, ib_v, ie, iv_ctr, il_e, ib_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: slev,elev 
    TYPE(t_subset_range), POINTER :: all_edges
!!$    CHARACTER(len=max_char_length), PARAMETER :: &
!!$      &      routine = 'mo_oce_state: set_lateral_boundary_values'

!---------------------------------------------------------------
    all_edges => p_patch%edges%all

! blocking
  slev         = 1
  elev         = n_zlev

  DO jb = all_edges%start_block, all_edges%end_block
    CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
    DO jk = slev, elev
      DO je= i_startidx_e, i_endidx_e
        IF ( v_base%lsm_oce_e(je,jk,jb) >= BOUNDARY ) THEN
          vn(je,jk,jb) = 0.0_wp
        ENDIF

!         il_v = p_patch%edges%vertex_idx(je,jb,1)
!         ib_v = p_patch%edges%vertex_blk(je,jb,1)
!         iv_ctr = 0
! 
!         DO ie = 1, p_patch%verts%num_edges(il_v,ib_v)
!           il_e = p_patch%verts%edge_idx(il_v,ib_v,ie)
!           ib_e = p_patch%verts%edge_blk(il_v,ib_v,ie)
! 
!           IF ( v_base%lsm_oce_e(il_e,jk,ib_e) /=sea ) THEN
!             iv_ctr = iv_ctr+1
!           ENDIF
!         END DO
!         IF(iv_ctr==1)THEN
!           DO ie = 1, p_patch%verts%num_edges(il_v,ib_v)
!             il_e = p_patch%verts%edge_idx(il_v,ib_v,ie)
!             ib_e = p_patch%verts%edge_blk(il_v,ib_v,ie)
!             vn(il_e,jk,ib_e) = 0.0_wp
!           END DO
!         ENDIF
      END DO
    END DO
  END DO

  END SUBROUTINE  set_lateral_boundary_values

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Initializes 3-dimensional structures for the ocean state variable.
  !!
  !! Calculates the vertical grid levels in z (meter) using the thickness of the elemental
  !! prisms (del_zlev_i) that are read from the ocean namelist.
  !!
  !! The 3-dimensional land-sea-mask is filled with values for interieur ocean, boundary
  !! ocean, and land, where parameter values from mo_impl_constants are used. The three
  !! dimensions are two for the nproma-blocking  and the middle one for the vertical levels.
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-02-19)
  !! Modified by Stephan Lorenz,        MPI-M (2010-08)
  !!  - fill dolic and sea_boundary as well
  !! Modified by Stephan Lorenz,        MPI-M (2011-05)
  !!  - level below surface receives same land-sea-mask
  !! Modified by Stephan Lorenz,        MPI-M (2011-06)
  !! - all 3-dim structures moved from patch_oce to type  t_hydro_ocean_base
  !!
  !!  mpi parallelized 
  SUBROUTINE init_ho_base( p_patch, p_ext_data, v_base )

    TYPE(t_patch),  TARGET,   INTENT(INOUT)    :: p_patch
    TYPE(t_external_data),    INTENT(INOUT)    :: p_ext_data
    TYPE(t_hydro_ocean_base), INTENT(INOUT)    :: v_base

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &        routine = 'mo_oce_state:init_ho_base'

    INTEGER :: jb, jc, je, jk, ji, nblks_c, nblks_e, npromz_c, npromz_e
    INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: noct1_c, noct1_e, inolsm, nowet_c
    INTEGER :: nolnd_c(n_zlev), nosea_c(n_zlev), nogllnd_c, noglsea_c
    INTEGER :: nolnd_e(n_zlev), nosea_e(n_zlev), nogllnd_e, noglsea_e
    INTEGER :: nobnd_e(n_zlev), nosbd_c(n_zlev), nolbd_c(n_zlev)
    INTEGER :: noglbnd_e, noglsbd_c, nogllbd_c
    INTEGER :: iic1, ibc1, iic2, ibc2, idxe, ible, idxn, ibln
    INTEGER :: n_zlvp, n_zlvm
    INTEGER :: jiter, niter, ctr, ctr_jk, ctr_glb
    INTEGER :: lsm_c   (nproma,p_patch%nblks_c)
    REAL(wp):: perc_lnd_c(n_zlev), perc_gllnd_c
    REAL(wp):: perc_lnd_e(n_zlev), perc_gllnd_e

    REAL(wp):: z_sync_c(nproma,p_patch%nblks_c)
    REAL(wp):: z_sync_e(nproma,p_patch%nblks_e)
    REAL(wp):: z_lat, z_lat_deg, z_north, z_south
    
    TYPE(t_subset_range), POINTER :: owned_cells, all_cells
    TYPE(t_subset_range), POINTER :: owned_edges
    INTEGER :: all_nobnd_e, all_nosbd_c, all_nolbd_c

    LOGICAL :: LIMITED_AREA = .FALSE.
    LOGICAL :: is_p_test_run

    !-----------------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    
    owned_cells => p_patch%cells%owned
    all_cells   => p_patch%cells%all
    owned_edges => p_patch%edges%owned
    
    is_p_test_run = p_test_run

    z_sync_c(:,:) = 0.0_wp
    z_sync_e(:,:) = 0.0_wp
    lsm_c   (:,:) = 0

    z_lat     = 0.0_wp 
    z_lat_deg = 0.0_wp 
    z_north   = 0.0_wp 
    z_south   = 0.0_wp
    !-----------------------------
    !
    ! Basic z-level configuration:
    !
    ! n_zlev    : number of z-coordinate surfaces - module global_variables, read in from namelist
    ! del_zlev_m: thickness of elemental prism - read in from namelist
    ! zlev_m    : position of coordinate surfaces in meters below zero surface
    ! zlev_i    : surface at the top of the respective z-coordinate surface (intermediate level)
    ! del_zlev_i: distance between two z-coordinate surfaces
    !
    !-----------------------------

    ! values for the blocking
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    !nblks_v = p_patch%nblks_v
    npromz_c = p_patch%npromz_c
    npromz_e = p_patch%npromz_e

    ! number of vertical levels from the namelist in mo_global_variables
    n_zlvp = n_zlev + 1
    n_zlvm = n_zlev - 1


    CALL set_del_zlev(n_zlev, dzlev_m,                  &
      &           v_base%del_zlev_i, v_base%del_zlev_m, &
      &           v_base%zlev_i    , v_base%zlev_m)

    !-----------------------------
    !
    ! Fill the 3-dim land-sea-mask and number of deepest ocean layer in column 'dolic'
    !
    !-----------------------------

    ! surface level: as read in ext_data:
    v_base%lsm_oce_c(:,1,:) = p_ext_data%oce%lsm_ctr_c(:,:)
    v_base%lsm_oce_e(:,1,:) = p_ext_data%oce%lsm_ctr_e(:,:)

    nogllnd_c = 0
    noglsea_c = 0
    nogllnd_e = 0
    noglsea_e = 0
    noglbnd_e = 0
    noglsbd_c = 0
    nogllbd_c = 0
    noct1_c = 0
    noct1_e = 0

    !  surface level and second level of lsm_oce_c defined by gridgenerator, not the current bathymetry
    v_base%lsm_oce_c(:,1,:) = p_ext_data%oce%lsm_ctr_c(:,:)

   IF(n_zlev>=2)&
   & v_base%lsm_oce_c(:,2,:) = p_ext_data%oce%lsm_ctr_c(:,:)

    !  first and second level of dolic_c defined by gridgenerator
    WHERE (p_ext_data%oce%lsm_ctr_c(:,:) <= SEA_BOUNDARY) v_base%dolic_c(:,:) = 2
    

    ! Coordinate surfaces - n_zlev z-levels:
    ! First vertical level loop to set wet cells below surface (and second layer) only

    INIT_ZLOOP: DO jk = 3, n_zlev

      !-----------------------------
      ! set dolic and wet grid points on cells:
      !  - if bathymetry is deeper than or equal to the coordinate surface (zlev_m)
      !    then grid point is wet; dolic is in that level
      !  - values for BOUNDARY set below

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx

          IF (p_ext_data%oce%bathymetry_c(jc,jb) <= -v_base%zlev_m(jk)) THEN
            v_base%lsm_oce_c(jc,jk,jb) = SEA
            v_base%dolic_c(jc,jb)      = jk
          ELSE IF (p_ext_data%oce%bathymetry_c(jc,jb)>-v_base%zlev_m(jk)) THEN
            v_base%lsm_oce_c(jc,jk,jb) = LAND
          END IF

        END DO
      END DO

      ! synchronize lsm on cells
      ! LL: this is done on all cells consistently on al procs
      ! so no sync is required

    END DO INIT_ZLOOP


    !-------------------------------------------------
    IF(LIMITED_AREA)THEN
    !  #slo# must be parallelized accordingly

      DO jk = 1, n_zlev
        z_south=-80.0_wp
        
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
          DO jc = i_startidx, i_endidx
    
              !get latitude of actual cell
              z_lat = p_patch%cells%center(jc,jb)%lat
              z_lat_deg = z_lat*rad2deg
    
              !If latitude of cell is above 80 N or below 80 S set triangle to land
              IF(z_lat_deg>z_north.OR.z_lat_deg<z_south)THEN
                v_base%lsm_oce_c(jc,:,jb)          = LAND
                p_ext_data%oce%bathymetry_c(jc,jb) = 100.0_wp
                v_base%dolic_c(jc,jb)              = 0
                v_base%wet_c(jc,:,jb)              = 0.0_wp
                !Set also all 3 edges to land
                DO ji = 1, 3
                  ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
                  idxe                                   = p_patch%cells%edge_idx(jc,jb,ji)
                  ible                                   = p_patch%cells%edge_blk(jc,jb,ji)
                  v_base%lsm_oce_e(idxe,:,ible)          = LAND
                  p_ext_data%oce%bathymetry_e(idxe,ible) = 100.0_wp
                  v_base%dolic_e(idxe,ible)              = 0
                  v_base%wet_e(idxe,:,ible)              = 0.0_wp
                END DO
             ENDIF
          END DO
        END DO
      END DO
    ENDIF



    !---------------------------------------------------------------------------------------------
    ! Correction loop for cells in all levels, similar to surface done in grid generator
    !  - through all levels each wet cell has at most one dry cell as neighbour
    !  - otherwise it is set to dry grid cell

    rl_start = 1           
    rl_end = min_rlcell_int   !  same as if not owner_mask?
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)

    niter=30

    DO jk=1,n_zlev
    
      ctr_jk = 0

      ! working on 2D lsm_c inside the loop
      lsm_c(:,:) = v_base%lsm_oce_c(:,jk,:)

      ! LL: disable checks here, the changes in halos will differ from seq run
      !     as the access patterns differ
      p_test_run = .false.
            
      DO jiter=1,niter
        !
        ctr = 0 ! no changes initially
        
        ! loop through owned patch cells    
        DO jb = owned_cells%start_block, owned_cells%end_block
          CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
        
            DO jc =  i_startidx, i_endidx
          
              nowet_c = 0

              ! LL: here we probably want to check if the above cell is land
              !     and change this into land accordingly
              
              IF (lsm_c(jc,jb) <= SEA_BOUNDARY) THEN
                DO ji = 1, 3
                  ! Get indices/blks of cells 1 to 3 adjacent to cell (jc,jb)
                  idxn = p_patch%cells%neighbor_idx(jc,jb,ji)
                  ibln = p_patch%cells%neighbor_blk(jc,jb,ji)
                  ! counts number of land-cells for all three neighbors
                  !  - only one land-point neighbor is allowed
                  IF ( lsm_c(idxn,ibln) > SEA_BOUNDARY) &
                    &  nowet_c=nowet_c + 1
                END DO

                ! More than 1 wet neighbor-cell then set cell to land
                !  - edges are set in the correction loop below
                IF ( nowet_c >= 2 ) THEN 
                  lsm_c(jc,jb)=LAND_BOUNDARY
                  ctr = ctr+1

                  IF (jk<3) THEN
                    WRITE(message_text,'(a,i2,a,i8)') &
                      &   'WARNING: Found 2 land neighbors at jc, jk=',jc,jk
                    CALL message(TRIM(routine), TRIM(message_text))
                  END IF

                END IF ! 2 land neighbors

              END IF ! lsm_c(jc,jb) <= SEA_BOUNDARY
              
          END DO  ! jc =  i_startidx, i_endidx
        END DO ! jb = owned_cells%start_block, owned_cells%end_block

        ! see what is the sum of changes of all procs
        ctr_glb = global_sum_array(ctr)
        
        WRITE(message_text,'(a,i2,a,i2,a,i8)') 'Level:', jk, &
          & ' Corrected wet cells with 2 land neighbors - iter=', &
          &                              jiter,' no of cor:',ctr_glb
        CALL message(TRIM(routine), TRIM(message_text))
        
        ! if no changes have been done, we are done with this level. Exit
        IF (ctr_glb == 0) EXIT

        ! we need to sync the halos here
        z_sync_c(:,:) =  REAL(lsm_c(:,:),wp)
        CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
        lsm_c(:,:) = INT(z_sync_c(:,:))

      END DO   ! jiter

      ! get back into 3D the slm
      v_base%lsm_oce_c(:,jk,:) = lsm_c(:,:)

      IF (is_p_test_run) THEN
        ! check if we have the correct slm
        p_test_run = is_p_test_run
        z_sync_c(:,:) =  REAL(lsm_c(:,:),wp)
        CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
        lsm_c(:,:) = INT(z_sync_c(:,:))               
      ENDIF

    END DO  ! jk=1,n_zlev

    ! restore p_test_run
    p_test_run = is_p_test_run

    !---------------------------------------------------------------------------------------------
    ! Now run through whole zlevel_loop after correction of cells for calculation of boundaries

    nogllnd_c = 0
    noglsea_c = 0
    nogllnd_e = 0
    noglsea_e = 0
    noglbnd_e = 0
    noglsbd_c = 0
    nogllbd_c = 0

    ! set dolic once more after jiter-correction
    v_base%dolic_c = 0
    v_base%dolic_e = 0

    ! Main loop for edges, dolic, boundaries, diagnosis and output
    !  - using lsm_oce_c after complete correction in jk>2 as input
    !  - (1) set land and sea values at cells <0 (sea) and >0 (land) - no boundaries
    !  - (2) set land and sea values at edges including boundaries
    !  - (3) set land and sea boundary values (-1 = SEA_BOUNDARY, 1=LAND_BOUNDARY)
    !
    ZLEVEL_LOOP: DO jk = 1, n_zlev

      !-----------------------------
      ! (1) set wet grid points and dolic at cells:
      !  - values for BOUNDARY set below

      nolnd_c(jk)=0
      nosea_c(jk)=0

      
      lsm_c(:,:) =  v_base%lsm_oce_c(:,jk,:)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx

        ! IF (.NOT.p_patch%cells%owner_mask(jc,jb)) CYCLE  ! access inner domain only
          IF (v_base%lsm_oce_c(jc,jk,jb) <= SEA_BOUNDARY) THEN
            nosea_c(jk)=nosea_c(jk)+1
            v_base%dolic_c(jc,jb) = jk
          ELSE
            ! -after correction: all other grid points are set to dry
            lsm_c(jc,jb) = LAND
            nolnd_c(jk)=nolnd_c(jk)+1
          END IF

          ! counting surface conditions as read from bathymetry - all wet cells
          IF (jk == 1) THEN
            IF (p_ext_data%oce%bathymetry_c(jc,jb) <= -v_base%zlev_m(jk)) &
              &   noct1_c = noct1_c+1
          ENDIF

        END DO

      END DO

      ! synchronize dolic_c
      z_sync_c(:,:) =  REAL(v_base%dolic_c(:,:),wp)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      v_base%dolic_c(:,:) = INT(z_sync_c(:,:))

      ! now synchronize auxiliary lsm_c and set to lsm_oce_c
      z_sync_c(:,:) =  REAL(lsm_c(:,:),wp)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      lsm_c(:,:) = INT(z_sync_c(:,:))
      v_base%lsm_oce_c(:,jk,:) = lsm_c(:,:)

      !  percentage of land area per level and global value 
      !   - here: nosea/nolnd include boundaries
      ctr         = global_sum_array(nolnd_c(jk))
      nolnd_c(jk) = ctr
      ctr         = global_sum_array(nosea_c(jk))
      nosea_c(jk) = ctr
      inolsm      = nolnd_c(jk) + nosea_c(jk)
      IF (inolsm == 0 ) THEN
        IF (jk < 3 ) CALL message (TRIM(routine), 'WARNING - number of cell points is zero?')
        perc_lnd_c(jk) = 0.0_wp
      ELSE
        perc_lnd_c(jk) = REAL(nolnd_c(jk),wp)/REAL(nosea_c(jk)+nolnd_c(jk),wp)*100.0_wp
        nogllnd_c = nogllnd_c + nolnd_c(jk)
        noglsea_c = noglsea_c + nosea_c(jk)
      END IF


      !-----------------------------
      ! (2) set wet grid points and dolic at edges (get values of neighbouring cells)
      !  - if the two corresponding cells are differing in sign then edge is BOUNDARY
      !  - if the two corresponding cells are <0 then edge is SEA
      !  - if the two corresponding cells are >0 then edge is LAND

      nolnd_e(jk)=0
      nosea_e(jk)=0
      nobnd_e(jk)=0

      rl_start = 1           !  #slo# - cannot run with holes on land in grid
      rl_end = min_rledge

      ! values for the blocking
      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,1)
      !
      ! loop through owned patch edges
      DO jb = owned_edges%start_block, owned_edges%end_block
        CALL get_index_range(owned_edges, jb, i_startidx, i_endidx)
        DO je = i_startidx, i_endidx

!           IF (.NOT.p_patch%edges%owner_mask(je,jb)) CYCLE  ! access inner domain only

          ! get indices/blks of cells 1 and 2 adjacent to edge (je,jb)
          iic1 = p_patch%edges%cell_idx(je,jb,1)
          ibc1 = p_patch%edges%cell_blk(je,jb,1)
          iic2 = p_patch%edges%cell_idx(je,jb,2)
          ibc2 = p_patch%edges%cell_blk(je,jb,2)
          !

          ! set land/sea for all edges
          IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) < 0)  .and.   &
            &  (v_base%lsm_oce_c(iic2,jk,ibc2) < 0) )        &
            &   v_base%lsm_oce_e(je,jk,jb) = SEA
          IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) > 0)  .and.   &
            &  (v_base%lsm_oce_c(iic2,jk,ibc2) > 0) )        &
            &   v_base%lsm_oce_e(je,jk,jb) = LAND

          ! set boundary values at edges
          IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) < 0)  .and.   &
            &  (v_base%lsm_oce_c(iic2,jk,ibc2) > 0) )        &
            &   v_base%lsm_oce_e(je,jk,jb) = BOUNDARY
          IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) > 0)  .and.   &
            &  (v_base%lsm_oce_c(iic2,jk,ibc2) < 0) )        &
            &   v_base%lsm_oce_e(je,jk,jb) = BOUNDARY

          ! count land/sea/boundary values (sum of nosea_e no_lnd_e nobnd_e is global value)
          IF ( v_base%lsm_oce_e(je,jk,jb) <  BOUNDARY )      &
            &  nosea_e(jk)=nosea_e(jk)+1
          IF ( v_base%lsm_oce_e(je,jk,jb) >  BOUNDARY )      &
            &  nolnd_e(jk)=nolnd_e(jk)+1
          IF ( v_base%lsm_oce_e(je,jk,jb) == BOUNDARY )      &
            &  nobnd_e(jk)=nobnd_e(jk)+1

          ! set dolic to jk if lsm_oce_e is wet or boundary (maximum depth)
          IF ( v_base%lsm_oce_e(je,jk,jb) <= BOUNDARY )      &
            &  v_base%dolic_e(je,jb) = jk

          ! counting surface conditions as read from bathymetry - all wet edges, boundary edges
          IF (jk == 1) THEN
            IF (p_ext_data%oce%bathymetry_e(je,jb) <= -v_base%zlev_m(jk)) noct1_e = noct1_e+1
          ENDIF

        END DO

      END DO

      ! synchronize lsm on edges
      z_sync_e(:,:) =  REAL(v_base%lsm_oce_e(:,jk,:),wp)
      CALL sync_patch_array(SYNC_E, p_patch, z_sync_e(:,:))
      v_base%lsm_oce_e(:,jk,:) = INT(z_sync_e(:,:))

      ! synchronize dolic_e
      z_sync_e(:,:) =  REAL(v_base%dolic_e(:,:),wp)
      CALL sync_patch_array(SYNC_e, p_patch, z_sync_e(:,:))
      v_base%dolic_e(:,:) = INT(z_sync_e(:,:))

      !  percentage of land area per level and global value
      ctr         = global_sum_array(nolnd_e(jk))
      nolnd_e(jk) = ctr
      ctr         = global_sum_array(nosea_e(jk))
      nosea_e(jk) = ctr
      inolsm = nolnd_e(jk) + nosea_e(jk)
      IF (inolsm == 0 ) THEN
        IF (jk == 1 ) CALL message (TRIM(routine), 'WARNING - number of edge points is zero?')
        perc_lnd_e(jk) = 0.0_wp
      ELSE
        perc_lnd_e(jk) = REAL(nolnd_e(jk),wp)/REAL(nosea_e(jk)+nolnd_e(jk),wp)*100.0_wp
        nogllnd_e = nogllnd_e + nolnd_e(jk)
        noglsea_e = noglsea_e + nosea_e(jk)
      END IF

      !-----------------------------
      ! (3) set values for LAND_BOUNDARY and SEA_BOUNDARY at cells
      !  - get values of neighbouring edges
      !  - if 1 of 3 edges of a sea-cell is BOUNDARY then cell is SEA_BOUNDARY
      !  - if 1 (or 2) of 3 edges of a land-cell is BOUNDARY then cell is LAND_BOUNDARY

      nosbd_c(jk)=0
      nolbd_c(jk)=0

      rl_start = 1           !  #slo# - cannot run with holes on land in grid
      rl_end = min_rlcell

      ! values for the blocking
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)
      !
      ! loop through all patch cells
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
        DO jc =  i_startidx, i_endidx

!           IF (.NOT.p_patch%cells%owner_mask(jc,jb)) CYCLE  ! access inner domain only

          ! sea points
          IF (v_base%lsm_oce_c(jc,jk,jb) < 0) THEN

            DO ji = 1, 3
              ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
              idxe = p_patch%cells%edge_idx(jc,jb,ji)
              ible = p_patch%cells%edge_blk(jc,jb,ji)
              ! if one of lsm_e is boundary then lsm_c is sea_boundary
              IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY ) &
                &  v_base%lsm_oce_c(jc,jk,jb) = SEA_BOUNDARY
            END DO

            ! count sea boundary for all levels
            IF ( v_base%lsm_oce_c(jc,jk,jb) == SEA_BOUNDARY )  &
              &  nosbd_c(jk)=nosbd_c(jk)+1
          END IF  !  lsm_c < 0

          ! land points
          IF (v_base%lsm_oce_c(jc,jk,jb) > BOUNDARY) THEN

            DO ji = 1, 3
              ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
              idxe = p_patch%cells%edge_idx(jc,jb,ji)
              ible = p_patch%cells%edge_blk(jc,jb,ji)
              ! if one of lsm_e is boundary then lsm_c is land_boundary
              IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY ) &
                &  v_base%lsm_oce_c(jc,jk,jb) = LAND_BOUNDARY
            END DO

            IF ( v_base%lsm_oce_c(jc,jk,jb) == LAND_BOUNDARY )   &
              &  nolbd_c(jk)=nolbd_c(jk)+1

          END IF  !  lsm_c > 0

        END DO

      END DO
      
      all_nobnd_e = global_sum_array( nobnd_e(jk))
      all_nosbd_c = global_sum_array( nosbd_c(jk))
      all_nolbd_c = global_sum_array( nolbd_c(jk))

      nobnd_e(jk) = all_nobnd_e
      nosbd_c(jk) = all_nosbd_c
      nolbd_c(jk) = all_nolbd_c
      
      noglbnd_e = noglbnd_e + all_nobnd_e
      noglsbd_c = noglsbd_c + all_nolbd_c
      nogllbd_c = nogllbd_c + all_nolbd_c

      ! synchronize lsm on cells
      ! LL: not necessary since we go throu all cells
!       z_sync_c(:,:) =  REAL(v_base%lsm_oce_c(:,jk,:),wp)
!       CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
!       v_base%lsm_oce_c(:,jk,:) = INT(z_sync_c(:,:))

    END DO ZLEVEL_LOOP

!---------------------------------------------------------------------------------------------
    ! Output the levels
    WRITE(message_text,'(a,a)') &
    &     'LEVEL   zlev_m  Thickness   zlev_i  Distance ', &
    &     '    SEA_c    LAND_c  PERC_LND     SEA_e    LAND_e   BND_e   SEA_B. LAND_B.'
    CALL message('', TRIM(message_text))
    DO jk = 1, n_zlev
      WRITE(message_text,'(a,i3,4f10.2,2i10,f10.2,2i10,3i8)') '.',  &
        &   jk, v_base%zlev_m(jk), v_base%del_zlev_m(jk), &
        &       v_base%zlev_i(jk), v_base%del_zlev_i(jk), &
        &   nosea_c(jk), nolnd_c(jk), perc_lnd_c(jk), nosea_e(jk), nolnd_e(jk),     &
        &   nobnd_e(jk), nosbd_c(jk), nolbd_c(jk)
      CALL message('', message_text)
    END DO

    ! Output last level
    inolsm = nogllnd_c + noglsea_c
    IF ( inolsm == 0 ) THEN
      CALL message (TRIM(routine), 'WARNING - number of global cell points is zero?')
      perc_gllnd_c = 0.0_wp
    ELSE
      perc_gllnd_c = REAL(nogllnd_c,wp)/REAL(noglsea_c+nogllnd_c,wp)*100.0_wp
    END IF
    inolsm = nogllnd_e + noglsea_e
    IF ( inolsm == 0 ) THEN
      CALL message (TRIM(routine), 'WARNING - number of global edge points is zero?')
      perc_gllnd_e = 0.0_wp
    ELSE
      perc_gllnd_e = REAL(nogllnd_e,wp)/REAL(noglsea_e+nogllnd_e,wp)*100.0_wp
    END IF
    n_zlvp = n_zlev + 1
    WRITE(message_text,'(a,f20.2,a,i9,i10,f10.2,2i10,3i8)') &
    &     'Bottom Level: ', v_base%zlev_i(n_zlvp), &
    &     '    GLOBAL:',     noglsea_c, nogllnd_c, perc_gllnd_c, &
    &                        noglsea_e, nogllnd_e, noglbnd_e, noglsbd_c, nogllbd_c
    CALL message('', TRIM(message_text))

    ! Warnings occur if create_ocean_grid parameter mindepth is not half the
    ! depth of the surface level dzlev_m(1) in namelist ocean_ctl (must not be an error)
    IF ( nosea_c(1) /= noct1_c ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        &   'WARNING - surface sea-cells read = ',nosea_c(1), &
        &   ' - calculated from bathymetry = ',noct1_c
      CALL message(routine, TRIM(message_text))
    END IF
    IF ( nosea_e(1) /= noct1_e ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        &   'WARNING - surface sea-edges read = ',nosea_e(1), &
        &   ' - calculated from bathymetry = ',noct1_e
      CALL message(routine, TRIM(message_text))
    END IF

    !-----------------------------
    ! real bathymetry should not be used since individual bottom layer thickness is not implemented
    ! set values of bathymetry to new non-individual dolic values

  ! DO jb = 1, nblks_c
  !   i_endidx=nproma
  !   IF (jb==nblks_c) i_endidx=npromz_c
  !     DO jc = 1, i_endidx
  !       v_base%bathymetry_c(jc,jb) = &
  !         &  -v_base%zlev_i(v_base%dolic_c(jc,jb)+1)
  !   ENDDO
  ! ENDDO

  ! DO jb = 1, nblks_e
  !   i_endidx=nproma
  !   IF (jb==nblks_e) i_endidx=npromz_e
  !     DO je = 1, i_endidx
  !       v_base%bathymetry_e(je,jb) = &
  !         &  -v_base%zlev_i(v_base%dolic_e(je,jb)+1)
  !   ENDDO
  ! ENDDO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_ho_base

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Initializes 2-dimensional definitions of basins for the ocean state variable.
  !!
  !! The 2-dimensional land-sea-mask is used to define basins for calculation of
  !! meridional overturning circulation (MOC) and to define areas of certain interest.
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2012-02)
  !! Modified by Stephan Lorenz,        MPI-M (2012-02)
  !! 
  !!  no-mpi parallelized
  SUBROUTINE init_ho_basins( p_patch, v_base )

    TYPE(t_patch), TARGET, INTENT(IN)          :: p_patch
    TYPE(t_hydro_ocean_base), INTENT(INOUT)    :: v_base

    REAL(wp) :: z_sync_c(nproma,p_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,p_patch%nblks_e)
    INTEGER  :: ibase   (nproma,p_patch%nblks_c)
    INTEGER  :: iarea   (nproma,p_patch%nblks_c)

    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: jb, jc, jk, i_endidx, nblks_c, npromz_c, i
    INTEGER  :: no_cor, g_cor, no_glb, jiter, iter
    INTEGER  :: n_idx(3), n_blk(3)
    REAL(wp) :: z60n, z30n, z30s, z85s, z10n, z100w
    REAL(wp) :: z_lat_deg, z_lon_deg
    REAL(wp) :: z_lon_pta, z_lon_ati, z_lon_itp, z_lon_ind, z_lon_nam, z_lon_med

    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain
    LOGICAL :: is_area_5, is_area_8, p_test_run_bac

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &        routine = 'mo_oce_state:init_ho_basins'

    !-----------------------------------------------------------------------------
    all_cells => p_patch%cells%all
    cells_in_domain => p_patch%cells%in_domain
    CALL message (TRIM(routine), 'start')

    z_sync_c(:,:) = 0.0_wp
    z_sync_e(:,:) = 0.0_wp

    ! values for the blocking
    nblks_c = p_patch%nblks_c
    !nblks_e = p_patch%nblks_e
    !nblks_v = p_patch%nblks_v
    npromz_c = p_patch%npromz_c
    !npromz_e = p_patch%npromz_e

    !-----------------------------
    !
    ! Ocean basins:
    !  1: Atlantic; 2: Indian; 3: Pacific; 4: Southern Ocean (for global)
    !
    !-----------------------------
  
    !-----------------------------
    !
    ! Ocean areas/regions:
    !  0 = land point
    !  1 = Greenland-Iceland-Norwegian Sea
    !  2 = Arctic Ocean
    !  3 = Labrador Sea
    !  4 = North Atlantic Ocean
    !  5 = Tropical Atlantic Ocean
    !  6 = Southern Ocean
    !  7 = Indian Ocean
    !  8 = Tropical Pacific Ocean
    !  9 = North Pacific Ocean
    !
    !-----------------------------

    !-----------------------------
    ! Define borders of region:
    !  two problematic regions remain that can be accessed via space filling curves
    !   - Caribbian Sea is partly divided in Pacific/Atlantic (border are land points)
    !     iterative loop is used as first guess, see below
    !   - Indonesian Region is both in Pacific/Indian Ocean 
    !     there is no clear land border since the Indonesian Throughflow(s) exist
    !     here a slanted geographic line can be implemented

    z60n     =  61.0_wp
    z30n     =  30.0_wp
    z30s     = -30.0_wp
    z85s     = -85.0_wp
    z_lon_pta = -70.0_wp   !  Pac/Atl - Drake Passage
    z_lon_ati =  25.0_wp   !  Atl/Ind - South Africa
    z_lon_itp = 115.0_wp   !  Ind/Pac - Australia-Indonesia
    z_lon_ind = 100.0_wp   !  Ind/Pac - Indonesia (north of Equator)
    z_lon_nam = -90.0_wp   !  Pac/Atl - North America
    z_lon_med =  50.0_wp   !  Atl/Ind - Mediterranean

    ! for Caribbean: 5N and 105W are crucial
    z10n     =    5.0_wp
    z100w    = -105.0_wp

    !-----------------------------
    ! Fill ocean areas:

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
  
      DO jc = i_startidx_c, i_endidx_c
  
         ! get lat/lon of actual cell
         z_lat_deg = rad2deg * p_patch%cells%center(jc,jb)%lat
         z_lon_deg = rad2deg * p_patch%cells%center(jc,jb)%lon
         IF (z_lon_deg >  180.0_wp) z_lon_deg = z_lon_deg-360.0_wp
         IF (z_lon_deg < -180.0_wp) z_lon_deg = z_lon_deg+360.0_wp

         ! Arctic Ocean: default
         iarea(jc,jb) = 2

         ! GIN Sea (not yet)

         ! Labrador Sea (not yet)

         ! 4 = North Atlantic
         IF (                                                           &
           & (z_lat_deg >= z30n      .AND. z_lat_deg < z60n)     .AND.  &
           & (z_lon_deg >= z_lon_nam .AND. z_lon_deg < z_lon_med)       &
           & ) iarea(jc,jb) = 4

         ! 9 = North Pacific
         IF (                                                           &
           & (z_lat_deg >= z30n      .AND. z_lat_deg < z60n)     .AND.  &
           & (z_lon_deg >= z_lon_itp .OR.  z_lon_deg < z_lon_nam)       &
           & ) iarea(jc,jb) = 9

         ! 5 = Tropical Atlantic (without Caribbean - yet)
         IF (                                                           &
           & (z_lat_deg >= z30s      .AND. z_lat_deg < z30n)     .AND.  &
           & (z_lon_deg >= z_lon_pta .AND. z_lon_deg < z_lon_med)       &
           & ) iarea(jc,jb) = 5

         ! 6 = Southern Ocean
         IF (z_lat_deg < z30s) iarea(jc,jb) = 6

         ! 7 = Indian (including Indonesian Pacific - yet)
         IF (                                                           &
           & (z_lat_deg >= z30s      .AND. z_lat_deg < z30n)     .AND.  &
           & (z_lon_deg >= z_lon_ati .AND. z_lon_deg < z_lon_itp)       &
           & ) iarea(jc,jb) = 7

         ! 8 = Tropical Pacific
         IF (                                                           &
           & (z_lat_deg >= z30s      .AND. z_lat_deg < z30n)     .AND.  &
           & (z_lon_deg >= z_lon_itp .OR.  z_lon_deg < z_lon_pta)       &
           & ) iarea(jc,jb) = 8

         ! Crucial Region: Caribbean undefined (-33)
         IF (                                                           &
           & (z_lat_deg >= z10n      .AND. z_lat_deg < z30n)     .AND.  &
           & (z_lon_deg >= z100w     .AND. z_lon_deg < z_lon_pta)       &
           & ) iarea(jc,jb) = -33

         ! Land points
         IF (v_base%lsm_oce_c(jc,1,jb) >= BOUNDARY) iarea(jc,jb) = 0
  
      END DO
    END DO

    !-----------------------------
    ! Fill crucial areas:

    ! Caribbean: border is land point
    !  - border is at least one line of land points (not valid for Indonesian)
    !  - no cell has Pacific and Atlantic neighbor
    !  - iterative procedure (time consuming in higher resolution)

    ! disable p_test_run since iterations will be different
    p_test_run_bac = p_test_run
    p_test_run = .false.

    g_cor=0
!   Do jiter=1,1000 ! should be adjusted to higher resolution
    Do jiter=1,100
!   Do jiter=1,10
!   Do jiter=1,1

      no_cor = 0
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
           is_area_5 = .false.
           is_area_8 = .false.
           IF (iarea(jc,jb) == -33) THEN
             DO i = 1, 3  !  no_of_edges
               ! coordinates of neighbouring cells
               n_blk(i) = p_patch%cells%neighbor_blk(jc,jb,i)
               n_idx(i) = p_patch%cells%neighbor_idx(jc,jb,i)
               IF (iarea(n_idx(i),n_blk(i)) == 5) THEN       ! neighbor is TropAtl
                 is_area_5 = .true.
               ELSE IF (iarea(n_idx(i),n_blk(i)) == 8) THEN  ! neighbor is TropPac
                 is_area_8 = .true.
               END IF
             END DO

             IF (is_area_5) THEN 
                iarea(jc,jb)=5
                no_cor =  no_cor + 1
             ENDIF
             IF (is_area_8) THEN
               iarea(jc,jb)=8
               no_cor =  no_cor + 1
             ENDIF
           END IF
     
        END DO
      END DO
      
      no_glb = global_sum_array(no_cor)
      no_cor = no_glb
      g_cor=g_cor+no_cor

      iter=jiter
      IF (no_cor == 0) exit

      WRITE(message_text,'(a,i4,a,i8)') 'Corrected Caribbean region - iter=', &
        &                              jiter,' no of cor:',no_cor
      CALL message(TRIM(routine), TRIM(message_text))

      ! do sync
      
      z_sync_c(:,:) =  REAL(iarea(:,:),wp)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      iarea(:,:) = INT(z_sync_c(:,:))

    END DO

    WRITE(message_text,'(a,i4,a,i8)') 'Corrected Caribbean region - iterations=', &
      &                              iter,' no of cor:',g_cor
    CALL message(TRIM(routine), TRIM(message_text))
    
    p_test_run = p_test_run_bac
    !chekc if iarea is the same
    IF (p_test_run) THEN
       z_sync_c(:,:) =  REAL(iarea(:,:),wp)
       CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
       iarea(:,:) = INT(z_sync_c(:,:))
    ENDIF
    !-----------------------------

    !-----------------------------
    ! Fill ocean basins using ocean areas:

    WHERE ( iarea(:,:) <= 5 )
      ibase(:,:) = 1
    ELSEWHERE ( iarea(:,:) == 6 )
      ibase(:,:) = 4
    ELSEWHERE ( iarea(:,:) == 7 )
      ibase(:,:) = 2
    ELSEWHERE ( iarea(:,:) >= 8 )
      ibase(:,:) = 3
    END WHERE

    WHERE ( iarea(:,:) == 0 )
      ibase(:,:) = 0
    END WHERE

    v_base%basin_c(:,:) = ibase(:,:)
    v_base%regio_c(:,:) = iarea(:,:)
    v_base%rbasin_c(:,:) = REAL(ibase(:,:),wp)
    v_base%rregio_c(:,:) = REAL(iarea(:,:),wp)

  ! write(66,*) 'IBASE:'
  ! write(66,'(100i1)') ibase(:,:)
  ! write(66,*) 'IAREA:'
  ! write(66,'(100i1)') iarea(:,:)

    !-----------------------------
    ! set wet_c and wet_e to 1 at sea points including boundaries

    ! cells
    WHERE ( v_base%lsm_oce_c(:,:,:) <= SEA_BOUNDARY )
      v_base%wet_c(:,:,:) = 1.0_wp
    END WHERE

    ! edges
    WHERE ( v_base%lsm_oce_e(:,:,:) <= SEA_BOUNDARY )
      v_base%wet_e(:,:,:) = 1.0_wp
    END WHERE

    ! #slo# for test:
    !v_base%wet_c(:,:,:) = real(v_base%lsm_oce_c(:,:,:),wp)
    !v_base%wet_e(:,:,:) = real(v_base%lsm_oce_e(:,:,:),wp)

    ! intermediate levels: same as wet_c
    !WHERE ( v_base%lsm_oce_c(:,:,:) <= SEA_BOUNDARY )
    !  v_base%wet_i(:,:,:) = 1.0_wp
    !END WHERE

    ! synchronize all elements of v_base - not necessary
    z_sync_c(:,:) =  REAL(v_base%basin_c(:,:),wp)
    CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
    v_base%basin_c(:,:) = INT(z_sync_c(:,:))
    z_sync_c(:,:) =  REAL(v_base%regio_c(:,:),wp)
    CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
    v_base%regio_c(:,:) = INT(z_sync_c(:,:))

    z_sync_c(:,:) =  v_base%rbasin_c(:,:)
    CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
    v_base%rbasin_c(:,:) = z_sync_c(:,:)
    z_sync_c(:,:) =  v_base%rregio_c(:,:)
    CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
    v_base%rregio_c(:,:) = z_sync_c(:,:)

    DO jk = 1, n_zlev

      z_sync_c(:,:) =  v_base%wet_c(:,jk,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      v_base%wet_c(:,jk,:) = z_sync_c(:,:)

      z_sync_e(:,:) =  v_base%wet_e(:,jk,:)
      CALL sync_patch_array(SYNC_e, p_patch, z_sync_e(:,:))
      v_base%wet_e(:,jk,:) = z_sync_e(:,:)

    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_ho_basins

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Modifies the already calculated Coriolis force, if beta-, f-plane or the nonrotating case 
  !! is selected in the namelist. The tangent plane is associated to the center of the basin that is 
  !! specified in the namelist. An alternative would be to associate it to the nearest edge/vertex,
  !! but this is not implemented yet, and i expect it to have a minor effect.
  !!
  !! The Coriolis parameter is specified for edges (needed in RBF-discretization) and at vertices
  !! (needed in mimetic discreization).
  !! The land-sea masks are not taken into account here. This would require to extend the
  !! 2D-coriolis-structure to a 3D one
  !!
  !! @par Revision History
  !!  developed by Peter Korn, 2011
  !!  mpi parallelized LL (no sync required)
  !!
  SUBROUTINE init_coriolis_oce( ptr_patch )
    !
    IMPLICIT NONE
    !
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    !
    INTEGER :: jb, je, jv
    INTEGER :: i_startidx_e, i_endidx_e
    INTEGER :: i_startidx_v, i_endidx_v
    TYPE(t_geographical_coordinates) :: gc1,gc2 
    TYPE(t_cartesian_coordinates) :: xx1, xx2
    REAL(wp) :: z_y, z_lat_basin_center
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    & routine = ('mo_oce_state:init_coriolis_oce')
    TYPE(t_subset_range), POINTER :: all_verts, all_edges
    !-----------------------------------------------------------------------
    all_verts => ptr_patch%verts%all
    all_edges => ptr_patch%edges%all

    CALL message (TRIM(routine), 'start')

    SELECT CASE (CORIOLIS_TYPE)

    CASE(BETA_PLANE_CORIOLIS)

      CALL message (TRIM(routine), 'BETA_PLANE_CORIOLIS: set to linear approximation')

      z_lat_basin_center = basin_center_lat * deg2rad
      gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
      gc1%lon = 0.0_wp
      xx1=gc2cc(gc1)

      DO jb = all_verts%start_block, all_verts%end_block
        CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v
            !z_y = re*(ptr_patch%verts%vertex(jv,jb)%lat - z_lat_basin_center) 
            gc2%lat = ptr_patch%verts%vertex(jv,jb)%lat!*deg2rad
            gc2%lon = 0.0_wp
            xx2=gc2cc(gc2)        
            z_y = re*arc_length(xx2,xx1)
            ptr_patch%verts%f_v(jv,jb) = 2.0_wp*omega*( sin(z_lat_basin_center)     &
            &                          + (cos(z_lat_basin_center)/re)*z_y)
         !  write(*,*)'beta', jv,jb,z_beta_plane_vort,2.0_wp*omega*sin(z_lat_basin_center),&
         !  &2.0_wp*omega*((cos(z_lat_basin_center)/re)*z_y)
        END DO
      END DO

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e
          ! depends on basin_center_lat only - not dependent on center_lon, basin_width or height
            gc2%lat = ptr_patch%edges%center(je,jb)%lat!*deg2rad
            gc2%lon = 0.0_wp
            xx2=gc2cc(gc2)        
            z_y = re*arc_length(xx2,xx1)

            !z_y = ptr_patch%edges%center(je,jb)%lat - z_lat_basin_center
            ptr_patch%edges%f_e(je,jb) = 2.0_wp*omega*( sin(z_lat_basin_center)     &
            &                          + (cos(z_lat_basin_center)/re)*z_y)
        END DO
      END DO
    CASE(F_PLANE_CORIOLIS)

      CALL message (TRIM(routine), 'F_PLANE_CORIOLIS: set to constant value')
   
      z_lat_basin_center = basin_center_lat * deg2rad
   
      ptr_patch%edges%f_e  = 2.0_wp*omega*sin(z_lat_basin_center)
      ptr_patch%verts%f_v  = 2.0_wp*omega*sin(z_lat_basin_center)
   
    CASE(ZERO_CORIOLIS)
   
      CALL message (TRIM(routine), 'ZERO_CORIOLIS: set to zero')
      ptr_patch%verts%f_v = 0.0_wp
      ptr_patch%edges%f_e = 0.0_wp

    CASE(FULL_CORIOLIS)

      CALL message (TRIM(routine), 'FULL_CORIOLIS: Nothing to do, coriolis not modified')

    END SELECT

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_coriolis_oce
!-------------------------------------------------------------------------  
!
!!! Helper functions for computing the vertical layer structure  
!>
!!
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2011).
!!

  SUBROUTINE set_zlev(zlev_i, zlev_m)
    REAL(wp), INTENT(INOUT) :: zlev_i(n_zlev+1)    , zlev_m(n_zlev)
!--------------------------------------
    INTEGER :: jk

    zlev_m(1) = 0.5_wp * dzlev_m(1)
    zlev_i(1) = 0.0_wp
    ! zlev_i    : upper border surface of vertical cells
    DO jk = 2, n_zlev+1
      zlev_i(jk) = zlev_i(jk-1) + dzlev_m(jk-1)
    END DO

    ! zlev_m    : position of coordinate surfaces in meters below zero surface.
    DO jk = 2, n_zlev
      zlev_m(jk) = 0.5_wp * ( zlev_i(jk+1) + zlev_i(jk)  )
    END DO
  END SUBROUTINE set_zlev
!-------------------------------------------------------------------------  
!
!!Subroutine calculates vertical coordinates
!>
!!
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2011).
!!
  SUBROUTINE set_del_zlev(n_zlev, dzlev_m, del_zlev_i, del_zlev_m, zlev_i, zlev_m)
    INTEGER,  INTENT(IN)  :: n_zlev
    REAL(wp), INTENT(IN) :: dzlev_m(100)
    REAL(wp) :: del_zlev_i(n_zlev), del_zlev_m(n_zlev)
    REAL(wp) :: zlev_i(n_zlev+1)    , zlev_m(n_zlev)
    
    INTEGER :: jk
!!-------------------------------------
    CALL set_zlev(zlev_i, zlev_m)
    ! del_zlev_i: distance between two z-coordinate surfaces.
    !             The first is the distance from the ocean surface = zlev_m(1)
    del_zlev_i(1) = zlev_m(1)
    DO jk = 2, n_zlev
      del_zlev_i(jk) = zlev_m(jk) -  zlev_m(jk-1)
    END DO
!TODO    del_zlev_i(n_zlev+1) = 0.5*dzlev_m(n_zlev)
    del_zlev_m(:) = dzlev_m(1:n_zlev)
  END SUBROUTINE set_del_zlev
!-------------------------------------------------------------------------  
!
!!Subroutine 
!>
!!
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2011).
!!
  SUBROUTINE set_oce_tracer_info(max_oce_tracer,&
  &                              oce_tracer_names,&
  &                              oce_tracer_longnames,&
  &                              oce_tracer_codes,&
  &                              oce_tracer_units,&
  &                              suffix)

    INTEGER, INTENT(IN)            :: max_oce_tracer
    CHARACTER(len=max_char_length) :: oce_tracer_names(max_oce_tracer),&
      &                               oce_tracer_units(max_oce_tracer),&
      &                               oce_tracer_longnames(max_oce_tracer)
    INTEGER                        :: oce_tracer_codes(max_oce_tracer)
    CHARACTER(len=max_char_length), OPTIONAL :: suffix

    IF (max_oce_tracer < no_tracer) THEN
      CALL finish('set_oce_tracer_info','Too many tracers! Please provide trace info')
    ENDIF
    IF (PRESENT(suffix)) THEN
!     write(0,*)'suffix:',suffix
    END IF
    oce_tracer_names(1)     = 'T'
    IF (PRESENT(suffix)) THEN
      oce_tracer_names(1) = 'T'//TRIM(suffix)
    END IF
    oce_tracer_longnames(1) = 'potential temperature'
    oce_tracer_units(1)     = 'deg C'
    oce_tracer_codes(1)     = 200

    oce_tracer_names(2)     = 'S'
    IF (PRESENT(suffix)) THEN
      oce_tracer_names(2) = 'S'//TRIM(suffix)
    END IF
    oce_tracer_longnames(2) = 'salinity'
    oce_tracer_units(2)     = 'psu'
    oce_tracer_codes(2)     = 201

  END SUBROUTINE
!-------------------------------------------------------------------------  

  SUBROUTINE init_oce_config()
    oce_config%tracer_names(1)     = 'T'
    oce_config%tracer_longnames(1) = 'potential temperature'
    oce_config%tracer_units(1)     = 'deg C'
    oce_config%tracer_codes(1)     = 200
    oce_config%tracer_tags(1)      = '_'//TRIM(oce_config%tracer_names(1))

    oce_config%tracer_names(2)     = 'S'
    oce_config%tracer_longnames(2) = 'salinity'
    oce_config%tracer_units(2)     = 'psu'
    oce_config%tracer_codes(2)     = 201
    oce_config%tracer_tags(2)      = '_'//TRIM(oce_config%tracer_names(2))
  END SUBROUTINE
  FUNCTION is_initial_timestep(timestep)
    INTEGER :: timestep
    LOGICAL is_initial_timestep

    IF (timestep == 1 .AND. .NOT. is_restart_run()) THEN
      is_initial_timestep = .TRUE.
    ELSE
      is_initial_timestep = .FALSE.
    END IF
  END FUNCTION is_initial_timestep

END MODULE mo_oce_state
