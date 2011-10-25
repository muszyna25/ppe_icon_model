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
  USE mo_parallel_config,     ONLY: nproma
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    &                               success, max_char_length, min_rledge, min_rlcell,  &
    &                               min_rlvert, min_rlvert_int,                        &
    &                               full_coriolis, beta_plane_coriolis,                &
    &                               f_plane_coriolis, zero_coriolis
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, t_ref, s_ref,          &
    &                               CORIOLIS_TYPE, basin_center_lat, basin_height_deg
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: n_dom
  USE mo_ext_data,            ONLY: t_external_data
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, t_cartesian_coordinates,      &
    &                               t_geographical_coordinates, vector_product, &
    &                               arc_length
  USE mo_math_constants,      ONLY: pi, deg2rad
  USE mo_physical_constants,  ONLY: re, omega
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: add_var,                  &
    &                               new_var_list,             &
    &                               delete_var_list,          &
    &                               default_var_list_settings,&
    &                               add_ref
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants


  IMPLICIT NONE
  PRIVATE

! !VERSION CONTROL:
  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !public interface
  !
  ! subroutines
  PUBLIC :: init_ho_base
 !PUBLIC :: fill_ho_base
  PUBLIC :: construct_hydro_ocean_base
  PUBLIC :: destruct_hydro_ocean_base
  PUBLIC :: construct_hydro_ocean_state
  PUBLIC :: destruct_hydro_ocean_state
  PUBLIC :: set_lateral_boundary_values
  PUBLIC :: init_scalar_product_base
  PUBLIC :: init_geo_factors_base
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
  INTEGER, PRIVATE :: i_cell_type=3

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



    ! The following two arrays are required for the reconstruction process that 
    ! is used within the ocean model. Once the new version is implemented this 
    ! could eventually be shifted to the edge/vertex datatypes. It is currently 
    ! placed here to reduce interference with the atmospheric code (P.K.).
    !
    ! Vector pointing from cell circumcenter to edge midpoint. In the associated 
    ! cell2edge_weight-array the cell2edge_vec is multiplied by some other geometric 
    ! quantities (edge-length, cell area). The weight is used in the reconstruction 
    ! the vector is used in the transposed reconstruction.
    ! index=1,nproma, index2=1,nblks_c, index3=1,3
    ! other choice would be index2=1,nblks_e, index3=1,2
    ! Eventually switch to other second indexing if this is more appropriate

    ! Vector pointing from vertex (dual center) to midpoint of dual edge 
    ! (/= midpoint of primal edge).
    ! In the associated vertex2dualedge_mid_weight-array the vertex2dualedge_mid_vec 
    ! is multiplied by some other geometric quantities (dual edge-length, dual cell 
    ! area). The weight is used in the reconstruction the vector is used in the 
    ! transposed reconstruction.
    ! index=1,nproma, index2=1,nblks_v, index3=1,6
    ! other choice index2=1,nblks_e, index3=1,2
    ! Eventually switch to other second indexing if this is more appropriate
    ! new constructs for mimetic core:
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc(:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2cell_coeff_cc_t(:,:,:)


    !REAL(wp),                      ALLOCATABLE :: edge2vert_coeff(:,:,:,:)
    !TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_t(:,:,:)

    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc(:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_coeff_cc_t(:,:,:)
    TYPE(t_cartesian_coordinates), ALLOCATABLE :: edge2vert_vector_cc(:,:,:)

    REAL(wp), ALLOCATABLE :: fixed_vol_norm(:,:)
    REAL(wp), ALLOCATABLE :: variable_vol_norm(:,:,:)
    REAL(wp), ALLOCATABLE :: variable_dual_vol_norm(:,:,:)


    ! Location of midpoint of dual edge
!!$    TYPE(t_geographical_coordinates), ALLOCATABLE :: mid_dual_edge(:,:)
    ! Cartesian distance from vertex1 to vertex2 via dual edge midpoint
    REAL(wp), ALLOCATABLE :: dist_cell2edge(:,:,:)


    ! Data Structures from mo_interpolation:

    ! factor for divergence (nproma,i_cell_type,nblks_c)
    REAL(wp), ALLOCATABLE :: geofac_div(:,:,:)
    ! factor for quad-cell divergence (nproma,4,nblks_e)
    REAL(wp), ALLOCATABLE :: geofac_qdiv(:,:,:)
    ! factor for divergence (nproma,9-i_cell_type,nblks_v)
    REAL(wp), ALLOCATABLE :: geofac_rot(:,:,:)
    ! factor for nabla2-scalar (nproma,i_cell_type+1,nblks_c)
    REAL(wp), ALLOCATABLE :: geofac_n2s(:,:,:)
    ! factor for Green-Gauss gradient (nproma,4,nblks_c,2)
    !REAL(wp), ALLOCATABLE :: geofac_grg(:,:,:,:) !  not used


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
      &  ptp_vn(:,:,:)         ,& ! normal velocity after mapping P^T P
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_pred(:,:,:)        ,& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_impl_vert_diff(:,:,:),& ! predicted normal velocity vector at edges.
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
      CALL new_var_list(ocean_var_list, listname)
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
    INTEGER :: nblks_c, nblks_e, nblks_v, n_zlvp, n_zlvm, ie
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

    ! 
    ! arrays that are required for #slo OLD# reconstruction
    !
    ALLOCATE(v_base%dist_cell2edge(nproma,nblks_e,2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating dist_cell2edge failed')
    ENDIF

    !
    ! arrays that are required for setting up the scalar product
    !
    !coefficients for edge to cell mapping, one half of the scalar product.
    !Dimension: nproma,nblks_c encode number of cells, 1:3 corresponds to number
    !of edges per cell, 1:2 is for u and v component of cell vector
!     ALLOCATE(v_base%edge2cell_coeff(nproma,nblks_c,1:3, 1:2),STAT=ist)
!     IF (ist /= SUCCESS) THEN
!       CALL finish (routine,'allocating edge2cell_coeff failed')
!     ENDIF
    ALLOCATE(v_base%edge2cell_coeff_cc(nproma,nblks_c,1:3),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating edge2cell_coeff_cc failed')
    ENDIF

    !coefficients for transposed of edge to cell mapping, second half of the scalar product.
    !Dimension: nproma,nblks_e encode number of edges, 1:2 is for cell neighbors of an edge
    ALLOCATE(v_base%edge2cell_coeff_cc_t(nproma,nblks_e,1:2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating transposed edge2cell_coeff failed')
    ENDIF

    !
    !coefficients for edge to vertex mapping.
    !
    !Dimension: nproma,nblks_v encode number of vertices, 
    !1:6 is number of edges of a vertex,
    !1:2 is for u and v component of vertex vector
    ALLOCATE(v_base%edge2vert_coeff_cc(nproma,nblks_v,1:6),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(v_base%edge2vert_coeff_cc_t(nproma,nblks_e,1:2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(v_base%edge2vert_vector_cc(nproma,nblks_v,1:6),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating edge2vert_vector failed')
    ENDIF

    !
    !normalizing factors for edge to cell mapping.
    !
    !Either by fixed volume or by variable one taking the surface elevation
    !into account. The later one depends on time and space.
    ALLOCATE(v_base%fixed_vol_norm(nproma,nblks_c),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating fixed_vol_norm failed')
    ENDIF
    ALLOCATE(v_base%variable_vol_norm(nproma,nblks_c,1:3),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating variable_vol_norm failed')
    ENDIF

    ALLOCATE(v_base%variable_dual_vol_norm(nproma,nblks_v,1:6),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating variable_dual_vol_norm failed')
    ENDIF

    !
    ! Init geometrical factors
    !
    ! This part is copied from the interpolation modules.
    ! Interpolation module is no longer used in ocean model
    ! but some of the grid related information is stored in interpolation type.

    ALLOCATE (v_base%geofac_div(nproma, i_cell_type, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation for geofac_div failed')
    ENDIF

    ALLOCATE (v_base%geofac_qdiv(nproma, 4, nblks_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation for geofac_qdiv failed')
    ENDIF

    SELECT CASE (i_cell_type)
    CASE(3)
      ALLOCATE (v_base%geofac_rot(nproma, 6, nblks_v), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine,'allocation for geofac_rot failed')
      ENDIF
    CASE(6)
      ALLOCATE (v_base%geofac_rot(nproma, 4, nblks_e), STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish (routine,'allocation for geofac_rot failed')
      ENDIF
    END SELECT
    ALLOCATE (v_base%geofac_n2s(nproma, i_cell_type+1, nblks_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocation for geofac_n2s failed')
    ENDIF

!   ALLOCATE (v_base%geofac_grg(nproma, i_cell_type+1, nblks_c, 2), STAT=ist )
!   IF (ist /= SUCCESS) THEN
!     CALL finish (routine,'allocation for geofac_grg failed')
!   ENDIF

    v_base%del_zlev_m = 0._wp
    v_base%del_zlev_i = 0._wp
    v_base%zlev_m     = 0._wp
    v_base%zlev_i     = 0._wp

    v_base%lsm_oce_c = 0
    v_base%lsm_oce_e = 0
    v_base%dolic_c = 0
    v_base%dolic_e = 0

    v_base%wet_c = 0.0_wp
    v_base%wet_e = 0.0_wp

    DO ie = 1,3
      v_base%edge2cell_coeff_cc%x(ie)   = 0._wp
      v_base%edge2cell_coeff_cc_t%x(ie) = 0._wp
      v_base%edge2vert_coeff_cc%x(ie)   = 0._wp
      v_base%edge2vert_coeff_cc_t%x(ie) = 0._wp
      v_base%edge2vert_vector_cc%x(ie)  = 0._wp
    END DO

    v_base%fixed_vol_norm         = 0._wp
    v_base%variable_vol_norm      = 0._wp
    v_base%variable_dual_vol_norm = 0._wp

    v_base%dist_cell2edge = 0._wp

    v_base%geofac_div   = 0._wp
    v_base%geofac_qdiv  = 0._wp
    v_base%geofac_rot   = 0._wp
    v_base%geofac_n2s   = 0._wp

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

    !reconstrcuted velocity in cartesian coordinates
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
!
!
!>
!! Sbr set boundary values for velocity field.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!
!
  SUBROUTINE  set_lateral_boundary_values( p_patch, vn )

    TYPE(t_patch), INTENT(in) :: p_patch
    REAL(wp)                  :: vn(:,:,:)

    ! local variables
    INTEGER :: jb, je, jk !, il_v, ib_v, ie, iv_ctr, il_e, ib_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: slev,elev 
!!$    CHARACTER(len=max_char_length), PARAMETER :: &
!!$      &      routine = 'mo_oce_state: set_lateral_boundary_values'

!---------------------------------------------------------------

! blocking
i_startblk_e = p_patch%edges%start_blk(1,1)
i_endblk_e   = p_patch%edges%end_blk(min_rledge,1)
rl_start_e   = 1
rl_end_e     = min_rledge
slev         = 1
elev         = n_zlev

DO jb=i_startblk_e, i_endblk_e
  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
    DO jk = slev, elev
      DO je= i_startidx_e, i_endidx_e
        IF ( v_base%lsm_oce_e(je,jk,jb) >= boundary ) THEN
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
  SUBROUTINE init_ho_base( p_patch, p_ext_data, v_base )

    TYPE(t_patch),            INTENT(IN)       :: p_patch
    TYPE(t_external_data),    INTENT(IN)       :: p_ext_data
    TYPE(t_hydro_ocean_base), INTENT(INOUT)    :: v_base

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &        routine = 'mo_oce_state:init_ho_base'

    INTEGER :: jb, jc, je, jk, ji, nblks_c, nblks_e, npromz_c, npromz_e
    INTEGER :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: noct1_c, noct1_e, noctb_e, nocsb_c, noclb_c, inolsm, nowet_c
    INTEGER :: nolnd_c(n_zlev), nosea_c(n_zlev), nogllnd_c, noglsea_c
    INTEGER :: nolnd_e(n_zlev), nosea_e(n_zlev), nogllnd_e, noglsea_e
    INTEGER :: nobnd_e(n_zlev), nosbd_c(n_zlev), nolbd_c(n_zlev)
    INTEGER :: noglbnd_e, noglsbd_c, nogllbd_c
    INTEGER :: iic1, ibc1, iic2, ibc2, idxe, ible
    INTEGER :: n_zlvp, n_zlvm
    INTEGER :: ctr, jiter, niter
    REAL(wp):: perc_lnd_c(n_zlev), perc_gllnd_c
    REAL(wp):: perc_lnd_e(n_zlev), perc_gllnd_e

    REAL(wp) :: z_sync_c(nproma,p_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,p_patch%nblks_e)
    !REAL(wp) :: z_sync_v(nproma,p_patch%nblks_v)


    !-----------------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    z_sync_c(:,:) = 0.0_wp
    z_sync_e(:,:) = 0.0_wp

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
    noctb_e = 0
    nocsb_c = 0
    noclb_c = 0

    ! coordinate surfaces - n_zlev z-levels:
    ZLEVEL_LOOP: DO jk = 1, n_zlev

      !-----------------------------
      ! cells
      !  - values for BOUNDARY set below

      nolnd_c(jk)=0
      nosea_c(jk)=0

      !i_startblk = p_patch%cells%start_blk(1,1)
      !DO jb = i_startblk, nblks_c

      DO jb = 1, nblks_c

        !CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
        !  &                i_startidx, i_endidx, 1)

        i_endidx=nproma
        IF (jb==nblks_c) i_endidx=npromz_c

        !-----------------------------
        ! set dolic and wet grid points:
        !  - if bathymetry is deeper than or equal to the coordinate surface (zlev_m)
        !    then grid point is wet; dolic is in that level

        DO je = 1, i_endidx

          !  surface level of lsm defined by gridgenerator, not the current bathymetry
          !  read in from ext_data
          IF (jk == 1) THEN

            ! counts sea cells from lsm
            IF (p_ext_data%oce%lsm_ctr_c(je,jb) <= -1) THEN
              nosea_c(jk) = nosea_c(jk)+1
              v_base%dolic_c(je,jb) = jk
            ELSE IF (p_ext_data%oce%lsm_ctr_c(je,jb) >=  1) THEN
              nolnd_c(jk) = nolnd_c(jk)+1
            ELSE ! 0 not defined
              STOP ' lsm_ctr_c = 0'
            END IF

            ! counts sea points from bathymetry
            IF (p_ext_data%oce%bathymetry_c(je,jb) <= -v_base%zlev_m(jk)) &
              &   noct1_c = noct1_c+1

          !  second level of lsm and dolic defined by surface level, not the current bathymetry
          ELSE IF (jk == 2) THEN

            ! dependent on jk-1:
            v_base%lsm_oce_c(je,jk,jb) = v_base%lsm_oce_c(je,jk-1,jb)
            IF (v_base%lsm_oce_c(je,jk,jb) <= -1) THEN
              nosea_c(jk) = nosea_c(jk)+1
              v_base%dolic_c(je,jb) = jk
            ELSE IF (v_base%lsm_oce_c(je,jk-1,jb) >=  1) THEN
              nolnd_c(jk) = nolnd_c(jk)+1
            ELSE ! 0 not defined
              STOP ' lsm_oce_c = 0'
            END IF

          ELSE  ! jk>2

            IF (p_ext_data%oce%bathymetry_c(je,jb) <= -v_base%zlev_m(jk)) THEN
              nosea_c(jk)=nosea_c(jk)+1
              v_base%lsm_oce_c(je,jk,jb) = SEA
              v_base%dolic_c(je,jb) = jk
            ELSE IF (p_ext_data%oce%bathymetry_c(je,jb)>-v_base%zlev_m(jk)) THEN
              nolnd_c(jk)=nolnd_c(jk)+1
              v_base%lsm_oce_c(je,jk,jb) = LAND
            END IF

          END IF

        END DO

      END DO

      !  percentage of land area per level and global value
      inolsm = nolnd_c(jk) + nosea_c(jk)
      IF (inolsm == 0 ) THEN
        IF (jk == 1 ) CALL message (TRIM(routine), 'WARNING - number of cell points is zero?')
        perc_lnd_c(jk) = 0.0_wp
      ELSE
        perc_lnd_c(jk) = REAL(nolnd_c(jk),wp)/REAL(nosea_c(jk)+nolnd_c(jk),wp)*100.0_wp
        nogllnd_c = nogllnd_c + nolnd_c(jk)
        noglsea_c = noglsea_c + nosea_c(jk)
      END IF


      !-----------------------------
      ! edges
      !  - values for BOUNDARY set below, LAND, SEA only

      nolnd_e(jk)=0
      nosea_e(jk)=0

      DO jb = 1, nblks_e

        i_endidx=nproma
        IF (jb==nblks_e) i_endidx=npromz_e

        DO je = 1, i_endidx

          !  surface level of lsm and dolic defined by gridgenerator, not the current bathymetry
          IF (jk == 1) THEN

            ! count and define sea edges from lsm - boundary edges are counted as land
            IF (v_base%lsm_oce_e(je,jk,jb) == -2 ) THEN
              nosea_e(jk)=nosea_e(jk)+1
              v_base%dolic_e(je,jb) = jk
            ELSE
              nolnd_e(jk)=nolnd_e(jk)+1
            END IF

            ! counts sea points from bathymetry
            IF (p_ext_data%oce%bathymetry_e(je,jb) <= -v_base%zlev_m(jk)) THEN
              noct1_e = noct1_e+1
            ENDIF

          !  second level of lsm and dolic defined by surface level, not the current bathymetry
          ELSE IF (jk == 2) THEN

            ! dependent on jk-1:
            v_base%lsm_oce_e(je,jk,jb) = v_base%lsm_oce_e(je,jk-1,jb)
            IF (v_base%lsm_oce_e(je,jk,jb) == -2 ) THEN
              nosea_e(jk)=nosea_e(jk)+1
              v_base%dolic_e(je,jb) = jk
            ELSE
              nolnd_e(jk)=nolnd_e(jk)+1
            END IF

          ELSE  ! jk>2

            IF (p_ext_data%oce%bathymetry_e(je,jb) <= -v_base%zlev_m(jk)) THEN
              nosea_e(jk)=nosea_e(jk)+1
              v_base%lsm_oce_e(je,jk,jb) = SEA
              v_base%dolic_e(je,jb) =jk
            ELSE IF (p_ext_data%oce%bathymetry_e(je,jb)>-v_base%zlev_m(jk)) THEN
              nolnd_e(jk)=nolnd_e(jk)+1
              v_base%lsm_oce_e(je,jk,jb) = LAND
            END IF

          END IF

        END DO

      END DO

      !  percentage of land area per level and global value
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
      ! set values for BOUNDARY at edges (get values of neighbouring cells)
      !  - if the two corresponding cells are differing then edge is BOUNDARY
      !    (they are not both LAND or SEA)
      !  - done for jk>2 only, checks for read lsm in jk=1

      nobnd_e(jk)=0

      rl_start = 1           !  #slo# - cannot run with holes on land in grid
      rl_end = min_rledge

      ! values for the blocking
      i_startblk = p_patch%edges%start_blk(rl_start,1)
      i_endblk   = p_patch%edges%end_blk(rl_end,1)
      !
      ! loop through all patch edges
      BLK_LOOP_E: DO jb = i_startblk, i_endblk

        CALL get_indices_e  &
          &  (p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        IDX_LOOP_E: DO je =  i_startidx, i_endidx

          ! Get indices/blks of cells 1 and 2 adjacent to edge (je,jb)
          ! #slo# 2011-05-17:
          !  - set all layers except for the two surface layers which are set by gridgen
          !  - count numbers for all layers
          !  - check number for surface layer
          iic1 = p_patch%edges%cell_idx(je,jb,1)
          ibc1 = p_patch%edges%cell_blk(je,jb,1)
          iic2 = p_patch%edges%cell_idx(je,jb,2)
          ibc2 = p_patch%edges%cell_blk(je,jb,2)
          !
          ! cells may have -2, -1, 1, 2 for sea, sea_boundary, land_boundary, land:
          IF (jk == 1) THEN
            ! counts number of boundaries for jk=1
            IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) < 0)  .and.   &
              &  (v_base%lsm_oce_c(iic2,jk,ibc2) > 0) )        &
              &   noctb_e=noctb_e + 1
            IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) > 0)  .and.   &
              &  (v_base%lsm_oce_c(iic2,jk,ibc2) < 0) )        &
              &   noctb_e=noctb_e + 1
          END IF  !  jk = 1
          IF (jk > 2) THEN
            IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) < 0)  .and.   &
              &  (v_base%lsm_oce_c(iic2,jk,ibc2) > 0) )        &
              &   v_base%lsm_oce_e(je,jk,jb) = BOUNDARY
            IF ( (v_base%lsm_oce_c(iic1,jk,ibc1) > 0)  .and.   &
              &  (v_base%lsm_oce_c(iic2,jk,ibc2) < 0) )        &
              &   v_base%lsm_oce_e(je,jk,jb) = BOUNDARY
          END IF  !  jk > 2
          IF ( v_base%lsm_oce_e(je,jk,jb) == BOUNDARY )      &
            &  nobnd_e(jk)=nobnd_e(jk)+1


        END DO IDX_LOOP_E

      END DO BLK_LOOP_E

      !-----------------------------
      ! set values for LAND_BOUNDARY and SEA_BOUNDARY at cells
      !  - get values of neighbouring edges
      !  - if one of 3 edges of a sea-cell is BOUNDARY then cell is SEA_BOUNDARY
      !  - if one of 3 edges of a land-cell is BOUNDARY then cell is LAND_BOUNDARY
      !  - done for jk>2 only, checks for read lsm in jk=1

      nosbd_c(jk)=0
      nolbd_c(jk)=0

      rl_start = 1           !  #slo# - cannot run with holes on land in grid
      rl_end = min_rlcell

      ! values for the blocking
      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,1)
      !
      ! loop through all patch cells
      BLK_LOOP_C: DO jb = i_startblk, i_endblk

        CALL get_indices_c  &
          &  (p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        IDX_LOOP_C: DO jc =  i_startidx, i_endidx

          ! #slo# 2011-05-17
          !  - set all layers except for the two surface layers which are set by gridgen
          !  - count numbers for all layers
          !  - check number for surface layer

          ! sea points
          IF (v_base%lsm_oce_c(jc,jk,jb) < 0) THEN
            DO ji = 1, 3
              ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
              idxe = p_patch%cells%edge_idx(jc,jb,ji)
              ible = p_patch%cells%edge_blk(jc,jb,ji)
              ! if one of lsm_e is boundary then lsm_c is sea_boundary
              IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY .AND. jk  > 2) &
                &  v_base%lsm_oce_c(jc,jk,jb) = SEA_BOUNDARY
              ! counts number of sea-boundaries for jk=1 - only one boundary is allowed
              IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY .AND. jk == 1) &
                &  nocsb_c=nocsb_c + 1

            END DO
            IF ( v_base%lsm_oce_c(jc,jk,jb) == SEA_BOUNDARY )  &
              &  nosbd_c(jk)=nosbd_c(jk)+1
          END IF  !  lsm_c < 0

          ! land points
          IF (v_base%lsm_oce_c(jc,jk,jb) > 0) THEN

            DO ji = 1, 3
              ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
              idxe = p_patch%cells%edge_idx(jc,jb,ji)
              ible = p_patch%cells%edge_blk(jc,jb,ji)
              ! if one of lsm_e is boundary then lsm_c is land_boundary
              IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY .AND. jk  > 2 ) &
                &  v_base%lsm_oce_c(jc,jk,jb) = LAND_BOUNDARY
              ! counts number of land-boundaries for jk=1 - one land cell may have 2 boundaries
              IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY .AND. jk == 1 ) THEN
                   noclb_c=noclb_c + 1
                   EXIT
              END IF
            END DO

            IF ( v_base%lsm_oce_c(jc,jk,jb) == LAND_BOUNDARY )   &
              &  nolbd_c(jk)=nolbd_c(jk)+1

          END IF  !  lsm_c > 0

        END DO  IDX_LOOP_C

      END DO  BLK_LOOP_C
      noglbnd_e = noglbnd_e + nobnd_e(jk)
      noglsbd_c = noglsbd_c + nosbd_c(jk)
      nogllbd_c = nogllbd_c + nolbd_c(jk)

    END DO ZLEVEL_LOOP

!---------------------------------------------------------------------------------------------
! preliminary - through all levels each wet point has at most one dry point as
! neighbour

    rl_start = 1           
    rl_end = min_rlcell
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)

    niter=10
    DO jiter=1,niter

      DO jk=1,n_zlev
        !
        ! loop through all patch cells 
        ctr=0
        DO jb = i_startblk, i_endblk
          CALL get_indices_c  &
            &  (p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          DO jc =  i_startidx, i_endidx
            nowet_c = 0
            IF (v_base%lsm_oce_c(jc,jk,jb) < 0) THEN
              DO ji = 1, 3
                ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
                idxe = p_patch%cells%edge_idx(jc,jb,ji)
                ible = p_patch%cells%edge_blk(jc,jb,ji)
                ! if one of lsm_e is boundary then lsm_c is sea_boundary
                ! counts number of sea-boundaries for jk=1 - only one boundary is allowed
                IF ( v_base%lsm_oce_e(idxe,jk,ible) == BOUNDARY) &
                  &  nowet_c=nowet_c + 1
              END DO
              !More than 1 wet edge -> set cell to land
              IF ( nowet_c >= 2 ) THEN 
                v_base%lsm_oce_c(jc,jk,jb)=LAND
                ctr = ctr+1
                DO ji = 1, 3
                  ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
                  idxe = p_patch%cells%edge_idx(jc,jb,ji)
                  ible = p_patch%cells%edge_blk(jc,jb,ji)
                  !set all edges to boundary
                  v_base%lsm_oce_e(idxe,jk,ible) = BOUNDARY
                END DO
!                 !get adjacent triangles
!                iic1 = p_patch%edges%cell_idx(je,jb,1)
!                ibc1 = p_patch%edges%cell_blk(je,jb,1)
!                iic2 = p_patch%edges%cell_idx(je,jb,2)
!                ibc2 = p_patch%edges%cell_blk(je,jb,2)
!                !triangle 1 is primary triangle
!                IF(iic1==jc.AND.ibc1==jb)THEN
!                  v_base%lsm_oce_c(iic1,jk,ibc1)=LAND_c
!                !triangle 1 is primary triangle
!                ELSEIF(iic1==jc.AND.ibc1==jb)THEN
!                  v_base%lsm_oce_c(iic2,jk,ibc2)=LAND_c
!                ENDIF
              ENDIF
            END IF
          END DO
        END DO
        write(*,*)'triangles with 2 land edges present at jiter=',jiter,jk,ctr 
      END DO
    END DO  ! jiter

!---------------------------------------------------------------------------------------------
! preliminary - run through whole zlevel_loop once more - after correction in
! all levels

    nogllnd_c = 0
    noglsea_c = 0
    nogllnd_e = 0
    noglsea_e = 0
    noglbnd_e = 0
    noglsbd_c = 0
    nogllbd_c = 0

    ! set once more after jiter-correction
    v_base%dolic_c = 0
    v_base%dolic_e = 0

    ! 2011-10-24: second loop for edges, dolic, boundaries, diagnosis and output
    !  - using lsm_oce_c after jiter-correction as input
    !  - (1) set land and sea values at cells <0 (sea) and >0 (land) - no boundaries
    !  - (2) set land and sea values at edges including boundaries
    !  - (3) set land and sea boundary values (-1 = SEA_BOUNDARY, 1=LAND_BOUNDARY)
    !
    ZLEVEL_LOOP_cor: DO jk = 1, n_zlev

      !-----------------------------
      ! (1) set wet grid points and dolic at cells:
      !  - values for BOUNDARY set below

      nolnd_c(jk)=0
      nosea_c(jk)=0

      DO jb = 1, nblks_c

        !CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
        !  &                i_startidx, i_endidx, 1)

        i_endidx=nproma
        IF (jb==nblks_c) i_endidx=npromz_c

        DO je = 1, i_endidx
          IF (v_base%lsm_oce_c(je,jk,jb) <= SEA_BOUNDARY) THEN
            nosea_c(jk)=nosea_c(jk)+1
            v_base%dolic_c(je,jb) = jk
          ELSE
            ! -after correction: all other grid points are set to dry
            v_base%lsm_oce_c(je,jk,jb) = LAND
            nolnd_c(jk)=nolnd_c(jk)+1
          END IF
        END DO

      END DO

      !  percentage of land area per level and global value 
      !   - here: nosea/nolnd include boundaries
      inolsm = nolnd_c(jk) + nosea_c(jk)
      IF (inolsm == 0 ) THEN
        IF (jk == 1 ) CALL message (TRIM(routine), 'WARNING - number of cell points is zero?')
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
      ! loop through all patch edges
      DO jb = i_startblk, i_endblk

        CALL get_indices_e  &
          &  (p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO je =  i_startidx, i_endidx

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

        END DO

      END DO

      !  percentage of land area per level and global value
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
      DO jb = i_startblk, i_endblk

        CALL get_indices_c  &
          &  (p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc =  i_startidx, i_endidx

          ! #slo# 2011-10-25
          !  - set and count all layers
          !  - count of surface layer from grid-generator done in first zlevel_loop

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
          IF (v_base%lsm_oce_c(jc,jk,jb) > 0) THEN

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
      noglbnd_e = noglbnd_e + nobnd_e(jk)
      noglsbd_c = noglsbd_c + nosbd_c(jk)
      nogllbd_c = nogllbd_c + nolbd_c(jk)

    END DO ZLEVEL_LOOP_cor

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
    ! depth of the surface level dzlev_m(1) in ocean_ctl (must not be an error)
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

    ! Boundary warnings in case of inconsistency in read land-sea-mask
    IF ( nobnd_e(1) /= noctb_e ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        &   'WARNING - surface boundary edges read = ',nobnd_e(1), &
        &   ' - calculated from bathymetry = ',noctb_e
      CALL message(routine, TRIM(message_text))
    END IF
    IF ( nosbd_c(1) /= nocsb_c ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        &   'WARNING - surface sea-boundary cells read = ',nosbd_c(1), &
        &   ' - calculated from bathymetry = ',nocsb_c
      CALL message(routine, TRIM(message_text))
    END IF
    IF ( nolbd_c(1) /= noclb_c ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        &   'WARNING - surface land-boundary cells read = ',nolbd_c(1), &
        &   ' - calculated from bathymetry = ',noclb_c
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

    ! synchronize all elements of v_base:

    DO jk = 1, n_zlev

      z_sync_c(:,:) =  REAL(v_base%lsm_oce_c(:,jk,:),wp)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      v_base%lsm_oce_c(:,jk,:) = INT(z_sync_c(:,:))


      z_sync_c(:,:) =  v_base%wet_c(:,jk,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      v_base%wet_c(:,jk,:) = z_sync_c(:,:)

      z_sync_e(:,:) =  REAL(v_base%lsm_oce_e(:,jk,:),wp)
      CALL sync_patch_array(SYNC_e, p_patch, z_sync_e(:,:))
      v_base%lsm_oce_e(:,jk,:) = INT(z_sync_e(:,:))

      z_sync_e(:,:) =  v_base%wet_e(:,jk,:)
      CALL sync_patch_array(SYNC_e, p_patch, z_sync_e(:,:))
      v_base%wet_e(:,jk,:) = z_sync_e(:,:)

    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_ho_base

  !-------------------------------------------------------------------------
  !
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !! @par Revision History
  !!  developed by Peter Korn, MPI-M  2010-09
  !!  Modification by Stephan Lorenz, 2010-11
  !!
  SUBROUTINE init_scalar_product_base( p_patch, v_base)

    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout)    :: p_patch
    TYPE(t_hydro_ocean_base), INTENT(INOUT) :: v_base

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    & routine = ('mo_oce_state:init_scalar_product_base')

    INTEGER, PARAMETER :: no_cell_edges = 3
    INTEGER, PARAMETER :: no_vert_edges = 6
    INTEGER :: jb, je, jv, ie, ie_1, ie_2, icc
    INTEGER :: il_e,ib_e,k  
    INTEGER :: il_c1, ib_c1, il_c2, ib_c2
    INTEGER :: il_v1, il_v2, ib_v1, ib_v2
    !INTEGER :: jc, ile,ibe
    !INTEGER :: jl_v1, jl_v2!, jb_v1, jb_v2

    INTEGER :: iil_c1(no_cell_edges), iil_c2(no_cell_edges)
    INTEGER :: iib_c1(no_cell_edges), iib_c2(no_cell_edges)

    !INTEGER :: il_c1_e1, ib_c1_e1, il_c1_e2, ib_c1_e2!, il_c1_e3, ib_c1_e3
    !INTEGER :: il_c2_e1, il_c2_e2, il_c2_e3!, ib_c2_e3!ib_c2_e1, ib_c2_e2,

    INTEGER :: jil_c1, jib_c1,jil_c2, jib_c2
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx

    !REAL(wp) :: z_lon, z_lat, z_rlong, z_rlat, z_long_c, z_lat_c
    REAL(wp) :: z_tmp!z_twopi, z_longmax, z_longmin
    !REAL(wp) :: cell_edge_dist(no_cell_edges,2)
    !REAL(wp) :: cell_edge_distance
    REAL(wp) :: norm_c1_c2, norm_v1_v2, norm!, norm_v1_e0, norm_v2_e0
    !REAL(wp) :: dual_edge_length_v1(no_vert_edges)
    !REAL(wp) :: dual_edge_length_v2(no_vert_edges)
    REAL(wp) :: dual_edge_length(no_vert_edges)
!    REAL(wp) :: edge_length(no_cell_edges)
!     REAL(wp) :: vert_edge_dist_v1(no_vert_edges,2)
!     REAL(wp) :: vert_edge_dist_v2(no_vert_edges,2)
    REAL(wp) :: vert_edge_dist(no_vert_edges,2)!, new_lon,new_lat
!     REAL(wp) :: vert_dual_mid_dist_v1(no_vert_edges,2)
!     REAL(wp) :: vert_dual_mid_dist_v2(no_vert_edges,2)
    REAL(wp) :: vert_dual_mid_dist(no_vert_edges,2)
    !REAL(wp) :: vert2vert_dist
    REAL(wp) :: vert_edge_distance, vert_dual_mid_distance
    TYPE(t_geographical_coordinates) :: gc_mid_dual_edge(no_vert_edges)!gc_edge(no_cell_edges) 
    !TYPE(t_geographical_coordinates) :: gc_v0,gc_e0,gc_c0, gc_dual_edge(no_vert_edges)!gc_v1, gc_v2,
    !TYPE(t_geographical_coordinates) :: ll1, ll2!, ll3
    TYPE(t_cartesian_coordinates)    :: cc_dual_edge(no_vert_edges), cc_edge(no_cell_edges)
    TYPE(t_cartesian_coordinates)    :: xx1,xx2!,xx3
    !TYPE(t_cartesian_coordinates)    :: vert2vert_cc(nproma,p_patch%nblks_v,no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: vert1_midedge_cc(nproma,p_patch%nblks_v,no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: vert2_midedge_cc(nproma,p_patch%nblks_v,no_vert_edges)
    !TYPE(t_cartesian_coordinates)    :: normal_cc(no_cell_edges)!, norm_c0_e
    TYPE(t_cartesian_coordinates)    :: cell2cell_cc
    TYPE(t_cartesian_coordinates)    :: cc_e0, cc_c1,cc_c2,cc_v0!, cc_c0, cc_v1, cc_v2
    TYPE(t_cartesian_coordinates)    :: cv_c1_e0, cv_c2_e0, cv_c1_c2
    !TYPE(t_cartesian_coordinates)    :: cv_v1_e0, cv_v2_e0!, cv_v1_v2
    TYPE(t_cartesian_coordinates)    :: cc_mid_dual_edge(no_vert_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc
    TYPE(t_geographical_coordinates) :: gc1,gc2!recon_vec_gc(no_vert_edges), gc_tmp
    !TYPE(t_geographical_coordinates) :: recon_vec_gc(no_vert_edges), gc_tmp
    !TYPE(t_geographical_coordinates) :: normal_gc(no_cell_edges)
    !TYPE(t_cartesian_coordinates)    :: recon_vec_cc1(no_vert_edges), recon_vec_cc2(no_vert_edges)

    !TYPE(t_cartesian_coordinates)    :: cc_c1_e1, cc_c1_e2, cc_c1_e3
    !TYPE(t_cartesian_coordinates)    :: cc_c2_e1, cc_c2_e2, cc_c2_e3
    !TYPE(t_cartesian_coordinates)    :: cv_c1_e1, cv_c1_e2, cv_c1_e3
    !TYPE(t_cartesian_coordinates)    :: cv_c2_e1, cv_c2_e2, cv_c2_e3
    TYPE(t_cartesian_coordinates)    :: z_vec_c1(no_cell_edges),z_vec_c2(no_cell_edges)
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v1(no_vert_edges) 
    TYPE(t_cartesian_coordinates)    :: recon_vec_cc_v2(no_vert_edges)
    !TYPE(t_cartesian_coordinates)    :: vec_mid_dual_edge_2v1(no_vert_edges)
    !TYPE(t_cartesian_coordinates)    :: vec_mid_dual_edge_2v2(no_vert_edges)
    !REAL(wp) :: length
    REAL(wp) :: z_edge_length(no_cell_edges)!, z_e_length!, z_ce_dist
    REAL(wp) :: z_cell_edge_dist_c1(no_cell_edges,2),z_cell_edge_dist_c2(no_cell_edges,2)
    !REAL(wp) :: z_cell_edge_dist(no_cell_edges,2)
    REAL(wp) :: z_y!z_cell_edge_dist(no_cell_edges,2)

    !REAL(wp) :: p_c(3), r1(3,3), r2(3,3), rot_p_c(3), barlon, barlat
    !REAL(wp) :: edge_length_c1_e1, edge_length_c1_e2,edge_length_c1_e3
    !REAL(wp) :: edge_length_c2_e1, edge_length_c2_e2,edge_length_c2_e3

    REAL(wp) :: z_sync_c(nproma,p_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,p_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,p_patch%nblks_v)
    !REAL(wp) :: z_sync_v(nproma,p_patch%verts%end_blk(min_rlvert_int,1))
    !REAL(wp) :: z_sync_e(nproma,p_patch%edges%end_blk(min_rledge,1))  ! #slo# min_rledge_int ??
    !REAL(wp) :: z_sync_v(nproma,p_patch%nblks_int_v)
    !REAL(wp) :: z_sync_e(nproma,p_patch%nblks_int_e)

    LOGICAL, PARAMETER :: MID_POINT_DUAL_EDGE = .TRUE. !Please do not change this unless
                                                       !you are sure, you know what you do.
    LOGICAL, PARAMETER :: LARC_LENGTH = .FALSE.
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')

    rl_start = 1 
    rl_end = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)

     rl_start_e = 1
     rl_end_e = min_rledge

     i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
     i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

     !STEP 1: edge2cell and cell2edge coefficients
     EDGE_BLK_LOOP_PRIMAL: DO jb = i_startblk_e, i_endblk_e

       CALL get_indices_e(p_patch, jb,&
                        & i_startblk_e, i_endblk_e,&
                        & i_startidx_e, i_endidx_e,&
                        & rl_start_e, rl_end_e)

       EDGE_IDX_LOOP_PRIMAL: DO je =  i_startidx_e, i_endidx_e

         !Get indices of two adjacent triangles
         il_c1 = p_patch%edges%cell_idx(je,jb,1)
         ib_c1 = p_patch%edges%cell_blk(je,jb,1)
         il_c2 = p_patch%edges%cell_idx(je,jb,2)
         ib_c2 = p_patch%edges%cell_blk(je,jb,2)

         cc_c1 = gc2cc(p_patch%cells%center(il_c1, ib_c1))
         cc_c2 = gc2cc(p_patch%cells%center(il_c2, ib_c2))

         z_cell_edge_dist_c1 = 0.0_wp
         z_cell_edge_dist_c2 = 0.0_wp

         !normals in cell 1
         DO ie = 1, no_cell_edges

           !actual edges of cell c1
           iil_c1(ie) = p_patch%cells%edge_idx(il_c1,ib_c1,ie)
           iib_c1(ie) = p_patch%cells%edge_blk(il_c1,ib_c1,ie)

           cc_edge(ie)   = gc2cc(p_patch%edges%center(iil_c1(ie),iib_c1(ie)))

           !calculate edge length
           !get vertex indices adjacent to actual edge
           il_v1 = p_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),1)
           ib_v1 = p_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),1)
           il_v2 = p_patch%edges%vertex_idx(iil_c1(ie),iib_c1(ie),2)
           ib_v2 = p_patch%edges%vertex_blk(iil_c1(ie),iib_c1(ie),2)

           !get vertex positions
           xx1 = gc2cc(p_patch%verts%vertex(il_v1,ib_v1))
           xx2 = gc2cc(p_patch%verts%vertex(il_v2,ib_v2))
 
           IF(LARC_LENGTH)THEN
             norm=SQRT(SUM(xx1%x*xx1%x))
             xx1%x= xx1%x/norm
             norm=SQRT(SUM(xx2%x*xx2%x))
             xx2%x= xx2%x/norm
             z_edge_length(ie) = arc_length(xx2,xx1)
             !z_edge_length(ie) = p_patch%edges%primal_edge_length(iil_c1(ie),iib_c1(ie))/re
           ELSE
             z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
            !write(*,*)'length 3:', z_edge_length(ie),arc_length(xx1,xx2) 
           ENDIF

           !calculate cell-edge distance as half of cell-cell distance
           !get cell indices adjacent to actual edge
           jil_c1 = p_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),1)
           jib_c1 = p_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),1)
           jil_c2 = p_patch%edges%cell_idx(iil_c1(ie),iib_c1(ie),2)
           jib_c2 = p_patch%edges%cell_blk(iil_c1(ie),iib_c1(ie),2)

           !get cell positions
           xx1 = gc2cc(p_patch%cells%center(jil_c1,jib_c1))
           xx2 = gc2cc(p_patch%cells%center(jil_c2,jib_c2))

           IF(jil_c1==il_c1.AND.jib_c1==ib_c1)THEN
            k=1
           ELSEIF(jil_c2==il_c1.AND.jib_c2==ib_c1)THEN
            k=2
           ENDIF

           IF(LARC_LENGTH)THEN
             norm=SQRT(SUM(xx1%x*xx1%x))
             xx1%x= xx1%x/norm
             norm=SQRT(SUM(xx2%x*xx2%x))
             xx2%x= xx2%x/norm
             norm          = SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
             cc_edge(ie)%x = cc_edge(ie)%x/norm
             z_cell_edge_dist_c1(ie,1) = arc_length(cc_edge(ie),xx1)
             z_cell_edge_dist_c1(ie,2) = arc_length(cc_edge(ie),xx2)
             !z_cell_edge_dist_c1(ie,1) = p_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re 
             !z_cell_edge_dist_c1(ie,2) = p_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re 
           ELSE
             z_cell_edge_dist_c1(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
             z_cell_edge_dist_c1(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))  
             !write(*,*)'length 4',z_cell_edge_dist_c1(ie,1), z_cell_edge_dist_c1(ie,1),&
             !&p_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),1)/re 
           ENDIF
           v_base%dist_cell2edge(iil_c1(ie),iib_c1(ie),1) = z_cell_edge_dist_c1(ie,1)
           v_base%dist_cell2edge(iil_c1(ie),iib_c1(ie),2) = z_cell_edge_dist_c1(ie,2)

           z_vec_c1(ie)%x =  cc_edge(ie)%x - cc_c1%x     !p_patch%edges%primal_cart_normal(iil_c1(ie),iib_c1(ie))
           norm = SQRT(SUM( z_vec_c1(ie)%x* z_vec_c1(ie)%x))

           v_base%edge2cell_coeff_cc(il_c1,ib_c1,ie)%x&
           & = z_vec_c1(ie)%x*p_patch%cells%edge_orientation(il_c1,ib_c1,ie)*z_edge_length(ie)

!test
!z_vec_c1(ie)%x =z_vec_c1(ie)%x/norm
!write(*,*)'vec:normal:',z_vec_c1(ie)%x,p_patch%edges%primal_cart_normal(iil_c1(ie),iib_c1(ie))%x 
!---------
                 v_base%fixed_vol_norm(il_c1,ib_c1) = &
             &   v_base%fixed_vol_norm(il_c1,ib_c1) + 0.5_wp*norm*z_edge_length(ie)
           v_base%variable_vol_norm(il_c1,ib_c1,ie) = 0.5_wp*norm*z_edge_length(ie)

           !write(*,*)'edge length   :',z_edge_length(ie),p_patch%edges%primal_edge_length(iil_c1(ie),iib_c1(ie))/re
           !write(*,*)'cell-edge dist:', z_cell_edge_dist_c1(ie,k),p_patch%edges%edge_cell_length(iil_c1(ie),iib_c1(ie),k)/re
         END DO

         !normals in cell 2
         DO ie = 1, no_cell_edges

           !actual edges of cell c2
           iil_c2(ie) = p_patch%cells%edge_idx(il_c2,ib_c2,ie)
           iib_c2(ie) = p_patch%cells%edge_blk(il_c2,ib_c2,ie)

           cc_edge(ie) = gc2cc(p_patch%edges%center(iil_c2(ie),iib_c2(ie)))

           !calculate edge length
           !get vertex indices adjacent to actual edge
           il_v1 = p_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),1)
           ib_v1 = p_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),1)
           il_v2 = p_patch%edges%vertex_idx(iil_c2(ie),iib_c2(ie),2)
           ib_v2 = p_patch%edges%vertex_blk(iil_c2(ie),iib_c2(ie),2)

           !get vertex positions
           xx1 = gc2cc(p_patch%verts%vertex(il_v1,ib_v1))
           xx2 = gc2cc(p_patch%verts%vertex(il_v2,ib_v2))

           IF(LARC_LENGTH)THEN
             norm=SQRT(SUM(xx1%x*xx1%x))
             xx1%x= xx1%x/norm

             norm=SQRT(SUM(xx2%x*xx2%x))
             xx2%x= xx2%x/norm

             z_edge_length(ie) = arc_length(xx2,xx1)
             !z_edge_length(ie) = p_patch%edges%primal_edge_length(iil_c2(ie),iib_c2(ie))/re
              !write(*,*)'arc length',arc_length(xx2,xx1),z_edge_length(ie),SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
           ELSE
             z_edge_length(ie) = SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
           ENDIF
           !calculate cell-edge distance as half of cell-cell distance
           !get cell indices adjacent to actual edge
           jil_c1 = p_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),1)
           jib_c1 = p_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),1)
           jil_c2 = p_patch%edges%cell_idx(iil_c2(ie),iib_c2(ie),2)
           jib_c2 = p_patch%edges%cell_blk(iil_c2(ie),iib_c2(ie),2)

           !get cell positions
           xx1 = gc2cc(p_patch%cells%center(jil_c1,jib_c1))  
           xx2 = gc2cc(p_patch%cells%center(jil_c2,jib_c2))

           IF(jil_c1==il_c2.AND.jib_c1==ib_c2)THEN
             k=1
           ELSEIF(jil_c2==il_c2.AND.jib_c2==ib_c2)THEN
             k=2
           ENDIF  

           IF(LARC_LENGTH)THEN
             norm=SQRT(SUM(xx1%x*xx1%x))
             xx1%x= xx1%x/norm
             norm=SQRT(SUM(xx2%x*xx2%x))
             xx2%x= xx2%x/norm
             norm=SQRT(SUM(cc_edge(ie)%x*cc_edge(ie)%x))
             cc_edge(ie)%x =  cc_edge(ie)%x/norm
             z_cell_edge_dist_c2(ie,1) = arc_length(cc_edge(ie),xx1)
             z_cell_edge_dist_c2(ie,2) = arc_length(cc_edge(ie),xx2)
             !z_cell_edge_dist_c2(ie,1) = p_patch%edges%edge_cell_length(iil_c2(ie),iib_c2(ie),1)/re
             !z_cell_edge_dist_c2(ie,2) = p_patch%edges%edge_cell_length(iil_c2(ie),iib_c2(ie),2)/re
             !write(*,*)'arc length',0.5_wp*arc_length(xx2,xx1),z_cell_edge_dist_c2(ie,k),0.5_wp*SQRT(SUM((xx2%x-xx1%x)*(xx2%x-xx1%x)))
           ELSE 
             z_cell_edge_dist_c2(ie,1) = SQRT(SUM((cc_edge(ie)%x-xx1%x)*(cc_edge(ie)%x-xx1%x)))
             z_cell_edge_dist_c2(ie,2) = SQRT(SUM((cc_edge(ie)%x-xx2%x)*(cc_edge(ie)%x-xx2%x)))
           ENDIF
           v_base%dist_cell2edge(iil_c2(ie),iib_c2(ie),1) = z_cell_edge_dist_c1(ie,1)
           v_base%dist_cell2edge(iil_c2(ie),iib_c2(ie),2) = z_cell_edge_dist_c1(ie,2)

           z_vec_c2(ie)%x =  cc_edge(ie)%x - cc_c2%x  !p_patch%edges%primal_cart_normal(iil_c2(ie),iib_c2(ie))
           norm = SQRT(SUM( z_vec_c2(ie)%x* z_vec_c2(ie)%x))

           v_base%edge2cell_coeff_cc(il_c2,ib_c2,ie)%x&
             & = z_vec_c2(ie)%x*p_patch%cells%edge_orientation(il_c2,ib_c2,ie)*z_edge_length(ie)

           v_base%fixed_vol_norm(il_c2,ib_c2) &
             & = v_base%fixed_vol_norm(il_c2,ib_c2) + 0.5_wp*norm*z_edge_length(ie)

           v_base%variable_vol_norm(il_c2,ib_c2,ie) = 0.5_wp*norm*z_edge_length(ie)


         END DO
       END DO EDGE_IDX_LOOP_PRIMAL
     END DO EDGE_BLK_LOOP_PRIMAL
     !In he edge loop above each triangle is visisted three times. Since the "fixed_vol_norm" is
     !accumulated we correct its value here:
     v_base%fixed_vol_norm = v_base%fixed_vol_norm/3.0_wp

    rl_start = 1 
    rl_end = min_rledge

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk(rl_end,1)

    EDGE_BLK_LOOP_0: DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
      EDGE_IDX_LOOP_0: DO je =  i_startidx, i_endidx

        !Get indices of two adjacent triangles
        il_c1 = p_patch%edges%cell_idx(je,jb,1)
        ib_c1 = p_patch%edges%cell_blk(je,jb,1)
        il_c2 = p_patch%edges%cell_idx(je,jb,2)
        ib_c2 = p_patch%edges%cell_blk(je,jb,2)

        !cartesian coordinates of edge and neighbor cells on 1-sphere
        cc_e0 = gc2cc(p_patch%edges%center(je,jb))
        cc_c1 = gc2cc(p_patch%cells%center(il_c1,ib_c1))
        cc_c2 = gc2cc(p_patch%cells%center(il_c2,ib_c2))

        !cartesian vectors from:
        !cell 2 to cell 1, cell 1 to edge je and cell 2 to edge je
        cv_c1_c2%x = cc_c1%x - cc_c2%x
        cv_c1_e0%x = cc_e0%x - cc_c1%x
        cv_c2_e0%x = cc_e0%x - cc_c2%x

        IF(LARC_LENGTH)THEN
          norm=SQRT(SUM(cc_e0%x*cc_e0%x))
          cc_e0%x= cc_e0%x/norm
          norm=SQRT(SUM(cc_c1%x*cc_c1%x))
          cc_c1%x= cc_c1%x/norm
          norm=SQRT(SUM(cc_c2%x*cc_c2%x))
          cc_c2%x= cc_c2%x/norm
          norm_c1_c2 = arc_length(cc_c1, cc_c2)
          !norm_c1_c2 = p_patch%edges%dual_edge_length(je,jb)/re 
                       !p_patch%edges%edge_cell_length(je,jb,1)/re&
                      !&+p_patch%edges%edge_cell_length(je,jb,2)/re
        ELSE
          norm_c1_c2 = SQRT(SUM(cv_c1_e0%x*cv_c1_e0%x))+SQRT(SUM(cv_c2_e0%x*cv_c2_e0%x)) !SQRT(SUM(cv_c1_c2%x*cv_c1_c2%x))!
        ENDIF

        !Determine which edge of both of the two adjacent cells corresponds to the
        !actual edge "je". This information is used below for the edge-orientation.
        DO ie = 1, no_cell_edges
          IF (p_patch%cells%edge_idx(il_c1,ib_c1,ie) == je.AND.&
            & p_patch%cells%edge_blk(il_c1,ib_c1,ie) == jb) THEN
            ie_1 = ie
          END IF
          IF (p_patch%cells%edge_idx(il_c2,ib_c2,ie) == je.AND.&
            & p_patch%cells%edge_blk(il_c2,ib_c2,ie) == jb) THEN
            ie_2 = ie
          END IF
        END DO

        v_base%edge2cell_coeff_cc_t(je,jb,1)%x&
        & = cv_c1_e0%x * p_patch%cells%edge_orientation(il_c1,ib_c1,ie_1)/norm_c1_c2

        v_base%edge2cell_coeff_cc_t(je,jb,2)%x&
        & = cv_c2_e0%x * p_patch%cells%edge_orientation(il_c2,ib_c2,ie_2)/norm_c1_c2

      END DO EDGE_IDX_LOOP_0
    END DO EDGE_BLK_LOOP_0

    !------------------------------------------------------------------------------
    !STEP 2: edge2vert coefficients for dual grid
    !------------------------------------------------------------------------------

    rl_start = 1
 !  rl_end   = min_rlvert
    rl_end   = min_rlvert_int  ! inner part of decomposition only - no halo (!!)

    i_startblk = p_patch%verts%start_blk(rl_start,1)
    i_endblk   = p_patch%verts%end_blk(rl_end,1)

    VERT_BLK_LOOP: DO jb = i_startblk, i_endblk
      CALL get_indices_v(p_patch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
!     write(*,*) 'jblock, strt/end_idx:', jb, i_startidx, i_endidx
      VERT_IDX_LOOP: DO jv =  i_startidx, i_endidx

        ! current number of edges around vertex (5 or 6)
        cc_v0        = gc2cc(p_patch%verts%vertex(jv,jb))
 
        DO ie = 1, no_vert_edges  ! #slo# it_vertedges ??

          il_e = p_patch%verts%edge_idx(jv,jb,ie)
          ib_e = p_patch%verts%edge_blk(jv,jb,ie)

          ! #slo# - I assume all geographical coordinates are already synchronized

          cc_dual_edge(ie) = gc2cc(p_patch%edges%center(il_e,ib_e))
          !Parts of this code parrallels the implementation in the grid-generator
          !module "mo_geometry".
          !
          !1) determine normal vector from adjacent cell to adjacent cell
          !   in cartesian coordinate for moved dual cell
          !Get indices of two adjacent triangles
          il_c1 = p_patch%edges%cell_idx(il_e,ib_e,1)
          ib_c1 = p_patch%edges%cell_blk(il_e,ib_e,1)
          il_c2 = p_patch%edges%cell_idx(il_e,ib_e,2)
          ib_c2 = p_patch%edges%cell_blk(il_e,ib_e,2)

          xx1 = gc2cc(p_patch%cells%center(il_c1,ib_c1))
          norm=SQRT(SUM(xx1%x*xx1%x))
          xx1%x= xx1%x/norm

          xx2 = gc2cc(p_patch%cells%center(il_c2,ib_c2))
          norm=SQRT(SUM(xx2%x*xx2%x))
          xx2%x= xx2%x/norm

          cell2cell_cc%x       = xx2%x - xx1%x
          IF(LARC_LENGTH)THEN
            norm_c1_c2 = arc_length(xx1,xx2)
          ELSE
            norm_c1_c2 = SQRT(SUM(cell2cell_cc%x*cell2cell_cc%x)) 
          ENDIF 
          dual_edge_length(ie) = norm_c1_c2
!          cell2cell_cc%x       = cell2cell_cc%x/norm_c1_c2 

          IF(MID_POINT_DUAL_EDGE)THEN
            cc_mid_dual_edge(ie)%x = 0.5_wp*(xx2%x+xx1%x)
            gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie)) 

            IF(CORIOLIS_TYPE==full_coriolis)THEN
              p_patch%edges%f_e(il_e, ib_e) = 2._wp*omega*SIN(gc_mid_dual_edge(ie)%lat) 
            ELSEIF(CORIOLIS_TYPE==BETA_PLANE_CORIOLIS)THEN
              gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
              gc1%lon = 0.0_wp
              xx1=gc2cc(gc1)

              gc2%lat = gc_mid_dual_edge(ie)%lat!*deg2rad
              gc2%lon = 0.0_wp
              xx2=gc2cc(gc2)        
              z_y = re*arc_length(xx2,xx1)
 
              !z_y = ptr_patch%edges%center(je,jb)%lat - z_lat_basin_center
              p_patch%edges%f_e(il_e, ib_e) &
              &= 2.0_wp*omega*( sin(basin_center_lat * deg2rad)     &
              &   + (cos(basin_center_lat * deg2rad)/re)*z_y)
            ENDIF
          ELSE
            cc_mid_dual_edge(ie)%x = cc_dual_edge(ie)%x
            gc_mid_dual_edge(ie)   = cc2gc(cc_mid_dual_edge(ie)) 
          ENDIF

          !2) determine vector from adjacent vertex to adjacent vertex
          !   in cartesian coordinate for moved dual cell
          !Get indices of two adjacent vertices
          il_v1 = p_patch%edges%vertex_idx(il_e,ib_e,1)
          ib_v1 = p_patch%edges%vertex_blk(il_e,ib_e,1)
          il_v2 = p_patch%edges%vertex_idx(il_e,ib_e,2)
          ib_v2 = p_patch%edges%vertex_blk(il_e,ib_e,2)

          xx1 = gc2cc(p_patch%verts%vertex(il_v1,ib_v1))
          norm=SQRT(SUM(xx1%x*xx1%x))
          xx1%x= xx1%x/norm

          xx2 = gc2cc(p_patch%verts%vertex(il_v2,ib_v2))
          norm=SQRT(SUM(xx2%x*xx2%x))
          xx2%x= xx2%x/norm

          vert1_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx1%x
          vert2_midedge_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx2%x
          !vert2vert_cc(jv, jb, ie)%x = cc_mid_dual_edge(ie)%x - xx1%x !xx2%x - xx1%x
!           norm                       = SQRT(SUM(vert2vert_cc(jv,jb,ie)%x*vert2vert_cc(jv,jb,ie)%x))
!           vert2vert_cc(jv,jb,ie)%x   = vert2vert_cc(jv, jb, ie)%x/norm
          norm = SQRT(SUM(vert1_midedge_cc(jv,jb,ie)%x*vert1_midedge_cc(jv,jb,ie)%x))
          vert1_midedge_cc(jv, jb, ie)%x = vert1_midedge_cc(jv, jb, ie)%x/norm

          norm = SQRT(SUM(vert2_midedge_cc(jv,jb,ie)%x*vert2_midedge_cc(jv,jb,ie)%x))
          vert2_midedge_cc(jv, jb, ie)%x = vert2_midedge_cc(jv, jb, ie)%x/norm


          !calculate vertex edge distance 
          IF(LARC_LENGTH)THEN
            vert_edge_dist(ie,1) = arc_length (cc_dual_edge(ie), xx1) 
            vert_edge_dist(ie,2) = arc_length (cc_dual_edge(ie), xx2) 
            vert_dual_mid_dist(ie,1)= arc_length (cc_mid_dual_edge(ie), xx1) 
            vert_dual_mid_dist(ie,2)= arc_length (cc_mid_dual_edge(ie), xx2)
          ELSE
            vert_edge_dist(ie,1)&
            & = SQRT(SUM((cc_dual_edge(ie)%x - xx1%x)*(cc_dual_edge(ie)%x - xx1%x)))
            vert_edge_dist(ie,2)&
            & = SQRT(SUM((cc_dual_edge(ie)%x - xx2%x)*(cc_dual_edge(ie)%x - xx2%x)))
            vert_dual_mid_dist(ie,1)&
            & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx1%x)*(cc_mid_dual_edge(ie)%x - xx1%x)))
            vert_dual_mid_dist(ie,2)&
            & = SQRT(SUM((cc_mid_dual_edge(ie)%x - xx2%x)*(cc_mid_dual_edge(ie)%x - xx2%x)))
          ENDIF

          !calculate normal vector that is perpendicular to vertex-vertex- and edge position vector
          !If one uses the edge position vector this results in the moved primal normal. Later
          !edge position vector has to be replaced by the midpoint of the dual edge.
          recon_vec_cc_v1(ie) = vector_product(vert1_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
          norm                = SQRT(SUM(recon_vec_cc_v1(ie)%x*recon_vec_cc_v1(ie)%x))
          recon_vec_cc_v1(ie)%x = recon_vec_cc_v1(ie)%x/norm

          recon_vec_cc_v2(ie)   = vector_product(vert2_midedge_cc(jv,jb,ie), cc_mid_dual_edge(ie))
          norm                  = SQRT(SUM(recon_vec_cc_v2(ie)%x*recon_vec_cc_v2(ie)%x))
          recon_vec_cc_v2(ie)%x = recon_vec_cc_v2(ie)%x/norm

          !Fix orientation
          z_tmp = DOT_PRODUCT(recon_vec_cc_v1(ie)%x, p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
          IF (z_tmp <0._wp) recon_vec_cc_v1(ie)%x = -1._wp * recon_vec_cc_v1(ie)%x

          z_tmp = DOT_PRODUCT(recon_vec_cc_v2(ie)%x, p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
          IF (z_tmp <0._wp) recon_vec_cc_v2(ie)%x = -1._wp * recon_vec_cc_v2(ie)%x


          !write(*,*)'recon vec:primal 1:', il_e,ib_e,ie,recon_vec_cc%x, p_patch%edges%primal_cart_normal(il_e,ib_e)%x
          IF      ( (p_patch%edges%vertex_idx(il_e,ib_e,1) == jv) .and. &
                    (p_patch%edges%vertex_blk(il_e,ib_e,1) == jb) ) THEN

            vert_edge_distance     = vert_edge_dist(ie,1)
            vert_dual_mid_distance = vert_dual_mid_dist(ie,1)
            recon_vec_cc           = recon_vec_cc_v1(ie)

            v_base%edge2vert_vector_cc(jv,jb,ie)=&
            &vector_product(vert1_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
            z_tmp = DOT_PRODUCT(v_base%edge2vert_vector_cc(jv,jb,ie)%x,&
            &p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
            IF (z_tmp <0._wp) v_base%edge2vert_vector_cc(jv,jb,ie)%x&
            & = -1._wp * v_base%edge2vert_vector_cc(jv,jb,ie)%x


          ELSE IF ( (p_patch%edges%vertex_idx(il_e,ib_e,2) == jv) .and. &
                    (p_patch%edges%vertex_blk(il_e,ib_e,2) == jb) ) THEN

            vert_edge_distance     = vert_edge_dist(ie,2)
            vert_dual_mid_distance = vert_dual_mid_dist(ie,2)
            recon_vec_cc           = recon_vec_cc_v2(ie)

            v_base%edge2vert_vector_cc(jv,jb,ie)=&
            &vector_product(vert2_midedge_cc(jv, jb, ie), cc_mid_dual_edge(ie))
            z_tmp = DOT_PRODUCT(v_base%edge2vert_vector_cc(jv,jb,ie)%x,&
            & p_patch%edges%primal_cart_normal(il_e,ib_e)%x)
            IF (z_tmp <0._wp) v_base%edge2vert_vector_cc(jv,jb,ie)%x =&
            & -1._wp * v_base%edge2vert_vector_cc(jv,jb,ie)%x
          ELSE
            CALL message (TRIM(routine), 'WARNING - vert_edge_distance not found')
            write(*,'(a,7i5)') 'jv, jb, edge, p_patch%edges%vertex_idx/blk(il_e,ib_e,1-2)=', &
              &         jv, jb, ie, &
              &         p_patch%edges%vertex_idx(il_e,ib_e,1), &
              &         p_patch%edges%vertex_blk(il_e,ib_e,1), &
              &         p_patch%edges%vertex_idx(il_e,ib_e,2), &
              &         p_patch%edges%vertex_blk(il_e,ib_e,2)

          END IF

          v_base%variable_dual_vol_norm(jv,jb,ie) = &
            &                           0.5_wp*dual_edge_length(ie)*vert_dual_mid_distance
                                               !vert_edge_distance*dual_edge_length(ie)!

          v_base%edge2vert_coeff_cc(jv,jb,ie)%x   = &
            &                      recon_vec_cc%x*dual_edge_length(ie)*vert_dual_mid_distance

          norm_v1_v2 = SQRT(SUM(vert1_midedge_cc(jv, jb, ie)%x*vert1_midedge_cc(jv, jb, ie)%x))&
                    &+ SQRT(SUM(vert2_midedge_cc(jv, jb, ie)%x*vert2_midedge_cc(jv, jb, ie)%x))
          
          v_base%edge2vert_coeff_cc_t(il_e,ib_e,1)%x = vert1_midedge_cc(jv, jb, ie)%x * &
            &    ( p_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )
          
          v_base%edge2vert_coeff_cc_t(il_e,ib_e,2)%x = vert2_midedge_cc(jv, jb, ie)%x * &
            &    ( p_patch%edges%system_orientation(il_e,ib_e)/norm_v1_v2 )

        END DO

      ENDDO VERT_IDX_LOOP
    END DO VERT_BLK_LOOP

    ! synchronize all elements of v_base:

    ! synchronize elements on cells
    DO ie = 1, no_cell_edges

      DO icc = 1, 3

        z_sync_c(:,:) =  v_base%edge2cell_coeff_cc(:,:,ie)%x(icc)
        CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
        v_base%edge2cell_coeff_cc(:,:,ie)%x(icc) = z_sync_c(:,:)

      END DO

      z_sync_c(:,:) = v_base%variable_vol_norm(:,:,ie)
      CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
      v_base%variable_vol_norm(:,:,ie) = z_sync_c(:,:)
   
    END DO

    z_sync_c(:,:) = v_base%fixed_vol_norm(:,:)
    CALL sync_patch_array(SYNC_C, p_patch, z_sync_c(:,:))
    v_base%fixed_vol_norm(:,:) = z_sync_c(:,:)

    ! synchronize elements on edges
    DO ie = 1, 2

      DO icc = 1, 3

        z_sync_e(:,:) =  v_base%edge2vert_coeff_cc_t(:,:,ie)%x(icc)
        CALL sync_patch_array(SYNC_E, p_patch, z_sync_e(:,:))
        v_base%edge2vert_coeff_cc_t(:,:,ie)%x(icc) = z_sync_e(:,:)

        z_sync_e(:,:) =  v_base%edge2cell_coeff_cc_t(:,:,ie)%x(icc)
        CALL sync_patch_array(SYNC_E, p_patch, z_sync_e(:,:))
        v_base%edge2cell_coeff_cc_t(:,:,ie)%x(icc) = z_sync_e(:,:)

      END DO

    END DO

    DO ie = 1, no_vert_edges

      ! synchronize cartesian coordinates on vertices:
      DO icc = 1, 3

        z_sync_v(:,:) =  v_base%edge2vert_vector_cc(:,:,ie)%x(icc)
        CALL sync_patch_array(SYNC_V, p_patch, z_sync_v(:,:))
        v_base%edge2vert_vector_cc(:,:,ie)%x(icc) = z_sync_v(:,:)

        z_sync_v(:,:) = v_base%edge2vert_coeff_cc(:,:,ie)%x(icc)
        CALL sync_patch_array(SYNC_V, p_patch, z_sync_v(:,:))
        v_base%edge2vert_coeff_cc(:,:,ie)%x(icc) = z_sync_v(:,:)

      END DO

      z_sync_v(:,:) = v_base%variable_dual_vol_norm(:,:,ie)
      CALL sync_patch_array(SYNC_V, p_patch, z_sync_v(:,:))
      v_base%variable_dual_vol_norm(:,:,ie) = z_sync_v(:,:)

    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_scalar_product_base


!-------------------------------------------------------------------------
  !
  !
  !>
  !! Precomputes the geometrical factors used in the divergence, rotation.
  !!
  !! @par Revision History
  !!  developed by Guenther Zaengl, 2009-03-17
  !!  Modification by Almut Gassmann, 2009-12-19
  !!  - Vorticity is computed on quads in case of the hexagonal grid
  !!  Modification by Almut Gassmann, 2010-02-05
  !!  - Added feature for poor men's 3rd order advection, where a directional
  !!    laplace is needed at the edges.
  !!  Modification by Stephan Lorenz, 2010-06-02
  !!  - Storage moved from int_state into patch_oce since it is static
  !!    geometric information used in the ocean model
  !!  Modification by Peter Korn, 2010-11
  !!  - Calculation of cell area changed to achieve compatibility with
  !!    sw-model (cell area and consequently divergence different)
  !!  Modification by Stephan Lorenz, 2011-07
  !!   - 3-dim structures moved from patch_oce to hydro_ocean_base for parallelization
  !!
  SUBROUTINE init_geo_factors_base( ptr_patch, v_base )
    !
    IMPLICIT NONE
    !
    !  patch on which computation is performed
    !
    TYPE(t_patch), TARGET, INTENT(inout)    :: ptr_patch
    TYPE(t_hydro_ocean_base), INTENT(INOUT) :: v_base
    !

    INTEGER :: jc, jb, je, jv, je1, ie
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

    INTEGER :: ile, ibe, ilc1, ibc1, ilc2, ibc2, ifac, ic, ilnc, ibnc
    INTEGER :: ile1, ibe1,ile2,ibe2,ile3,ibe3
    !TYPE(cartesian_coordinates)::z_pn_k,z_pn_j
    !REAL(wp) :: z_lon, z_lat, z_nu, z_nv, z_proj
    REAL(wp) :: cell_area

    REAL(wp) :: z_sync_c(nproma,ptr_patch%nblks_c)
    REAL(wp) :: z_sync_e(nproma,ptr_patch%nblks_e)
    REAL(wp) :: z_sync_v(nproma,ptr_patch%nblks_v)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    & routine = ('mo_oce_state:init_geo_factors_base')

    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    i_nchdom   = MAX(1,ptr_patch%n_childdom)


!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk,ifac)
    ! a) Geometrical factor for divergence
    rl_start = 1
    rl_end = min_rlcell

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch cells (and blocks)
    !
!$OMP DO PRIVATE(jb,je,jc,i_startidx,i_endidx,ile,ibe)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk,      &
        &                i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
        ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
        ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
        ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
        ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
        ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

        cell_area =  0.25_wp&
  & *( ptr_patch%edges%primal_edge_length(ile1,ibe1)*ptr_patch%edges%dual_edge_length(ile1,ibe1)&
  &   +ptr_patch%edges%primal_edge_length(ile2,ibe2)*ptr_patch%edges%dual_edge_length(ile2,ibe2)&
  &   +ptr_patch%edges%primal_edge_length(ile3,ibe3)*ptr_patch%edges%dual_edge_length(ile3,ibe3))


       DO je = 1, i_cell_type

          IF (je > ptr_patch%cells%num_edges(jc,jb)) CYCLE ! relevant for hexagons

           ile = ptr_patch%cells%edge_idx(jc,jb,je)
           ibe = ptr_patch%cells%edge_blk(jc,jb,je)

           v_base%geofac_div(jc,je,jb) =                &
         &   ptr_patch%edges%primal_edge_length(ile,ibe) * &
         &   ptr_patch%cells%edge_orientation(jc,jb,je)  / &
         &   ptr_patch%cells%area(jc,jb)

           v_base%geofac_div(jc,je,jb) =                &
         &   ptr_patch%edges%primal_edge_length(ile,ibe) * &
         &   ptr_patch%cells%edge_orientation(jc,jb,je)  / &
         &   ptr_patch%cells%area(jc,jb)!cell_area

        ENDDO !edge loop

      ENDDO !idx loop

    END DO !block loop
!$OMP END DO

    ! b) Geometrical factor for rotation
    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rlvert

    ! Vorticity should have the right sign
      ifac = 0
      SELECT CASE (i_cell_type)
      CASE (3)
        ifac = 1
      CASE (6)
        ifac = -1
      END SELECT
      ! values for the blocking
      i_startblk = ptr_patch%verts%start_blk(rl_start,1)
      i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)
      !
      ! loop through all patch cells (and blocks)
      !
!$OMP DO PRIVATE(jb,je,jv,i_startidx,i_endidx,ile,ibe)
      DO jb = i_startblk, i_endblk

        CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)

        DO je = 1, 9-i_cell_type
          DO jv = i_startidx, i_endidx

            IF (je > ptr_patch%verts%num_edges(jv,jb)) CYCLE   ! relevant for hexagons

            ile = ptr_patch%verts%edge_idx(jv,jb,je)
            ibe = ptr_patch%verts%edge_blk(jv,jb,je)

            v_base%geofac_rot(jv,je,jb) =              &
         &    ptr_patch%edges%dual_edge_length(ile,ibe) * &
         &    ptr_patch%verts%edge_orientation(jv,jb,je)/ &
         &    ptr_patch%verts%dual_area(jv,jb) * REAL(ifac,wp)

          ENDDO !vertex loop
        ENDDO

      END DO !block loop
!$OMP END DO

      ! c) Geometrical factor for nabla2_scalar
      rl_start = 1  ! #slo# changed to 1 - 2010-12-07
      rl_end = min_rlcell

      ! values for the blocking
      i_startblk = ptr_patch%cells%start_blk(rl_start,1)
      i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
      !
      ! loop through all patch cells (and blocks)
      !
!$OMP DO PRIVATE(jb,je,jc,ic,i_startidx,i_endidx,ile,ibe,ilc1,ibc1,&
!$OMP    ilc2,ibc2,ilnc,ibnc)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO je = 1, i_cell_type
          DO jc = i_startidx, i_endidx

            ile = ptr_patch%cells%edge_idx(jc,jb,je)
            ibe = ptr_patch%cells%edge_blk(jc,jb,je)

            ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
            ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
            ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
            ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

            IF (jc == ilc1 .AND. jb == ibc1) THEN
              IF (i_cell_type == 3) THEN
                v_base%geofac_n2s(jc,1,jb)     =  &
                &  v_base%geofac_n2s(jc,1,jb)  -  &
                &  v_base%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)
            ELSE IF (i_cell_type == 6) THEN
              v_base%geofac_n2s(jc,1,jb)       =  &
                &  v_base%geofac_n2s(jc,1,jb)  -  &
                &  v_base%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)*  &
                &  ptr_patch%edges%system_orientation(ile,ibe)
            ENDIF
          ELSE IF (jc == ilc2 .AND. jb == ibc2) THEN
            IF (i_cell_type == 3) THEN
              v_base%geofac_n2s(jc,1,jb)       =  &
                &  v_base%geofac_n2s(jc,1,jb)  +  &
                &  v_base%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)
            ELSE IF (i_cell_type == 6) THEN
              v_base%geofac_n2s(jc,1,jb)       =  &
                &  v_base%geofac_n2s(jc,1,jb)  +  &
                &  v_base%geofac_div(jc,je,jb) /  &
                &  ptr_patch%edges%dual_edge_length(ile,ibe)*  &
                &  ptr_patch%edges%system_orientation(ile,ibe)
            ENDIF
          ENDIF
          DO ic = 1, i_cell_type
            ilnc = ptr_patch%cells%neighbor_idx(jc,jb,ic)
            ibnc = ptr_patch%cells%neighbor_blk(jc,jb,ic)
            IF (ilnc == ilc1 .AND. ibnc == ibc1) THEN
              IF (i_cell_type == 3) THEN
                v_base%geofac_n2s(jc,ic+1,jb)     = &
                  &  v_base%geofac_n2s(jc,ic+1,jb)- &
                  &  v_base%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe)
              ELSE IF (i_cell_type == 6) THEN
                v_base%geofac_n2s(jc,ic+1,jb)     = &
                  &  v_base%geofac_n2s(jc,ic+1,jb)- &
                  &  v_base%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe) * &
                  &  ptr_patch%edges%system_orientation(ile,ibe)
              ENDIF
            ELSE IF (ilnc == ilc2 .AND. ibnc == ibc2) THEN
              IF (i_cell_type == 3) THEN
                v_base%geofac_n2s(jc,ic+1,jb)     = &
                  &  v_base%geofac_n2s(jc,ic+1,jb)+ &
                  &  v_base%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe)
              ELSE IF (i_cell_type == 6) THEN
                v_base%geofac_n2s(jc,ic+1,jb)     = &
                  &  v_base%geofac_n2s(jc,ic+1,jb)+ &
                  &  v_base%geofac_div(jc,je,jb)  / &
                  &  ptr_patch%edges%dual_edge_length(ile,ibe) * &
                  &  ptr_patch%edges%system_orientation(ile,ibe)
              ENDIF
            ENDIF
          ENDDO

          ! To ensure that dummy edges have a factor of 0:
          IF (je > ptr_patch%cells%num_edges(jc,jb)) THEN
            v_base%geofac_n2s(jc,je+1,jb) = 0._wp
          ENDIF

        ENDDO !cell loop
      ENDDO

    END DO !block loop
!$OMP END DO

    ! d) Geometrical factor for quad-cell divergence (triangles only)
    IF (i_cell_type == 3) THEN

      rl_start = 1  ! #slo# changed to 1 - 2010-12-07
      rl_end = min_rledge

      ! values for the blocking
      i_startblk = ptr_patch%edges%start_blk(rl_start,1)
      i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,je,je1,i_startidx,i_endidx,ile,ibe)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO je1 = 1, 4
          DO je = i_startidx, i_endidx

          ile = ptr_patch%edges%quad_idx(je,jb,je1)
          ibe = ptr_patch%edges%quad_blk(je,jb,je1)

          v_base%geofac_qdiv(je,je1,jb) =               &
            ptr_patch%edges%primal_edge_length(ile,ibe) *  &
            ptr_patch%edges%quad_orientation(je,jb,je1) /  &
            ptr_patch%edges%quad_area(je,jb)

          ENDDO !edge loop
        ENDDO

      END DO !block loop
!$OMP END DO

    ENDIF

    ! f) compute inverse dual edge length (used in math_operators for the ocean)

    rl_start = 1  ! #slo# changed to 1 - 2010-12-07
    rl_end = min_rledge

    ! Second step: computed projected orientation vectors and related information
    i_startblk = ptr_patch%edges%start_blk(rl_start,1)
    i_endblk   = ptr_patch%edges%end_blk(rl_end,i_nchdom)
    !
    ! loop through all patch edges
    !
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(ptr_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO je =  i_startidx, i_endidx

        ! compute inverse dual edge length (undefined for refin_ctrl=1)

        ptr_patch%edges%inv_dual_edge_length(je,jb) = &
          1._wp/ptr_patch%edges%dual_edge_length(je,jb)

      ENDDO

    END DO !block loop
!$OMP END DO

!$OMP END PARALLEL

    ! synchronize all elements of v_base:

    DO ie = 1, i_cell_type

        z_sync_c(:,:) = v_base%geofac_div(:,ie,:)
        CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:))
        v_base%geofac_div(:,ie,:) = z_sync_c(:,:)

        z_sync_c(:,:) = v_base%geofac_n2s(:,ie,:)
        CALL sync_patch_array(SYNC_C, ptr_patch, z_sync_c(:,:))
        v_base%geofac_n2s(:,ie,:) = z_sync_c(:,:)

    END DO

    DO ie = 1, 4

        z_sync_e(:,:) = v_base%geofac_qdiv(:,ie,:)
        CALL sync_patch_array(SYNC_e, ptr_patch, z_sync_e(:,:))
        v_base%geofac_qdiv(:,ie,:) = z_sync_e(:,:)

    END DO

    CALL sync_patch_array(SYNC_E, ptr_patch, ptr_patch%edges%inv_dual_edge_length(:,:))

    DO ie = 1, 9-i_cell_type

        z_sync_v(:,:) = v_base%geofac_rot(:,ie,:)
        CALL sync_patch_array(SYNC_v, ptr_patch, z_sync_v(:,:))
        v_base%geofac_rot(:,ie,:) = z_sync_v(:,:)

    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_geo_factors_base
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
  !!
  SUBROUTINE init_coriolis_oce( ptr_patch )
    !
    IMPLICIT NONE
    !
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: ptr_patch
    !
    INTEGER :: jb, je, jv
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start_v, rl_end_v
    INTEGER :: i_startblk_v, i_endblk_v, i_startidx_v, i_endidx_v
    TYPE(t_geographical_coordinates) :: gc1,gc2 
    TYPE(t_cartesian_coordinates) :: xx1, xx2
    REAL(wp) :: z_y, z_lat_basin_center
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    & routine = ('mo_oce_state:init_coriolis_oce')
    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    rl_start_e = 1
    rl_end_e   = min_rledge
    rl_start_v = 1
    rl_end_v   = min_rlvert

    i_startblk_e = ptr_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = ptr_patch%edges%end_blk(rl_end_e,1)
    i_startblk_v = ptr_patch%verts%start_blk(rl_start_v,1)
    i_endblk_v   = ptr_patch%verts%end_blk(rl_end_v,1)

    SELECT CASE (CORIOLIS_TYPE)

    CASE(BETA_PLANE_CORIOLIS)

      CALL message (TRIM(routine), 'BETA_PLANE_CORIOLIS: set to linear approximation')

      z_lat_basin_center = basin_center_lat * deg2rad
      gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
      gc1%lon = 0.0_wp
      xx1=gc2cc(gc1)

      DO jb = i_startblk_v, i_endblk_v
        CALL get_indices_v(ptr_patch, jb, i_startblk_v, i_endblk_v, &
                          i_startidx_v, i_endidx_v, rl_start_v, rl_end_v)
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

      DO jb = i_startblk_e, i_endblk_e
        CALL get_indices_e(ptr_patch, jb, i_startblk_e, i_endblk_e, &
                          i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
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
    REAL(wp), INTENT(OUT) :: zlev_i(n_zlev+1)    , zlev_m(n_zlev)
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
