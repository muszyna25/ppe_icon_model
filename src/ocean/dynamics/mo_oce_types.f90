MODULE  mo_oce_types
  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    &                               success, max_char_length, MIN_DOLIC,               &
    &                               full_coriolis, beta_plane_coriolis,                &
    &                               f_plane_coriolis, zero_coriolis, halo_levels_ceiling
  USE mo_math_utilities,      ONLY: gc2cc,t_cartesian_coordinates,cvec2gvec,      &
    &                               t_geographical_coordinates, &!vector_product, &
    &                               arc_length,set_del_zlev

  PUBLIC :: t_hydro_ocean_base
  PUBLIC :: t_hydro_ocean_state
  PUBLIC :: t_hydro_ocean_prog
  PUBLIC :: t_hydro_ocean_diag
  PUBLIC :: t_hydro_ocean_aux
  PUBLIC :: t_hydro_ocean_acc
  PUBLIC :: t_ptr3d
  PUBLIC :: t_oce_config
  PUBLIC :: t_ocean_tracer

  PUBLIC :: t_ocean_regions
  PUBLIC :: t_ocean_region_volumes
  PUBLIC :: t_ocean_region_areas
  PUBLIC :: t_ocean_basins
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
    ! index1=1,nproma, index2=1,n_zlev, index3=1,alloc_cell_blocks
    INTEGER, ALLOCATABLE :: lsm_c(:,:,:)
    ! land-sea-mask for cell edges
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e
    INTEGER, ALLOCATABLE :: lsm_e(:,:,:)
    ! land-sea-mask for cell vertices
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_v
    ! INTEGER, ALLOCATABLE :: lsm_v(:,:,:)


    ! To simplify the acess to the required information within these loops
    ! we store an cell and edge based version of the deepest ocean layer
    ! in column. dolic_e(edge1) and dolic_c(cell1) are identical if 'edge1'
    ! is one of the edges of 'cell1'.
    ! If the ocean bottom is flat dolic_c and dolic_e are identical and equal
    ! to the number of z-coodinate surfaces.

    ! index1=1,nproma, index2=1,alloc_cell_blocks
    INTEGER, ALLOCATABLE :: dolic_c(:,:)
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: dolic_e(:,:)

    ! For diagnosis like stream functions and area calculations we add surface arrays
    ! index1=1,nproma, index2=1,alloc_cell_blocks
    INTEGER,  ALLOCATABLE :: basin_c(:,:)  ! basin information Atlantic/Indian/Pacific
    INTEGER,  ALLOCATABLE :: regio_c(:,:)  ! area information like tropical Atlantic etc.

    ! To simply set land points to zero we store additional 3-dim wet points
    ! dimensions as in lsm_oce:
    REAL(wp), ALLOCATABLE :: wet_c(:,:,:)  ! cell centers
    REAL(wp), ALLOCATABLE :: wet_e(:,:,:)  ! cell edges
    !REAL(wp), ALLOCATABLE :: wet_i(:,:,:)  ! vertical velocity points
    !                                       ! on intermediate levels


!!$    ! Arrays that describe vertical connectivity of the triangles.
!!$    ! The indices of triangles are stored in a whole vertical column.
!!$    ! index1=1,nproma, index2=1,alloc_cell_blocks, index3=1,n_zlev
!!$    INTEGER, ALLOCATABLE :: neighbor_c(:,:,:)
!!$    ! index1=1,nproma, index2=1,nblks_e, index3=1,n_zlev
!!$    INTEGER, ALLOCATABLE :: neighbor_e(:,:,:)
  END TYPE t_hydro_ocean_base

  TYPE t_ptr3d
    REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
  END TYPE t_ptr3d

  TYPE t_ocean_tracer
    REAL(wp),POINTER :: concentration(:,:,:)
    REAL(wp),POINTER :: concentration_x_height(:,:,:)
  END TYPE t_ocean_tracer
!
!! prognostic variables
!
  TYPE t_hydro_ocean_prog

    REAL(wp), POINTER ::    &
      &  h(:,:)                ,& ! height of the free surface. Unit: [m]
                                  ! dimension:(nproma, alloc_cell_blocks)
      &  vn(:,:,:)             ,& ! velocity component normal to cell edge. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  tracer(:,:,:,:)          ! tracer concentration.
                                  ! dimension: (nproma, n_zlev, alloc_cell_blocks, no_tracer)
                                  ! Ordering of tracers:
                                  !   1) pot_temp:= potential temperature, Unit: [deg C]
                                  !   2) salinity:= salinity, Unit [psu]

    TYPE(t_ocean_tracer), ALLOCATABLE :: ocean_tracers(:)

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
                                  ! dimension: (nproma,n_zlev, alloc_cell_blocks)
      &  rhopot(:,:,:)         ,& ! potential density. Unit: [kg/m^3]
                                  ! dimension: (nproma,n_zlev, alloc_cell_blocks)
      &  zgrad_rho(:,:,:)      ,& ! vertical density gradient. Unit: [kg/m^2]
                                  ! dimension: (nproma,n_zlev, alloc_cell_blocks)
      &  h_e(:,:)              ,& ! surface height at cell edges. Unit [m].
                                  ! dimension: (nproma, nblks_e)
      &  thick_c(:,:)          ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma, alloc_cell_blocks)
      &  thick_e(:,:)          ,& ! individual fluid column thickness at edges. Unit [m].
                                  ! dimension: (nproma, nblks_e)
      &  mass_flx_e(:,:,:)     ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  div_mass_flx_c(:,:,:) ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma,n_zlev, alloc_cell_blocks)
      &  w(:,:,:)              ,& ! vertical velocity. Unit [m/s].
                                  ! dimension: (nproma, n_zlev+1, alloc_cell_blocks)
      &  w_old(:,:,:)          ,& ! vertical velocity from previous timestep. Unit [m/s].
                                  ! dimension: (nproma, n_zlev+1, alloc_cell_blocks)
      &  w_e(:,:,:)            ,& ! vertical velocity at edges. Unit [m/s]
                                  ! dimension: (nproma, n_zlev+1, nblks_e)
      &  w_prev(:,:,:)         ,& ! vertical velocity at cells, from previous timestep. Unit [m/s]
                                  ! dimension: (nproma, n_zlev+1, alloc_cell_blocks)
      &  u(:,:,:)              ,& ! reconstructed zonal velocity component. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, alloc_cell_blocks)
      &  v(:,:,:)              ,& ! reconstructed meridional velocity component. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, alloc_cell_blocks)
      &  u_vint(:,:)           ,& ! barotropic zonal velocity. Unit [m*m/s]
                                  ! dimension: (nproma, alloc_cell_blocks)
      &  ptp_vn(:,:,:)         ,& ! normal velocity after mapping P^T P
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_pred(:,:,:)        ,& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_impl_vert_diff(:,:,:),& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_time_weighted(:,:,:),&  ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  w_time_weighted(:,:,:),& ! predicted normal velocity vector at cells
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  vort(:,:,:)           ,& ! vorticity at triangle vertices. Unit [1/s]
                                  ! dimension: (nproma, n_zlev, nblks_v)
      &  vort_e(:,:,:)         ,& ! vorticity interpolated to triangle edges. Unit [1/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  kin(:,:,:)            ,& ! kinetic energy. Unit [m/s].
                                  ! (nproma, n_zlev, alloc_cell_blocks)
      &  mld(:,:)              ,& ! mixed layer depth [m].
                                  ! (nproma,  alloc_cell_blocks)
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
                                  ! dimension: (nproma, n_zlev, alloc_cell_blocks)
      &  press_hyd(:,:,:)      ,& ! hydrostatic pressure. Unit [m]
                                  ! dimension: (nproma, n_zlev, alloc_cell_blocks)
      &  press_grad(:,:,:)     ,& ! hydrostatic pressure gradient term. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  temp_insitu(:,:,:)

    INTEGER, POINTER :: &
      & condep(:,:)               ! convection depth index
                                  ! (nproma,  alloc_cell_blocks)

    TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_vn(:,:,:)              ! reconstructed velocity at cell center in cartesian coordinates
                                  ! dimension: (nproma, n_zlev, alloc_cell_blocks)
    TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_vn_dual(:,:,:)         ! reconstructed velocity at vertex in cartesian coordinates
                                  ! dimension: (nproma, n_zlev, nblks_v)
    TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_vn_mean(:,:,:)         ! reconstructed velocity at vertex in cartesian coordinates
                                  ! dimension: (nproma, n_zlev, nblks_v). For mimetic miura scheme

   TYPE(t_cartesian_coordinates), POINTER :: &
      &  p_mass_flux_sfc_cc(:,:)  ! mass flux at surface in cartesian coordinates
                                  ! dimension: (nproma, alloc_cell_blocks).
   !-----------------------------------------------------------------------------------
   ! dummy pointers for prognostic variables:
   REAL(wp), POINTER :: h(:,:),vn(:,:,:),t(:,:,:),s(:,:,:) ! dummy pointer for output variabless
   !-----------------------------------------------------------------------------------
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
    REAL(wp), POINTER ::       &
      &  bc_top_vn(:,:)       ,& ! normal velocity boundary condition at surface
                                 ! dimension: (nproma,nblks_e)
      &  bc_bot_vn(:,:)       ,& ! normal velocity boundary condition at bottom
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_top_u(:,:)        ,& ! zonal velocity boundary condition at surface
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_top_v(:,:)        ,& ! meridional velocity boundary condition at surface
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_bot_u(:,:)        ,& ! zonal velocity boundary condition at bottom
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_bot_v(:,:)        ,& ! meridional velocity boundary condition at bottom
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_top_w(:,:)        ,& ! vertical velocity boundary condition at surface
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_bot_w(:,:)        ,& ! vertical velocity boundary condition at bottom
      &  bc_top_tracer(:,:,:) ,& ! vertical velocity boundary condition at surface
                                 ! dimension: (nproma,alloc_cell_blocks)
      &  bc_bot_tracer(:,:,:) ,& ! vertical velocity boundary condition at bottom
      &  p_rhs_sfc_eq(:,:)!,   & ! right hand side of surface equation
                                         ! dimension: (nproma,alloc_cell_blocks)
     TYPE(t_cartesian_coordinates), POINTER :: bc_top_veloc_cc(:,:), &
                                  &                bc_bot_veloc_cc(:,:)
     TYPE(t_ptr3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer

    ! Variables for 3-dim tracer relaxation:
    REAL(wp), POINTER ::         &
      &  relax_3d_data_T(:,:,:), & ! 3-dim temperature relaxation data (T*)
                                   ! dimension: (nproma,n_zlev,alloc_cell_blocks)
      &  relax_3d_forc_T(:,:,:), & ! 3-dim temperature relaxation forcing (1/tau*(T-T*))
                                   ! dimension: (nproma,n_zlev,alloc_cell_blocks)
      &  relax_3d_data_S(:,:,:), & ! 3-dim salinity relaxation data (T*)
                                   ! dimension: (nproma,n_zlev,alloc_cell_blocks)
      &  relax_3d_forc_S(:,:,:)    ! 3-dim salinity relaxation forcing (1/tau*(T-T*))
                                   ! dimension: (nproma,n_zlev,alloc_cell_blocks)

  END TYPE t_hydro_ocean_aux

  ! variables to be accumulated
  TYPE t_hydro_ocean_acc
    REAL(wp), POINTER :: &
      & h(:,:)                  ,&
      & u(:,:,:)                ,&
      & v(:,:,:)                ,&
      & w(:,:,:)                ,& ! vertical velocity. Unit [m/s].
      & vt(:,:,:)               ,& ! tangential velocity component at edges. Unit [m/s].
      & rho(:,:,:)              ,& ! density. Unit: [kg/m^3]
      & rhopot(:,:,:)           ,& ! potential density. Unit: [kg/m^3]
      & mass_flx_e(:,:,:)       ,& ! mass flux at edges. Unit [?].
      & div_mass_flx_c(:,:,:)   ,& ! divergence of mass flux at cells. Unit [?].
      & u_vint(:,:)             ,& ! barotropic zonal velocity. Unit [m*m/s]
      & ptp_vn(:,:,:)           ,& ! normal velocity after mapping P^T P
      & vn_pred(:,:,:)          ,& ! predicted normal velocity vector at edges.
      & vn_impl_vert_diff(:,:,:),& ! predicted normal velocity vector at edges.
      & vn_time_weighted(:,:,:) ,&  ! predicted normal velocity vector at edges.
      & w_time_weighted(:,:,:)  ,& ! predicted normal velocity vector at cells.
      & vort(:,:,:)             ,& ! vorticity at triangle vertices. Unit [1/s]
      & kin(:,:,:)              ,& ! kinetic energy. Unit [m/s].
      & veloc_adv_horz(:,:,:)   ,& ! horizontal velocity advection
      & veloc_adv_vert(:,:,:)   ,& ! vertical velocity advection
      & laplacian_horz(:,:,:)   ,& ! horizontal diffusion of horizontal velocity
      & laplacian_vert(:,:,:)   ,& ! vertical diffusion of horizontal velocity
      & grad(:,:,:)             ,& ! gradient of kinetic energy. Unit [m/s]
      & div(:,:,:)              ,& ! divergence. Unit [m/s]
      & press_hyd(:,:,:)        ,& ! hydrostatic pressure. Unit [m]
      & press_grad(:,:,:)       ,& ! hydrostatic pressure gradient term. Unit [m/s]
      & temp_insitu(:,:,:)      ,&
      & tracer(:,:,:,:)
    TYPE(t_ptr3d),ALLOCATABLE :: tracer_ptr(:)  !< pointer array: one pointer for each tracer
  END TYPE

!! array of states
!
  TYPE t_hydro_ocean_state

    TYPE(t_hydro_ocean_prog), POINTER :: p_prog(:)    ! time array of prognostic states at different time levels
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_aux)  :: p_aux
    TYPE(t_hydro_ocean_acc)  :: p_acc

  END TYPE t_hydro_ocean_state

  INTEGER, PARAMETER               :: max_tracers = 2
  TYPE t_oce_config
    CHARACTER(len=max_char_length) :: tracer_names(max_tracers)
    CHARACTER(len=max_char_length) :: tracer_longnames(max_tracers)
    CHARACTER(len=max_char_length) :: tracer_units(max_tracers)
    CHARACTER(len=max_char_length) :: tracer_tags(max_tracers)
    INTEGER                        :: tracer_codes(max_tracers)
  END TYPE t_oce_config

  !----------------------------------------------------------------------------
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
  TYPE t_ocean_regions
    INTEGER            :: &
      & land                            = 0,&
      & greenland_iceland_norwegian_sea = 1,&
      & arctic_ocean                    = 2,&
      & labrador_sea                    = 3,&
      & north_atlantic                  = 4,&
      & tropical_atlantic               = 5,&
      & southern_ocean                  = 6,&
      & indian_ocean                    = 7,&
      & tropical_pacific                = 8,&
      & north_pacific                   = 9,&
      & caribbean                       = -33
  END TYPE t_ocean_regions
  TYPE t_ocean_region_volumes
    REAL(wp)            :: &
      & land                            = 0.0_wp,&
      & greenland_iceland_norwegian_sea = 0.0_wp,&
      & arctic_ocean                    = 0.0_wp,&
      & labrador_sea                    = 0.0_wp,&
      & north_atlantic                  = 0.0_wp,&
      & tropical_atlantic               = 0.0_wp,&
      & southern_ocean                  = 0.0_wp,&
      & indian_ocean                    = 0.0_wp,&
      & tropical_pacific                = 0.0_wp,&
      & north_pacific                   = 0.0_wp,&
      & caribbean                       = 0.0_wp,&
      & total                           = 0.0_wp
  END TYPE t_ocean_region_volumes
  TYPE t_ocean_region_areas
    REAL(wp)            :: &
      & land                            = 0.0_wp,&
      & greenland_iceland_norwegian_sea = 0.0_wp,&
      & arctic_ocean                    = 0.0_wp,&
      & labrador_sea                    = 0.0_wp,&
      & north_atlantic                  = 0.0_wp,&
      & tropical_atlantic               = 0.0_wp,&
      & southern_ocean                  = 0.0_wp,&
      & indian_ocean                    = 0.0_wp,&
      & tropical_pacific                = 0.0_wp,&
      & north_pacific                   = 0.0_wp,&
      & caribbean                       = 0.0_wp,&
      & total                           = 0.0_wp
  END TYPE t_ocean_region_areas
  !-----------------------------
  !
  ! Ocean basins:
  !  1: Atlantic; 3: Pacific, for Indean and pacific the area values ara used
  !
  !-----------------------------
  TYPE t_ocean_basins
    INTEGER            :: &
      & atlantic = 1, pacific = 3
  END TYPE t_ocean_basins
END MODULE mo_oce_types

