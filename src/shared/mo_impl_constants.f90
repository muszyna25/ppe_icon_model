!>
!!  Module contains some constants relevant for implementational issues.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn (2005)
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Peter Korn (2006-08):
!!  - added identifier for edges, cells, vertices (used for memory allocation)
!! @par
!!  Modification by Peter Korn (2006-12):
!!   - added identifier for variables (used for land-sea mask application)
!!  Modification by Hui Wan (2007-02):
!!  - added "nproma" and some comments
!!  Modification by Guenther Zaengl (2008-10-23):
!!  - added parameters defining the range of the refin_ctrl flags
!!    (these are also set in the patch generator and must match each other)
!! @par
!!  Modification by Hui Wan (2009-08-07)
!!  - added identifiers for time stepping methods.
!!  Modifications by Daniel Reinert (2010-10-06)
!!  - added identifiers for NWP physics time control
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
MODULE mo_impl_constants
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  USE mo_kind,               ONLY: wp
  USE mtime,                 ONLY: MAX_TIMEDELTA_STR_LEN
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e

  IMPLICIT NONE

  PUBLIC

  ! define model name and version
  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  INTEGER, PARAMETER :: MAX_CHAR_LENGTH     = 1024

  INTEGER, PARAMETER :: SUCCESS             = 0
  INTEGER, PARAMETER :: CELLS               = 123
  INTEGER, PARAMETER :: EDGES               = 345
  INTEGER, PARAMETER :: VERTS               = 678
  
  INTEGER, PARAMETER :: ON_CELLS            = 1
  INTEGER, PARAMETER :: ON_EDGES            = 2
  INTEGER, PARAMETER :: ON_VERTICES         = 3
  INTEGER, PARAMETER :: HALO_LEVELS_CEILING = 256 ! should be greater than the max level
                                         ! of halo levels
!-------------------------------------------------------------------------------
! Comments by Hui:
! According to Luis' explanation, the declarations above related to the blocking
! are not correct. There should be a single NPROMA, and different NBLKS values
! for edges, cells and vertices.
!
! Take the triangular cells for example. Considering that the number of cells is
! different from patch to patch, it may be a good idea NOT to declare "nblks_c" as
! a 2D array HERE, but as a component of the "patch" type. Since we want to do
! calculations ONLY for the internal cells, we also need to know how many blocks
! (e.g. "nblks_i_c") we have for the internal cells. Furthermore, it is also
! necessary to define a a parameter "npromz_i_c" for each patch, which stores the
! number of valid items in the last block of the internal cells.
!   At the very beginning of the execution of the program, "nproma" can be read
! from the namelist file. Then, after reading in the number of cells and external halo
! cells of a certain patch, we calculate the value of "nblks_c", "nblks_i_c"
! and "npromz_i_c" for that patch by:
!
!   n_patch_cell_all     = ptr_patch%ncells + ptr_patch%n_e_halo_cells
!   ptr_patch%nblks_c    = ( n_patch_cell_all - 1 )/nproma + 1
!
!   ptr_patch%nblks_i_c  = ( ptr_patch%ncell - 1 )/nproma + 1
!   ptr_patch%npromz_i_c = ptr_patch%ncell - (ptr_patch%nblks_i_c - 1)*nproma
!
!   (cf. echam5, mo_decomposition.f90,
!        the last few lines of SUBROUTINE decompose_grid)
!
! These calculations can be done in {\it mo\_model\_domain\_import} .
!-------------------------------------------------------------------------------

  !identifiers for prognostic and diagnostic variables
  ! required for mask application
  !
  !prognostic
  INTEGER, PARAMETER :: VELOCITY_NORMAL     = 1
  INTEGER, PARAMETER :: VELOCITY_TANGENTIAL = 2
  INTEGER, PARAMETER :: HEIGHT              = 3
  INTEGER, PARAMETER :: TRACER              = 4
  !
  !diagnostic
  INTEGER, PARAMETER :: VORTICITY           = 5
  INTEGER, PARAMETER :: DIVERGENCE          = 6
  INTEGER, PARAMETER :: KINETIC_ENERGY      = 7
  INTEGER, PARAMETER :: THICK_ED            = 8

  ! external
  ! #slo# changed 2010-08-03
  INTEGER,PARAMETER :: LAND          =  2  !  inner land
  INTEGER,PARAMETER :: LAND_BOUNDARY =  1  !  e.g. land cell with neighbouring wet cell
  INTEGER,PARAMETER :: BOUNDARY      =  0  !  edge where cells are differing
  INTEGER,PARAMETER :: SEA_BOUNDARY  = -1  !  e.g. wet cell with neighbouring land cell
  INTEGER,PARAMETER :: SEA           = -2  !  inner sea

  INTEGER,PARAMETER :: ZERO_CORIOLIS       = 0
  INTEGER,PARAMETER :: FULL_CORIOLIS       = 1
  INTEGER,PARAMETER :: BETA_PLANE_CORIOLIS = 2
  INTEGER,PARAMETER :: F_PLANE_CORIOLIS    = 3


  ! dimensions of index list fields and index ranges for which grid points are reordered
  ! according to the refin_ctrl values
  ! Specifically:
  ! max_rl* indicates the number of cell/edge/vertex rows along the lateral boundary of nested
  !   domains for which grid points are reordered, i.e. moved to the beginning of the index lists;
  !   the number of cell rows for which the refin_ctrl flag is set is determined by the variable
  !   bdy_indexing_depth in prepare_gridref; it is in general larger than max_rlcell
  !  (the refin_ctrl flag here counts the distance from the lateral boundary in units of cell rows)
  ! ABS(min_rl*_int)-1 indicates the number of cell/edge/vertex rows overlapping with the lateral boundary
  !   of a nested domain for which grid points are reordered, i.e. moved to the end of the index lists;
  !   min_rl*_int refers to grid points overlapping with interior points of nested domains
  ! Finally, the indices between min_rl*_int-1 and min_rl*_int are reserved for halo points emerging
  !   from the MPI domain decomposition; these parts of the index lists are empty on exit of the
  !   grid generator and are filled only on mo_subdivision. However, the index list fields are always
  !   dimensioned with (min_rl*:max_rl*). The values set below are sufficient for a halo
  !   width of two full cell rows; normally we use one, but stencils for high-order schemes may
  !   sometime require a halo width of two full rows
!
!   ------------------------------------------
!   LL: copied from the icon_flowcontrol as described by Guenther:
!
!   The ordering of the halo points is as follows:
!   min_rlcell_int- 1: halo cells having a prognostic cell as neighbor
!   min_rlcell_int- 2: halo cells in the first cell row having no prognostic cell as neighbor
!   and analogously for the second halo cell row if present. For n_ghost_rows = 1, the index segments
!   corresponding to min_rlcell_int - 3 and min_rlcell_int - 4 (= min_rlcell) are empty.
! 
!   For edges and vertices, one needs to be aware of the fact that outer boundary edges/vertices of a prognostic
!   cell may not be owned by the current PE because the PE of the neighboring cell has the ownership (otherwise
!   there would be double-counting). There are, however, operations for which even such edges/vertices can be 
!   excluded from prognostic computation because a halo synchronization follows immediately afterwards (and
!   has to be there anyway). Thus, the following ordering is applied:
!   min_rledge_int - 1: outer boundary edges of a prognostic cell not owned by the current PE\\
!   min_rledge_int - 2: edges connecting halo cells of the first row
!   min_rledge_int - 3: outer boundary edges of the first halo cells row, or edges connecting cells
!   of the first halo cell row with cells of the second halo cell row.
!   For n_ghost_rows = 2, an analogous setting applies to min_rledge_int - 4 and
!   min_rledge_int - 5 (= min_rledge). For vertices, we have
!   min_rlvert_int - 1: outer boundary vertices of a prognostic cell not owned by the current PE
!   min_rlvert_int - 2: outer boundary vertices of the first halo cells row, or vertices connecting cells
!   of the first halo cell row with cells of the second halo cell row.
!   For n_ghost_rows = 2, an analogous setting applies to min_rlvert_int - 3  (= min_rlvert).
  !---------------------------------------------
  !
  ! Ordering Scheme:
  !
  ! Following is the order of the grid entities (an all the associated variables)
  ! in ascending order.
  !
  ! A. The indices from 1 to max_rl
  !    Mark the lateral boundaries of the patch
  !    start_idx(1) = start of the first boundary level. It is always 1
  !    end_idx(1)   = end of the first boundary level end_idx(1)
  !    start_idx(2) = start of the second boundary level, it is always end_idx(1)+1
  !    ..... etc until end_idx(max_rl) which the the end of the lateral boundary levels
  !
  ! B. The index 0
  !    Marks the internal entities, that do not overlap with child patches
  !    start_idx(0) = end_idx(max_rl) + 1
  !    end_idx(0)   = start_idx(-1) -1
  !
  ! C. The indices from -1 to min_rl_int
  !    Mark the internal entities that overlap with child patches
  !    (they are defined for each child patch)
  !    start_idx(-1) = start of the internal entities overlapping with the first (two) levels
  !                    of the lateral boundaries of the child patch
  !    end_idx(-1)   = end of the internal entities overlapping with the first (two) levels
  !                    of the lateral boundaries of the child patch
  !    start_idx(-2) = start of the internal entities overlapping with the next (two) levels
  !                    of the lateral boundaries of the child patch = end_idx(-1) + 1
  !    ... etc
  !    end_idx(minrl_int) = end of all the internal entities overlapping with the the child patch
  !
  ! D. The indices from min_rl_int-1 to min_rl
  !    Mark the halo entities, when they do not overlap with a child patch
  !    Note: See above
  !
  !---------------------------------------------
  !
  ! Examples:
  !  A. Get all entities in the grid:    start_idx(1) -- end_idx(min_rl)
  !     This is the default range for most operators
  !  B. Get all owned enitities: start_idx(1) -- end_idx(min_rl_int)
  !     Note that this may still contain halo entities if they ovelap with child patches
  !  C. Get all enitities, except halos: for cells: start_idx(1) -- end_idx(min_rl_int)
  !        For verts/edges: start_idx(1) -- end_idx(min_rl_int - 1)
  !
  !---------------------------------------------

  INTEGER, PARAMETER :: max_hw         = 2                         ! maximum halo width (n_ghost_rows)
  !
  INTEGER, PARAMETER :: min_rlcell_int = -4                        ! previously -6
  INTEGER, PARAMETER :: min_rlcell     = min_rlcell_int - 2*max_hw ! = -8
  INTEGER, PARAMETER :: max_rlcell     = 5                         ! previously 8
  INTEGER, PARAMETER :: min_rlvert_int = min_rlcell_int
  INTEGER, PARAMETER :: min_rlvert     = min_rlvert_int - (max_hw+1)
  INTEGER, PARAMETER :: max_rlvert     = max_rlcell
  INTEGER, PARAMETER :: min_rledge_int = 2*min_rlcell_int          ! -8
  INTEGER, PARAMETER :: min_rledge     = min_rledge_int - (2*max_hw+1)  ! -13
  INTEGER, PARAMETER :: max_rledge     = 2*max_rlcell              ! 10

  ! Aliases for loop control
  !
  ! start index level for prognostic cells (excluding a possible nest boundary interpolation zone)
  INTEGER, PARAMETER :: start_prog_cells = grf_bdywidth_c+1
  ! start index level for prognostic edges (excluding a possible nest boundary interpolation zone)
  INTEGER, PARAMETER :: start_prog_edges = grf_bdywidth_e+1
  ! remark: the corresponding parameter for vertices would be unused
  !
  ! end index level for computations excluding all halo cells
  INTEGER, PARAMETER :: end_prog_cells = min_rlcell_int
  ! end index level for computations including halo level 1 cells (direct neighbors, i.e. halo cells 
  ! sharing an edge with a prognostic cell)
  INTEGER, PARAMETER :: end_halo_lev1_cells = min_rlcell_int - 1
  ! end index for computations including all halo cells (this adds indirect neighbors, i.e.
  ! halo cells sharing only a vertex with a prognostic cell)
  ! remark: n_ghost_rows=2, i.e. using a full second row of halo cells, is currently not foreseen in the code,
  ! so an end index level of min_rlcell_int - 2 is equivalent to min_rlcell
  INTEGER, PARAMETER :: end_all_cells = min_rlcell
  !
  ! end index level for computations excluding all halo edges
  INTEGER, PARAMETER :: end_prog_edges = min_rledge_int
  ! end index level for edges of prognostic cells (including those not owned by the current PE)
  INTEGER, PARAMETER :: end_edges_of_prog_cells = min_rledge_int - 1
  ! end index level for computations including all edges of halo level 1 cells
  INTEGER, PARAMETER :: end_halo_lev1_edges = min_rledge_int - 2
  ! end index for computations including all halo edges (remark: consistent with what has been said above,
  ! min_rledge is equivalent to min_rledge_int - 3)
  INTEGER, PARAMETER :: end_all_edges = min_rledge
  !
  ! end index level for computations excluding all halo verices
  INTEGER, PARAMETER :: end_prog_verts = min_rlvert_int
  ! end index level for vertices of prognostic cells (including those not owned by the current PE)
  INTEGER, PARAMETER :: end_verts_of_prog_cells = min_rlvert_int - 1
  ! end index for computations including all halo vertices (remark: consistent with what has been said above,
  ! min_rledge is equivalent to min_rlvert_int - 2)
  INTEGER, PARAMETER :: end_all_verts = min_rlvert


  ! maximum allowed number of model domains (10 should be enough for the time being)
  INTEGER, PARAMETER :: max_dom = 10

  ! Maximum allowed number of physical model domains
  INTEGER, PARAMETER :: max_phys_dom = 30

  ! Maximum number of extra model levels added to the reduced radiation grid in case of vertical nesting
  INTEGER, PARAMETER :: nexlevs_rrg_vnest = 8

  ! maximum allowed number of tracers (20 should be enough for the time being) ! DRIEG: For ART, more than 20 tracers are needed
  INTEGER, PARAMETER :: max_ntracer = 200

  ! identifiers for model initialization
  INTEGER, PARAMETER :: ianalytic      =  0 ! - from analytical functions
  INTEGER, PARAMETER :: irestart       =  1 ! - from restart file

  ! identifiers for atm time stepping schemes
  INTEGER,PARAMETER :: TRACER_ONLY   = 1 ! pure tracer advection

  ! the hydrostatic model
  INTEGER,PARAMETER :: TWO_TL_SI     = 12 ! semi-implicit two time level
  INTEGER,PARAMETER :: LEAPFROG_EXPL = 13 ! explicit leapfrog
  INTEGER,PARAMETER :: LEAPFROG_SI   = 14 ! semi-implicit leapfrog
  INTEGER,PARAMETER :: RK4           = 15 ! standard 4th-order Runge-Kutta method
  INTEGER,PARAMETER :: SSPRK54       = 16 ! SSP RK(5,4)

  ! Scheme for the "slow" component in the TWO_TL_SI time stepping
  INTEGER,PARAMETER :: EULER_FORWARD = 1
  INTEGER,PARAMETER :: AB2           = 2

  ! Rayleigh damping identifiers
  INTEGER,PARAMETER :: RAYLEIGH_CLASSIC = 1  ! classical Rayleigh damping, which makes use of 
                                             ! a reference state.
  INTEGER,PARAMETER :: RAYLEIGH_KLEMP   = 2  ! Klemp (2008) type Rayleigh damping

  ! identifiers for physical processes (NWP)
  INTEGER, PARAMETER :: itconv   =  1
  INTEGER, PARAMETER :: itccov   =  2
  INTEGER, PARAMETER :: itrad    =  3
  INTEGER, PARAMETER :: itsso    =  4
  INTEGER, PARAMETER :: itgwd    =  5
  INTEGER, PARAMETER :: itsatad  =  6
  INTEGER, PARAMETER :: itturb   =  7
  INTEGER, PARAMETER :: itgscp   =  8
  INTEGER, PARAMETER :: itsfc    =  9
  INTEGER, PARAMETER :: itradheat=  10 !calculation of radiative heating rates from radiative
                                       !fluxes with updated solar zenith angle

  INTEGER, PARAMETER :: itfastphy=  6  !index for accessing fast physics

  INTEGER, PARAMETER :: iphysproc = 10! for NWP:
                                      ! number of physical processes:
                                      ! convection, cloud cover, radiation, radheat, sso,
                                      ! microphysics, saturation adjustment, gwd, 
                                      ! turbulence, surface

  INTEGER, PARAMETER :: iphysproc_short = 6 ! for NWP:
                                            ! number of physical processes:
                                            ! convection, cloud cover, radiation,
                                            ! sso, gwd, fastphysics
                                            ! i.e. fastphysics processes are treated 
                                            ! as a combined process


  ! identifiers for different NWP turbulent schemes
  INTEGER, PARAMETER :: icosmo  =  1
  INTEGER, PARAMETER :: igme    =  2
  INTEGER, PARAMETER :: iedmf   =  3
  INTEGER, PARAMETER :: ismag   =  5

  ! identifiers for aerosol classes of Tegen climatology 
  INTEGER, PARAMETER :: iss   =  1
  INTEGER, PARAMETER :: iorg  =  2
  INTEGER, PARAMETER :: ibc   =  3
  INTEGER, PARAMETER :: iso4  =  4
  INTEGER, PARAMETER :: idu   =  5
  INTEGER, PARAMETER :: nclass_aero = 5

  ! external parameter for radiation

  INTEGER, PARAMETER :: io3_interact =  1
  INTEGER, PARAMETER :: io3_clim     =  2
  INTEGER, PARAMETER :: io3_ape      =  4
  INTEGER, PARAMETER :: io3_amip     =  8
  INTEGER, PARAMETER :: iaero_kinne  =  3
  INTEGER, PARAMETER :: io3_art      =  10

  ! identifier for landcover classification
  INTEGER, PARAMETER :: GLOBCOVER2009 =  1
  INTEGER, PARAMETER :: GLC2000       =  2

  !
  ! transport identifiers
  !
  ! identifier for horizontal transport scheme
  INTEGER, PARAMETER :: NO_HADV = 0
  INTEGER, PARAMETER :: UP      = 1
  INTEGER, PARAMETER :: MIURA   = 2
  INTEGER, PARAMETER :: MIURA3  = 3
  INTEGER, PARAMETER :: FFSL    = 4
  INTEGER, PARAMETER :: FFSL_HYB= 5
  INTEGER, PARAMETER :: UP3     = 6
  INTEGER, PARAMETER :: MCYCL   = 20
  INTEGER, PARAMETER :: MIURA_MCYCL  = 22
  INTEGER, PARAMETER :: MIURA3_MCYCL = 32
  INTEGER, PARAMETER :: FFSL_MCYCL = 42
  INTEGER, PARAMETER :: FFSL_HYB_MCYCL = 52

  ! identifier for vertical transport scheme
  INTEGER, PARAMETER :: NO_VADV     = 0
  INTEGER, PARAMETER :: iup_v       = 1
  INTEGER, PARAMETER :: ippm_vcfl   = 3
  INTEGER, PARAMETER :: ippm_v      = 30

  ! identifier for horizontal limiter
  INTEGER, PARAMETER :: inol       = 0
  INTEGER, PARAMETER :: ifluxl_m   = 3
  INTEGER, PARAMETER :: ifluxl_sm  = 4

  ! identifier for vertical limiter
  INTEGER, PARAMETER :: inol_v      = 0
  INTEGER, PARAMETER :: islopel_vsm = 1
  INTEGER, PARAMETER :: islopel_vm  = 2
  INTEGER, PARAMETER :: ifluxl_vpd  = 4

  ! identifier for upper boundary condition (ubc)
  INTEGER, PARAMETER :: ino_flx     = 0
  INTEGER, PARAMETER :: izero_grad  = 1
  INTEGER, PARAMETER :: iparent_flx = 2


  ! equations to be solved
  INTEGER, PARAMETER :: ihs_atm_temp   =  1 ! - hydrostatic atmosphere, T as progn. var.
  INTEGER, PARAMETER :: ihs_atm_theta  =  2 ! - hydrostatic atmosphere, Theta as progn. var.
  INTEGER, PARAMETER :: inh_atmosphere =  3 ! - non-hydrost.atm.
  INTEGER, PARAMETER :: ishallow_water =  0 ! - shallow water model
  INTEGER, PARAMETER :: ihs_ocean      = -1 ! - hydrostatic ocean

  ! cell geometry
  INTEGER, PARAMETER :: itri           =  3 ! - triangles
  INTEGER, PARAMETER :: ihex           =  6 ! - hexagons/pentagons

  ! parameterized forcing (right hand side) of dynamics
  INTEGER, PARAMETER :: inoforcing     =  0 ! - no forcing
                                            ! - atmosphere
  INTEGER, PARAMETER :: iheldsuarez    =  1 !   - Held-Suarez test
  INTEGER, PARAMETER :: iecham         =  2 !   - ECHAM physics
  INTEGER, PARAMETER :: inwp           =  3 !   - NWP physics
  INTEGER, PARAMETER :: ildf_dry       =  4 !   - local diabatic forcing test without physics
  INTEGER, PARAMETER :: ildf_echam     =  5 !   - local diabatic forcing test with physics
                                            ! - ocean
  INTEGER, PARAMETER :: impiom         = -1 !   - MPIOM physics

  

  ! NWP SST-ICE modes
  INTEGER, PARAMETER :: SSTICE_ANA         = 1     ! SST and sea ice read from analysis and kept constant
  INTEGER, PARAMETER :: SSTICE_ANA_CLINC   = 2     ! SST and sea ice read from analysis. SST is updated 
                                                   ! by climatological increments on a daily basis
  INTEGER, PARAMETER :: SSTICE_CLIM        = 3     ! SST and sea ice based on climatology (monthly fields)
  INTEGER, PARAMETER :: SSTICE_AVG_MONTHLY = 4     ! SST and sea ice based on monthly averages
  INTEGER, PARAMETER :: SSTICE_AVG_DAILY   = 5     ! SST and sea ice based on daily averages


  !---------------------!
  !        LAND         !
  !---------------------!

  ! full level heights [m]
  REAL(wp), PARAMETER, DIMENSION(8)::                               &
    & zml_soil=(/ 0.005_wp,0.02_wp,0.06_wp,0.18_wp,0.54_wp,1.62_wp, &
    & 4.86_wp,14.58_wp /)

  ! Soil layer thicknesses [m]
  REAL(wp), PARAMETER, DIMENSION(8)::                               &
    & dzsoil=(/ 0.01_wp,0.02_wp,0.06_wp,0.18_wp,0.54_wp,1.62_wp,    &
    & 4.86_wp,14.58_wp/)

  ! identifier for MODIS albedo
  INTEGER, PARAMETER :: MODIS   = 2


  !---------------------!
  !        OCEAN        !
  !---------------------!

  ! identifier for parameterized forcing of the ocean model (iforc_oce)
  INTEGER, PARAMETER :: analyt_stat    = 11   ! stationary harmonic wind forcing
  INTEGER, PARAMETER :: core_forc      = 12   ! forcing from CORE database
  INTEGER, PARAMETER :: core_annwind   = 13   ! annual mean CORE winds
  INTEGER, PARAMETER :: full_forc      = 14   ! mpiom-type forcing

  ! identifier for ocean model test cases (itestcase_oce)
  ! (should probably be moved to some testcase-module)
  INTEGER, PARAMETER :: oce_testcase_zero  =  0   ! no or zero forcing
  INTEGER, PARAMETER :: oce_testcase_init  = 21   ! simply defined test case
  INTEGER, PARAMETER :: oce_testcase_file  = 22   ! test case read from file

  ! number of tracers used in ocean state
  INTEGER, PARAMETER :: ntrac_oce = 2
  ! ocean surface level
  INTEGER, PARAMETER :: toplev  = 1
  INTEGER, PARAMETER :: MIN_DOLIC=2  !mimal number of vertical layers that have to be present in 3D
                                     !ocean. Not relevant for shallow-water.

  !---------------------!
  !  PARALLELIZATION    !
  !---------------------!

  ! Division method for area subdivision
  INTEGER, PARAMETER :: div_geometric = 1  ! Geometric subdivision

  !-----  horizontal interpolation: type of interpolation
  CHARACTER(len=32), PARAMETER :: STR_HINTP_TYPE(4) = &
    (/ "HINTP_TYPE_NONE       ",  &
    &  "HINTP_TYPE_LONLAT_RBF ",  &
    &  "HINTP_TYPE_LONLAT_NNB ",  &
    &  "HINTP_TYPE_LONLAT_BCTR" /)
  INTEGER, PARAMETER :: HINTP_TYPE_NONE        = 1
  INTEGER, PARAMETER :: HINTP_TYPE_LONLAT_RBF  = 2
  INTEGER, PARAMETER :: HINTP_TYPE_LONLAT_NNB  = 3
  INTEGER, PARAMETER :: HINTP_TYPE_LONLAT_BCTR = 4

  !-----  vertical interpolation algorithms
  INTEGER, PARAMETER :: VINTP_METHOD_VN    = 1
  INTEGER, PARAMETER :: VINTP_METHOD_LIN   = 2
  INTEGER, PARAMETER :: VINTP_METHOD_QV    = 3
  INTEGER, PARAMETER :: VINTP_METHOD_PRES  = 4
  INTEGER, PARAMETER :: VINTP_METHOD_LIN_NLEVP1 = 5

  !----- RBF lon-lat interpolation: shape parameter mode
  INTEGER, PARAMETER :: SCALE_MODE_TABLE   = 1
  INTEGER, PARAMETER :: SCALE_MODE_AUTO    = 2
  INTEGER, PARAMETER :: SCALE_MODE_PRESET  = 3

  !----- init ICON operation modes -----
  INTEGER, PARAMETER :: MODE_DWDANA      = 1
  INTEGER, PARAMETER :: MODE_IFSANA      = 2
  INTEGER, PARAMETER :: MODE_COMBINED    = 3
  INTEGER, PARAMETER :: MODE_COSMO       = 4
  INTEGER, PARAMETER :: MODE_IAU         = 5
  INTEGER, PARAMETER :: MODE_IAU_OLD     = 6
  INTEGER, PARAMETER :: MODE_ICONVREMAP  = 7

  !----- MPI parallelization -----
  INTEGER, PARAMETER :: MAX_NUM_IO_PROCS = 100      !< max. number of output ranks

  !----------------------!
  !  VARIABLE DATA TYPES !
  !----------------------!

  INTEGER, PARAMETER :: REAL_T   = 1
  INTEGER, PARAMETER :: SINGLE_T = 2
  INTEGER, PARAMETER :: BOOL_T   = 3
  INTEGER, PARAMETER :: INT_T    = 4

  !----------------!
  !  MODEL OUTPUT  !
  !----------------!

  ! maximum string length for variable names
  INTEGER, PARAMETER :: VARNAME_LEN = 256 

  INTEGER, PARAMETER :: &
    max_var_lists  = 256, & ! max number of output var_lists
    MAX_NVARS      = 999, & ! maximum number of output variables (total)
    max_var_ml     = 999, & ! maximum number of output model-level variables
    max_var_pl     = 150, & ! maximum number of pressure-level variables
    max_var_hl     = 150, & ! maximum number of height-level variables
    max_var_il     = 150, & ! maximum number of variables on isentropes
    vname_len      = VARNAME_LEN ! variable name length in I/O namelists

  INTEGER, PARAMETER :: &
    MAX_TIME_INTERVALS = 10 ! maximum number of time intervals specified in "output_nml"

  ! Method for computation of mean sea level pressure:
  INTEGER, PARAMETER :: &
    PRES_MSL_METHOD_GME = 1,  &   ! GME-type extrapolation
    PRES_MSL_METHOD_SAI = 2,  &   ! stepwise analytical integration 
    PRES_MSL_METHOD_IFS = 3,  &   ! current IFS method
    PRES_MSL_METHOD_IFS_CORR = 4,&! modified IFS method that is consistent with 
                                  ! geopotential computation 
    PRES_MSL_METHOD_DWD = 5       ! mixture between GME and IFS method (elevation-dependent departure 
                                  ! level for downward extraplation)

  ! Method for computation of relative humidity:
  INTEGER, PARAMETER :: &
    RH_METHOD_WMO      = 1,  &   ! WMO-type: e_s=e_s_water (water only)
    RH_METHOD_IFS      = 2,  &   ! IFS-type: (mixed phases, water and ice)
    RH_METHOD_IFS_CLIP = 3       ! IFS-type + clipping to rh<=100%

  ! Max number of time levels:
  INTEGER, PARAMETER :: MAX_TIME_LEVELS = 5

  INTEGER, PARAMETER :: MAX_NPLEVS = 100 !< max. no. of pressure levels
  INTEGER, PARAMETER :: MAX_NZLEVS = 100 !< max. no. of height levels
  INTEGER, PARAMETER :: MAX_NILEVS = 100 !< max. no. of isentropic levels

  !-----------------------------------!
  !  POST PROCESSING SCHEDULER TASKS  !
  !-----------------------------------!

  INTEGER, PARAMETER, PUBLIC :: TASK_NONE              = 0 
  !------ setup tasks (coefficients,...)
  INTEGER, PARAMETER, PUBLIC :: TASK_INIT_VER_Z        = 1  !< task: setup z-interpolation
  INTEGER, PARAMETER, PUBLIC :: TASK_INIT_VER_P        = 2  !< task: setup p-interpolation
  INTEGER, PARAMETER, PUBLIC :: TASK_INIT_VER_I        = 3  !< task: setup i-interpolation
  INTEGER, PARAMETER, PUBLIC :: TASK_FINALIZE_IPZ      = 5  !< task: deallocate ipz-interpolation
  !------ interpolation tasks:
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_HOR_LONLAT   = 6  !< task: lon-lat
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_VER_PLEV     = 7  !< task: vertical p-levels
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_VER_ZLEV     = 8  !< task: vertical z-levels
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_VER_ILEV     = 9  !< task: vertical isentropic levels
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_SYNC         = 10 !< task: synchronizes halo regions
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_MSL          = 11 !< task: intp. to mean sea level
  INTEGER, PARAMETER, PUBLIC :: TASK_INTP_EDGE2CELL    = 12 !< task: intp. from edge midpoints to cell centers
  !------ computation of optional diagnostic fields
  INTEGER, PARAMETER, PUBLIC :: TASK_COMPUTE_RH        = 13 !< task: compute relative humidity
  INTEGER, PARAMETER, PUBLIC :: TASK_COMPUTE_OMEGA     = 14 !< task: compute vertical velocity
  INTEGER, PARAMETER, PUBLIC :: TASK_COMPUTE_PV        = 15 !< task: compute potential vorticity
  INTEGER, PARAMETER, PUBLIC :: TASK_COMPUTE_SMI       = 16 !< task: compute soil moisture index

  !--------------------------------------------------------------------!
  !  VARIABLE TIMELEVEL SPECIFICATION (FOR POST-PROCESSING SCHEDULER)  !
  !--------------------------------------------------------------------!

  INTEGER, PARAMETER, PUBLIC :: UNDEF_TIMELEVEL        = -2 !< means "undefined time level"
  INTEGER, PARAMETER, PUBLIC :: ALL_TIMELEVELS         = -1 !< means "active at all time levels"

  !-----------------------------!
  !  PROFILING OUTPUT VERBOSITY !
  !-----------------------------!

  INTEGER, PARAMETER, PUBLIC :: TIMER_MODE_AGGREGATED  = 1
  INTEGER, PARAMETER, PUBLIC :: TIMER_MODE_DETAILED    = 2
  INTEGER, PARAMETER, PUBLIC :: TIMER_MODE_WRITE_FILES = 3

  !-------------------------------------------------!
  !  TIME LEVEL SOURCE (WHICH "NNOW"/"NNEW" TO USE) !
  !-------------------------------------------------!

  INTEGER, PARAMETER, PUBLIC :: TLEV_NNOW     = 0
  INTEGER, PARAMETER, PUBLIC :: TLEV_NNOW_RCF = 1
  INTEGER, PARAMETER, PUBLIC :: TLEV_NNEW     = 2
  INTEGER, PARAMETER, PUBLIC :: TLEV_NNEW_RCF = 3

  !-------------------------!
  !  RTTOV FIELD CATEGORIES !
  !-------------------------!

  INTEGER, PARAMETER, PUBLIC :: RTTOV_BT_CL  = 1
  INTEGER, PARAMETER, PUBLIC :: RTTOV_BT_CS  = 2
  INTEGER, PARAMETER, PUBLIC :: RTTOV_RAD_CL = 3
  INTEGER, PARAMETER, PUBLIC :: RTTOV_RAD_CS = 4

  !------------------------!
  !  CALENDAR TYPES        !
  !------------------------!

  INTEGER,  PARAMETER :: julian_gregorian    = 0 !< historic Julian / Gregorian
  INTEGER,  PARAMETER :: proleptic_gregorian = 1 !< proleptic Gregorian
  INTEGER,  PARAMETER :: cly360              = 2 !< constant 30 dy/mo and 360 dy/yr

  !------------------------------------------------!
  !  MISSING VALUE FOR BOUNDARY INTERPOLATION ZONE !
  !------------------------------------------------!

  REAL(WP), PARAMETER, PUBLIC :: BOUNDARY_MISSVAL = -999.e-10

  !------------------------------------------------!
  !  MISSING VALUE FOR NWP SEAICE ALBEDO           !
  !------------------------------------------------!

  REAL(wp), PARAMETER, PUBLIC :: ALB_SI_MISSVAL = -1._wp

  ! The lon-lat parameterization of the torus is 
  !    (lon,lat) = [0, 2*pi] x [-max_lat, max_lat]
  ! where max_lat := pi/180 = 10 degrees 
  ! (hard-coded in the torus grid generator)
  REAL(wp), PARAMETER :: TORUS_MAX_LAT = 4._wp / 18._wp * ATAN(1._wp)

!--------------------------------------------------------------------
END MODULE mo_impl_constants
