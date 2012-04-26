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
  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  ! define model name and version
  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  INTEGER, PARAMETER :: MAX_CHAR_LENGTH     = 1024

  INTEGER, PARAMETER :: SUCCESS             = 0
  INTEGER, PARAMETER :: CELLS               = 123
  INTEGER, PARAMETER :: EDGES               = 345
  INTEGER, PARAMETER :: VERTS               = 678
  
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
  ! A. The indexes from 1 to max_rl
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
  ! C. The indexes from -1 to min_rl_int
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
  ! D. The indexes from min_rl_int-1 to min_rl
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
  INTEGER, PARAMETER :: min_rledge_int = 2*min_rlcell_int
  INTEGER, PARAMETER :: min_rledge     = min_rledge_int - (2*max_hw+1)
  INTEGER, PARAMETER :: max_rledge     = 2*max_rlcell

  ! maximum allowed number of model domains (10 should be enough for the time being)
  INTEGER, PARAMETER :: max_dom = 10

  ! Maximum allowed number of physical model domains
  INTEGER, PARAMETER :: max_phys_dom = 30

  ! maximum allowed number of tracers (20 should be enough for the time being)
  INTEGER, PARAMETER :: max_ntracer = 20

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

  ! the non-hydrostatic model
  INTEGER,PARAMETER :: MATSUNO_DEF  = 4 !31 future enumeration to make the belonging clear,  Matsuno scheme
  INTEGER,PARAMETER :: MATSUNO_COR  = 3 !32 Matsuno, comp of velocity tendencies on corretor step only
  INTEGER,PARAMETER :: MATSUNO_UNK  = 5 !34 Matsuno, variation unknown
  INTEGER,PARAMETER :: MATSUNO_AVE  = 6 !33 Matsuno with velocitiy tendendcies averaged over 2 time steps


  ! Scheme for the "slow" component in the TWO_TL_SI time stepping
  INTEGER,PARAMETER :: EULER_FORWARD = 1
  INTEGER,PARAMETER :: AB2           = 2

  ! identifiers for NWP time control variables lcall_phy, dt_phy, t_elapsed_phy
  INTEGER, PARAMETER :: itconv   =  1
  INTEGER, PARAMETER :: itccov   =  2
  INTEGER, PARAMETER :: itrad    =  3
  INTEGER, PARAMETER :: itsso    =  4
  INTEGER, PARAMETER :: itgwd    =  5
  INTEGER, PARAMETER :: itupdate =  6
  INTEGER, PARAMETER :: itsatad  =  7
  INTEGER, PARAMETER :: itturb   =  8
  INTEGER, PARAMETER :: itgscp   =  9
  INTEGER, PARAMETER :: itsfc    =  10
  INTEGER, PARAMETER :: itradheat=  11 !calculation of radiative heating rates from radiative
                                       !fluxes with updated solar zenith angle
  INTEGER, PARAMETER :: itfastphy=  6

  INTEGER, PARAMETER :: iphysproc = 11! for NWP:
                                      ! number of physical processes:
                                      ! convection, cloud cover, radiation, radheat, sso,
                                      ! microphysics, saturation adjustment, tracerupdate, 
                                      ! gwd, turbulence, surface

  INTEGER, PARAMETER :: iphysproc_short = 6 ! for NWP:
                                            ! number of physical processes:
                                            ! convection, cloud cover, radiation,
                                            ! sso, gwd, fastphysics
                                            ! i.e. fastphysics processes are treated 
                                            ! as a combined process



  ! external parameter for radiation

  INTEGER, PARAMETER :: io3_clim     =  2
  INTEGER, PARAMETER :: io3_ape      =  4
  INTEGER, PARAMETER :: iaero_kinne  =  3

  !
  ! transport identifiers
  !
  ! identifier for horizontal transport scheme
  INTEGER, PARAMETER :: NO_HADV = 0
  INTEGER, PARAMETER :: UP      = 1
  INTEGER, PARAMETER :: MIURA   = 2
  INTEGER, PARAMETER :: MIURA3  = 3
  INTEGER, PARAMETER :: FFSL    = 4
  INTEGER, PARAMETER :: UP3     = 5
  INTEGER, PARAMETER :: MCYCL   = 20
  INTEGER, PARAMETER :: MIURA_MCYCL  = 22
  INTEGER, PARAMETER :: MIURA3_MCYCL = 32

  ! identifier for vertical transport scheme
  INTEGER, PARAMETER :: ino_vadv    = 0
  INTEGER, PARAMETER :: iup_v       = 1
  INTEGER, PARAMETER :: imuscl_vcfl = 2
  INTEGER, PARAMETER :: imuscl_v    = 20
  INTEGER, PARAMETER :: ippm_vcfl   = 3
  INTEGER, PARAMETER :: ippm_v      = 30

  ! identifier for horizontal limiter
  INTEGER, PARAMETER :: inol       = 0
  INTEGER, PARAMETER :: islopel_sm = 1
  INTEGER, PARAMETER :: islopel_m  = 2
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

  
  ! auxiliary parameter to access single field of the 4D array prm_diag%tot_cld
  INTEGER, PARAMETER :: icc = 4    !! diagnostic cloud fraction in prm_diag%tot_cld


  !---------------------!
  !        LAND         !
  !---------------------!

  ! full level heights [m]
  REAL(wp), PARAMETER, DIMENSION(8)::                               &
    & zml_soil=(/ 0.005_wp,0.02_wp,0.06_wp,0.18_wp,0.54_wp,1.62_wp, &
    & 4.86_wp,14.58_wp /)

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
  INTEGER, PARAMETER :: div_from_file = 0  ! Read from file
  INTEGER, PARAMETER :: div_geometric = 1  ! Geometric subdivision
  INTEGER, PARAMETER :: div_metis     = 2  ! Use Metis

  !--------------------------!
  !  VERTICAL INTERPOLATION  !
  !--------------------------!

  !-----  vertical interpolation: type of interpolation
  INTEGER, PARAMETER :: VINTP_TYPE_NONE    = 0
  INTEGER, PARAMETER :: VINTP_TYPE_Z       = 1
  INTEGER, PARAMETER :: VINTP_TYPE_P_OR_Z  = 2
  !-----  horizontal interpolation: type of interpolation
  INTEGER, PARAMETER :: HINTP_TYPE_NONE    = 0
  INTEGER, PARAMETER :: HINTP_TYPE_LONLAT  = 1
  !-----  vertical interpolation algorithms
  INTEGER, PARAMETER :: VINTP_METHOD_UV    = 1
  INTEGER, PARAMETER :: VINTP_METHOD_LIN   = 2
  INTEGER, PARAMETER :: VINTP_METHOD_QV    = 3
  INTEGER, PARAMETER :: VINTP_METHOD_PRES  = 4
  INTEGER, PARAMETER :: VINTP_METHOD_LIN_NLEVP1 = 5

  !----------------!
  !  MODEL OUTPUT  !
  !----------------!

  INTEGER, PARAMETER :: &
    max_var_lists  = 256, & ! max number of output var_lists
    MAX_NVARS      = 999, & ! maximum number of output variables (total)
    max_var_ml     = 400, & ! maximum number of output model-level variables
    max_var_pl     = 100, & ! maximum number of pressure-level variables
    max_var_hl     = 100, & ! maximum number of height-level variables
    max_bounds     = 100, & ! maximum number of output_bounds
    max_levels     = 100, & ! maximum number of pressure/height levels
    vname_len      =  32    ! variable name length in I/O namelists


!--------------------------------------------------------------------
END MODULE mo_impl_constants
