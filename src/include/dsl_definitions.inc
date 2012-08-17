
!--------------------- Notes on parsing ---------------------------------
! define a rule:
! !define_rule <name>: {current order} -> {new order}
!
! note: ID_RULE means the indexes will be the same: {current order} -> {current order}
!
! after every #define of an array structure follows the rule to be used for handling the structure
!
! An example of parsing is the files example_old.f90 -> example_new.f90
!
! Notes:
!   1. structures inhereted by "use" or "types" should be treated with the same rules !
!   2. the #include dimensions.inc line should be omitted when parsing
!
!--------------------- End of Notes --------------------------------------


!--------------------- define the rules --------------------------
!define_rule ID_RULE: 
!define_rule 3D:       {1,2,3}   ->   {2,1,3}
!define_rule CONNECTIVITY:  {1,2,3}   ->   {3,1,2}
!define_rule 3D_BLOCK:    {1,2}     ->   {2,1}
!define_rule 3D_TRACE: {1,2,3,4} ->   {2,1,3,4}
!define_rule 2D_COEFFS_IN_2: {1,2,3} -> {2,1,3}
!define_rule 2D_2COEFFS_IN_2_4: {1,2,3,4} -> {2, 4, 1, 3}
!define_rule 3D_COEFFS_IN_1:    {1,2,3,4} -> {1,3,2,4}
!define_rule 3D_COEFFS_IN_2:    {1,2,3,4} -> {2,3,1,4}
!define_rule 3D_COEFFS_IN_4:    {1,2,3,4} -> {4,2,1,3}
!define_rule 2D_2COEFFS_IN_3_4: {1,2,3,4} -> {4,3,1,2}
!define_rule 4D_TIME:           {1,2,3,4} -> {4,2,1,3}

!--------------------- positions --------------------------
#define NPROMA_POSITION 1
#define LEVELS_POSITION 2
#define BLOCKS_POSITION 3

!--------------------- general --------------------------
! general 2D structure, maybe cells, verts or edges
#define GENERAL_2D_SIZE(nproma, blocks) DIMENSION(nproma, blocks)
!use_rule ID_RULE

! general 3D structure, maybe cells, verts or edges
#define GENERAL_3D_SIZE(nproma, levels, blocks) DIMENSION(nproma, levels, blocks)
!use_rule 3D

! general 3D structure, maybe cells, verts or edges, currenlty (nproma, levels, blocks)
#define GENERAL_3D DIMENSION(:,:,:)
!use_rule 3D

! general connectivity, maybe cells, verts or edges, currenlty (nproma, blocks, neigbors)
#define GENERAL_CONNECTIVITY DIMENSION(:,:,:)
!use_rule CONNECTIVITY

! general one 3D_BLOCK, currenlty defined as (nproma, levels)
#define BLOCK_3D_SIZE(nproma, levels) DIMENSION(nproma, levels)
!use_rule 3D_BLOCK

! the general 3D with tracers, currenlty defined as (nproma, levels, blocks, tracers)
#define GENERAL_3D_TRACE DIMENSION(:,:,:,:)
!use_rule 3D_TRACE

!--------------------- vertical (z_axis) --------------------------
 
! the 1D vertical dimension,  currenlty defined as (levels)
#define VERTICAL_1D DIMENSION(:)
!use_rule ID_RULE
    
!--------------------- cells --------------------------
! the 2D cells, currenlty defined as (nproma, blocks)
#define CELLS_2D DIMENSION(:,:)
!use_rule ID_RULE

! the 2D cells coefficients, currenlty defined as (nproma, connectivity, blocks), should be changed
#define CELLS_2D_COEFFS_IN_2 DIMENSION(:,:,:)
!use_rule 2D_COEFFS_IN_2

! the 2D cells 2 coefficients, currenlty defined as (nproma, 1st coeff, blocks, 2nd coeff),
!  should be (1st coeff,2nd coeff, nproma, blocks)
#define CELLS_2D_2COEFFS_IN_2_4 DIMENSION(:,:,:,:)
!use_rule 2D_2COEFFS_IN_2_4

! the 3D cells, currenlty defined as (nproma, levels, blocks)
#define CELLS_3D DIMENSION(:,:,:)
!use_rule 3D

! the 4D cells, (space+time) currenlty defined as (nproma, levels, blocks, timestep 1-3)
#define CELLS_4D_TIMESTEPS DIMENSION(:,:,:,:)
!use_rule 4D_TIME

! the 3D cells with tracers, currenlty defined as (nproma, levels, blocks, tracers)
#define CELLS_3D_TRACE DIMENSION(:,:,:,:)
!use_rule 3D_TRACE

! the 3D cells with 4 variables calculated at the same time.currenlty defined as nproma, vars, levels, blocks)
! vars should be in the 1st place
#define CELLS_3D_4VARS_IN_2  DIMENSION(:,:,:,:)
!use_rule 3D_COEFFS_IN_2

! the 2D cells with defined size
#define CELLS_2D_SIZE(nproma, blocks) DIMENSION(nproma, blocks)
!use_rule ID_RULE

! the 3D cells with defined size
#define CELLS_3D_SIZE(nproma, levels, blocks) DIMENSION(nproma, levels, blocks)
!use_rule 3D

! the 3D cells with tracers (or other valibales that are independent) with defined size,
! currenlty defined as (nproma, levels, blocks, tracers)
#define CELLS_3D_TRACE_SIZE(nproma, levels, blocks, tracers) DIMENSION(nproma, levels, blocks, tracers)
!use_rule 3D_TRACE

! the 3D cells and ORDER dimension, with defined size
#define CELLS_3D_ORDER_SIZE(order, nproma, levels, blocks) DIMENSION(order, nproma, levels, blocks)
!use_rule 3D_COEFFS_IN_1

! the 3D cells with 4 variables calculated at the same time, with defined size. vars should be in the 1st place
#define CELLS_3D_4VARS_IN_2_SIZE(nproma, vars, levels, blocks) DIMENSION(nproma, vars, levels, blocks)
!use_rule 3D_COEFFS_IN_2

! the 2D cells quadratic coefficients, currenlty defined as (nproma, 3, blocks)
#define CELLS_2D_QUADRATIC  DIMENSION(:,:,:)
!use_rule 2D_COEFFS_IN_2

! predefined shapes 
#define SHAPE_CELLS_2D(nproma, blocks) (/nproma, blocks /)
!use_rule ID_RULE


#define SHAPE_CELLS_2D_QUADRATIC(nproma, coeff, blocks) (/nproma, coeff, blocks /)
!use_rule 2D_COEFFS_IN_2

#define SHAPE_CELLS_3D(nproma, levels, blocks) (/nproma,  levels,  blocks /)
!use_rule 3D

#define SHAPE_CELLS_2D_TRACE(nproma, blocks, vars) (/nproma, blocks, vars /)
!use_rule ID_RULE

#define SHAPE_CELLS_3D_TRACE(nproma, levels, blocks, vars) (/nproma,  levels,  blocks, vars /)
!use_rule 3D_TRACE

! connectivity
! cells to neigboring cells, currenlty (nproma, blocks, neigbor 1-3)
#define CELLS_CONNECT_TO_CELLS   DIMENSION(:,:,:)
!use_rule CONNECTIVITY

! cells to cell edges, currenlty (nproma, blocks, neigbor 1-3)
#define CELLS_CONNECT_TO_EDGES   DIMENSION(:,:,:)
!use_rule CONNECTIVITY

! cells to cell verts, currenlty (nproma, blocks, neigbor 1-3)
#define CELLS_CONNECT_TO_VERTS   DIMENSION(:,:,:)
!use_rule CONNECTIVITY



!--------------------- edges --------------------------
! the 2D edges, currenlty defined as (nproma, blocks)
#define EDGES_2D DIMENSION(:,:)
!use_rule ID_RULE

! the 2D edges coefficients, currenlty defined as (nproma, connectivity, blocks), should be changed
#define EDGES_2D_COEFFS_IN_2 DIMENSION(:,:,:)
!use_rule 2D_COEFFS_IN_2

! 2D edges 4 coefficients, currenlty defined as (nproma, blocks, 1st, 2nd),
! should be ( 2nd,  1st, nproma, blocks)
#define EDGES_2D_2COEFFS_IN_3_4 DIMENSION(:,:,:,:)
!use_rule 2D_2COEFFS_IN_3_4

! the 3D edges, currenlty defined as (nproma, levels, blocks)
#define EDGES_3D DIMENSION(:,:,:)
!use_rule 3D

! the 4D edges, (space+time) currenlty defined as (nproma, levels, blocks, timestep 1-3)
#define EDGES_4D_TIMESTEPS DIMENSION(:,:,:,:)
!use_rule 4D_TIME

! the 2D edges quadratic coefficients, currenlty defined as (nproma, 3, blocks)
#define EDGES_2D_QUADRATIC  DIMENSION(:,:,:)
!use_rule 2D_COEFFS_IN_2

! the 3D edges with tracers the 4th dimension: currenlty defined as (nproma, levels, blocks, tracers)
#define EDGES_3D_TRACE DIMENSION(:,:,:,:)
!use_rule 3D_TRACE

! the 3D edges with 2-connectivity in the 4th dimension: currenlty defined as (nproma, levels, blocks, connectivity)
#define EDGES_3D_2VARS_IN_4 DIMENSION(:,:,:,:)
!use_rule 3D_COEFFS_IN_4

! the 3D edges with 2-connectivity in the 1st dimension: currenlty defined as (connectivity, nproma, levels, blocks)
#define EDGES_3D_COEFF_IN_1 DIMENSION(:,:,:,:)
!use_rule 3D_COEFFS_IN_1

! the 2D edges with defined size
#define EDGES_2D_SIZE(nproma, blocks) DIMENSION(nproma, blocks)
!use_rule ID_RULE

! the 3D edges with defined size
#define EDGES_3D_SIZE(nproma, levels, blocks) DIMENSION(nproma, levels, blocks)
!use_rule 3D

! the 3D edges in one 3D_BLOCK, currenlty defined as (nproma, levels)
#define EDGES_BLOCK_3D_SIZE(nproma, levels) DIMENSION(nproma, levels)
!use_rule 3D_BLOCK

! the 3D edges with 2-connectivity in the 4th dimension: currenlty defined as (nproma, levels, blocks, connectivity)
#define EDGES_3D_2VARS_IN_4_SIZE(nproma, levels, blocks, neigbors) DIMENSION(nproma, levels, blocks, neigbors)
!use_rule 3D_COEFFS_IN_4

! predefined shapes 
#define SHAPE_EDGES_2D(nproma, blocks) (/nproma, blocks /)
!use_rule ID_RULE

#define SHAPE_EDGES_2D_QUADRATIC(nproma, coeff, blocks) (/nproma, coeff, blocks /)
!use_rule 2D_COEFFS_IN_2

#define SHAPE_EDGES_3D(nproma, levels, blocks) (/nproma,  levels,  blocks /)
!use_rule 3D

#define SHAPE_EDGES_2D_TRACE(nproma, blocks, vars) (/nproma, blocks, vars /)
!use_rule ID_RULE

#define SHAPE_EDGES_3D_TRACE(nproma, levels, blocks, vars) (/nproma,  levels,  blocks, vars /)
!use_rule 3D_TRACE

#define SHAPE_EDGES_3D_COEFF_IN_1(coeff, nproma, levels, blocks)(/coeff, nproma,  levels,  blocks /)
!use_rule 3D_COEFFS_IN_1

#define SHAPE_EDGES_4D_TIMESTEP(nproma, levels, blocks, time) (/nproma,  levels,  blocks, time /)
!use_rule 4D_TIME

! connectivity
! edges to neigboring cells, currenlty (nproma, blocks, neigbor 1-2)
#define EDGES_CONNECT_TO_CELLS   DIMENSION(:,:,:)
!use_rule CONNECTIVITY

! edges to its verts, currenlty (nproma, blocks, neigbor 1-2)
#define EDGES_CONNECT_TO_VERTS   DIMENSION(:,:,:)
!use_rule CONNECTIVITY


! edges to the containing quads, currenlty (nproma, blocks, 1-4)
#define EDGES_CONNECT_TO_QUADS  DIMENSION(:,:,:)
!use_rule CONNECTIVITY

!--------------------- verts --------------------------
! the 2D verts, currenlty defined as (nproma, blocks)
#define VERTS_2D DIMENSION(:,:)
!use_rule ID_RULE

! the 3D verts, currenlty defined as (nproma, levels, blocks)
#define VERTS_3D DIMENSION(:,:,:)
!use_rule 3D

! the 2D verts coefficients, currenlty defined as (nproma, connectivity, blocks), should be changed
#define VERTS_2D_COEFFS_IN_2 DIMENSION(:,:,:)
!use_rule 2D_COEFFS_IN_2

! the 3D verts with defined size
#define VERTS_3D_SIZE(nproma, levels, blocks) DIMENSION(nproma, levels, blocks)
!use_rule 3D

! predefined shapes 
#define SHAPE_VERTS_3D(nproma, levels, blocks) (/nproma,  levels,  blocks /)
!use_rule 3D

! connectivity
! verts to neigboring cells, currenlty (nproma, blocks, neigbor 1-6)
#define VERTS_CONNECT_TO_CELLS   DIMENSION(:,:,:)
!use_rule CONNECTIVITY

! verts to neigboring edges, currenlty (nproma, blocks, neigbor 1-6)
#define VERTS_CONNECT_TO_EDGES   DIMENSION(:,:,:)
!use_rule CONNECTIVITY

! verts to neigboring verts, currenlty (nproma, blocks, neigbor 1-6)
#define VERTS_CONNECT_TO_VERTS   DIMENSION(:,:,:)
!use_rule CONNECTIVITY


!-----------------------------------------------

