
!-----------------------------------------------------------------------
! access definitions for arrays
#define get_cell_index(i,j)     a_cell_index(j,i)
#define get_edge_index(i,j)     a_edge_index(j,i)
#define get_vertex_index(i,j)   a_vertex_index(j,i)
#define get_neighbor_index(i,j) a_neighbor_index(j,i) 
#define get_edge_orient(i,j)    a_edge_orientation(j,i)
#define get_edge_cell_length(i,j) a_edge_cell_length(j,i)
#define get_edge_vert_length(i,j) a_edge_vert_length(j,i)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Basic geometry functions definitions
#define d_norma_3d(v) SQRT(DOT_PRODUCT(v%x,v%x))
#define d_normalize(v) v%x=v%x/d_norma_3d(v)
#define d_arc_of_hord_normalsphere(l) (2.0_wp*ASIN((l)*0.5_wp))

!-----------------------------------------------------------------------
