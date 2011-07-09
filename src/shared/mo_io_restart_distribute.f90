#define NOMPI

MODULE mo_io_restart_distribute
  !
  USE mo_kind,          ONLY: wp
#ifndef NOMPI
  USE mo_mpi,           ONLY: p_io, p_bcast, p_comm_work
  USE mo_model_domain,  ONLY: t_patch
  USE mo_communication, ONLY: idx_no, blk_no, exchange_data
#endif
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: gather_cells
  PUBLIC :: gather_edges
  PUBLIC :: gather_vertices
  !
  PUBLIC :: scatter_cells
  PUBLIC :: scatter_edges
  PUBLIC :: scatter_vertices
  !
  !------------------------------------------------------------------------------------------------
  !
  INTERFACE gather_cells
    MODULE PROCEDURE gather_cells_r2d
    MODULE PROCEDURE gather_cells_r3d
    MODULE PROCEDURE gather_cells_i2d
    MODULE PROCEDURE gather_cells_i3d
    MODULE PROCEDURE gather_cells_l2d
    MODULE PROCEDURE gather_cells_l3d
  END INTERFACE gather_cells
  !
  INTERFACE gather_edges
    MODULE PROCEDURE gather_edges_r2d
    MODULE PROCEDURE gather_edges_r3d
    MODULE PROCEDURE gather_edges_i2d
    MODULE PROCEDURE gather_edges_i3d
    MODULE PROCEDURE gather_edges_l2d
    MODULE PROCEDURE gather_edges_l3d
  END INTERFACE gather_edges
  !
  INTERFACE gather_vertices
    MODULE PROCEDURE gather_vertices_r2d
    MODULE PROCEDURE gather_vertices_r3d
    MODULE PROCEDURE gather_vertices_i2d
    MODULE PROCEDURE gather_vertices_i3d
    MODULE PROCEDURE gather_vertices_l2d
    MODULE PROCEDURE gather_vertices_l3d
  END INTERFACE gather_vertices
  !
  INTERFACE scatter_cells
    MODULE PROCEDURE scatter_cells_r2d
    MODULE PROCEDURE scatter_cells_r3d
    MODULE PROCEDURE scatter_cells_i2d
    MODULE PROCEDURE scatter_cells_i3d
    MODULE PROCEDURE scatter_cells_l2d
    MODULE PROCEDURE scatter_cells_l3d
  END INTERFACE scatter_cells
  !
  INTERFACE scatter_edges
    MODULE PROCEDURE scatter_edges_r2d
    MODULE PROCEDURE scatter_edges_r3d
    MODULE PROCEDURE scatter_edges_i2d
    MODULE PROCEDURE scatter_edges_i3d
    MODULE PROCEDURE scatter_edges_l2d
    MODULE PROCEDURE scatter_edges_l3d
  END INTERFACE scatter_edges
  !
  INTERFACE scatter_vertices
    MODULE PROCEDURE scatter_vertices_r2d
    MODULE PROCEDURE scatter_vertices_r3d
    MODULE PROCEDURE scatter_vertices_i2d
    MODULE PROCEDURE scatter_vertices_i3d
    MODULE PROCEDURE scatter_vertices_l2d
    MODULE PROCEDURE scatter_vertices_l3d
  END INTERFACE scatter_vertices
  !------------------------------------------------------------------------------------------------
  !
  INTERFACE reorder
    MODULE PROCEDURE reorder_foreward_r2d
    MODULE PROCEDURE reorder_foreward_r3d
    MODULE PROCEDURE reorder_backward_r2d
    MODULE PROCEDURE reorder_backward_r3d
    MODULE PROCEDURE reorder_foreward_i2d
    MODULE PROCEDURE reorder_foreward_i3d
    MODULE PROCEDURE reorder_backward_i2d
    MODULE PROCEDURE reorder_backward_i3d
    MODULE PROCEDURE reorder_foreward_l2d
    MODULE PROCEDURE reorder_foreward_l3d
    MODULE PROCEDURE reorder_backward_l2d
    MODULE PROCEDURE reorder_backward_l3d
  END INTERFACE reorder
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE gather_cells_r2d(in_array, out_array, p_patch)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r1d(:)
    !
#ifndef NOMPI
    REAL(wp), ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_cells_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_c, RECV=tmp_array, SEND=in_array) 
    r1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, r1d)
    !
    DEALLOCATE(tmp_array)
#else
    r1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, r1d)
#endif
  END SUBROUTINE gather_cells_r2d
  !
  SUBROUTINE gather_cells_r3d(in_array, out_array, p_patch)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r2d(:,:)
    !
#ifndef NOMPI
    REAL(wp), ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_cells_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_c, RECV=tmp_array, SEND=in_array) 
    r2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, r2d)
    !
    DEALLOCATE(tmp_array)
#else
    r2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, r2d)
#endif
  END SUBROUTINE gather_cells_r3d
  !
  SUBROUTINE gather_edges_r2d(in_array, out_array, p_patch)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r1d(:)
    !
#ifndef NOMPI
    REAL(wp), ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_edges_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_e, RECV=tmp_array, SEND=in_array) 
    r1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, r1d)
    !
    DEALLOCATE(tmp_array)
#else
    r1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, r1d)
#endif
  END SUBROUTINE gather_edges_r2d
  !
  SUBROUTINE gather_edges_r3d(in_array, out_array, p_patch)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r2d(:,:)
    !
#ifndef NOMPI
    REAL(wp), ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_edges_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_e, RECV=tmp_array, SEND=in_array) 
    r2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, r2d)
    !
    DEALLOCATE(tmp_array)
#else
    r2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, r2d)
#endif
  END SUBROUTINE gather_edges_r3d
  !
  SUBROUTINE gather_vertices_r2d(in_array, out_array, p_patch)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r1d(:)
    !
#ifndef NOMPI
    REAL(wp), ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_verts_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_v, RECV=tmp_array, SEND=in_array) 
    r1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, r1d)
    !
    DEALLOCATE(tmp_array)
#else
    r1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, r1d)
#endif
  END SUBROUTINE gather_vertices_r2d
  !
  SUBROUTINE gather_vertices_r3d(in_array, out_array, p_patch)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r2d(:,:)
    !
#ifndef NOMPI
    REAL(wp), ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_verts_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_v, RECV=tmp_array, SEND=in_array) 
    r2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, r2d)
    !
    DEALLOCATE(tmp_array)
#else
    r2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, r2d)
#endif
  END SUBROUTINE gather_vertices_r3d
  !
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE gather_cells_i2d(in_array, out_array, p_patch)
    INTEGER,                   INTENT(in) :: in_array(:,:)
    INTEGER, POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i1d(:)
    !
#ifndef NOMPI
    INTEGER, ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_cells_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_c, RECV=tmp_array, SEND=in_array) 
    i1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, i1d)
    !
    DEALLOCATE(tmp_array)
#else
    i1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, i1d)
#endif
  END SUBROUTINE gather_cells_i2d
  !
  SUBROUTINE gather_cells_i3d(in_array, out_array, p_patch)
    INTEGER,                   INTENT(in) :: in_array(:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i2d(:,:)
    !
#ifndef NOMPI
    INTEGER, ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_cells_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_c, RECV=tmp_array, SEND=in_array) 
    i2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, i2d)
    !
    DEALLOCATE(tmp_array)
#else
    i2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, i2d)
#endif
  END SUBROUTINE gather_cells_i3d
  !
  SUBROUTINE gather_edges_i2d(in_array, out_array, p_patch)
    INTEGER,                   INTENT(in) :: in_array(:,:)
    INTEGER, POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i1d(:)
    !
#ifndef NOMPI
    INTEGER, ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_edges_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_e, RECV=tmp_array, SEND=in_array) 
    i1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, i1d)
    !
    DEALLOCATE(tmp_array)
#else
    i1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, i1d)
#endif
  END SUBROUTINE gather_edges_i2d
  !
  SUBROUTINE gather_edges_i3d(in_array, out_array, p_patch)
    INTEGER,                   INTENT(in) :: in_array(:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i2d(:,:)
    !
#ifndef NOMPI
    INTEGER, ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_edges_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_e, RECV=tmp_array, SEND=in_array) 
    i2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, i2d)
    !
    DEALLOCATE(tmp_array)
#else
    i2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, i2d)
#endif
  END SUBROUTINE gather_edges_i3d
  !
  SUBROUTINE gather_vertices_i2d(in_array, out_array, p_patch)
    INTEGER,                   INTENT(in) :: in_array(:,:)
    INTEGER, POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i1d(:)
    !
#ifndef NOMPI
    INTEGER, ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_verts_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_v, RECV=tmp_array, SEND=in_array) 
    i1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, i1d)
    !
    DEALLOCATE(tmp_array)
#else
    i1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, i1d)
#endif
  END SUBROUTINE gather_vertices_i2d
  !
  SUBROUTINE gather_vertices_i3d(in_array, out_array, p_patch)
    INTEGER,                   INTENT(in) :: in_array(:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i2d(:,:)
    !
#ifndef NOMPI
    INTEGER, ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_verts_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_v, RECV=tmp_array, SEND=in_array) 
    i2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, i2d)
    !
    DEALLOCATE(tmp_array)
#else
    i2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, i2d)
#endif
  END SUBROUTINE gather_vertices_i3d
  !
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE gather_cells_l2d(in_array, out_array, p_patch)
    LOGICAL,                   INTENT(in) :: in_array(:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l1d(:)
    !
#ifndef NOMPI
    LOGICAL, ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_cells_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_c, RECV=tmp_array, SEND=in_array) 
    l1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, l1d)
    !
    DEALLOCATE(tmp_array)
#else
    l1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, l1d)
#endif
  END SUBROUTINE gather_cells_l2d
  !
  SUBROUTINE gather_cells_l3d(in_array, out_array, p_patch)
    LOGICAL,                   INTENT(in) :: in_array(:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l2d(:,:)
    !
#ifndef NOMPI
    LOGICAL, ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_cells_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_c, RECV=tmp_array, SEND=in_array) 
    l2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, l2d)
    !
    DEALLOCATE(tmp_array)
#else
    l2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, l2d)
#endif
  END SUBROUTINE gather_cells_l3d
  !
  SUBROUTINE gather_edges_l2d(in_array, out_array, p_patch)
    LOGICAL,                   INTENT(in) :: in_array(:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l1d(:)
    !
#ifndef NOMPI
    LOGICAL, ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_edges_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_e, RECV=tmp_array, SEND=in_array) 
    l1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, l1d)
    !
    DEALLOCATE(tmp_array)
#else
    l1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, l1d)
#endif
  END SUBROUTINE gather_edges_l2d
  !
  SUBROUTINE gather_edges_l3d(in_array, out_array, p_patch)
    LOGICAL,                   INTENT(in) :: in_array(:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l2d(:,:)
    !
#ifndef NOMPI
    LOGICAL, ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_edges_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_e, RECV=tmp_array, SEND=in_array) 
    l2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, l2d)
    !
    DEALLOCATE(tmp_array)
#else
    l2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, l2d)
#endif
  END SUBROUTINE gather_edges_l3d
  !
  SUBROUTINE gather_vertices_l2d(in_array, out_array, p_patch)
    LOGICAL,                   INTENT(in) :: in_array(:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:,:,:) 
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l1d(:)
    !
#ifndef NOMPI
    LOGICAL, ALLOCATABLE :: tmp_array(:,:) 
    INTEGER :: n1, n2
    !
    n1 = SIZE(in_array, 1)
    n2 = (p_patch%n_patch_verts_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2))
    !
    CALL exchange_data(p_patch%comm_pat_gather_v, RECV=tmp_array, SEND=in_array) 
    l1d => out_array(:,1,1,1,1)
    CALL reorder(tmp_array, l1d)
    !
    DEALLOCATE(tmp_array)
#else
    l1d => out_array(:,1,1,1,1)
    CALL reorder(in_array, l1d)
#endif
  END SUBROUTINE gather_vertices_l2d
  !
  SUBROUTINE gather_vertices_l3d(in_array, out_array, p_patch)
    LOGICAL,                   INTENT(in) :: in_array(:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l2d(:,:)
    !
#ifndef NOMPI
    LOGICAL, ALLOCATABLE :: tmp_array(:,:,:) 
    INTEGER :: n1, n2, n3
    !
    n1 = SIZE(in_array, 1)
    n2 = SIZE(in_array, 2)
    n3 = (p_patch%n_patch_verts_g-1)/n1+1
    ALLOCATE(tmp_array(n1, n2, n3))
    !
    CALL exchange_data(p_patch%comm_pat_gather_v, RECV=tmp_array, SEND=in_array) 
    l2d => out_array(:,:,1,1,1)
    CALL reorder(tmp_array, l2d)
    !
    DEALLOCATE(tmp_array)
#else
    l2d => out_array(:,:,1,1,1)
    CALL reorder(in_array, l2d)
#endif
  END SUBROUTINE gather_vertices_l3d
  !
  !------------------------------------------------------------------------------------------------
  !
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_cells_r2d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r1d(:)
    r1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(r1d, out_array)
#else
    CALL scatter_array_r2d(r1d, out_array, p_patch%cells%glb_index)
#endif
  END SUBROUTINE scatter_cells_r2d
  !
  SUBROUTINE scatter_cells_r3d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r2d(:,:)
    r2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(r2d, out_array)
#else
    CALL scatter_array_r3d(r2d, out_array, p_patch%cells%glb_index)
#endif
  END SUBROUTINE scatter_cells_r3d
  !
  SUBROUTINE scatter_edges_r2d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r1d(:)
    r1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(r1d, out_array)
#else
    CALL scatter_array_r2d(r1d, out_array, p_patch%edges%glb_index)
#endif
  END SUBROUTINE scatter_edges_r2d
  !
  SUBROUTINE scatter_edges_r3d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r2d(:,:)
    r2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(r2d, out_array)
#else
    CALL scatter_array_r3d(r2d, out_array, p_patch%edges%glb_index)
#endif
  END SUBROUTINE scatter_edges_r3d
  !
  SUBROUTINE scatter_vertices_r2d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r1d(:)
    r1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(r1d, out_array)
#else
    CALL scatter_array_r2d(r1d, out_array, p_patch%verts%glb_index)
#endif
  END SUBROUTINE scatter_vertices_r2d
  !
  SUBROUTINE scatter_vertices_r3d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    REAL(wp), POINTER :: r2d(:,:)
    r2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(r2d, out_array)
#else
    CALL scatter_array_r3d(r2d, out_array, p_patch%verts%glb_index)
#endif
  END SUBROUTINE scatter_vertices_r3d
  !
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_cells_i2d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i1d(:)
    i1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(i1d, out_array)
#else
    CALL scatter_array_i2d(i1d, out_array, p_patch%cells%glb_index)
#endif
  END SUBROUTINE scatter_cells_i2d
  !
  SUBROUTINE scatter_cells_i3d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i2d(:,:)
    i2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(i2d, out_array)
#else
    CALL scatter_array_i3d(i2d, out_array, p_patch%cells%glb_index)
#endif
  END SUBROUTINE scatter_cells_i3d
  !
  SUBROUTINE scatter_edges_i2d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i1d(:)
    i1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(i1d, out_array)
#else
    CALL scatter_array_i2d(i1d, out_array, p_patch%edges%glb_index)
#endif
  END SUBROUTINE scatter_edges_i2d
  !
  SUBROUTINE scatter_edges_i3d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i2d(:,:)
    i2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(i2d, out_array)
#else
    CALL scatter_array_i3d(i2d, out_array, p_patch%edges%glb_index)
#endif
  END SUBROUTINE scatter_edges_i3d
  !
  SUBROUTINE scatter_vertices_i2d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i1d(:)
    i1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(i1d, out_array)
#else
    CALL scatter_array_i2d(i1d, out_array, p_patch%verts%glb_index)
#endif
  END SUBROUTINE scatter_vertices_i2d
  !
  SUBROUTINE scatter_vertices_i3d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    INTEGER, POINTER :: i2d(:,:)
    i2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(i2d, out_array)
#else
    CALL scatter_array_i3d(i2d, out_array, p_patch%verts%glb_index)
#endif
  END SUBROUTINE scatter_vertices_i3d
  !
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_cells_l2d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l1d(:)
    l1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(l1d, out_array)
#else
    CALL scatter_array_l2d(l1d, out_array, p_patch%cells%glb_index)
#endif
  END SUBROUTINE scatter_cells_l2d
  !
  SUBROUTINE scatter_cells_l3d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l2d(:,:)
    l2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(l2d, out_array)
#else
    CALL scatter_array_l3d(l2d, out_array, p_patch%cells%glb_index)
#endif
  END SUBROUTINE scatter_cells_l3d
  !
  SUBROUTINE scatter_edges_l2d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l1d(:)
    l1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(l1d, out_array)
#else
    CALL scatter_array_l2d(l1d, out_array, p_patch%edges%glb_index)
#endif
  END SUBROUTINE scatter_edges_l2d
  !
  SUBROUTINE scatter_edges_l3d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l2d(:,:)
    l2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(l2d, out_array)
#else
    CALL scatter_array_l3d(l2d, out_array, p_patch%edges%glb_index)
#endif
  END SUBROUTINE scatter_edges_l3d
  !
  SUBROUTINE scatter_vertices_l2d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l1d(:)
    l1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(l1d, out_array)
#else
    CALL scatter_array_l2d(l1d, out_array, p_patch%verts%glb_index)
#endif
  END SUBROUTINE scatter_vertices_l2d
  !
  SUBROUTINE scatter_vertices_l3d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:)
#ifndef NOMPI
    TYPE(t_patch),    OPTIONAL, INTENT(in) :: p_patch
#else
    INTEGER,          OPTIONAL, INTENT(in) :: p_patch
#endif
    LOGICAL, POINTER :: l2d(:,:)
    l2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(l2d, out_array)
#else
    CALL scatter_array_l3d(l2d, out_array, p_patch%verts%glb_index)
#endif
  END SUBROUTINE scatter_vertices_l3d
  !
  !------------------------------------------------------------------------------------------------
  !
#ifndef NOMPI
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_array_r2d (in_array, out_array, global_index)
    REAL(wp), INTENT(inout) :: in_array(:)
    REAL(wp), INTENT(out)   :: out_array(:,:)
    INTEGER,  INTENT(in)    :: global_index(:)
    !
    INTEGER :: j, jl, jb
    !
    CALL p_bcast(in_array, p_io, p_comm_work)
    !
    out_array(:,:) = 0.0_wp
    !
    DO j = 1, SIZE(out_array)
      jb = blk_no(j)
      jl = idx_no(j)
      out_array(jl,jb) = in_array(global_index(j))
    ENDDO
    !
  END SUBROUTINE scatter_array_r2d
  !
  SUBROUTINE scatter_array_r3d (in_array, out_array, global_index)
    REAL(wp), INTENT(inout) :: in_array(:,:)
    REAL(wp), INTENT(out)   :: out_array(:,:,:)
    INTEGER,  INTENT(in)    :: global_index(:)
    !
    INTEGER :: j, jl, jb, jk
    !
    CALL p_bcast(in_array, p_io, p_comm_work)
    !
    out_array(:,:,:) = 0.0_wp
    !
    DO jk = 1, SIZE(out_array,2)
      DO j = 1, SIZE(out_array,1)*SIZE(out_array,3)
        jb = blk_no(j)
        jl = idx_no(j)
        out_array(jl,jk,jb) = in_array(global_index(j),jk)
      ENDDO
    ENDDO
    !
  END SUBROUTINE scatter_array_r3d
  !
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_array_i2d (in_array, out_array, global_index)
    INTEGER, INTENT(inout) :: in_array(:)
    INTEGER, INTENT(out)   :: out_array(:,:)
    INTEGER, INTENT(in)    :: global_index(:)
    !
    INTEGER :: j, jl, jb
    !
    CALL p_bcast(in_array, p_io, p_comm_work)
    !
    out_array(:,:) = 0
    !
    DO j = 1, SIZE(out_array)
      jb = blk_no(j)
      jl = idx_no(j)
      out_array(jl,jb) = in_array(global_index(j))
    ENDDO
    !
  END SUBROUTINE scatter_array_i2d
  !
  SUBROUTINE scatter_array_i3d (in_array, out_array, global_index)
    INTEGER, INTENT(inout) :: in_array(:,:)
    INTEGER, INTENT(out)   :: out_array(:,:,:)
    INTEGER, INTENT(in)    :: global_index(:)
    !
    INTEGER :: j, jl, jb, jk
    !
    CALL p_bcast(in_array, p_io, p_comm_work)
    !
    out_array(:,:,:) = 0
    !
    DO jk = 1, SIZE(out_array,2)
      DO j = 1, SIZE(out_array,1)*SIZE(out_array,3)
        jb = blk_no(j)
        jl = idx_no(j)
        out_array(jl,jk,jb) = in_array(global_index(j),jk)
      ENDDO
    ENDDO
    !
  END SUBROUTINE scatter_array_i3d
  !
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_array_l2d (in_array, out_array, global_index)
    LOGICAL, INTENT(inout) :: in_array(:)
    LOGICAL, INTENT(out)   :: out_array(:,:)
    INTEGER, INTENT(in)    :: global_index(:)
    !
    INTEGER :: j, jl, jb
    !
    CALL p_bcast(in_array, p_io, p_comm_work)
    !
    out_array(:,:) = .FALSE.
    !
    DO j = 1, SIZE(out_array)
      jb = blk_no(j)
      jl = idx_no(j)
      out_array(jl,jb) = in_array(global_index(j))
    ENDDO
    !
  END SUBROUTINE scatter_array_l2d
  !
  SUBROUTINE scatter_array_l3d (in_array, out_array, global_index)
    LOGICAL, INTENT(inout) :: in_array(:,:)
    LOGICAL, INTENT(out)   :: out_array(:,:,:)
    INTEGER, INTENT(in)    :: global_index(:)
    !
    INTEGER :: j, jl, jb, jk
    !
    CALL p_bcast(in_array, p_io, p_comm_work)
    !
    out_array(:,:,:) = .FALSE.
    !
    DO jk = 1, SIZE(out_array,2)
      DO j = 1, SIZE(out_array,1)*SIZE(out_array,3)
        jb = blk_no(j)
        jl = idx_no(j)
        out_array(jl,jk,jb) = in_array(global_index(j),jk)
      ENDDO
    ENDDO
    !
  END SUBROUTINE scatter_array_l3d
  !
#endif
  !------------------------------------------------------------------------------------------------
  !
  !================================================================================================ 
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE reorder_backward_r2d(in, out)
    REAL(wp), INTENT(in)    :: in(:,:)
    REAL(wp), INTENT(inout) :: out(:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::  isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in  = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_in-isize_out
    !
    IF(idiscrep == 0 )THEN
      out = RESHAPE(in,(/ isize_out /))
    ELSE
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
      out = PACK(RESHAPE(in,(/isize_in/)),lmask)
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_backward_r2d
  !
  SUBROUTINE reorder_backward_r3d(in, out)
    REAL(wp), INTENT(in)    :: in(:,:,:)
    REAL(wp), INTENT(inout) :: out(:,:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in  = SIZE(in,1)*SIZE(in,3)
    isize_out = SIZE(out,1)
    isize_lev = SIZE(in,2)
    idiscrep = isize_in-isize_out
    !
    IF (idiscrep /= 0 )THEN
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep /= 0 )THEN
        out(:,k) = PACK(RESHAPE(in(:,k,:),(/isize_in/)),lmask)
      ELSE
        out(:,k) =      RESHAPE(in(:,k,:),(/isize_out/))
      ENDIF
    ENDDO
    !   
    IF (idiscrep /= 0 )THEN
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_backward_r3d
  !
  SUBROUTINE reorder_foreward_r2d(in, out)
    REAL(wp), INTENT(in)    :: in(:)
    REAL(wp), INTENT(inout) :: out(:,:)
    !
    REAL(wp), ALLOCATABLE :: rpad(:)
    INTEGER :: isize_nproma, isize_nblks
    INTEGER :: isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_out-isize_in
    !
    isize_nproma = SIZE(out,1)
    isize_nblks = SIZE(out,2)
    !
    IF (idiscrep == 0) THEN
      out = RESHAPE(in,(/isize_nproma,isize_nblks/))
    ELSE
      ALLOCATE(rpad(idiscrep))
      rpad = 0.0_wp
      out = RESHAPE(in,(/isize_nproma,isize_nblks/),rpad)
      DEALLOCATE(rpad)
    ENDIF
    !
  END SUBROUTINE reorder_foreward_r2d
  !
  SUBROUTINE reorder_foreward_r3d(in, out)
    REAL(wp), INTENT(in)    :: in(:,:)
    REAL(wp), INTENT(inout) :: out(:,:,:)
    !
    !
    REAL(wp), ALLOCATABLE :: rpad(:)
    INTEGER :: isize_nproma, isize_nblks
    INTEGER :: isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in = SIZE(in,1)
    isize_out = SIZE(out,1)*SIZE(out,3)
    isize_lev = SIZE(in,2)
    idiscrep = isize_out-isize_in
    !
    isize_nproma = SIZE(out,1)
    isize_nblks = SIZE(out,3)
    !
    IF (idiscrep /= 0) THEN
      ALLOCATE(rpad(idiscrep))
      rpad = 0.0_wp
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep == 0) THEN
        out(:,k,:) = RESHAPE(in(:,k),(/isize_nproma,isize_nblks/))
      ELSE
        out(:,k,:) = RESHAPE(in(:,k),(/isize_nproma,isize_nblks/), rpad)
      ENDIF
    ENDDO
    !
    IF (idiscrep /= 0) THEN
      DEALLOCATE(rpad)
    ENDIF
    !
  END SUBROUTINE reorder_foreward_r3d
  !
  !================================================================================================ 
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE reorder_backward_i2d(in, out)
    INTEGER, INTENT(in)    :: in(:,:)
    INTEGER, INTENT(inout) :: out(:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::  isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in  = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_in-isize_out
    !
    IF(idiscrep == 0 )THEN
      out = RESHAPE(in,(/ isize_out /))
    ELSE
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
      out = PACK(RESHAPE(in,(/isize_in/)),lmask)
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_backward_i2d
  !
  SUBROUTINE reorder_backward_i3d(in, out)
    INTEGER, INTENT(in)    :: in(:,:,:)
    INTEGER, INTENT(inout) :: out(:,:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in  = SIZE(in,1)*SIZE(in,3)
    isize_out = SIZE(out,1)
    isize_lev = SIZE(in,2)
    idiscrep = isize_in-isize_out
    !
    IF (idiscrep /= 0 )THEN
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep /= 0 )THEN
        out(:,k) = PACK(RESHAPE(in(:,k,:),(/isize_in/)),lmask)
      ELSE
        out(:,k) =      RESHAPE(in(:,k,:),(/isize_out/))
      ENDIF
    ENDDO
    !   
    IF (idiscrep /= 0 )THEN
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_backward_i3d
  !
  SUBROUTINE reorder_foreward_i2d(in, out)
    INTEGER, INTENT(in)    :: in(:)
    INTEGER, INTENT(inout) :: out(:,:)
    !
    INTEGER, ALLOCATABLE :: ipad(:)
    INTEGER :: isize_nproma, isize_nblks
    INTEGER :: isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_out-isize_in
    !
    isize_nproma = SIZE(out,1)
    isize_nblks = SIZE(out,2)
    !
    IF (idiscrep == 0) THEN
      out = RESHAPE(in,(/isize_nproma,isize_nblks/))
    ELSE
      ALLOCATE(ipad(idiscrep))
      ipad = 0
      out = RESHAPE(in,(/isize_nproma,isize_nblks/), ipad)
      DEALLOCATE(ipad)
    ENDIF
    !
  END SUBROUTINE reorder_foreward_i2d
  !
  SUBROUTINE reorder_foreward_i3d(in, out)
    INTEGER, INTENT(in)    :: in(:,:)
    INTEGER, INTENT(inout) :: out(:,:,:)
    !
    !
    INTEGER, ALLOCATABLE :: ipad(:)
    INTEGER :: isize_nproma, isize_nblks
    INTEGER :: isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in = SIZE(in,1)
    isize_out = SIZE(out,1)*SIZE(out,3)
    isize_lev = SIZE(in,2)
    idiscrep = isize_out-isize_in
    !
    isize_nproma = SIZE(out,1)
    isize_nblks = SIZE(out,3)
    !
    IF (idiscrep /= 0) THEN
      ALLOCATE(ipad(idiscrep))
      ipad = 0
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep == 0) THEN
        out(:,k,:) = RESHAPE(in(:,k),(/isize_nproma,isize_nblks/))
      ELSE
        out(:,k,:) = RESHAPE(in(:,k),(/isize_nproma,isize_nblks/), ipad)
      ENDIF
    ENDDO
    !
    IF (idiscrep /= 0) THEN
      DEALLOCATE(ipad)
    ENDIF
    !
  END SUBROUTINE reorder_foreward_i3d
  !
  !================================================================================================ 
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE reorder_backward_l2d(in, out)
    LOGICAL, INTENT(in)    :: in(:,:)
    LOGICAL, INTENT(inout) :: out(:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::  isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in  = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_in-isize_out
    !
    IF(idiscrep == 0 )THEN
      out = RESHAPE(in,(/ isize_out /))
    ELSE
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
      out = PACK(RESHAPE(in,(/isize_in/)),lmask)
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_backward_l2d
  !
  SUBROUTINE reorder_backward_l3d(in, out)
    LOGICAL, INTENT(in)    :: in(:,:,:)
    LOGICAL, INTENT(inout) :: out(:,:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in  = SIZE(in,1)*SIZE(in,3)
    isize_out = SIZE(out,1)
    isize_lev = SIZE(in,2)
    idiscrep = isize_in-isize_out
    !
    IF (idiscrep /= 0 )THEN
      ALLOCATE (lmask(isize_in))
      lmask(1:isize_out) = .TRUE.
      lmask(isize_out+1:isize_in) = .FALSE.
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep /= 0 )THEN
        out(:,k) = PACK(RESHAPE(in(:,k,:),(/isize_in/)),lmask)
      ELSE
        out(:,k) =      RESHAPE(in(:,k,:),(/isize_out/))
      ENDIF
    ENDDO
    !   
    IF (idiscrep /= 0 )THEN
      DEALLOCATE (lmask)
    ENDIF
    !
  END SUBROUTINE reorder_backward_l3d
  !
  SUBROUTINE reorder_foreward_l2d(in, out)
    LOGICAL, INTENT(in)    :: in(:)
    LOGICAL, INTENT(inout) :: out(:,:)
    !
    LOGICAL, ALLOCATABLE :: lpad(:)
    INTEGER :: isize_nproma, isize_nblks
    INTEGER :: isize_in, isize_out
    INTEGER :: idiscrep
    !
    isize_in = SIZE(in)
    isize_out = SIZE(out)
    idiscrep = isize_out-isize_in
    !
    isize_nproma = SIZE(out,1)
    isize_nblks = SIZE(out,2)
    !
    IF (idiscrep == 0) THEN
      out = RESHAPE(in,(/isize_nproma,isize_nblks/))
    ELSE
      ALLOCATE(lpad(idiscrep))
      lpad = .FALSE.
      out = RESHAPE(in,(/isize_nproma,isize_nblks/), lpad)
      DEALLOCATE(lpad)
    ENDIF
    !
  END SUBROUTINE reorder_foreward_l2d
  !
  SUBROUTINE reorder_foreward_l3d(in, out)
    LOGICAL, INTENT(in)    :: in(:,:)
    LOGICAL, INTENT(inout) :: out(:,:,:)
    !
    !
    LOGICAL, ALLOCATABLE :: lpad(:)
    INTEGER :: isize_nproma, isize_nblks
    INTEGER :: isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in = SIZE(in,1)
    isize_out = SIZE(out,1)*SIZE(out,3)
    isize_lev = SIZE(in,2)
    idiscrep = isize_out-isize_in
    !
    isize_nproma = SIZE(out,1)
    isize_nblks = SIZE(out,3)
    !
    IF (idiscrep /= 0) THEN
      ALLOCATE(lpad(idiscrep))
      lpad = .FALSE.
    ENDIF
    !
    DO k = 1, isize_lev
      IF (idiscrep == 0) THEN
        out(:,k,:) = RESHAPE(in(:,k),(/isize_nproma,isize_nblks/))
      ELSE
        out(:,k,:) = RESHAPE(in(:,k),(/isize_nproma,isize_nblks/), lpad)
      ENDIF
    ENDDO
    !
    IF (idiscrep /= 0) THEN
      DEALLOCATE(lpad)
    ENDIF
    !
  END SUBROUTINE reorder_foreward_l3d
  !
  !------------------------------------------------------------------------------------------------
  !
END MODULE mo_io_restart_distribute
