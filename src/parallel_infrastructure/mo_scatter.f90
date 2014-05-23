#ifdef __xlC__
@PROCESS smp=noopt
@PROCESS noopt
#endif
!>
!! This module provides basic MPI scatter
!!
!!
!! @author Luis Kornblueh, MPI
!! @author Leonidas Linardakis, MPI
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_scatter
  !
  USE mo_kind,          ONLY: wp
  USE mo_communication, ONLY: idx_no, blk_no
  USE mo_model_domain,  ONLY: t_patch
  USE mo_mpi,           ONLY: process_mpi_root_id, p_bcast, p_comm_work, &
                              my_process_is_mpi_seq, my_process_is_mpi_workroot
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: broadcast_array
  !
  PUBLIC :: scatter_cells
  PUBLIC :: scatter_cells_2D
  PUBLIC :: scatter_cells_2D_time
  PUBLIC :: scatter_cells_3D_time
  PUBLIC :: scatter_edges
  PUBLIC :: scatter_vertices
  !
  !------------------------------------------------------------------------------------------------
  !

  INTERFACE broadcast_array
    MODULE PROCEDURE broadcast_real_array_1D
    MODULE PROCEDURE broadcast_real_array_2D
    MODULE PROCEDURE broadcast_real_array_3D
    MODULE PROCEDURE broadcast_real_array_4D
    MODULE PROCEDURE broadcast_int_array_1D
  END INTERFACE broadcast_array

  INTERFACE scatter_cells_2D
    MODULE PROCEDURE scatter_real_cells_2D_noblocks_2blocks
  END INTERFACE scatter_cells_2D

  INTERFACE scatter_cells_2D_time
    MODULE PROCEDURE scatter_real_cells_2D_time_noblocks_2blocks
  END INTERFACE scatter_cells_2D_time

  INTERFACE scatter_cells_3D_time
    MODULE PROCEDURE scatter_real_cells_3D_time_noblocks_2blocks
  END INTERFACE scatter_cells_3D_time

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
  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE broadcast_int_array_1D(in_array)
    INTEGER, TARGET  :: in_array(:)

    IF (my_process_is_mpi_seq()) RETURN
#ifndef NOMPI
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
#endif
  END SUBROUTINE broadcast_int_array_1D
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE scatter_real_cells_2D_time_noblocks_2blocks(in_array, out_array, p_patch)
    REAL(wp), POINTER                   :: in_array(:,:)
    REAL(wp), POINTER                   :: out_array(:,:,:)
    TYPE(t_patch), INTENT(in)           :: p_patch

!    WRITE(0,*) "LBOUND(in_array,2)= ", LBOUND(in_array,2)
!    WRITE(0,*) "UBOUND(in_array,2)= ", UBOUND(in_array,2)
!    WRITE(0,*) "LBOUND(out_array,3)= ", LBOUND(out_array,3)
!    WRITE(0,*) "UBOUND(out_array,3)= ", UBOUND(out_array,3)
    CALL scatter_array_r2d_time(in_array, out_array, p_patch%cells%decomp_info%glb_index)

  END SUBROUTINE scatter_real_cells_2D_time_noblocks_2blocks
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE scatter_real_cells_3D_time_noblocks_2blocks(in_array, out_array, p_patch)
    REAL(wp), POINTER                   :: in_array(:,:,:)
    REAL(wp), POINTER                   :: out_array(:,:,:,:)
    TYPE(t_patch), INTENT(in)           :: p_patch

!    WRITE(0,*) "LBOUND(in_array,3)= ", LBOUND(in_array,3)
!    WRITE(0,*) "UBOUND(in_array,3)= ", UBOUND(in_array,3)
!    WRITE(0,*) "LBOUND(out_array,4)= ", LBOUND(out_array,4)
!    WRITE(0,*) "UBOUND(out_array,4)= ", UBOUND(out_array,4)
    CALL scatter_array_r4d(in_array, out_array, p_patch%cells%decomp_info%glb_index)

  END SUBROUTINE scatter_real_cells_3D_time_noblocks_2blocks
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE scatter_real_cells_2D_noblocks_2blocks(in_array, out_array, p_patch)
    REAL(wp), POINTER                   :: in_array(:)
    REAL(wp), POINTER                   :: out_array(:,:)
    TYPE(t_patch), INTENT(in)           :: p_patch

    IF (my_process_is_mpi_seq()) THEN
      CALL reorder(in_array, out_array)
#ifndef NOMPI
    ELSE
      CALL scatter_array_r2d(in_array, out_array, p_patch%cells%decomp_info%glb_index)
#endif
    ENDIF

  END SUBROUTINE scatter_real_cells_2D_noblocks_2blocks
  !
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE scatter_array_r3d (in_array, out_array, global_index)
    REAL(wp), POINTER :: in_array(:,:)
    REAL(wp), POINTER :: out_array(:,:,:)
    INTEGER           :: global_index(:)

    INTEGER :: j, jl, jb, jk, jk1

    out_array(:,:,:) = 0.0_wp

    DO jk = 1, SIZE(out_array,2)
      IF (my_process_is_mpi_workroot()) THEN
        CALL p_bcast(in_array(:,jk), process_mpi_root_id, p_comm_work)
        jk1 = jk
      ELSE
        CALL p_bcast(in_array(:,1), process_mpi_root_id, p_comm_work)
        jk1 = 1
      ENDIF
      DO j = 1, SIZE(global_index)
        jb = blk_no(j)
        jl = idx_no(j)
        out_array(jl,jk,jb) = in_array(global_index(j),jk1)
      ENDDO
    ENDDO
    !
  END SUBROUTINE scatter_array_r3d
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  ! standard input array shape: (cells, vertical_levels, file_time_steps)
  ! standard output array shape: (nproma, vertical_levels, blocks, file_time_steps)
  SUBROUTINE scatter_array_r2d_time (in_array, out_array, global_index)
    REAL(wp), POINTER :: in_array(:,:)
    REAL(wp), POINTER :: out_array(:,:,:)
    INTEGER           :: global_index(:)

    INTEGER :: j, jl, jb, time_step

    ! this has internal check if sequential
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)

    out_array(:,:,:) = 0.0_wp

    DO time_step=LBOUND(in_array,2), UBOUND(in_array,2)
      DO j = 1,   SIZE(global_index)
        jb = blk_no(j)
        jl = idx_no(j)
        out_array(jl, jb, time_step) = in_array(global_index(j), time_step)
      ENDDO
    ENDDO

  END SUBROUTINE scatter_array_r2d_time
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  ! standard input array shape: (cells, vertical_levels, file_time_steps)
  ! standard output array shape: (nproma, vertical_levels, blocks, file_time_steps)
  SUBROUTINE scatter_array_r4d (in_array, out_array, global_index)
    REAL(wp), POINTER :: in_array(:,:,:)
    REAL(wp), POINTER :: out_array(:,:,:,:)
    INTEGER            :: global_index(:)

    INTEGER :: j, jl, jb, jk, time_step

    ! this has internal check if sequential
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)

    out_array(:,:,:,:) = 0.0_wp

    DO time_step=LBOUND(in_array,3), UBOUND(in_array,3)
      DO jk = 1,    SIZE(out_array,2)
        DO j = 1,   SIZE(global_index)
          jb = blk_no(j)
          jl = idx_no(j)
          out_array(jl, jk, jb, time_step) = in_array(global_index(j), jk, time_step)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE scatter_array_r4d

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE broadcast_real_array_1D(in_array)
    REAL(wp), TARGET  :: in_array(:)

    IF (my_process_is_mpi_seq()) RETURN
#ifndef NOMPI
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
#endif
  END SUBROUTINE broadcast_real_array_1D
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE broadcast_real_array_2D(in_array)
    REAL(wp), TARGET  :: in_array(:,:)

    IF (my_process_is_mpi_seq()) RETURN
#ifndef NOMPI
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
#endif
  END SUBROUTINE broadcast_real_array_2D
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE broadcast_real_array_3D(in_array)
    REAL(wp), TARGET  :: in_array(:,:,:)

    IF (my_process_is_mpi_seq()) RETURN
#ifndef NOMPI
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
#endif
  END SUBROUTINE broadcast_real_array_3D
  !--------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE broadcast_real_array_4D(in_array)
    REAL(wp), TARGET  :: in_array(:,:,:,:)

    IF (my_process_is_mpi_seq()) RETURN
#ifndef NOMPI
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
#endif
  END SUBROUTINE broadcast_real_array_4D
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  !>
  SUBROUTINE scatter_cells_r2d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
    TYPE(t_patch),              INTENT(in) :: p_patch
    REAL(wp), POINTER :: r1d(:)
    r1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(r1d, out_array)
#else
    CALL scatter_array_r2d(r1d, out_array, p_patch%cells%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_cells_r2d
  !
  SUBROUTINE scatter_cells_r3d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),              INTENT(in) :: p_patch
    REAL(wp), POINTER :: r2d(:,:)
    r2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(r2d, out_array)
#else
    CALL scatter_array_r3d(r2d, out_array, p_patch%cells%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_cells_r3d
  !
  SUBROUTINE scatter_edges_r2d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
    TYPE(t_patch),              INTENT(in) :: p_patch
    REAL(wp), POINTER :: r1d(:)
    r1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(r1d, out_array)
#else
    CALL scatter_array_r2d(r1d, out_array, p_patch%edges%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_edges_r2d
  !
  SUBROUTINE scatter_edges_r3d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),              INTENT(in) :: p_patch
    REAL(wp), POINTER :: r2d(:,:)
    r2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(r2d, out_array)
#else
    CALL scatter_array_r3d(r2d, out_array, p_patch%edges%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_edges_r3d
  !
  SUBROUTINE scatter_vertices_r2d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
    TYPE(t_patch),              INTENT(in) :: p_patch
    REAL(wp), POINTER :: r1d(:)
    r1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(r1d, out_array)
#else
    CALL scatter_array_r2d(r1d, out_array, p_patch%verts%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_vertices_r2d
  !
  SUBROUTINE scatter_vertices_r3d(in_array, out_array, p_patch)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),              INTENT(in) :: p_patch
    REAL(wp), POINTER :: r2d(:,:)
    r2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(r2d, out_array)
#else
    CALL scatter_array_r3d(r2d, out_array, p_patch%verts%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_vertices_r3d
  !
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_cells_i2d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    INTEGER, POINTER :: i1d(:)
    i1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(i1d, out_array)
#else
    CALL scatter_array_i2d(i1d, out_array, p_patch%cells%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_cells_i2d
  !
  SUBROUTINE scatter_cells_i3d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    INTEGER, POINTER :: i2d(:,:)
    i2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(i2d, out_array)
#else
    CALL scatter_array_i3d(i2d, out_array, p_patch%cells%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_cells_i3d
  !
  SUBROUTINE scatter_edges_i2d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    INTEGER, POINTER :: i1d(:)
    i1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(i1d, out_array)
#else
    CALL scatter_array_i2d(i1d, out_array, p_patch%edges%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_edges_i2d
  !
  SUBROUTINE scatter_edges_i3d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    INTEGER, POINTER :: i2d(:,:)
    i2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(i2d, out_array)
#else
    CALL scatter_array_i3d(i2d, out_array, p_patch%edges%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_edges_i3d
  !
  SUBROUTINE scatter_vertices_i2d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    INTEGER, POINTER :: i1d(:)
    i1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(i1d, out_array)
#else
    CALL scatter_array_i2d(i1d, out_array, p_patch%verts%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_vertices_i2d
  !
  SUBROUTINE scatter_vertices_i3d(in_array, out_array, p_patch)
    INTEGER, POINTER                      :: in_array(:,:,:,:,:)
    INTEGER, POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    INTEGER, POINTER :: i2d(:,:)
    i2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(i2d, out_array)
#else
    CALL scatter_array_i3d(i2d, out_array, p_patch%verts%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_vertices_i3d
  !
  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_cells_l2d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    LOGICAL, POINTER :: l1d(:)
    l1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(l1d, out_array)
#else
    CALL scatter_array_l2d(l1d, out_array, p_patch%cells%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_cells_l2d
  !
  SUBROUTINE scatter_cells_l3d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    LOGICAL, POINTER :: l2d(:,:)
    l2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(l2d, out_array)
#else
    CALL scatter_array_l3d(l2d, out_array, p_patch%cells%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_cells_l3d
  !
  SUBROUTINE scatter_edges_l2d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    LOGICAL, POINTER :: l1d(:)
    l1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(l1d, out_array)
#else
    CALL scatter_array_l2d(l1d, out_array, p_patch%edges%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_edges_l2d
  !
  SUBROUTINE scatter_edges_l3d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    LOGICAL, POINTER :: l2d(:,:)
    l2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(l2d, out_array)
#else
    CALL scatter_array_l3d(l2d, out_array, p_patch%edges%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_edges_l3d
  !
  SUBROUTINE scatter_vertices_l2d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    LOGICAL, POINTER :: l1d(:)
    l1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(l1d, out_array)
#else
    CALL scatter_array_l2d(l1d, out_array, p_patch%verts%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_vertices_l2d
  !
  SUBROUTINE scatter_vertices_l3d(in_array, out_array, p_patch)
    LOGICAL, POINTER                      :: in_array(:,:,:,:,:)
    LOGICAL, POINTER                      :: out_array(:,:,:)
    TYPE(t_patch),             INTENT(in) :: p_patch
    LOGICAL, POINTER :: l2d(:,:)
    l2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(l2d, out_array)
#else
    CALL scatter_array_l3d(l2d, out_array, p_patch%verts%decomp_info%glb_index)
#endif
  END SUBROUTINE scatter_vertices_l3d
  !--------------------------------------------------------------------------------------
  !-----------------------------------------------------------------
  !>
  SUBROUTINE reorder_backward_r3d(in, out)
    REAL(wp), TARGET  :: in(:,:,:)
    REAL(wp), TARGET  :: out(:,:)
    !
    LOGICAL, ALLOCATABLE    :: lmask(:)
    INTEGER ::isize_in, isize_out, isize_lev
    INTEGER :: idiscrep, k
    !
    isize_in  = SIZE(in,1)*SIZE(in,3)
    isize_out = SIZE(out,1)
    isize_lev = SIZE(in,2)
    idiscrep = isize_in - isize_out
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
  !-----------------------------------------------------------------


  !------------------------------------------------------------------------------------------------
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
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
    !
    out_array(:,:) = 0.0_wp
    !
    DO j = 1, SIZE(global_index)
      jb = blk_no(j)
      jl = idx_no(j)
      out_array(jl,jb) = in_array(global_index(j))
    ENDDO
    !
  END SUBROUTINE scatter_array_r2d
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
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
    !
    out_array(:,:) = 0
    !
    DO j = 1, SIZE(global_index)
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
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
    !
    out_array(:,:,:) = 0
    !
    DO jk = 1, SIZE(out_array,2)
      DO j = 1, SIZE(global_index)
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
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
    !
    out_array(:,:) = .FALSE.
    !
    DO j = 1, SIZE(global_index)
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
    CALL p_bcast(in_array, process_mpi_root_id, p_comm_work)
    !
    out_array(:,:,:) = .FALSE.
    !
    DO jk = 1, SIZE(out_array,2)
      DO j = 1, SIZE(global_index)
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
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !>
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


  !------------------------------------------------------------------------------------------------
END MODULE mo_scatter
