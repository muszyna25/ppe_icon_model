#define NOMPI
MODULE mo_io_restart_distribute
  !
  USE mo_kind, ONLY: wp
#ifndef NOMPI
  USE mo_mpi,  ONLY: p_parallel_io
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
#ifdef NOMPI
  LOGICAL :: p_parallel_io = .TRUE.
#endif
  !------------------------------------------------------------------------------------------------
  !
  INTERFACE gather_cells
    MODULE PROCEDURE gather_cells_2d
    MODULE PROCEDURE gather_cells_3d
  END INTERFACE gather_cells
  !
  INTERFACE gather_edges
    MODULE PROCEDURE gather_edges_2d
    MODULE PROCEDURE gather_edges_3d
  END INTERFACE gather_edges
  !
  INTERFACE gather_vertices
    MODULE PROCEDURE gather_vertices_2d
    MODULE PROCEDURE gather_vertices_3d
  END INTERFACE gather_vertices
  !
  INTERFACE scatter_cells
    MODULE PROCEDURE scatter_cells_2d
    MODULE PROCEDURE scatter_cells_3d
  END INTERFACE scatter_cells
  !
  INTERFACE scatter_edges
    MODULE PROCEDURE scatter_edges_2d
    MODULE PROCEDURE scatter_edges_3d
  END INTERFACE scatter_edges
  !
  INTERFACE scatter_vertices
    MODULE PROCEDURE scatter_vertices_2d
    MODULE PROCEDURE scatter_vertices_3d
  END INTERFACE scatter_vertices
  !------------------------------------------------------------------------------------------------
  !
  INTERFACE reorder
    MODULE PROCEDURE reorder_foreward_2d
    MODULE PROCEDURE reorder_foreward_3d
    MODULE PROCEDURE reorder_backward_2d
    MODULE PROCEDURE reorder_backward_3d
  END INTERFACE reorder
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE gather_cells_2d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => out_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(in_array, z1d)
#endif
  END SUBROUTINE gather_cells_2d
  !
  SUBROUTINE gather_cells_3d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => out_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(in_array, z2d)
#endif
  END SUBROUTINE gather_cells_3d
  !
  SUBROUTINE gather_vertices_2d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => out_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(in_array, z1d)
#endif
  END SUBROUTINE gather_vertices_2d
  !
  SUBROUTINE gather_vertices_3d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => out_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(in_array, z2d)
#endif
  END SUBROUTINE gather_vertices_3d
  !
  SUBROUTINE gather_edges_2d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => out_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(in_array, z1d)
#endif
  END SUBROUTINE gather_edges_2d
  !
  SUBROUTINE gather_edges_3d(in_array, out_array, name)
    REAL(wp),                   INTENT(in) :: in_array(:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => out_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(in_array, z2d)
#endif
  END SUBROUTINE gather_edges_3d
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE scatter_cells_2d(in_array, out_array, name)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(z1d, out_array)
#endif
  END SUBROUTINE scatter_cells_2d
  !
  SUBROUTINE scatter_cells_3d(in_array, out_array, name)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(z2d, out_array)
#endif
  END SUBROUTINE scatter_cells_3d
  !
  SUBROUTINE scatter_vertices_2d(in_array, out_array, name)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(z1d, out_array)
#endif
  END SUBROUTINE scatter_vertices_2d
  !
  SUBROUTINE scatter_vertices_3d(in_array, out_array, name)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(z2d, out_array)
#endif
  END SUBROUTINE scatter_vertices_3d
  !
  SUBROUTINE scatter_edges_2d(in_array, out_array, name)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z1d(:)
    z1d => in_array(:,1,1,1,1)
#ifdef NOMPI
    CALL reorder(z1d, out_array)
#endif
  END SUBROUTINE scatter_edges_2d
  !
  SUBROUTINE scatter_edges_3d(in_array, out_array, name)
    REAL(wp), POINTER                      :: in_array(:,:,:,:,:)
    REAL(wp), POINTER                      :: out_array(:,:,:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
    REAL(wp), POINTER :: z2d(:,:)
    z2d => in_array(:,:,1,1,1)
#ifdef NOMPI
    CALL reorder(z2d, out_array)
#endif
  END SUBROUTINE scatter_edges_3d
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE reorder_backward_2d(in, out)
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
  END SUBROUTINE reorder_backward_2d
  !
  SUBROUTINE reorder_backward_3d(in, out)
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
  END SUBROUTINE reorder_backward_3d
  !
  SUBROUTINE reorder_foreward_2d(in, out)
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
  END SUBROUTINE reorder_foreward_2d
  !
  SUBROUTINE reorder_foreward_3d(in, out)
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
  END SUBROUTINE reorder_foreward_3d
!!$  !------------------------------------------------------------------------------------------------
!!$  !
!!$  SUBROUTINE gather_array_2d(typ, in_field, out_field, name)
!!$    !
!!$    INTEGER,                    INTENT(in)  :: typ
!!$    REAL(wp), POINTER                       :: in_field(:,:)
!!$    REAL(wp), POINTER                       :: out_field(:,:,:,:,:)
!!$    CHARACTER(len=*), OPTIONAL, INTENT(in)  :: name
!!$    !
!!$    REAL(wp), ALLOCATABLE :: out_field2(:,:)
!!$    !
!!$    IF (p_parallel_io) THEN
!!$      ALLOCATE(out_field2(UBOUND(out_field,1),1))
!!$    ELSE
!!$      ALLOCATE(out_field2(0,0))
!!$    ENDIF
!!$
!!$    CALL gather_array_3d(typ,p_patch,&
!!$         &               RESHAPE(in_field,(/SIZE(in_field,1),1,SIZE(in_field,2)/)), &
!!$         &               out_field2,name=name)
!!$
!!$    IF(p_parallel_io) THEN
!!$      out_field(:) = out_field2(:,1)
!!$    ENDIF
!!$
!!$    DEALLOCATE(out_field2)
!!$
!!$  END SUBROUTINE gather_array_2d
!!$  !
!!$  SUBROUTINE gather_array_3d(typ, in_field, out_field, name)
!!$    !
!!$    INTEGER,                    INTENT(in) :: typ
!!$    REAL(wp),    INTENT(in)         :: in_field(:,:,:)
!!$    REAL(wp),    INTENT(inout)      :: out_field(:,:)
!!$    CHARACTER(len=*), OPTIONAL, INTENT(in) :: name
!!$    !
!!$    REAL(wp), ALLOCATABLE :: tmp_field(:,:,:)
!!$    !
!!$    INTEGER :: isize_out, isize_lev 
!!$    !
!!$    INTEGER :: nblks, nproma, jb, jl, jk, jend
!!$    !
!!$#ifndef NOMPI
!!$    TYPE(t_comm_pattern), POINTER :: p_comm_pat
!!$#endif
!!$    !
!!$    !-----------------------------------------------------------------------
!!$
!!$    IF (SIZE(in_field,1) /= nproma) THEN
!!$      CALL finish('gather_array_3d','Illegal 1st array dimension')
!!$    ENDIF
!!$    !
!!$#ifndef NOMPI
!!$    IF (p_io /= p_test_pe .AND. p_io /= p_work_pe0) THEN ! Safety check only
!!$      CALL finish('gather_array_3d','Illegal I/O PE number for this routine')
!!$    ENDIF
!!$#endif
!!$
!!$    SELECT CASE (type)
!!$    CASE (GRID_UNSTRUCTURED_CELL)
!!$
!!$      IF (SIZE(in_field,3) /= p_patch%nblks_c) &
!!$           CALL finish('gather_array_3d','Illegal 3rd array dimension')
!!$      
!!$      ALLOCATE(tmp_field(nproma,SIZE(in_field,2),(p_patch%n_patch_cells_g-1)/nproma+1))
!!$
!!$      p_comm_pat => p_patch%comm_pat_gather_c
!!$      nblks      =  p_patch%nblks_c
!!$      nproma     =  p_patch%nproma_c
!!$
!!$    CASE (GRID_UNSTRUCTURED_EDGE)
!!$
!!$      IF (SIZE(in_field,3) /= p_patch%nblks_e) &
!!$           CALL finish('gather_array_3d','Illegal 3rd array dimension')
!!$      
!!$      ALLOCATE(tmp_field(nproma,SIZE(in_field,2),(p_patch%n_patch_edges_g-1)/nproma+1))
!!$      
!!$      p_comm_pat => p_patch%comm_pat_gather_e
!!$      nblks      =  p_patch%nblks_e
!!$      nproma     =  p_patch%nproma_e
!!$      
!!$    CASE (GRID_UNSTRUCTURED_VERT)
!!$      
!!$      IF (SIZE(in_field,3) /= p_patch%nblks_v) &
!!$           CALL finish('gather_array_3d','Illegal 3rd array dimension')
!!$
!!$      ALLOCATE(tmp_field(nproma,SIZE(in_field,2),(p_patch%n_patch_verts_g-1)/nproma+1))
!!$
!!$      p_comm_pat => p_patch%comm_pat_gather_v
!!$      nblks      =  p_patch%nblks_v
!!$      nproma     =  p_patch%nproma_v
!!$      
!!$
!!$    CASE DEFAULT
!!$      
!!$      CALL finish('gather_array_3d','Illegal typ parameter')
!!$
!!$    END SELECT
!!$
!!$    tmp_field(:,:,:) = 0.0_wp
!!$
!!$    IF (p_test_run) THEN
!!$      IF (p_pe /= p_test_pe) THEN
!!$        ! Gather all data on p_work_pe0 and send it to p_test_pe for verification
!!$        CALL exchange_data(p_comm_pat, RECV=tmp_field, SEND=in_field)
!!$        IF (p_pe_work == 0) CALL p_send(tmp_field, p_test_pe, 1)
!!$      ELSE
!!$        ! Receive result from parallel worker PEs and check for correctness
!!$        CALL p_recv(tmp_field, p_work_pe0, 1)
!!$        DO jb = 1, nblks
!!$          jend = nproma
!!$          IF (jb == nblks) jend = nproma
!!$          DO jl = 1, jend
!!$            IF (ANY(tmp_field(jl,:,jb) /= in_field(jl,:,jb))) THEN
!!$              IF (PRESENT(name)) THEN
!!$                WRITE(message_text,'(a,a,i5,i5)') 'Error ', name, jl, jb !,tmp_field(jl,:,jb), in_field(jl,:,jb)
!!$              ELSE
!!$                WRITE(message_text,'(a,i5,i5)') 'Error ', jl, jb !,tmp_field(jl,:,jb), in_field(jl,:,jb)
!!$              ENDIF
!!$              CALL message('gather_array_3d','Sync error test PE/worker PEs')
!!$            ENDIF
!!$          ENDDO
!!$        ENDDO
!!$      ENDIF
!!$    ELSE
!!$      IF (num_work_procs == 1) THEN
!!$        ! We are running on 1 PE, just copy in_field
!!$        DO jb = 1, nblks
!!$          jend = nproma
!!$          IF (jb == nblks) jend = nproma
!!$          DO jl = 1, jend
!!$            tmp_field(jl,:,jb) = in_field(jl,:,jb)
!!$          ENDDO
!!$        ENDDO
!!$      ELSE
!!$        ! Gather all data on p_work_pe0
!!$        CALL exchange_data(p_comm_pat, RECV=tmp_field, SEND=in_field)
!!$      ENDIF
!!$    ENDIF
!!$
!!$    IF (p_pe == p_io) THEN
!!$      isize_out = SIZE(out_field,1)
!!$      isize_lev = SIZE(in_field,2)
!!$
!!$      DO jk = 1, isize_lev
!!$        out_field(:,jk) = RESHAPE(tmp_field(:,jk,:),(/isize_out/))
!!$      ENDDO
!!$    ENDIF
!!$    
!!$    DEALLOCATE(tmp_field)
!!$    
!!$  END SUBROUTINE gather_array_3d
!!$
!!$!  SUBROUTINE scatter_array_2d (ncid, varname, glb_arr_len, loc_arr_len, glb_index, var_out)
!!$
!!$    CHARACTER(len=*), INTENT(IN)  ::   varname
!!$
!!$    INTEGER, INTENT(IN) :: ncid       
!!$    INTEGER, INTENT(IN) :: glb_arr_len
!!$    INTEGER, INTENT(IN) :: loc_arr_len
!!$    INTEGER, INTENT(IN) :: glb_index(:)
!!$
!!$    REAL(wp), INTENT(INOUT) ::  var_out(:,:)
!!$
!!$    INTEGER :: varid, mpi_comm, j, jl, jb
!!$    REAL(wp):: z_dummy_array(glb_arr_len)
!!$
!!$
!!$    ! Get var ID
!!$    IF(p_pe==p_io) CALL nf(nf_inq_varid(ncid, TRIM(varname), varid))
!!$
!!$    IF(p_test_run) THEN
!!$      mpi_comm = p_comm_work_test
!!$    ELSE
!!$      mpi_comm = p_comm_work
!!$    ENDIF
!!$
!!$    ! I/O PE reads and broadcasts data
!!$
!!$    IF(p_pe==p_io) CALL nf(nf_get_var_double(ncid, varid, z_dummy_array(:)))
!!$    CALL p_bcast(z_dummy_array, p_io, mpi_comm)
!!$
!!$    var_out(:,:) = 0._wp
!!$
!!$    ! Set var_out from global data
!!$    DO j = 1, loc_arr_len
!!$
!!$      jb = blk_no(j) ! Block index in distributed patch
!!$      jl = idx_no(j) ! Line  index in distributed patch
!!$
!!$      var_out(jl,jb) = z_dummy_array(glb_index(j))
!!$    ENDDO
!!$
!!$  END SUBROUTINE scatter_array_2d

END MODULE mo_io_restart_distribute

