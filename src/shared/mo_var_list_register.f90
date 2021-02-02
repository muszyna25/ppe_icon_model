! @par Copyright and Li/cense
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_register
  USE mo_var_metadata_types, ONLY: var_metadata_fromBinary, var_metadata_toBinary, var_metadata_get_size
  USE mo_var,              ONLY: t_var
  USE mo_var_list,         ONLY: t_var_list_ptr
  USE mo_exception,        ONLY: message, finish
  USE mo_io_config,        ONLY: restart_file_type
  USE mo_packed_message,   ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_key_value_store,  ONLY: t_key_value_store
  USE mo_impl_constants,   ONLY: vlname_len

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vl_register_iter
  PUBLIC :: vlr_add, vlr_get, vlr_del, vlr_packer

  TYPE t_vl_register_iter
    INTEGER, PRIVATE :: cur_idx = -1
    TYPE(t_var_list_ptr), PUBLIC :: cur
    LOGICAL, PRIVATE :: anew = .TRUE.
  CONTAINS
    PROCEDURE, PUBLIC :: next => iter_next
  END TYPE t_vl_register_iter

  INTEGER, SAVE :: last_id = 0
  TYPE(t_key_value_store) :: map
  TYPE(t_var_list_ptr), ALLOCATABLE :: storage(:)
  CHARACTER(*), PARAMETER :: modname = 'mo_var_list_store'

CONTAINS
  !------------------------------------------------------------------------------------------------
  ! Create a new memory buffer / output var_list
  ! Get a pointer to the new var_list
  SUBROUTINE vlr_add(list, vlname, output_type, restart_type,    &
    & post_suf, rest_suf, init_suf, loutput, lrestart, linitial, patch_id,          &
    & vlevel_type, model_type, filename, compression_type)
    TYPE(t_var_list_ptr), INTENT(OUT) :: list    ! anchor
    CHARACTER(*), INTENT(IN) :: vlname         ! name of output var_list
    INTEGER, INTENT(IN), OPTIONAL :: output_type, restart_type, patch_id, vlevel_type, compression_type
    CHARACTER(*), INTENT(IN), OPTIONAL :: post_suf, rest_suf, init_suf, model_type, filename
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, linitial  ! in standard output/restart/initial file
    INTEGER :: vln_len, ivl
    TYPE(t_var_list_ptr), ALLOCATABLE :: tmp_stor(:)

    vln_len = LEN_TRIM(vlname)
    CALL message('','')
    CALL message('','adding new var_list '//vlname(1:vln_len))
    IF (last_id .EQ. 0) THEN
      ALLOCATE(storage(8))
      CALL map%init(.false.)
    END IF
    CALL vlr_get(list, vlname)
    IF (ASSOCIATED(list%p)) CALL finish('new_var_list', ' >'//vlname(1:vln_len)//'< already in use.')
    ALLOCATE(list%p)
    last_id = last_id + 1
    CALL map%put(vlname, last_id)
    ! set default list characteristics
    list%p%vlname   = vlname(1:vln_len)
    list%p%post_suf = '_'//vlname(1:vln_len)
    list%p%rest_suf = list%p%post_suf
    list%p%init_suf = list%p%post_suf
    ! set non-default list characteristics
    list%p%restart_type = restart_file_type
    IF (PRESENT(output_type))      list%p%output_type      = output_type
    IF (PRESENT(restart_type))     list%p%restart_type     = restart_type
    IF (PRESENT(post_suf))         list%p%post_suf         = post_suf
    IF (PRESENT(rest_suf))         list%p%rest_suf         = rest_suf
    IF (PRESENT(init_suf))         list%p%init_suf         = init_suf
    IF (PRESENT(loutput))          list%p%loutput          = loutput
    IF (PRESENT(lrestart))         list%p%lrestart         = lrestart
    IF (PRESENT(linitial))         list%p%linitial         = linitial
    IF (PRESENT(patch_id))         list%p%patch_id         = patch_id
    IF (PRESENT(vlevel_type))      list%p%vlevel_type      = vlevel_type
    IF (PRESENT(filename))         list%p%filename         = filename
    IF (PRESENT(compression_type)) list%p%compression_type = compression_type
    IF (PRESENT(model_type))       list%p%model_type       = model_type
    IF (SIZE(storage) .LT. last_id) THEN
      CALL MOVE_ALLOC(storage, tmp_stor)
      ALLOCATE(storage(SIZE(tmp_stor) + 8))
      DO ivl = 1, SIZE(tmp_stor)
        IF (ASSOCIATED(tmp_stor(ivl)%p)) storage(ivl)%p => tmp_stor(ivl)%p
      END DO
    END IF
    storage(last_id)%p => list%p
  END SUBROUTINE vlr_add

  ! Get a reference to a memory buffer/output var_list
  SUBROUTINE vlr_get(list, vlname)
    TYPE(t_var_list_ptr), INTENT(OUT) :: list ! pointer
    CHARACTER(*), INTENT(IN) :: vlname
    INTEGER :: ivl, ierr

    CALL map%get(vlname, ivl, ierr)
    NULLIFY(list%p)
    IF (ierr .EQ. 0) THEN
      IF (ASSOCIATED(storage(ivl)%p)) &
        list%p => storage(ivl)%p
    END IF
  END SUBROUTINE vlr_get

  LOGICAL FUNCTION iter_next(this) RESULT(valid)
    CLASS(t_vl_register_iter), INTENT(INOUT) :: this

    valid = .false.
    IF (last_id .GT. 0) THEN
      IF (this%anew) THEN
        this%cur_idx = 0
        this%anew = .false.
      END IF
      DO WHILE(this%cur_idx .LT. SIZE(storage) .AND. .NOT.valid)
        this%cur_idx = this%cur_idx + 1
        IF (ASSOCIATED(storage(this%cur_idx)%p)) THEN
          valid = .true.
          this%cur%p => storage(this%cur_idx)%p
        END IF
      END DO
      IF (.NOT.valid) THEN
        NULLIFY(this%cur%p)
        this%anew = .true.
      END IF
    END IF
  END FUNCTION iter_next

  ! Delete an output var_list, nullify the associated pointer
  SUBROUTINE vlr_del(list)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: list
    INTEGER :: ivl, ierr

    IF (ASSOCIATED(list%p)) THEN
      CALL map%get(list%p%vlname, ivl, ierr)
      IF (ierr .NE. 0) CALL finish(modname//":delete_var_list", &
        & "var_list <" // TRIM(list%p%vlname) // "> not registered")
      CALL list%delete()
      DEALLOCATE(storage(ivl)%p)
      NULLIFY(list%p)
    END IF
  END SUBROUTINE vlr_del

  ! (Un)pack the var_lists to a t_packed_message
  SUBROUTINE vlr_packer(op, pmsg, nv_all)
    INTEGER, INTENT(IN) :: op
    TYPE(t_PackedMessage), INTENT(INOUT) :: pmsg
    INTEGER, INTENT(OUT), OPTIONAL :: nv_all
    INTEGER :: ivl, nvl, iv, nv, nv_al, patch_id, r_type, vl_type, ierr, infosize
    INTEGER, ALLOCATABLE :: info_buf(:)
    TYPE(t_var), POINTER :: elem
    CHARACTER(LEN=vlname_len) :: vl_name
    CHARACTER(LEN=32) :: m_type
    LOGICAL :: lre, lout, l_end
    TYPE(t_var_list_ptr) :: vlp
    TYPE(t_vl_register_iter) :: iter
    CHARACTER(*), PARAMETER :: routine = modname//':varlistPacker'

    IF (op .EQ. kUnpackOp .AND. last_id .NE. 0) &
      & CALL finish(routine, "var_list_register must be empty if receiver")
    nvl = last_id
    CALL pmsg%packer(op, nvl)
    nv_al = 0
    l_end = .false.
    infosize = var_metadata_get_size()
    DO ivl = 1, nvl
      IF(op .EQ. kPackOp) THEN
        IF (.NOT.iter_next(iter)) THEN
          l_end = .true.
          EXIT
        END IF
        vlp%p => iter%cur%p
        ! copy the values needed for the new_var_list() CALL to local variables
        lre      = vlp%p%lrestart
        lout     = vlp%p%loutput
        vl_name  = vlp%p%vlname
        m_type   = vlp%p%model_type
        patch_id = vlp%p%patch_id
        r_type   = vlp%p%restart_type
        vl_type  = vlp%p%vlevel_type
        nv = vlp%p%nvars
        nv_al = nv_al + nv
      END IF
      CALL pmsg%packer(op, l_end)
      IF (l_end) EXIT
      CALL pmsg%packer(op, lre)
      CALL pmsg%packer(op, lout)
      CALL pmsg%packer(op, vl_name)
      CALL pmsg%packer(op, m_type)
      CALL pmsg%packer(op, patch_id)
      CALL pmsg%packer(op, r_type)
      CALL pmsg%packer(op, vl_type)
      CALL pmsg%packer(op, nv)
      IF (nv .EQ. 0) CYCLE ! check if there are valid restart fields
      IF (op .EQ. kPackOp) THEN
        DO iv = 1, vlp%p%nvars
          elem => vlp%p%vl(iv)%p
          info_buf = var_metadata_toBinary(elem%info, infosize)
          CALL pmsg%pack(info_buf)
          CALL pmsg%pack(vlp%p%tl(iv))
          CALL pmsg%pack(vlp%p%hgrid(iv))
          CALL pmsg%pack(vlp%p%key(iv))
          CALL pmsg%pack(vlp%p%key_notl(iv))
          CALL pmsg%pack(vlp%p%lout(iv))
        END DO
      ELSE
        ! create var list
        CALL vlr_add(vlp, vl_name, patch_id=patch_id, model_type=m_type, &
          & restart_type=r_type, vlevel_type=vl_type, lrestart=lre, loutput=lout)
        vlp%p%nvars = nv
        ! insert elements into var list
        ALLOCATE(vlp%p%vl(nv), vlp%p%tl(nv), vlp%p%hgrid(nv), &
          & vlp%p%key(nv), vlp%p%key_notl(nv), vlp%p%lout(nv))
        DO iv = 1, nv
          ALLOCATE(vlp%p%vl(iv)%p, STAT=ierr)
          IF (ierr .NE. 0) CALL finish(routine, "memory allocation failure")
          elem => vlp%p%vl(iv)%p
          NULLIFY(elem%r_ptr, elem%s_ptr, elem%i_ptr, elem%l_ptr)
          elem%var_base_size = 0 ! Unknown here
          CALL pmsg%unpack(info_buf)
          elem%info = var_metadata_fromBinary(info_buf, infosize)
          CALL pmsg%unpack(vlp%p%tl(iv))
          CALL pmsg%unpack(vlp%p%hgrid(iv))
          CALL pmsg%unpack(vlp%p%key(iv))
          CALL pmsg%unpack(vlp%p%key_notl(iv))
          CALL pmsg%unpack(vlp%p%lout(iv))
        END DO
      END IF
    END DO
    IF (op .EQ. kPackOp) THEN
      IF (iter_next(iter)) CALL finish(routine, "inconsistency -- second kind")
    END IF
    CALL pmsg%packer(op, nv_al)
    IF (PRESENT(nv_all)) nv_all = nv_al
  END SUBROUTINE vlr_packer

END MODULE mo_var_list_register
