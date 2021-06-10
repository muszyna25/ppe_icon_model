! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_gpu

  USE mo_impl_constants,      ONLY: vlname_len, REAL_T, SINGLE_T, INT_T, BOOL_T
  USE mo_var_metadata_types,  ONLY: t_var_metadata
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_var,                 ONLY: t_var
  USE mo_var_list_register,   ONLY: vlr_get

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: gpu_update_var_list

  CHARACTER(*), PARAMETER :: modname = 'mo_var_list_gpu'

CONTAINS

  !> Update data of variable list on device or host
  ! we expect "trimmed" character variables as input!
  SUBROUTINE gpu_update_var_list(vlname, to_device, domain, substr, timelev)
    CHARACTER(*), INTENT(IN) :: vlname    ! name of output var_list
    LOGICAL, INTENT(IN) :: to_device ! direction of update (true, if host->device)
    INTEGER, INTENT(IN), OPTIONAL :: domain, timelev  ! domain/timelev index to append
    CHARACTER(*), INTENT(IN), OPTIONAL :: substr  ! String after domain, before timelev
    TYPE(t_var_list_ptr) :: list
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_var), POINTER :: element
    CHARACTER(LEN=vlname_len) :: listname
    INTEGER :: ii,vln_pos, subs_len
    CHARACTER(LEN=2) :: i2a

    IF (PRESENT(domain)) THEN
      WRITE(listname, "(a,i2.2)") TRIM(vlname), domain
    ELSE
      listname = vlname
    END IF
    vln_pos = LEN_TRIM(listname)
    IF (PRESENT(substr)) THEN
      subs_len = LEN_TRIM(substr)
      listname(vln_pos+1:subs_len+vln_pos) = substr(1:subs_len)
      vln_pos = vln_pos + subs_len
    END IF
    IF (PRESENT(timelev)) THEN
      WRITE(listname(vln_pos+1:), "(i2.2)") timelev
      vln_pos = vln_pos + 2
    END IF
    CALL vlr_get(list, listname(1:vln_pos))
    IF (ASSOCIATED(list%p)) THEN
      DO ii = 1, list%p%nvars
        element => list%p%vl(ii)%p
        info    => element%info
        IF (to_device) THEN
          CALL upd_dev()
        ELSE
          CALL upd_host()
        END IF
      END DO
    END IF
  CONTAINS

    SUBROUTINE upd_dev()

      SELECT CASE(info%data_type)
      CASE (REAL_T)
        !$ACC UPDATE DEVICE(element%r_ptr) IF(info%lopenacc)
      CASE (SINGLE_T)
        !$ACC UPDATE DEVICE(element%s_ptr) IF(info%lopenacc)
      CASE (INT_T)
        !$ACC UPDATE DEVICE(element%i_ptr) IF(info%lopenacc)
      CASE (BOOL_T)
        !$ACC UPDATE DEVICE(element%l_ptr) IF(info%lopenacc)
      END SELECT
    END SUBROUTINE upd_dev

    SUBROUTINE upd_host()
      
      SELECT CASE(info%data_type)
      CASE (REAL_T)
        !$ACC UPDATE HOST(element%r_ptr) IF(info%lopenacc)
      CASE (SINGLE_T)
        !$ACC UPDATE HOST(element%s_ptr) IF(info%lopenacc)
      CASE (INT_T)
        !$ACC UPDATE HOST(element%i_ptr) IF(info%lopenacc)
      CASE (BOOL_T)
        !$ACC UPDATE HOST(element%l_ptr) IF(info%lopenacc)
      END SELECT
    END SUBROUTINE upd_host
  END SUBROUTINE gpu_update_var_list

END MODULE mo_var_list_gpu
