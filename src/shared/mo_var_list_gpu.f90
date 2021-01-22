! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_gpu

  USE mo_impl_constants,      ONLY: REAL_T, SINGLE_T, INT_T, BOOL_T
  USE mo_var_metadata_types,  ONLY: t_var_metadata
  USE mo_var_list,            ONLY: t_var_list_ptr
  USE mo_var,                 ONLY: t_var
  USE mo_var_list_register,   ONLY: vlr_get
!  USE mo_exception,           ONLY: message, message_text

  IMPLICIT NONE
  PRIVATE

  CHARACTER(*), PARAMETER :: modname = 'mo_var_list_gpu'

  PUBLIC :: gpu_update_var_list

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
    CHARACTER(:), ALLOCATABLE :: listname
!    CHARACTER(*), PARAMETER :: d2h = "dev => host", h2d = "host => dev"
    INTEGER :: ii
    CHARACTER(LEN=2) :: i2a

    listname = vlname
    IF (PRESENT(domain)) THEN
      WRITE(i2a, "(i2.2)") domain
      listname = listname//i2a
    END IF
    IF (PRESENT(substr))  listname = listname//substr
    IF (PRESENT(timelev)) THEN  
      WRITE(i2a, "(i2.2)") timelev
      listname = listname//i2a
    END IF
    CALL vlr_get(list, listname)
!    WRITE(message_text, "(a,l1)") MERGE(h2d, d2h, to_device)//" update <"//listname//"> found=", ASSOCIATED(list%p)
!    CALL message("", message_text)
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
