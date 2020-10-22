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
  USE mo_var_list,            ONLY: t_var_list_ptr, t_list_element
  USE mo_var_list_register,   ONLY: vlr_get
  USE mo_util_string,         ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_var_list_gpu'

  PUBLIC :: gpu_h2d_var_list
  PUBLIC :: gpu_d2h_var_list

CONTAINS

  !> Update device data of variable list
  !
  SUBROUTINE gpu_h2d_var_list(name, domain, substr, timelev)
    CHARACTER(len=*), INTENT(IN)            :: name    ! name of output var_list
    INTEGER, INTENT(IN), OPTIONAL           :: domain  ! domain index to append
    CHARACTER(len=*), INTENT(IN), OPTIONAL  :: substr  ! String after domain, before timelev
    INTEGER, INTENT(IN), OPTIONAL           :: timelev ! timelev index to append

    TYPE(t_var_list_ptr) :: list
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_list_element), POINTER :: element
    CHARACTER(:), ALLOCATABLE :: listname

    ! generate listname
    listname = TRIM(name) // &
      & MERGE(int2string(domain, "(i2.2)"), '', PRESENT(domain)) // &
      & MERGE(TRIM(substr),                 '', PRESENT(substr)) // &
      & MERGE(int2string(domain, "(i2.2)"), '', PRESENT(timelev))
    CALL vlr_get(list, listname)

    element => list%p%first_list_element
    for_all_list_elements: DO WHILE (ASSOCIATED(element))
      info    => element%field%info
      SELECT CASE(info%data_type)
      CASE (REAL_T)
        !$ACC UPDATE DEVICE(element%field%r_ptr) IF(info%lopenacc)
      CASE (SINGLE_T)
        !$ACC UPDATE DEVICE(element%field%s_ptr) IF(info%lopenacc)
      CASE (INT_T)
        !$ACC UPDATE DEVICE(element%field%i_ptr) IF(info%lopenacc)
      CASE (BOOL_T)
        !$ACC UPDATE DEVICE(element%field%l_ptr) IF(info%lopenacc)
      END SELECT
      element => element%next_list_element
    END DO for_all_list_elements

  END SUBROUTINE gpu_h2d_var_list

  !> Update host data of variable list
  !
  SUBROUTINE gpu_d2h_var_list( name, domain, substr, timelev )
    CHARACTER(len=*), INTENT(IN)            :: name    ! name of output var_list
    INTEGER, INTENT(IN), OPTIONAL           :: domain  ! domain index to append
    CHARACTER(len=*), INTENT(IN), OPTIONAL  :: substr  ! String after domain, before timelev
    INTEGER, INTENT(IN), OPTIONAL           :: timelev ! timelev index to append

    TYPE(t_var_list_ptr) :: list
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_list_element), POINTER :: element
    CHARACTER(:), ALLOCATABLE :: listname

    ! generate listname
    listname = TRIM(name) // &
      & MERGE(int2string(domain, "(i2.2)"), '', PRESENT(domain)) // &
      & MERGE(TRIM(substr),                 '', PRESENT(substr)) // &
      & MERGE(int2string(domain, "(i2.2)"), '', PRESENT(timelev))
    CALL vlr_get(list, listname)

    element => list%p%first_list_element
    for_all_list_elements: DO WHILE (ASSOCIATED(element))
      info    => element%field%info
      SELECT CASE(info%data_type)
      CASE (REAL_T)
        !$ACC UPDATE HOST(element%field%r_ptr) IF(info%lopenacc)
      CASE (SINGLE_T)
        !$ACC UPDATE HOST(element%field%s_ptr) IF(info%lopenacc)
      CASE (INT_T)
        !$ACC UPDATE HOST(element%field%i_ptr) IF(info%lopenacc)
      CASE (BOOL_T)
        !$ACC UPDATE HOST(element%field%l_ptr) IF(info%lopenacc)
      END SELECT
      element => element%next_list_element
    END DO for_all_list_elements

  END SUBROUTINE gpu_d2h_var_list

END MODULE mo_var_list_gpu
