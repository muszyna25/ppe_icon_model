!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
#include <icon_contiguous_defines.h>

MODULE mo_var_list_element

  USE mo_kind,               ONLY: dp, sp
  USE mo_var_metadata_types, ONLY: t_var_metadata, t_var_metadata_dynamic

  IMPLICIT NONE
  PRIVATE

  ! export index constants: model/pressure/height levels
  PUBLIC :: level_type_ml, level_type_pl, level_type_hl, level_type_il, lev_type_str
  PUBLIC :: t_var_list_element, t_p_var_list_element

  ! constants defining level type:
  INTEGER, PARAMETER             :: level_type_ml = 1
  INTEGER, PARAMETER             :: level_type_pl = 2
  INTEGER, PARAMETER             :: level_type_hl = 3
  INTEGER, PARAMETER             :: level_type_il = 4
  ! string defining level type:
  CHARACTER(LEN=2), PARAMETER    :: lev_type_str(4) = (/ 'ML', 'PL', 'HL', 'IL' /)

  TYPE t_var_list_element
    REAL(dp), CONTIGUOUS_POINTER :: r_ptr(:,:,:,:,:) => NULL()
    REAL(sp), CONTIGUOUS_POINTER :: s_ptr(:,:,:,:,:) => NULL()
    INTEGER,  CONTIGUOUS_POINTER :: i_ptr(:,:,:,:,:) => NULL()
    LOGICAL,  CONTIGUOUS_POINTER :: l_ptr(:,:,:,:,:) => NULL()
    INTEGER                      :: var_base_size      ! generic size in bytes of variable used
    TYPE(t_var_metadata)         :: info               ! meta data for this entry
    TYPE(t_var_metadata_dynamic) :: info_dyn           ! dynamic meta data for this entry (see type description)
  END TYPE t_var_list_element

  TYPE t_p_var_list_element
    TYPE(t_var_list_element), POINTER :: p => NULL()
  END TYPE t_p_var_list_element

END MODULE mo_var_list_element
