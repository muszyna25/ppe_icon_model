!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_var_list_element

  USE mo_kind,               ONLY: dp
  USE mo_var_metadata_types, ONLY: t_var_metadata, t_var_metadata_dynamic

  IMPLICIT NONE

  PRIVATE

  ! export index constants: model/pressure/height levels
  PUBLIC level_type_ml, level_type_pl, level_type_hl, level_type_il, lev_type_str

  ! constants defining level type:
  INTEGER, PARAMETER             :: level_type_ml = 1
  INTEGER, PARAMETER             :: level_type_pl = 2
  INTEGER, PARAMETER             :: level_type_hl = 3
  INTEGER, PARAMETER             :: level_type_il = 4
  ! string defining level type:
  CHARACTER(LEN=2), PARAMETER    :: lev_type_str(4) = (/ 'ML', 'PL', 'HL', 'IL' /)

  TYPE t_var_list_element
    REAL(dp), POINTER            :: r_ptr(:,:,:,:,:)   ! pointer to 4D-field
    INTEGER,  POINTER            :: i_ptr(:,:,:,:,:)   ! pointer to 4D-field
    LOGICAL,  POINTER            :: l_ptr(:,:,:,:,:)   ! pointer to 4D-field
    INTEGER                      :: var_base_size      ! generic size in bytes of variable used
    TYPE(t_var_metadata)         :: info               ! meta data for this entry
    TYPE(t_var_metadata_dynamic) :: info_dyn           ! dynamic meta data for this entry (see type description)
  END type t_var_list_element

  PUBLIC :: t_var_list_element

END MODULE mo_var_list_element
