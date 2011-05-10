MODULE mo_var_list_element

  USE mo_kind,         ONLY: dp
  USE mo_var_metadata, ONLY: t_var_metadata

  IMPLICIT NONE

  PRIVATE

  TYPE t_var_list_element
    REAL(dp), POINTER    :: r_ptr(:,:,:,:,:)   ! pointer to 4D-field
    INTEGER,  POINTER    :: i_ptr(:,:,:,:,:)   ! pointer to 4D-field
    LOGICAL,  POINTER    :: l_ptr(:,:,:,:,:)   ! pointer to 4D-field
    INTEGER              :: var_base_size      ! generic size in bytes of variable used
    TYPE(t_var_metadata) :: info               ! meta data for this entry
  END type t_var_list_element

  PUBLIC :: t_var_list_element

END MODULE mo_var_list_element
