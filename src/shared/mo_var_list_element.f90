MODULE mo_var_list_element

  USE mo_var_metadata

  IMPLICIT NONE

  TYPE var_list_element
    REAL(dp), POINTER   :: ptr(:,:,:,:)   ! pointer to 4D-field
    TYPE (var_metadata) :: info           ! meta data for this entry
  END type var_list_element

  PUBLIC :: var_list_element

END MODULE mo_var_list_element
