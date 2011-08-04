
MODULE test_coupler_mod

  USE mo_icon_cpl, ONLY : cplout, complist

  USE mo_kind, ONLY     : wp

  PUBLIC

  ! Loop count

  INTEGER                          :: i, n

  ! grid characteristic

  INTEGER, PARAMETER               :: grid_loc_dim = 100
  INTEGER                          :: grid_glo_dim
  INTEGER                          :: grid_cmp_dim

  ! process counting

  INTEGER                          :: nbr_ocean_procs
  INTEGER                          :: nbr_atmos_procs

  ! Arguments for cpl interfaces
  ! ----------------------------
  !
  ! Grid description

  INTEGER                          :: grid_id
  INTEGER                          :: grid_shape (2)

  ! Coupling field description

  CHARACTER(len=132)               :: field_name
  INTEGER                          :: field_id(4)
  INTEGER                          :: field_shape(3)

  ! Return code for error handling

  INTEGER                          :: info

END MODULE test_coupler_mod
