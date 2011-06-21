
MODULE test_coupler_mod

  USE mo_kind, ONLY     : wp
  USE mo_datetime

  PUBLIC

  ! Loop count

  INTEGER                          :: i, n

  ! grid characteristic

  INTEGER                          :: grid_loc_dim
  INTEGER                          :: grid_glo_dim
  INTEGER                          :: grid_cmp_dim
  INTEGER                          :: nbr_land_points

  ! process counting

  INTEGER                          :: n_pes, my_pe

  INTEGER                          :: nbr_ocean_procs
  INTEGER                          :: nbr_atmos_procs

  ! Arguments for cpl interfaces
  ! ----------------------------
  !
  ! Component description

  CHARACTER(len=32)                :: comp_name

  ! Grid description

  INTEGER                          :: grid_id
  INTEGER                          :: grid_shape (2)

  ! Coupling field description

  CHARACTER(len=132)               :: field_name
  INTEGER                          :: field_id(4)
  INTEGER                          :: field_shape(3)

  ! Return code for error handling

  CHARACTER(len=132)               :: err_string
  INTEGER                          :: info
  INTEGER                          :: len
  INTEGER                          :: ierror

  ! MPI Communicator handling

  INTEGER                          :: ICON_atmos_comm
  INTEGER                          :: n_atm_pes
  INTEGER                          :: atm_pe

  INTEGER                          :: ICON_ocean_comm
  INTEGER                          :: n_ocn_pes
  INTEGER                          :: ocn_pe

END MODULE test_coupler_mod
