!>
!! @brief Setup for synthetic radar data on the model grid
!!
!! configuration setup for synthetic radar data on the model grid
!!
!! @author Ulrich Blahak, DWD
!!
!!
!! @par Revision History
!! Initial revision by Ulrich Blahak, DWD (2021-11-29)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_synradar_config
  
  USE mo_io_units,                ONLY: filename_max
  USE radar_dbzcalc_params_type,  ONLY: dbzcalc_params

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Namelist parameters
  !--------------------------------------------------------------------------

  ! Meta data for reflectivity computations (DBZ, DBZ850, DBZ_CMAX, etc.) on the model grid by using advanced methods
  !  from EMVORADO (Mie-scattering, T-matrix):
  TYPE(dbzcalc_params)        :: synradar_meta
  CHARACTER(LEN=filename_max) :: ydir_mielookup_read
  CHARACTER(LEN=filename_max) :: ydir_mielookup_write


END MODULE mo_synradar_config
