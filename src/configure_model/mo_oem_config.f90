!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_oem_config

  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: filename_max
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Basic configuration setup for the online emission module
  !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! ghgctrl_nml:
    !------------------------------------------------------------------------
    CHARACTER(LEN=filename_max) :: &
      &  vertical_profile_nc,      & !< name of the oae vertical profile
      &  hour_of_day_nc,           & !< name of the oae hour of day file
      &  day_of_week_nc,           & !< name of the oae day of week file
      &  month_of_year_nc,         & !< name of the oae month of year file
      &  hour_of_year_nc,          & !< name of the oae hour of year file
      &  gridded_emissions_nc        !< name of the oae gridded emission file

    REAL(wp) ::                    restart_init_time


END MODULE mo_oem_config
