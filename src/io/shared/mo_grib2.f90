!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_grib2

  IMPLICIT NONE

  PRIVATE

  TYPE t_grib2_global
    INTEGER :: centre
    INTEGER :: subcentre
    INTEGER :: generating_process
  END TYPE t_grib2_global

  TYPE t_grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
  END TYPE t_grib2_var

  PUBLIC :: t_grib2_global
  PUBLIC :: t_grib2_var

END MODULE mo_grib2
