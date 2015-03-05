!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_grib1

  IMPLICIT NONE

  PRIVATE

  TYPE t_grib1_global
    INTEGER :: centre
    INTEGER :: subcentre
  END TYPE t_grib1_global

  TYPE t_grib1_var
    INTEGER :: table
    INTEGER :: parameter
    INTEGER :: bits
    INTEGER :: leveltype
  END type t_grib1_var

  PUBLIC :: t_grib1_global
  PUBLIC :: t_grib1_var

END MODULE mo_grib1
