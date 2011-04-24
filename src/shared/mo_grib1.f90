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
