MODULE mo_grib1

  IMPLICIT NONE

  PRIVATE

  TYPE grib1_global
    INTEGER :: centre
    INTEGER :: subcentre
  END TYPE grib1_global

  TYPE grib1_var
    INTEGER :: table
    INTEGER :: parameter
    INTEGER :: bits
    INTEGER :: leveltype
  END type grib1_var

  PUBLIC :: grib1_global
  PUBLIC :: grib1_var

END MODULE mo_grib1
