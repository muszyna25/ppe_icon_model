MODULE mo_grib2

  IMPLICIT NONE

  PRIVATE

  TYPE grib2_global
    INTEGER :: centre
    INTEGER :: subcentre
  END TYPE grib2_global

  TYPE grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: parameter
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
    INTEGER :: leveltype
  END type grib2_var

  PUBLIC :: grib2_global
  PUBLIC :: grib2_var

END MODULE mo_grib2
