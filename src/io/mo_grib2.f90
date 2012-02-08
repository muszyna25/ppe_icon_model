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
    INTEGER :: parameter
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
  END TYPE t_grib2_var

  PUBLIC :: t_grib2_global
  PUBLIC :: t_grib2_var

END MODULE mo_grib2
