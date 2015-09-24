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

  ! max. number of additional GRIB2 integer keys per  variable
  INTEGER, PARAMETER :: MAX_INT_KEYS = 5



  TYPE t_grib2_global
    INTEGER :: centre
    INTEGER :: subcentre
    INTEGER :: generating_process
  END TYPE t_grib2_global

  TYPE t_grib2_int_key
    CHARACTER(len=32) :: key
    INTEGER           :: val
  END TYPE t_grib2_int_key

  TYPE t_grib2_key_list
    INTEGER                :: nint_keys ! no. of integer keys
    TYPE (t_grib2_int_key) :: int_key(MAX_INT_KEYS)    
  END TYPE t_grib2_key_list

  TYPE t_grib2_var
    INTEGER :: discipline
    INTEGER :: category
    INTEGER :: number
    INTEGER :: bits
    INTEGER :: gridtype
    INTEGER :: subgridtype
    
    ! list of additional GRIB2 key/value pairs
    TYPE (t_grib2_key_list) :: additional_keys
  END TYPE t_grib2_var

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE grib2_key_list_plus
  END INTERFACE OPERATOR(+)

  PUBLIC :: t_grib2_global
  PUBLIC :: t_grib2_var
  PUBLIC :: t_grib2_int_key
  PUBLIC :: grib2_var
  PUBLIC :: OPERATOR(+)

CONTAINS

  ! constructor for GRIB2 derived data type
  FUNCTION grib2_var(discipline, category, number, bits, gridtype, subgridtype)
    INTEGER, INTENT(IN) :: discipline
    INTEGER, INTENT(IN) :: category
    INTEGER, INTENT(IN) :: number
    INTEGER, INTENT(IN) :: bits
    INTEGER, INTENT(IN) :: gridtype
    INTEGER, INTENT(IN) :: subgridtype
    TYPE(t_grib2_var) :: grib2_var

    grib2_var%discipline  = discipline 
    grib2_var%category    = category   
    grib2_var%number      = number     
    grib2_var%bits        = bits       
    grib2_var%gridtype    = gridtype   
    grib2_var%subgridtype = subgridtype

    grib2_var%additional_keys%nint_keys = 0
  END FUNCTION grib2_var

  FUNCTION grib2_key_list_plus(a, b)
    TYPE(t_grib2_var) :: grib2_key_list_plus
    TYPE(t_grib2_var),     INTENT(IN) :: a
    TYPE(t_grib2_int_key), INTENT(IN) :: b
    ! local variables
    INTEGER :: i

    grib2_key_list_plus = a
    IF (a%additional_keys%nint_keys < MAX_INT_KEYS) THEN
      i = grib2_key_list_plus%additional_keys%nint_keys
      grib2_key_list_plus%additional_keys%nint_keys = i + 1
      grib2_key_list_plus%additional_keys%int_key(i+1) = b
    END IF
  END FUNCTION grib2_key_list_plus

END MODULE mo_grib2
