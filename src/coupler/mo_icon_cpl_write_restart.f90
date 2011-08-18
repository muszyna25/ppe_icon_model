!>
!! Routine to check whether we need to provide a restart field
!! for the first coupling step
!!
!! <Description>
!! 
!! @author Rene Redler, MPI-M
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright
!! 2010-2011 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and WARRANTY conditions.
!! 
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! &ltol>
!! &ltli> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! &ltli> The code may not be re-distributed without the consent of the authors.
!! &ltli> The copyright notice and statement of authorship must appear in all
!!    copies.
!! &ltli> You accept the warranty conditions (see WARRANTY).
!! &ltli> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_icon_cpl_write_restart

  USE mo_kind, ONLY           : wp
  USE mo_io_units, ONLY       : find_next_free_unit
  USE mo_icon_cpl, ONLY       : t_cpl_field, cpl_fields

  IMPLICIT NONE

  PRIVATE

  TYPE(t_cpl_field), POINTER :: fptr

  INTEGER                :: len
  INTEGER                :: rest_unit
  INTEGER                :: year, month, day, hour, minute, second

  CHARACTER(len=132) :: file_name
  CHARACTER(len=132) :: field_name

  PUBLIC :: ICON_cpl_write_restart

CONTAINS

  ! ---------------------------------------------------------------------

  SUBROUTINE ICON_cpl_write_restart ( field_id, field_shape, coupling_field, count, ierror )

    INTEGER, INTENT(in)    :: field_id         !<  field id
    INTEGER, INTENT(in)    :: field_shape(3)   !<  shape of send field
    INTEGER, INTENT(in)    :: count            !<  number of accumulations

    REAL (wp), INTENT(in)  :: coupling_field (field_shape(1):field_shape(2),field_shape(3))

    INTEGER, INTENT(out)   :: ierror           !<  returned error code

    ierror = 0

    WRITE(file_name(1:7), '(A7)') 'restart'

    fptr => cpl_fields(field_id)

    len = LEN_TRIM(fptr%field_name) - 1
    WRITE(file_name(8:8+len), '(A7)') TRIM( fptr%field_name )

    rest_unit = find_next_free_unit(10,100)
    len       = len + 8

    OPEN   ( unit = rest_unit, file = file_name(1:len), status = 'UNKNOWN' )
    REWIND ( unit = rest_unit )
    WRITE  ( unit = rest_unit ) field_name
    WRITE  ( unit = rest_unit ) year, month, day, hour, minute, second
    WRITE  ( unit = rest_unit ) count
    WRITE  ( unit = rest_unit ) coupling_field
    CLOSE  ( unit = rest_unit )

  END SUBROUTINE ICON_cpl_write_restart

END MODULE mo_icon_cpl_write_restart
