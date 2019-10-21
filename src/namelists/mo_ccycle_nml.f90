!>
!! Read configuration parameters as Fortran namelist from an external file. 
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ccycle_nml

  USE mo_ccycle_config    ,ONLY: ccycle_config, init_ccycle_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_ccycle_nml

  NAMELIST /ccycle_nml/ ccycle_config

CONTAINS

    SUBROUTINE process_ccycle_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_ccycle_config
    !
    CALL process_nml(filename, 'ccycle_nml', nml_read, nml_write)
    !
  CONTAINS
    !
    SUBROUTINE nml_read(funit)
      INTEGER, INTENT(in) :: funit
      READ(funit, NML=ccycle_nml)
    END SUBROUTINE nml_read
    !
    SUBROUTINE nml_write(funit)
      INTEGER, INTENT(in) :: funit
      WRITE(funit, NML=ccycle_nml)
    END SUBROUTINE nml_write
    !
  END SUBROUTINE process_ccycle_nml

END MODULE mo_ccycle_nml
