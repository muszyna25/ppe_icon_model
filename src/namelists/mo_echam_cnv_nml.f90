!>
!! Read configuration parameters as Fortran namelist from an external file. 
!!
!! @author Marco, Giorgetta, MPI-M, 2017-12
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
MODULE mo_echam_cnv_nml

  USE mo_echam_cnv_config ,ONLY: echam_cnv_config, init_echam_cnv_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_echam_cnv_nml

  NAMELIST /echam_cnv_nml/ echam_cnv_config

CONTAINS

  SUBROUTINE process_echam_cnv_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_echam_cnv_config
    !
    CALL process_nml(filename, 'echam_cnv_nml', nml_read, nml_write)
    !
  CONTAINS
    !
    SUBROUTINE nml_read(funit)
      INTEGER, INTENT(in) :: funit
      READ(funit, NML=echam_cnv_nml)
    END SUBROUTINE nml_read
    !
    SUBROUTINE nml_write(funit)
      INTEGER, INTENT(in) :: funit
      WRITE(funit, NML=echam_cnv_nml)
    END SUBROUTINE nml_write
    !
  END SUBROUTINE process_echam_cnv_nml

END MODULE mo_echam_cnv_nml
