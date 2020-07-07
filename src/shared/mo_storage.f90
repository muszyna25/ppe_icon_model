!>
!! This module provides a simple key value storage based on string comparisons
!! using the generic hash tables provided by mo_hash_table. 
!! As method to calculate hash keys, the DJB algorithm is used.
!!
!!
!! @par Revision History
!! Initial release by Daniel Rieger, KIT (2016-12-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_storage
  USE mo_exception,       ONLY: message
  USE mo_key_value_store, ONLY: t_key_value_store

  IMPLICIT NONE
  PRIVATE

  TYPE, EXTENDS(t_key_value_store) :: t_storage
  CONTAINS
    PROCEDURE, PUBLIC :: dump => dump_storage
  END TYPE t_storage
  
  CHARACTER(*), PARAMETER :: modname = 'mo_storage'

  PUBLIC :: t_storage
  
CONTAINS

  SUBROUTINE dump_storage(me, opt_label_in)
    CLASS(t_storage),INTENT(IN), TARGET :: me
    CHARACTER(*), INTENT(IN), OPTIONAL :: opt_label_in
    CHARACTER(:), ALLOCATABLE :: opt_label
    CHARACTER(*), PARAMETER :: routine = modname // "dump_storage"
  
    opt_label = "noName"
    IF (PRESENT(opt_label_in)) opt_label = TRIM(opt_label_in)
    CALL message(routine,"START STORAGE DUMP (" // opt_label // ")")
    CALL me%output(label=opt_label_in)
    IF (me%getEntryCount() .EQ. 0) CALL message(routine,"structure is empty")
    CALL message(routine,"END OF STORAGE DUMP (" // opt_label // ")")
  END SUBROUTINE dump_storage
END MODULE mo_storage
