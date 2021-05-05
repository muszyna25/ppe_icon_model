!>
!! The base CLASS for the PUBLIC INTERFACE for restart writing.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

MODULE mo_restart_patch_data
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_var,                       ONLY: t_var_ptr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_RestartPatchData

  ! this type stores all the information that we need to know about a patch and its variables
  TYPE, ABSTRACT :: t_RestartPatchData
    TYPE(t_restart_patch_description) :: description
    TYPE(t_var_ptr), ALLOCATABLE :: varData(:)
    INTEGER :: restartType
  CONTAINS
    PROCEDURE(i_construct), DEFERRED :: construct
    PROCEDURE(i_writeData), DEFERRED :: writeData
    PROCEDURE(i_destruct), DEFERRED :: destruct
  END TYPE t_RestartPatchData

  ABSTRACT INTERFACE
  SUBROUTINE i_construct(me, modelType, jg)
    IMPORT t_RestartPatchData
    CLASS(t_RestartPatchData), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN) :: jg
  END SUBROUTINE i_construct

  SUBROUTINE i_writeData(me, file_handle)
    IMPORT t_RestartPatchData
    CLASS(t_RestartPatchData), INTENT(INOUT) :: me
    INTEGER, INTENT(IN) :: file_handle
  END SUBROUTINE i_writeData

  SUBROUTINE i_destruct(me)
    IMPORT t_RestartPatchData
    CLASS(t_RestartPatchData), INTENT(INOUT) :: me
  END SUBROUTINE i_destruct
  END INTERFACE

END MODULE mo_restart_patch_data
