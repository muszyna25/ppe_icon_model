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
  USE mo_mpi,                       ONLY: p_comm_work_2_restart
  USE mo_packed_message,            ONLY: t_PackedMessage, kPackOp, kUnpackOp
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_util,              ONLY: restartBcastRoot
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
    PROCEDURE :: transferToRestart => restartPatchData_transferToRestart
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

  CHARACTER(*), PARAMETER :: modname = "mo_restart_patch_data"

CONTAINS

  ! collective across restart AND worker PEs
  SUBROUTINE restartPatchData_transferToRestart(me)
    CLASS(t_RestartPatchData), INTENT(INOUT) :: me
#ifndef NOMPI
    TYPE(t_PackedMessage) :: packedMessage

    CALL me%description%packer(kPackOp, packedMessage)
    CALL packedMessage%bcast(restartBcastRoot(), p_comm_work_2_restart)   ! transfer data to restart PEs
    CALL me%description%packer(kUnpackOp, packedMessage)
#endif
  END SUBROUTINE restartPatchData_transferToRestart

END MODULE mo_restart_patch_data
