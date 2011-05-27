!>
!! @brief 
!!        
!! @par Revision History
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_master_nml

  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, POSITIONED

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC :: master_nml_setup

  !-------------------------------------------------------------------------
  ! Namelist variables 
  !-------------------------------------------------------------------------
  LOGICAL :: lrestart

  NAMELIST/master_nml/ lrestart

CONTAINS
  !>
  !! @brief 
  !!
  !! @par Revision History
  !!
  SUBROUTINE master_nml_setup
                                                
    INTEGER  :: ist
 
    ! Default values

    lrestart = .FALSE. 

    ! Read user's (new) specifications (Done so far by all MPI processes)

    CALL position_nml ('master_nml', status=ist)
    IF (ist == POSITIONED) THEN
      READ (nnml, master_nml)
    ENDIF

    ! Write the contents of the namelist to an ASCII file.

    IF(p_pe == p_io) WRITE(nnml_output,nml=master_nml)

 END SUBROUTINE master_nml_setup

END MODULE mo_master_nml
