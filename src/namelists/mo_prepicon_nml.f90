!>
!!        
!! @par Revision History
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
!!
MODULE mo_prepicon_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: max_char_length, max_dom, &
    &                              MODE_IFSANA, MODE_DWDANA, MODE_REMAP
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_prepicon_config,    ONLY: &
    & config_i_oper_mode        => i_oper_mode

  IMPLICIT NONE

  PUBLIC :: read_prepicon_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the prep_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: i_oper_mode   ! operation mode


  NAMELIST /prepicon_nml/ i_oper_mode

  
CONTAINS

!-------------------------------------------------------------------------
!
!
 !>
 !!  Initialization of the prepicon coordinate namelist
 !!
 !!
 !! @par Revision History
 !!  Initial version by Guenther Zaengl (2011-07-11)

 SUBROUTINE read_prepicon_namelist( filename )
    
  CHARACTER(LEN=*), INTENT(IN) :: filename

  !local variable
  INTEGER :: i_status

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_prepicon_nml: read_prepicon_namelist'

  !------------------------------------------------------------
  ! 2.0 set up the default values for prepicon
  !------------------------------------------------------------
  !
  !
  i_oper_mode = MODE_IFSANA


  !------------------------------------------------------------
  ! 3.0 Read the prepicon namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL open_nml(TRIM(filename))
  CALL position_nml ('prepicon_nml', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, prepicon_nml)
  END SELECT
  CALL close_nml

  !------------------------------------------------------------
  ! 4.0 Fill the configuration state
  !------------------------------------------------------------

  config_i_oper_mode       = i_oper_mode


  !------------------------------------------------------------
  ! 5.0 check the consistency of the parameters
  !------------------------------------------------------------


  ! write the contents of the namelist to an ASCII file

  IF(my_process_is_stdio()) WRITE(nnml_output,nml=prepicon_nml)

END SUBROUTINE read_prepicon_namelist

END MODULE mo_prepicon_nml
