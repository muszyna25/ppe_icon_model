!>
!! Contains the setup of the sleve coordinate
!!
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
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the prep_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: i_oper_mode   ! operation mode
  INTEGER  :: nlev_in       ! number of model levels of input data
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: l_w_in        ! Logical switch if w is provided as input
  LOGICAL  :: l_sfc_in      ! Logical switch if surface fields are provided as input
  LOGICAL  :: l_zp_out      ! Logical switch for diagnostic output on pressure and height levels

  NAMELIST /prepicon_nml/ i_oper_mode, nlev_in, zpbl1, zpbl2, &
                          l_w_in, l_zp_out, nlevsoil_in, l_sfc_in
  !
  !
  !

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

 SUBROUTINE prepicon_nml_setup( filename )
    
  CHARACTER(LEN=*), INTENT(IN) :: filename

  !local variable
  INTEGER :: i_status

  !------------------------------------------------------------
  ! 2.0 set up the default values for prepicon
  !------------------------------------------------------------
  !
  !
  i_oper_mode = 1           ! operation mode
  nlev_in     = 91          ! number of levels of input data
  nlevsoil_in = 4           ! number of soil levels of input data
  zpbl1       = 500._wp     ! AGL heights used for computing vertical 
  zpbl2       = 1000._wp    ! gradients
  l_w_in      = .FALSE.     ! true: w is provided as input
  l_sfc_in    = .TRUE.      ! true: surface fields are provided as input
  l_zp_out    = .FALSE.     ! true: diagnostic output on p and z levels
  !
  !
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
  !
  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  !currently no consistency check...

  ! write the contents of the namelist to an ASCII file

  IF(my_process_is_stdio()) WRITE(nnml_output,nml=prepicon_nml)

END SUBROUTINE prepicon_nml_setup

END MODULE mo_prepicon_nml
