!>
!! Subroutine setup_echam_phy calls various setup_xyz routines to
!! read the user's namelist specifications, and perform a
!! "cross check" to make sure the simulation configuration
!! is valid.
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, 2010-07-20
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
MODULE mo_echam_phy_setup

! USE mo_exception,          ONLY: message, finish, print_value
! USE mo_echam_phy_config,   ONLY: phy_config => echam_phy_config
! USE mo_echam_conv_config,  ONLY: config_echam_convection
! USE mo_run_nml,            ONLY: ltestcase, ntracer, io3
! USE mo_hydro_testcases,    ONLY: ctest_name

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: setup_echam_phy

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_nml'

CONTAINS
  !>
  !! Read and check namelist variables related to ECHAM physics
  !!
  SUBROUTINE setup_echam_phy

   !IF (phy_config%lconv) CALL config_echam_convection( nlev, vct_a, vct_b, ceta )

    ! Check whether echam_phy_nml is properly set for test cases;
    ! Check whether the process namelists are consistent with
    ! the configuration of dynamics and transport; 

  ! CALL cross_check

  END SUBROUTINE setup_echam_phy
  !-------------
  !>
  !!
 !SUBROUTINE cross_check

 !  IF (ltestcase) THEN

 !    SELECT CASE (TRIM(ctest_name))
 !    CASE('APE')
 !      CALL message('','Running hydrostatic atm model with ECHAM6 physics'// &
 !                  &' in aqua-planet mode.')

 !      llandsurf = .FALSE.
 !      lssodrag  = .FALSE.
 !      lice      = .FALSE.
 !      lmeltpond = .FALSE.
 !      lmlo      = .FALSE.
 !      lhd       = .FALSE.

 !    CASE('JWw-Moist')
 !      CALL message('','Running the Jablonowski-Williamson baroclinic test'// &
 !                  &' with the hydrostatic atm dynamical core'//              &
 !                  &' and ECHAM6 physics.'  )

 !      llandsurf = .FALSE.
 !      lssodrag  = .FALSE.
 !      lice      = .FALSE.
 !      lmeltpond = .FALSE.
 !      lmlo      = .FALSE.
 !      lhd       = .FALSE.

 !    CASE('LDF-Moist')
 !      CALL message('','Running the local diabatic forcing test'// &
 !                  &' with the hydrostatic atm dynamical core'//              &
 !                  &' and ECHAM6 physics.'  )

 !      llandsurf = .FALSE.
 !      lssodrag  = .FALSE.
 !      lice      = .FALSE.
 !      lmeltpond = .FALSE.
 !      lmlo      = .FALSE.
 !      lhd       = .FALSE.

 !    CASE DEFAULT
 !      CALL finish(TRIM(thismodule),'Invalid test case with ECHAM6 physics')
 !    END SELECT

 !  ENDIF

 !END SUBROUTINE cross_check
  !-------------

END MODULE mo_echam_phy_setup
