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

! USE mo_echam_phy_config,   ONLY: echam_phy_config, config_echam_phy
! USE mo_echam_conv_config,  ONLY: config_echam_convection
! USE mo_run_nml,            ONLY: ltestcase
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

   !IF (echam_phy_config%lconv) CALL config_echam_convection( nlev, vct_a, vct_b, ceta )
   !CALL config_echam_phy( ltestcase, ctest_name )

  END SUBROUTINE setup_echam_phy
  !-------------

END MODULE mo_echam_phy_setup
