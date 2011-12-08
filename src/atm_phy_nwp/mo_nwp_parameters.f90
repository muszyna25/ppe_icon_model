!>
!!  Declares parameters computed during the initialization of the physics
!!  parameterizations that have to be domain-dependent
!!
!! @par Revision History
!!  Guenther Zaengl, DWD, 2011-12-08
!!  - Restructuring the namelists
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
MODULE mo_nwp_parameters
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, finish
  USE mo_impl_constants,     ONLY: max_dom


  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  TYPE t_phy_params
    ! Level parameters for convection scheme
    INTEGER  :: kcon1, kcon2, kcon3, kcon4, kcon5
    ! resolution-dependent parameters for convection scheme
    REAL(wp) :: tau, mfcfl
    ! launch level for GWD scheme
    INTEGER  :: klaunch
  END TYPE t_phy_params

  TYPE (t_phy_params), ALLOCATABLE :: phy_params(:)

END MODULE mo_nwp_parameters
