!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_atm_dyn_config

  USE mo_impl_constants,     ONLY: max_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_atm_dyn_config, atm_dyn_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for atm dynamics
  !!--------------------------------------------------------------------------
  TYPE :: t_atm_dyn_config

    ! namelist variables

    INTEGER :: iequations      !< Choice of governing equation set
    INTEGER :: itime_scheme    !< Choice of time stepping scheme
    INTEGER :: idiv_method     !< Divergence operator
    INTEGER :: divavg_cntrwgt  !< Weight of central cell for divergence averaging
    LOGICAL :: ldry_dycore     !< if .TRUE., ignore the effact of water vapor,
                               !< cloud liquid and cloud ice on virtual temperature.

    ! derived variables

    LOGICAL :: ltwotime
    INTEGER,ALLOCATABLE :: nold(:)   !< variables denoting time levels
    INTEGER,ALLOCATABLE :: nnow(:)   !< variables denoting time levels
    INTEGER,ALLOCATABLE :: nnew(:)   !< variables denoting time levels

    INTEGER,ALLOCATABLE :: nsav1(:)  !< Extra 'time levels' of prognostic variables
    INTEGER,ALLOCATABLE :: nsav2(:)  !< needed to compute boundary tendencies and
                                     !< feedback increments

    INTEGER,ALLOCATABLE :: nnow_rcf(:)  !< Extra time levels for reduced
    INTEGER,ALLOCATABLE :: nnew_rcf(:)  !< calling frequency (rcf)

  END TYPE t_atm_dyn_config
  !>
  !!
  TYPE(t_atm_dyn_config) :: atm_dyn_config(MAX_DOM)

END MODULE mo_atm_dyn_config
