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
MODULE mo_dynamics_config

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: finish
  USE mo_impl_constants, ONLY: MAX_DOM, LEAPFROG_EXPL, LEAPFROG_SI, ISHALLOW_WATER
  USE mo_io_restart_attributes, ONLY: get_restart_attribute

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iequations, itime_scheme, idiv_method, divavg_cntrwgt
  PUBLIC :: sw_ref_height, lcoriolis, lshallow_water, ltwotime
  PUBLIC :: nold, nnow, nnew, nsav1, nsav2, nnow_rcf, nnew_rcf
  PUBLIC :: configure_dynamics

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic settings for the dynamical core 
  !--------------------------------------------------------------------------
  !TYPE :: t_dynamics_config

    ! namelist variables

    INTEGER  :: iequations      !< Choice of governing equation set
    INTEGER  :: itime_scheme    !< Choice of time stepping scheme
    INTEGER  :: idiv_method     !< Divergence operator
    REAL(wp) :: divavg_cntrwgt  !< Weight of central cell for divergence averaging
    REAL(wp) :: sw_ref_height   !< reference height to linearize around if using
                                !< lshallow_water and semi-implicit correction
    LOGICAL  :: lcoriolis       !< if .TRUE., Coriolis force is switched on     

    ! derived variables

    LOGICAL :: lshallow_water
    LOGICAL :: ltwotime

    INTEGER :: nold(MAX_DOM)      !< variables denoting time levels
    INTEGER :: nnow(MAX_DOM)      !< variables denoting time levels
    INTEGER :: nnew(MAX_DOM)      !< variables denoting time levels

    INTEGER :: nsav1(MAX_DOM)     !< Extra 'time levels' of prognostic variables
    INTEGER :: nsav2(MAX_DOM)     !< needed to compute boundary tendencies and
                                  !< feedback increments

    INTEGER :: nnow_rcf(MAX_DOM)  !< Extra time levels for reduced
    INTEGER :: nnew_rcf(MAX_DOM)  !< calling frequency (rcf)

  !END TYPE t_dynamics_config
  !>
  !!
  !TYPE(t_dynamics_config) :: dynamics_config(MAX_DOM)

CONTAINS
  !>
  !!
  SUBROUTINE configure_dynamics( lrestart,ndom )

    LOGICAL,INTENT(IN) :: lrestart
    INTEGER,INTENT(IN) :: ndom

    INTEGER :: jdom
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_config:setup_dynamics_config'

    !------------------------
    ! Set time level indices

    IF (lrestart) THEN
      ! Read time level indices from restart file.
      ! NOTE: this part will be modified later for a proper handling
      ! of multiple domains!!!

      IF (ndom>1) CALL finish(TRIM(routine), &
      'Restart functionality can not handle multiple domains (yet)')

      jdom = 1  ! only consider one domain at the moment
      !DO jdom = 1,ndom
        CALL get_restart_attribute( 'nold'    ,nold    (jdom) )
        CALL get_restart_attribute( 'nnow'    ,nnow    (jdom) )
        CALL get_restart_attribute( 'nnew'    ,nnew    (jdom) )
        CALL get_restart_attribute( 'nnow_rcf',nnow_rcf(jdom) )
        CALL get_restart_attribute( 'nnew_rcf',nnew_rcf(jdom) )
      !END DO

    ELSE ! not lrestart

      nnow(:) = 1
      nnew(:) = 2
      nold(:) = 3
      nnow_rcf(:) = 1
      nnew_rcf(:) = 2

    END IF

    nsav1(:) = 0
    nsav2(:) = 4

  END SUBROUTINE configure_dynamics

END MODULE mo_dynamics_config
