!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_dynamics_config

  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: MAX_DOM
  USE mo_restart_nml_and_att,   ONLY: getAttributesForRestarting
  USE mo_key_value_store,       ONLY: t_key_value_store

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iequations, idiv_method, divavg_cntrwgt
  PUBLIC :: sw_ref_height, lcoriolis, lshallow_water, ltwotime
  PUBLIC :: ldeepatmo
  PUBLIC :: nold, nnow, nnew, nsav1, nsav2, nnow_rcf, nnew_rcf
  PUBLIC :: configure_dynamics

  !--------------------------------------------------------------------------
  ! Basic settings for the dynamical core 
  !--------------------------------------------------------------------------
  !TYPE :: t_dynamics_config

    ! namelist variables

    INTEGER  :: iequations      !< Choice of governing equation set
    INTEGER  :: idiv_method     !< Divergence operator
    REAL(wp) :: divavg_cntrwgt  !< Weight of central cell for divergence averaging
    REAL(wp) :: sw_ref_height   !< reference height to linearize around if using
                                !< lshallow_water and semi-implicit correction
    LOGICAL  :: lcoriolis       !< if .TRUE., Coriolis force is switched on   
    LOGICAL  :: ldeepatmo       !< if .TRUE., dynamical core assumes a deep atmosphere
                                !< instead of a shallow atmosphere

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
  SUBROUTINE configure_dynamics( ndom )

    INTEGER,INTENT(IN) :: ndom

    INTEGER :: jdom
    CHARACTER(2) :: sdom
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_config:setup_dynamics_config'

    !------------------------
    ! Set time level indices

    CALL getAttributesForRestarting(restartAttributes)
    IF (restartAttributes%is_init) THEN
      ! Read time level indices from restart file.
      ! NOTE: this part will be modified later for a proper handling
      ! of multiple domains!!!

      DO jdom = 1,ndom
        WRITE (sdom, "(i2.2)") jdom
        CALL restartAttributes%get('nold_DOM'//sdom, nold(jdom))
        CALL restartAttributes%get('nnow_DOM'//sdom, nnow(jdom))
        CALL restartAttributes%get('nnew_DOM'//sdom, nnew(jdom))
        CALL restartAttributes%get('nnow_rcf_DOM'//sdom, nnow_rcf(jdom))
        CALL restartAttributes%get('nnew_rcf_DOM'//sdom, nnew_rcf(jdom))
      END DO

    ELSE ! not isRestart

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
