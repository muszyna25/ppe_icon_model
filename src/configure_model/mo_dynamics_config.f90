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
  USE mo_exception,             ONLY: finish, warning
  USE mo_impl_constants,        ONLY: MAX_DOM
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_master_control,        ONLY: is_restart_run
  USE mo_util_string,           ONLY: int2string

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iequations, idiv_method, divavg_cntrwgt
  PUBLIC :: sw_ref_height, lcoriolis, lshallow_water, ltwotime
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
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_config:setup_dynamics_config'

    !------------------------
    ! Set time level indices

    IF (is_restart_run()) THEN
      ! Read time level indices from restart file.
      ! NOTE: this part will be modified later for a proper handling
      ! of multiple domains!!!

      DO jdom = 1,ndom
        CALL get_restart_attribute( 'nold_DOM'//TRIM(int2string(jdom, "(i2.2)"))    ,nold    (jdom) )
        CALL get_restart_attribute( 'nnow_DOM'//TRIM(int2string(jdom, "(i2.2)"))    ,nnow    (jdom) )
        CALL get_restart_attribute( 'nnew_DOM'//TRIM(int2string(jdom, "(i2.2)"))    ,nnew    (jdom) )
        CALL get_restart_attribute( 'nnow_rcf_DOM'//TRIM(int2string(jdom, "(i2.2)")),nnow_rcf(jdom) )
        CALL get_restart_attribute( 'nnew_rcf_DOM'//TRIM(int2string(jdom, "(i2.2)")),nnew_rcf(jdom) )
      END DO

    ELSE ! not is_restart_run

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
