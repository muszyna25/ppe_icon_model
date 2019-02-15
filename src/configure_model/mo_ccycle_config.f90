!>
!! Configuration of the Carbon cycle settings
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!! First version by Monika Esch (2017-11)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ccycle_config

  USE mo_exception     ,ONLY: message, print_value
  USE mo_impl_constants,ONLY: max_dom
  USE mo_grid_config   ,ONLY: n_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                  name       !< name for this unit

  PUBLIC ::         ccycle_config
  PUBLIC ::    init_ccycle_config       !< initialize ccycle_config
  PUBLIC ::   print_ccycle_config       !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'ccycle'

  !>
  !! Configuration type containing switches for the configuration of the carbon cycle
  !!
  TYPE t_ccycle_config
     !
     ! configuration parameters
     ! ------------------------
     !
     INTEGER  :: iccy_co2conc  !< Type of co2 concentration
     INTEGER  :: iccy_co2flux  !< Type of co2 flux in the atmosphere
     !
  END TYPE t_ccycle_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_ccycle_config) :: ccycle_config(max_dom)
  
CONTAINS

  !----
  !>
  !! Initialize the configuration
  !!
  SUBROUTINE init_ccycle_config

    ! Carbon cycle configuration
    ! --------------------------

    ccycle_config(:)% iccy_co2conc  = 0 ! 0: constant co2 concentration, e.g. 284 ppmv
                                        ! 1: co2 tracer
                                        ! 2: co2 concentration idealized scenario (computed)
                                        ! 4: co2 concentration scenario from file
    ccycle_config(:)% iccy_co2flux  = 0 ! 0: constant, e.g. 0
                                        ! 1: idealized (computed)
                                        ! 2: from ocean and land (needs HAMOCC and JSBACH)

  END SUBROUTINE init_ccycle_config

  !>
  !! Print out the user controlled configuration
  !!
  SUBROUTINE print_ccycle_config

    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg

    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','Carbon cycle configuration')
    CALL message    ('','==========================')
    CALL message    ('','')

    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    ccycle_config('//TRIM(cg)//')% iccy_co2conc ',ccycle_config(jg)% iccy_co2conc         )
       CALL print_value('    ccycle_config('//TRIM(cg)//')% iccy_co2flux ',ccycle_config(jg)% iccy_co2flux         )
       CALL message    ('','')
       !
    END DO
    CALL message    ('','')
    CALL message    ('','------------------------------------------------------------------------')
    CALL message    ('','')

  END SUBROUTINE print_ccycle_config

END MODULE mo_ccycle_config
