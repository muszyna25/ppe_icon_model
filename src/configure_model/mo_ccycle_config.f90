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

  USE mo_exception     ,ONLY: message, print_value, finish
  USE mo_kind          ,ONLY: wp
  USE mo_impl_constants,ONLY: max_dom
  USE mo_grid_config   ,ONLY: n_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                  name       !< name for this unit

  ! configuration
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
     INTEGER  :: iccycle   !< c-cycle mode
     INTEGER  :: ico2conc  !< co2 concentration provided to land and ocean
     !
     REAL(wp) :: vmr_co2   !< co2 volume mixing ratio for c-cycle
     !
  END TYPE t_ccycle_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_ccycle_config) :: ccycle_config(max_dom)
  
CONTAINS

  !----
  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_ccycle_config
    !
    ! Carbon cycle configuration
    ! --------------------------
    !
    ccycle_config(:)% iccycle  = 0            ! 0: no c-cycle
    !                                           1: c-cycle with interactive atm. co2 concentration
    !                                           2: c-cycle with prescribed  atm. co2 concentration
    !
    ! For iccycle = 2:
    ccycle_config(:)% ico2conc = 2            ! 2: constant  co2 concentration vmr_co2
    !                                           4: transient co2 concentration scenario from file
    !
    ! For ico2conc = 2:
    ccycle_config(:)% vmr_co2  = 284.3e-06_wp ! co2 volume mixing ratio of 1850 (CMIP6)
    !
  END SUBROUTINE init_ccycle_config

  !----

  !>
  !! Print out the user controlled configuration
  !!
  SUBROUTINE print_ccycle_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','Carbon cycle configuration')
    CALL message    ('','==========================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    ccycle_config('//TRIM(cg)//')% iccycle  ',ccycle_config(jg)% iccycle )
       CALL print_value('    ccycle_config('//TRIM(cg)//')% ico2conc ',ccycle_config(jg)% ico2conc)
       CALL print_value('    ccycle_config('//TRIM(cg)//')% vmr_co2  ',ccycle_config(jg)% vmr_co2 )
       CALL message    ('','')
       !
       SELECT CASE(ccycle_config(jg)% iccycle)
       CASE (0)
          CALL message ('','C-cycle is switched off')
       CASE(1)
          CALL message ('','C-cycle is used with interactive atmospheric CO2 concentration')
       CASE (2)
          CALL message ('','C-cycle is used with prescribed atmospheric CO2 concentration')
          !
          SELECT CASE(ccycle_config(jg)% ico2conc)
          CASE (2)
             CALL print_value('    CO2 volume mixing ratio is constant (ppv)',ccycle_config(jg)% vmr_co2)
          CASE (4)
             CALL message('','CO2 volume mixing ratio is read from scenario file bc_greenhouse_gases.nc')
          CASE default
             CALL finish('print_ccycle_config','ccycle_config(jg)% ico2conc invalid, must be 2 or 4')
          END SELECT
          !
       CASE default
          CALL finish('print_ccycle_config','ccycle_config(jg)% iccycle invalid, must be 0, 1 or 2')
       END SELECT
       !
    END DO
    !
  END SUBROUTINE print_ccycle_config

  !----

END MODULE mo_ccycle_config
