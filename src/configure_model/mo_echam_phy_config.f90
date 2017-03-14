!>
!! Configuration of the ECHAM physics package. Includes main switches
!! for turning on/off parameterized processes.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_phy_config

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message, print_value

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_echam_phy_config, echam_phy_config   !< derived type and variable
  PUBLIC :: configure_echam_phy                    !< subroutine
  PUBLIC :: get_lrad, get_dt_rad, get_lvdiff,    & !< functions to retrieve values
    &       get_lconv, get_lcond,                & !<   of single parameters of the
    &       get_lgw_hines, get_lssodrag,         & !<   whole echam6 configuration
    &       get_lcariolle,                       & !<   state without 
    &       get_lmlo, get_lice, get_ljsbach,     & !<   USEing the
    &       get_lamip, get_lebudget                !<   whole state.

  !>
  !! Derived type containing main switches for configuring the echam physics package
  !!
  TYPE t_echam_phy_config

    INTEGER  :: idcphycpl   !<  determines the coupling between the dynamical core and the
                            !   echam phyiscs package
                            !   1: dynamics and physics update sequentially
                            !   2: dynamics uses physics forcing for updating
    LOGICAL  :: ldrymoist   !   .true.  physics assumes moist air dynamics
    LOGICAL  :: lrad        !<  .true. for radiation.
    REAL(wp) :: dt_rad      !<  [s] radiation time step
    LOGICAL  :: lvdiff      !<  .true. for vertical diffusion.
    LOGICAL  :: lconv       !<  .true. for moist convection
    LOGICAL  :: lcond       !<  .true. for large scale condensation
    LOGICAL  :: lgw_hines   !<  .true. for atmospheric gravity wave drag
    LOGICAL  :: lssodrag    !<  .true. for subgrid scale orographic drag,
                            !<         by blocking and gravity waves (lgwdrag in ECHAM6)
    LOGICAL  :: lcariolle   !<  .true. for Cariolle interactive ozone scheme
    LOGICAL  :: lmlo        !<  .true. for mixed layer ocean
    LOGICAL  :: lice        !<  .true. for sea-ice temperature calculation
    LOGICAL  :: ljsbach     !<  .true. for calculating the JSBACH land surface

    LOGICAL  :: lamip       !<  .true. for AMIP simulations with monthly transient boundary conditions   
    LOGICAL  :: lebudget    !<  .true. for echam physics energy budget calculation

  END TYPE t_echam_phy_config

  !>
  !! The configuration state (variable).
  !! So far we have not tried to use different configurations for different
  !! domains (grid levels) in experiments with nesting. Thus the variable
  !! is declared as a scalar. Later it might be changed into an array of
  !! shape (/n_dom/).
  !!
  TYPE(t_echam_phy_config)    :: echam_phy_config

CONTAINS
  !>
  !!
  SUBROUTINE configure_echam_phy

    CHARACTER(LEN=*),PARAMETER  :: method_name ='mo_echam_phy_config:configure_echam_phy'

    !------------------

#ifdef __NO_JSBACH__
    IF ( echam_phy_config% ljsbach ) THEN
      CALL finish(method_name, "This version was compiled without jsbach. Compile with __JSBACH__ to run this experiment")
    END IF
#endif

    CALL message('','')
    CALL message(method_name,'dynamics physics coupling:')
    CALL print_value('    idcphycpl  ',echam_phy_config% idcphycpl)
    CALL print_value('    ldrymoist  ',echam_phy_config% ldrymoist)
    CALL message('','')
    CALL message(method_name,'ECHAM6 physics configuration:')
    CALL print_value('    lrad       ',echam_phy_config% lrad     )
    CALL print_value('    dt_rad     ',echam_phy_config% dt_rad   )
    CALL print_value('    lvdiff     ',echam_phy_config% lvdiff   )
    CALL print_value('    lconv      ',echam_phy_config% lconv    )
    CALL print_value('    lcond      ',echam_phy_config% lcond    )
    CALL print_value('    lgw_hines  ',echam_phy_config% lgw_hines)
    CALL print_value('    lssodrag   ',echam_phy_config% lssodrag )
    CALL print_value('    lcariolle  ',echam_phy_config% lcariolle)
    CALL print_value('    lmlo       ',echam_phy_config% lmlo     )
    CALL print_value('    lice       ',echam_phy_config% lice     )
    CALL print_value('    ljsbach    ',echam_phy_config% ljsbach  )
    CALL message('','')
    CALL message(method_name,'ECHAM6 physics boundary conditions:')
    CALL print_value('    lamip      ',echam_phy_config% lamip    )
    CALL message('','')
    CALL message(method_name,'ECHAM6 physics energy budget control:')
    CALL print_value('    lebudget   ',echam_phy_config% lebudget )
    CALL message('','')


  END SUBROUTINE configure_echam_phy
  !------------------------------
  !>
  !!
  INTEGER FUNCTION get_idcphycpl()
    get_idcphycpl = echam_phy_config%idcphycpl
  END FUNCTION get_idcphycpl
  !>
  !!
  LOGICAL FUNCTION get_ldrymoist()
    get_ldrymoist = echam_phy_config%ldrymoist
  END FUNCTION get_ldrymoist
  !>
  !!
  LOGICAL FUNCTION get_lrad()
    get_lrad = echam_phy_config%lrad
  END FUNCTION get_lrad
  !>
  !!
  REAL(wp) FUNCTION get_dt_rad()
    get_dt_rad = echam_phy_config%dt_rad
  END FUNCTION get_dt_rad
  !>
  !!
  LOGICAL FUNCTION get_lvdiff()
    get_lvdiff = echam_phy_config%lvdiff
  END FUNCTION get_lvdiff
  !>
  !!
  LOGICAL FUNCTION get_lconv()
    get_lconv = echam_phy_config%lconv
  END FUNCTION get_lconv
  !>
  !!
  LOGICAL FUNCTION get_lcond()
    get_lcond = echam_phy_config%lcond
  END FUNCTION get_lcond
  !>
  !!
  LOGICAL FUNCTION get_lgw_hines()
    get_lgw_hines = echam_phy_config%lgw_hines
  END FUNCTION get_lgw_hines
  !>
  !!
  LOGICAL FUNCTION get_lssodrag()
    get_lssodrag = echam_phy_config%lssodrag
  END FUNCTION get_lssodrag
  !>
  !!
  LOGICAL FUNCTION get_lcariolle()
    get_lcariolle = echam_phy_config%lcariolle
  END FUNCTION get_lcariolle
  !>
  !!
  LOGICAL FUNCTION get_lmlo()
    get_lmlo = echam_phy_config%lmlo
  END FUNCTION get_lmlo
  !>
  !!
  LOGICAL FUNCTION get_lice()
    get_lice = echam_phy_config%lice
  END FUNCTION get_lice
  !>
  !!
  LOGICAL FUNCTION get_ljsbach()
    get_ljsbach = echam_phy_config%ljsbach
  END FUNCTION get_ljsbach
  !>
  !!
  LOGICAL FUNCTION get_lamip()
    get_lamip = echam_phy_config%lamip
  END FUNCTION get_lamip
  !>
  !!
  LOGICAL FUNCTION get_lebudget()
    get_lebudget = echam_phy_config%lebudget
  END FUNCTION get_lebudget

END MODULE mo_echam_phy_config
