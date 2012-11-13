!>
!! Configuration of the ECHAM physics package. Includes main switches
!! for turning on/off parameterized processes.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2010-07)
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
MODULE mo_echam_phy_config

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message, print_value

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_echam_phy_config, echam_phy_config   !< derived type and variable
  PUBLIC :: configure_echam_phy                       !< subroutine
  PUBLIC :: get_lrad, get_lcond, get_lcover        !< functions
  PUBLIC :: get_lconv, get_lvdiff, get_lgw_hines   !< functions
  PUBLIC :: get_ljsbach

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !>
  !! Derived type containing main swithes for configuring the echam physics package
  !!
  TYPE t_echam_phy_config

    LOGICAL :: lrad        !<  .true. for radiation.
    LOGICAL :: lvdiff      !<  .true. for vertical diffusion.
    LOGICAL :: lconv       !<  .true. for moist convection
    LOGICAL :: lcond       !<  .true. for large scale condensation
    LOGICAL :: lcover      !<  .true. for prognostic cloud cover scheme
    LOGICAL :: lgw_hines   !<  .true. for atmospheric gravity wave drag

    LOGICAL :: llandsurf   !<  .true. for surface exchanges. (lsurf in ECHAM6)
    LOGICAL :: lssodrag    !<  .true. for subgrid scale orographic drag,
                           !<   by blocking and gravity waves (lgwdrag in ECHAM6)
    LOGICAL :: lice        !<  .true. for sea-ice temperature calculation
    LOGICAL :: lmeltpond   !<  .true. for calculation of meltponds
    LOGICAL :: lmlo        !<  .true. for mixed layer ocean
    LOGICAL :: ljsbach     !<  .true. for calculating the JSBACH land surface
    LOGICAL :: lhd         !<  .true. for hydrologic discharge model
!!$    LOGICAL :: lmidatm     !<  .true. for middle atmosphere model version
    REAL(wp) :: dt_rad   !! "-"                     radiation
    

  END TYPE t_echam_phy_config

  !>
  !! The configuration state (variable).
  !! So far we have not tried to use different configurations for different
  !! domains (grid levels) in experiments with nesting. Thus the variable
  !! is declared as a scalar. Later it might be changed into an array of
  !! shape (/n_dom/).
  !!
  TYPE(t_echam_phy_config) :: echam_phy_config

CONTAINS
  !>
  !!
  SUBROUTINE configure_echam_phy( ltestcase, ctest_name )

    LOGICAL,         INTENT(IN) :: ltestcase
    CHARACTER(LEN=*),INTENT(IN) :: ctest_name
    CHARACTER(LEN=*),PARAMETER  :: &
             & routine ='mo_echam_phy_config:config_echam_phy'

    !------------------
    IF (ltestcase) THEN

      SELECT CASE (TRIM(ctest_name))
      CASE('APE','JWw-Moist','LDF-Moist')

        echam_phy_config% llandsurf = .FALSE.
        echam_phy_config% lssodrag  = .FALSE.
        echam_phy_config% lice      = .FALSE.
        echam_phy_config% lmeltpond = .FALSE.
        echam_phy_config% lmlo      = .FALSE.
        echam_phy_config% lhd       = .FALSE.

        CALL message('','')
        CALL message('','Running the hydrostatic atm model with ECHAM6 physics.')
        CALL message('','Testcase = '//TRIM(ctest_name))
        CALL message('','')

        CALL print_value('llandsurf ',echam_phy_config% llandsurf)
        CALL print_value('lssodrag  ',echam_phy_config% lssodrag )
        CALL print_value('lice      ',echam_phy_config% lice     )
        CALL print_value('lmeltpond ',echam_phy_config% lmeltpond)
        CALL print_value('lmlo      ',echam_phy_config% lmlo     )
        CALL print_value('lhd       ',echam_phy_config% lhd      )

        CALL message('','')

      CASE DEFAULT
        CALL message(TRIM(routine),'Testcase = '//TRIM(ctest_name))
        CALL finish (TRIM(routine),'Invalid test case with ECHAM6 physics')
      END SELECT
    ENDIF

  END SUBROUTINE configure_echam_phy
  !------------------------------
  !>
  !!
  LOGICAL FUNCTION get_lrad()
    get_lrad = echam_phy_config%lrad
  END FUNCTION get_lrad
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
  LOGICAL FUNCTION get_lcover()
    get_lcover = echam_phy_config%lcover
  END FUNCTION get_lcover
  !>
  !!
  LOGICAL FUNCTION get_lgw_hines()
    get_lgw_hines = echam_phy_config%lgw_hines
  END FUNCTION get_lgw_hines
  !>
  !!
  LOGICAL FUNCTION get_ljsbach()
    get_ljsbach = echam_phy_config%ljsbach
  END FUNCTION get_ljsbach

END MODULE mo_echam_phy_config
