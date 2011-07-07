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

  IMPLICIT NONE
  PUBLIC
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
    LOGICAL :: llandsurf   !<  .true. for surface exchanges. (lsurf in ECHAM6)
    LOGICAL :: lssodrag    !<  .true. for subgrid scale orographic drag,
                           !<   by blocking and gravity waves (lgwdrag in ECHAM6)
    LOGICAL :: lgw_hines   !<  .true. for atmospheric gravity wave drag
    LOGICAL :: lice        !<  .true. for sea-ice temperature calculation
    LOGICAL :: lmeltpond   !<  .true. for calculation of meltponds
    LOGICAL :: lmlo        !<  .true. for mixed layer ocean
    LOGICAL :: lhd         !<  .true. for hydrologic discharge model
    LOGICAL :: lmidatm     !<  .true. for middle atmosphere model version

  END TYPE t_echam_phy_config

  !>
  !! The configuration state (variable).
  !! So far we have not tried to use different configurations for different
  !! domains (grid levels) in experiments with nesting. Thus the variable
  !! is declared as a scalar. Later it might be changed into an array of
  !! shape (n_dom).
  !!
  TYPE(t_echam_phy_config) :: echam_phy_config

CONTAINS
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

END MODULE mo_echam_phy_config
