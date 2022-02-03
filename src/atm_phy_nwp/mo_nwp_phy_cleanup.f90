!>
!! Cleanup/Destruction wrapper for NWP physics suite
!!
!! This routine destructs all variable lists related to the NWP physics suite, 
!! and calls deallocation routines for specific parameterizations, 
!! if necessary.
!!
!! @author Daniel Reinert, DWD
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2021-04-08)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_phy_cleanup

  USE mo_nwp_phy_state,        ONLY: destruct_nwp_phy_state
  USE mo_nwp_lnd_state,        ONLY: destruct_nwp_lnd_state
  USE mo_nwp_reff_interface,   ONLY: reff_calc_dom
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_lnd_nwp_config,       ONLY: tile_list
  USE mo_grid_config,          ONLY: n_dom
  USE mo_aerosol_sources_types,ONLY: p_dust_source_const
  USE mo_aerosol_util,         ONLY: tegen_scal_factors

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: cleanup_nwp_phy


CONTAINS

  !>
  !! Wrapper routine for NWP physics cleanup
  !!
  !! Performs destruction of NWP-specific variable lists and calls 
  !! deallocation routines of individual parameterizations (if necessary).
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2021-04-08)
  !!
  SUBROUTINE cleanup_nwp_phy()

    ! local
    INTEGER :: jg
  !------------------------------------------------------

    ! destruct NWP physics and land variable lists
    !
    CALL destruct_nwp_phy_state()
    CALL destruct_nwp_lnd_state()

    ! destruct surface tile list
    CALL tile_list%destruct()

    DO jg = 1, n_dom
      IF ( atm_phy_nwp_config(jg)%icalc_reff > 0 ) CALL reff_calc_dom(jg)%destruct()
      !
      CALL atm_phy_nwp_config(jg)%finalize()
      
      IF ( iprog_aero > 0 ) CALL p_dust_source_const(jg)%finalize()
    ENDDO
    
    CALL tegen_scal_factors%finalize()

  END SUBROUTINE cleanup_nwp_phy

END MODULE mo_nwp_phy_cleanup

