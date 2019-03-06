!>
!! Some utilities which are specific to the transport algorithm.
!!
!! Module contains some functions and procedures which are specifically related
!! to the transport schemes. These subroutines or functions are needed at
!! various places within the transport scheme. Therefore outsourcing these
!! routines protects from possible circular dependencies.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-03-04)
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - implemented generalized Lax-Friedrich flux function
!!   laxfr_upflux_v, which allows to use the same transport
!!   code for pressure and height based vertical coordinate
!!   systems.
!! Modification by Daniel Reinert, DWD (2016-01-12)
!! - removed obsolete routine tupdate_tracer
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_utils

  USE mo_kind,                ONLY: wp

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: laxfr_upflux
  PUBLIC :: laxfr_upflux_v
  PUBLIC :: ptr_delp_mc_now
  PUBLIC :: ptr_delp_mc_new

  PUBLIC :: t_list2D

  
  TYPE t_list2D
    INTEGER, POINTER :: eidx(:,:)
    INTEGER, POINTER :: elev(:,:)
    INTEGER, POINTER :: len(:)
    INTEGER          :: npoints
  END TYPE t_list2D


  ! In order to avoid circular dependencies these two pointers
  ! have been moved from mo_advection_stepping to this module.
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_now(:,:,:) => NULL() !< pointer to old layer thickness
                                        !< at cell center
  REAL(wp), POINTER ::  &
    &  ptr_delp_mc_new(:,:,:) => NULL() !< pointer to new layer thickness
                                        !< at cell center

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux,.
  !!
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !!
  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
      &                   - ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to height based vertical coordinate systems. 
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized for p- and z-based vertical coordinate systems
  !! Modification by Daniel Reinert, DWD (2019-02-18)
  !! - revert generalization for p- and z-based vertical coordinate systems.
  !!   From now on only valid for z-based coordinate systems.
  !!
  FUNCTION laxfr_upflux_v( p_w, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_w
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_w  *( p_psi1 + p_psi2 )    &
      &                   + ABS( p_w )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux_v


END MODULE mo_advection_utils
