!>
!! Type definition for preparing tracer advection in ICONAM
!!
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD (2013-09-18)
!! - Moved here from mo_nh_stepping to avoid circular dependencies
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_nh_prepadv_types

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: t_prepare_adv, t_step_adv, prep_adv, jstep_adv

  ! for preparation of transport with optional reduced calling frequency
  TYPE :: t_prepare_adv
#ifdef _OPENACC
    REAL(wp), POINTER ::     &
#else
    REAL(wp), ALLOCATABLE :: &
#endif
    ! mass flux at full level edges (currently at N+1\2)
    mass_flx_me(:,:,:) ,     &
    !
    !mass flux at half level centers (currently at N+1\2)
    mass_flx_ic(:,:,:) ,     &
    !
    ! horizontal velocity at edges for computation of backward trajectories
    ! (currently at N+1/2)
    vn_traj(:,:,:) ,         &
    !
    ! vertical velocity at half level centers for computation of
    ! backward trajectories (currently at N+1/2)
    w_traj(:,:,:) ,          &
    !
    !< vertical tracer flux at domain top (time average; n+1/2)
    topflx_tra(:,:,:)

  END TYPE t_prepare_adv


  ! counter for Marchuk splitting
  TYPE :: t_step_adv
    ! Determines sequence of operations for Marchuk-splitting (for transport)
    INTEGER :: marchuk_order
  END TYPE t_step_adv


  TYPE(t_prepare_adv), ALLOCATABLE :: prep_adv(:)  ! n_dom

  TYPE(t_step_adv),    ALLOCATABLE :: jstep_adv(:) ! n_dom


END MODULE mo_nh_prepadv_types

