!>
!! Type definition for preparing tracer advection in ICONAM
!!
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD (2013-09-18)
!! - Moved here from mo_nh_stepping to avoid circular dependencies
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
!!
MODULE mo_nh_prepadv_types

  USE mo_kind,                 ONLY: wp, vp

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: t_prepare_adv, t_step_adv, prep_adv, jstep_adv

  ! for preparation of transport with optional reduced calling frequency
  TYPE :: t_prepare_adv
    ! mass flux at full level edges (currently at N+1\2)
    REAL(wp), ALLOCATABLE :: mass_flx_me(:,:,:)

    !mass flux at half level centers (currently at N+1\2)
    REAL(wp), ALLOCATABLE :: mass_flx_ic(:,:,:)

    ! horizontal velocity at edges for computation of backward trajectories
    ! (currently at N+1/2)
    REAL(wp), ALLOCATABLE :: vn_traj(:,:,:)

    ! vertical velocity at half level centers for computation of
    ! backward trajectories (currently at N+1/2)
    REAL(wp), ALLOCATABLE :: w_traj(:,:,:)

    ! density times layer thickness at cell center at time step N
    REAL(wp), ALLOCATABLE :: rhodz_mc_now(:,:,:)

    !< density times layer thickness at cell center at time step N+1
    REAL(wp), ALLOCATABLE :: rhodz_mc_new(:,:,:)

    !< average of rhodz_mc_now and rhodz_mc_new used for calculating the vertical Courant number
    REAL(wp), ALLOCATABLE :: rhodz_mc_avg(:,:,:)

    !< vertical tracer flux at domain top (time average; n+1/2)
    REAL(wp), ALLOCATABLE :: topflx_tra(:,:,:)

  END TYPE t_prepare_adv


  ! counter for 'reduced calling frequency' and Marchuk splitting
  TYPE :: t_step_adv
    ! Counts total number of dynamics time steps for each patch (necessary for
    ! generalization of rcf to arbitrary (even) number of iadv_rcf)
    INTEGER :: ntsteps

    ! Determines sequence of operations for Marchuk-splitting (for transport)
    INTEGER :: marchuk_order
  END TYPE t_step_adv


  TYPE(t_prepare_adv), ALLOCATABLE :: prep_adv(:)  ! n_dom

  TYPE(t_step_adv),    ALLOCATABLE :: jstep_adv(:) ! n_dom


END MODULE mo_nh_prepadv_types

