!>
!! Type definition for preparation of transport with optional reduced 
!! calling frequency.
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
MODULE mo_prepadv_types

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: t_prepare_adv
  PUBLIC :: t_step_adv


  ! for preparation of transport with optional reduced calling frequency
  !
  TYPE :: t_prepare_adv

    REAL(wp), POINTER        &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS             &
#endif
      ::                     &
      mass_flx_me(:,:,:) ,   & !< mass flux at full level edges   [kg/m^2/s]
                               !< averaged over dynamics substeps
      mass_flx_ic(:,:,:) ,   & !< mass flux at half level centers [kg/m^2/s]
                               !< averaged over dynamics substeps
      vn_traj    (:,:,:) ,   & !< horizontal velocity at edges for computation of backward trajectories [m/s]
                               !< averaged over dynamics substeps
      topflx_tra (:,:,:)       !< vertical tracer flux at domain top [kg/m^2/s]

  END TYPE t_prepare_adv


  ! counter for Marchuk splitting
  TYPE :: t_step_adv
    ! Determines sequence of operations for Marchuk-splitting (for transport)
    INTEGER :: marchuk_order
  END TYPE t_step_adv

END MODULE mo_prepadv_types

