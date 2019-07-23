!! @par Revision History
!! Initial revision by Klaus Stephan, DWD (2015-02-05)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_radar_data_types
  USE mo_kind,               ONLY: wp
  USE mo_linked_list,          ONLY: t_var_list
  USE mtime,                   ONLY: datetime

IMPLICIT NONE

PUBLIC t_radar_fields , t_lhn_diag, t_radar_td_fields, t_radar_ct_fields

TYPE t_radar_td_fields

  REAL(wp), POINTER      &
#ifdef _CRAYFTN
      , CONTIGUOUS             &
#endif
    & ::               &
    & obs(:,:,:)      ,& ! observations on model grid at six observation time levels
    & spqual(:,:,:)      ! spatial quality function on model grid at two observation time levels

  TYPE(datetime),POINTER ::  obs_date(:) ! reference date of observations

END TYPE t_radar_td_fields

TYPE t_radar_ct_fields

  REAL(wp), POINTER      &
#ifdef _CRAYFTN
      , CONTIGUOUS             &
#endif
    & ::               &
    & blacklist(:,:)  ,& ! blacklist for DX radar data
    & brightband(:,:) ,& ! bright band mask field
    & dxheight(:,:,:)    ! DX radar heights

END TYPE t_radar_ct_fields

TYPE t_radar_fields

  TYPE (t_radar_td_fields) :: radar_td
  TYPE (t_var_list) :: radar_td_list

  TYPE (t_radar_ct_fields) :: radar_ct
  TYPE (t_var_list) :: radar_ct_list

END TYPE t_radar_fields

TYPE t_lhn_diag
!
  REAL(wp), ALLOCATABLE      &
!#ifdef _CRAYFTN
!      , CONTIGUOUS             &
!#endif
    & ::                    &
    & ttend_lhn(:,:,:)       ,& ! temperature increment due to LHN
    & qvtend_lhn(:,:,:)      ,& ! moisture increment due to LHN
    & pr_obs_sum(:,:)       ,& ! cumulated precipitation (hourly)
    & pr_mod_sum(:,:)       ,& ! cumulated precipitation (hourly)
    & pr_ref_sum(:,:)          ! cumulated precipitation (hourly)

END TYPE t_lhn_diag

!===============================================================================

END MODULE mo_radar_data_types


