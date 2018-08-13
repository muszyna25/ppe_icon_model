!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_wmo

  USE mo_kind                ,ONLY: wp

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_wmo

  USE mo_echam_cld_config    ,ONLY: echam_cld_config

  USE mo_tropopause          ,ONLY: WMO_tropopause

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_wmo

CONTAINS

  SUBROUTINE interface_echam_wmo(jg,jb,jcs,jce,  &
       &                         nproma,nlev   )

    INTEGER                 ,INTENT(in) :: jg                  !< grid  index
    INTEGER                 ,INTENT(in) :: jb                  !< block index
    INTEGER                 ,INTENT(in) :: jcs, jce            !< start/end column index within this block
    INTEGER                 ,INTENT(in) :: nproma,nlev

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    !
    INTEGER  :: itrpwmo(nproma), itrpwmop1(nproma)         !< for submodel - dummy not in use yet
    LOGICAL  :: lresum                                     !< for WMO_tropopause
    ! Shortcuts to components of echam_cld_config
    !  for WMO_tropopause
    INTEGER, POINTER :: ncctop, nccbot

    ncctop => echam_cld_config(jg)% ncctop
    nccbot => echam_cld_config(jg)% nccbot

    IF (ltimer) call timer_start(timer_wmo)

    ! associate pointers
    field => prm_field(jg)

    !
    lresum=.FALSE.
    !
    CALL WMO_tropopause( jcs, jce, nproma, nlev,   &! in
                       & ncctop, nccbot, lresum,   &! in
                       & field% ta(:,:,jb),        &! in
                       & field% presm_old(:,:,jb), &! in
                       & field% ptp(:,jb),         &! inout for diagnostics
                       & itrpwmo, itrpwmop1        )! out for submodel
    !

    IF (ltimer) call timer_stop(timer_wmo)

  END SUBROUTINE interface_echam_wmo
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_wmo
