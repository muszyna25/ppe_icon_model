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

  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_wmo

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

    IF (ltimer) call timer_start(timer_wmo)

    ! associate pointers
    field => prm_field(jg)

    CALL WMO_tropopause( jg,                       &! in
                       & jcs, jce, nproma, nlev,   &! in
                       & field% ta(:,:,jb),        &! in
                       & field% pfull(:,:,jb),     &! in
                       & field% ptp(:,jb)          )! inout for diagnostics
    !

    IF (ltimer) call timer_stop(timer_wmo)

  END SUBROUTINE interface_echam_wmo
  !-------------------------------------------------------------------

END MODULE mo_interface_echam_wmo
