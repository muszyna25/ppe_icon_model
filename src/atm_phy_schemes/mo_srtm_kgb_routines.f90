!>
!! This module holds the subroutines srtm_kgbXY where XY stands for the sw bands from
!!   16 to 29. The modularization is necessary to avoid naming clashes with the externally linked
!!   ecRad radiation library.
!!
!! @author Daniel Rieger, Deutscher Wetterdienst, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2018-11-29)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_srtm_kgb_routines

! Dummy use statements are needed. Otherwise the dependency tree is not created correctly
!   for the use statements inside the #include-d files. On the other hand, the use statements
!   from the #include-d files can not be put on the module level, as they contain 
!   conflicting variable names. The ONLY statements are not strictly necessary, but
!   they avoid the usage of too many unnecessary variables. (So-called FAKE-USE)
  USE mo_kind,      ONLY: wp
  USE mo_yoesrta16, ONLY: dummy_ka16=>ka
  USE mo_yoesrta17, ONLY: dummy_ka17=>ka
  USE mo_yoesrta18, ONLY: dummy_ka18=>ka
  USE mo_yoesrta19, ONLY: dummy_ka19=>ka
  USE mo_yoesrta20, ONLY: dummy_ka20=>ka
  USE mo_yoesrta21, ONLY: dummy_ka21=>ka
  USE mo_yoesrta22, ONLY: dummy_ka22=>ka
  USE mo_yoesrta23, ONLY: dummy_ka23=>ka
  USE mo_yoesrta24, ONLY: dummy_ka24=>ka
  USE mo_yoesrta25, ONLY: dummy_ka25=>ka
  USE mo_yoesrta26, ONLY: dummy_rayl=>rayl
  USE mo_yoesrta27, ONLY: dummy_ka27=>ka
  USE mo_yoesrta28, ONLY: dummy_ka28=>ka
  USE mo_yoesrta29, ONLY: dummy_ka29=>ka

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: srtm_kgb16, srtm_kgb17, srtm_kgb18
  PUBLIC :: srtm_kgb19, srtm_kgb20, srtm_kgb21
  PUBLIC :: srtm_kgb22, srtm_kgb23, srtm_kgb24
  PUBLIC :: srtm_kgb25, srtm_kgb26, srtm_kgb27
  PUBLIC :: srtm_kgb28, srtm_kgb29

CONTAINS


#include "srtm_kgb16.inc"
#include "srtm_kgb17.inc"
#include "srtm_kgb18.inc"
#include "srtm_kgb19.inc"
#include "srtm_kgb20.inc"
#include "srtm_kgb21.inc"
#include "srtm_kgb22.inc"
#include "srtm_kgb23.inc"
#include "srtm_kgb24.inc"
#include "srtm_kgb25.inc"
#include "srtm_kgb26.inc"
#include "srtm_kgb27.inc"
#include "srtm_kgb28.inc"
#include "srtm_kgb29.inc"

END MODULE mo_srtm_kgb_routines
