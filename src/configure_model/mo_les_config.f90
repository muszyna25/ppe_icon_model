!>
!! Contains the setup of variables related to large eddy simulation setup
!!
!! @Anurag Dipankar, MPIM (2013-04)
!!
!!
!! @par Revision History
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
MODULE mo_les_config

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, print_value
  USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_les_config, les_config  !< derived type and variable

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic configuration setup for LES with or without TORUS grid
  !--------------------------------------------------------------------------
  TYPE t_les_config

    ! variables from namelist
    REAL(wp) :: sst        ! prescribed SST
    REAL(wp) :: shflx      ! prescribed sensible heat flux (W/m2)
    REAL(wp) :: lhflx      ! prescribed latent heat flux   (W/m2)
    INTEGER  :: isrfc_type ! 1=fixed sst, 2=fixed flux

    REAL(wp) :: ugeo(2)    ! ugeo(1)=constant, ugeo(2)=gradient
    REAL(wp) :: vgeo(2)    ! vgeo(1)=constant, vgeo(2)=gradient
    REAL(wp) :: umean(2)   ! umean(1)=constant, umean(2)=gradient
    REAL(wp) :: vmean(2)   ! vmean(1)=constant, vmean(2)=gradient
    REAL(wp) :: ufric      ! friction velocity
 
    LOGICAL  :: is_dry_cbl  !special case for CBL testcase
    LOGICAL  :: set_geowind !TRUE is geostrophic wind is set

    !Some parameters
    REAL(wp) :: karman_constant
    REAL(wp) :: rkarman_constant  !inverse karman constant
    REAL(wp) :: smag_constant
    REAL(wp) :: turb_prandtl 
    REAL(wp) :: rturb_prandtl     !inverse turbulent prandtl number

  END TYPE t_les_config
  !>
  !!
  TYPE(t_les_config),TARGET :: les_config(max_dom)


END MODULE mo_les_config
