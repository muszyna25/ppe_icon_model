!>
!! Implementation of physics utility routines, mostly domain
!! independent and elemental.
!!
!! @par Revision History
!!  Initial revision  :  F. Prill, DWD (2012-07-03)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_util_phys

  USE mo_kind,      ONLY: wp
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 cpd, p0ref, rd
  USE mo_satad,                 ONLY: sat_pres_water

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rel_hum

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

  !> POINTWISE computation of relative humidity as r=e/e_sat
  !!
  !! @par Revision History
  !! Initial revision  :  F. Prill, DWD (2012-07-03) 
  ELEMENTAL FUNCTION rel_hum(temp, qv, p_ex)
    REAL(wp) :: rel_hum
    REAL(wp), INTENT(IN) :: temp, &  ! temperature
      &                     qv,   &  ! spec. water vapor content
      &                     p_ex     ! exner pressure
    ! local variables
    REAL(wp) :: pres, e_s, e

    ! compute dynamic pressure from Exner pressure:
    pres = p0ref * EXP((cpd/rd)*LOG(p_ex))
    ! approx. saturation vapor pressure:
    e_s = sat_pres_water(temp)
    ! compute vapor pressure from formula for specific humidity:
    e   = pres*qv / (rdv + o_m_rdv*qv)

    rel_hum = 100._wp * e/e_s

  END FUNCTION rel_hum

END MODULE mo_util_phys
