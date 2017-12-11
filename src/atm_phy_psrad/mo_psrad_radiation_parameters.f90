!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!++mgs: new module 14.03.2010
!++mgs: added decl_sun_cur (for MOZ photolysis) 02.06.2010
!! @brief Module to provide parameters to radiation routines and avoid circular dependencies.
!!
!! @remarks
!!   This module contains the public parameters provided by the radiation module
!!   mo_radiation.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              extracted parameters from mo_radiation
!!         Dagmar Popke, MPI-M, Hamburg (2013-11-15):
!!              Implementation of RCE
!!
!! $ID: n/a$
!!
!! @par Origin
!!   see mo_radiation.f90
!!
!
MODULE mo_psrad_radiation_parameters

  USE mo_psrad_general, ONLY : wp

IMPLICIT NONE

  PRIVATE

  PUBLIC :: rad_perm, rad_undef
  PUBLIC :: i_overlap
  PUBLIC :: l_do_sep_clear_sky
  
  INTEGER :: rad_perm = 0                       ! Integer for perturbing random number seeds

  ! Radiation driver
  LOGICAL :: l_do_sep_clear_sky = .TRUE. ! Compute clear-sky fluxes by removing clouds
  INTEGER :: i_overlap = 1               ! 1 = max-ran, 2 = max, 3 = ran
  REAL, PARAMETER :: rad_undef = -999._wp

END MODULE mo_psrad_radiation_parameters
