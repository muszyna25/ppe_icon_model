!! mo_time_interpolation [module]
!!   routines for time interpolation of external data sets
!!
!! Authors;
!!   J.S. Rast, MPI February 2014     base version
!!  
!!----------------------------------------------------------------------
!!
!! This modules stores various time interpolation weights for 
!! interpolation to (i) the actual integration time step, and (ii)
!! the radiation time step (suffix _radt). The suffix _limm stands
!! for time weights and indices for linear interpolation 
!! of monthly means.
!!
!!----------------------------------------------------------------------
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
MODULE mo_time_interpolation_weights

  USE mo_kind,          ONLY: wp
  USE mo_datetime,      ONLY: t_datetime
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_wi_limm, wi_limm, wi_limm_radt

  TYPE t_wi_limm
    ! time and corresponding
    ! interpolation weights wgt1, wgt2 for monthly averages with
    ! inm1 and inm2 as indices [0,..,13], respectively
    TYPE(t_datetime) :: time
    LOGICAL          :: initialized=.FALSE.
    REAL(wp)         :: wgt1, wgt2
    INTEGER          :: inm1, inm2
  END TYPE t_wi_limm

  TYPE(t_wi_limm), SAVE   :: wi_limm, wi_limm_radt

  END MODULE mo_time_interpolation_weights
