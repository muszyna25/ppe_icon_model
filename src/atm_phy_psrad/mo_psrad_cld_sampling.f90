!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form. 
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide Monte Carlo samples based on cloud fraction
!!
!! @remarks
!!   This module draws Monte Carlo samples from profiles of the atmospheric
!!   state and on overlap assumption ( 1=max-ran, 2=maximum, 3=random). 
!!   Users provide profiles of cloud fraction and in-cloud optical thickness. 
!!   A single sample is drawn; cloud fraction is either 1 or 0 and optical 
!!   thickness is set to 0 where cloud fraction is 0. 
!!
!! @author Robert Pincus, U. Colorado, while visiting MPI-M, Hamburg (2010-08)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Written by Robert Pincus; simplifed from code originally written for the 
!!   GFDL atmospheric model AM2. 
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 

MODULE mo_psrad_cld_sampling 

  USE mo_psrad_general, ONLY : wp, finish
  USE mo_psrad_random, ONLY: get_random

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sample_cld_state

CONTAINS

  SUBROUTINE sample_cld_state(kproma, kbdim, klev, ksamps, &
    rnseeds, mode, fraction, is_cloudy)

    USE mo_psrad_general, ONLY : jTOA, jSFC, jINC

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: &
      kproma, kbdim, klev, &
      ksamps, & ! number of samples 
      mode !< 1=max-ran, 2=maximum, 3=random 
    INTEGER, INTENT(INOUT) :: rnseeds(:,:)
    REAL(wp), INTENT(IN) :: fraction(kbdim,klev) !cloud fraction
    LOGICAL, INTENT(INOUT) :: is_cloudy(kbdim,klev,ksamps)

    REAL(wp) :: rank(kbdim,klev,ksamps), one_minus(KBDIM,klev)
    INTEGER  :: jk, js 

    one_minus(1:kproma,:) = 1.0_wp - fraction(1:kproma,:)
    ! Here is_cloudy(:,:,1) indicates whether any cloud is present 
    is_cloudy(1:kproma,1:klev,1) = fraction(1:kproma,1:klev) > 0._wp
    SELECT CASE(mode) 
      ! Maximum-random overlap
      CASE(1) 
        DO js = 1, ksamps
          DO jk = jSFC, jTOA, jINC
            ! mask means we compute random numbers only when cloud is present 
            CALL get_random(kproma, kbdim, rnseeds, is_cloudy(:,jk,1), &
                rank(:,jk,js))
          END DO 
        END DO 
        ! There may be a better way to structure this calculation...
        DO js = 1, ksamps
          DO jk = jTOA-jINC, jSFC, -jINC
            rank(1:kproma,jk,js) = MERGE( &
              ! Max overlap:
              rank(1:kproma,jk+jINC,js), & 
              ! ... or random overlap in the clear sky portion:  
              rank(1:kproma,jk,js) * one_minus(1:kproma,jk+jINC), & 
              ! depending on whether or not you have cloud in the layer above 
              rank(1:kproma,jk+jINC,js) > one_minus(1:kproma,jk+jINC) ) 
          END DO
        END DO

      ! Max overlap means every cell in a column is identical 
      CASE(2) 
        DO js = 1, ksamps
          CALL get_random(kproma, kbdim, rnseeds, rank(:, 1, js))
          rank(1:kproma,2:klev,js) = SPREAD(rank(1:kproma,1,js), &
            DIM=2, NCOPIES=(klev-1))
        END DO 

      ! Random overlap means every cell is independent
      CASE(3) 
        DO js = 1, ksamps
          DO jk = jSFC, jTOA, jINC
            ! mask means we compute random numbers only when cloud is present 
            CALL get_random(kproma, kbdim, rnseeds, is_cloudy(:,jk,1), &
              rank(:,jk,js))
          END DO 
        END DO 
      CASE DEFAULT
        CALL finish('In sample_cld_state: unknown overlap assumption') 
    END SELECT
    ! Now is_cloudy indicates whether the sample (ks) is cloudy or not.    
    DO js = 1, ksamps
      is_cloudy(1:kproma,1:klev,js) = &
        rank(1:kproma,1:klev,js) > one_minus(1:kproma,1:klev)
    END DO 
  END SUBROUTINE sample_cld_state

END MODULE mo_psrad_cld_sampling
