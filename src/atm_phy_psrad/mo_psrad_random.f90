!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.

MODULE mo_psrad_random
  USE mo_psrad_general, ONLY: dp, i8
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER:: big_endian = (transfer(1_i8, 1) == 0)

  INTEGER, PARAMETER :: seed_size = 4
  
  PUBLIC :: get_random, seed_size, get_random_rank3
  INTERFACE get_random_rank3
    MODULE PROCEDURE kissvec_mask_rank3, kissvec_all_rank3
  END INTERFACE get_random_rank3
  INTERFACE get_random
    MODULE PROCEDURE kissvec_masked, kissvec_all 
  END INTERFACE get_random
  
CONTAINS

#define m(k,n) (ieor (k, ishft (k, n)))
#ifdef BIG_ENDIAN
#define low_byte(i) transfer(ishft(i,bit_size(1)),1)
#else
#define low_byte(i) transfer(i,1)
#endif

! George Masaglia's KISS random number generator
! adapted from public domain code available at http://www.fortran.com/kiss.f90

  SUBROUTINE kissvec_all(kproma, kbdim, seed, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM) 
    
    INTEGER(i8) :: kiss(kproma) 
    INTEGER :: jk 
    
    DO jk = 1, kproma
      seed(jk,1) = low_byte(69069_i8 * seed(jk,1) + 1327217885)
      seed(jk,2) = m (seed(jk,2), 13)
      seed(jk,2) = m (seed(jk,2), - 17)
      seed(jk,2) = m (seed(jk,2), 5)
      seed(jk,3) = 18000 * iand (seed(jk,3), 65535) + ishft (seed(jk,3), - 16)
      seed(jk,4) = 30903 * iand (seed(jk,4), 65535) + ishft (seed(jk,4), - 16)
      kiss(jk) = int(seed(jk,1), i8) + seed(jk,2) + ishft (seed(jk,3), 16) + seed(jk,4)
      harvest(jk) = low_byte(kiss(jk))*2.328306e-10_dp + 0.5_dp
    END DO
  END SUBROUTINE kissvec_all

  SUBROUTINE kissvec_all_rank3(kproma, kbdim, klev, ksamps, seed, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, ksamps
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM,klev,ksamps) 
    
    INTEGER(i8) :: kiss(kproma) 
    INTEGER :: i,j,k
    
    DO k = 1, ksamps
    DO j = klev, 1, -1
    DO i = 1, kproma
      seed(i,1) = low_byte(69069_i8 * seed(i,1) + 1327217885)
      seed(i,2) = m (seed(i,2), 13)
      seed(i,2) = m (seed(i,2), - 17)
      seed(i,2) = m (seed(i,2), 5)
      seed(i,3) = 18000 * iand (seed(i,3), 65535) + ishft (seed(i,3), - 16)
      seed(i,4) = 30903 * iand (seed(i,4), 65535) + ishft (seed(i,4), - 16)
      kiss(i) = int(seed(i,1), i8) + seed(i,2) + ishft (seed(i,3), 16) + seed(i,4)
      harvest(i,j,k) = low_byte(kiss(i))*2.328306e-10_dp + 0.5_dp
    END DO
    END DO
    END DO
  END SUBROUTINE kissvec_all_rank3

  SUBROUTINE kissvec_mask_rank3(kproma, kbdim, klev, ksamps, seed, &
    mask, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, ksamps
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    LOGICAL,  INTENT(IN) :: mask(KBDIM,klev)    
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM,klev,ksamps) 
    
    INTEGER(i8) :: kiss(kproma) 
    INTEGER :: i,j,k
    
    DO k = 1, ksamps
    DO j = klev, 1, -1
    DO i = 1, kproma
      IF (mask(i,j)) THEN
        seed(i,1) = low_byte(69069_i8 * seed(i,1) + 1327217885)
        seed(i,2) = m (seed(i,2), 13)
        seed(i,2) = m (seed(i,2), - 17)
        seed(i,2) = m (seed(i,2), 5)
        seed(i,3) = 18000 * iand (seed(i,3), 65535) + ishft (seed(i,3), - 16)
        seed(i,4) = 30903 * iand (seed(i,4), 65535) + ishft (seed(i,4), - 16)
        kiss(i) = int(seed(i,1), i8) + seed(i,2) + ishft (seed(i,3), 16) + seed(i,4)
        harvest(i,j,k) = low_byte(kiss(i))*2.328306e-10_dp + 0.5_dp
      ENDIF
    END DO
    END DO
    END DO
  END SUBROUTINE kissvec_mask_rank3


!NOTE: not invoked in runtime, see mo_psrad_cld_sampling.f90
  SUBROUTINE kissvec_masked(kproma, kbdim, seed, mask, harvest)
    INTEGER, INTENT(IN) :: kproma, kbdim
    INTEGER, INTENT(INOUT) :: seed(:,:) ! Dimension kbdim, seed_size or bigger
    LOGICAL,  INTENT(IN) :: mask(KBDIM)    
    REAL(DP), INTENT(INOUT) :: harvest(KBDIM) 
    
    INTEGER(i8) :: kiss(kproma) 
    INTEGER     :: jk 
    
    DO jk = 1, kproma
      IF(mask(jk)) THEN  
        seed(jk,1) = low_byte(69069_i8 * seed(jk,1) + 1327217885)
        seed(jk,2) = m (seed(jk,2), 13)
        seed(jk,2) = m (seed(jk,2), - 17)
        seed(jk,2) = m (seed(jk,2), 5)
        seed(jk,3) = 18000 * iand (seed(jk,3), 65535) + ishft (seed(jk,3), - 16)
        seed(jk,4) = 30903 * iand (seed(jk,4), 65535) + ishft (seed(jk,4), - 16)
        kiss(jk) = int(seed(jk,1), i8) + seed(jk,2) + ishft (seed(jk,3), 16) + seed(jk,4)
        harvest(jk) = low_byte(kiss(jk))*2.328306e-10_dp + 0.5_dp
      ELSE  
        harvest(jk) = 0._dp
      END IF 
    END DO
  END SUBROUTINE kissvec_masked
  
END MODULE mo_psrad_random
