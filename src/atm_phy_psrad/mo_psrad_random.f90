!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide random numbers
!!
!! @remarks
!!   This implements scalar and vector thread-safe random number generators 
!!   using the KISS (Keep It Simple, Stupid) algorithm of George Marsiglia
!!
!! @author Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2013-08)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   This module is based on code developed at NCAR, US, and released as part of CESM5, 
!!   where it was used in the mcica_subcol_gen_[lw]w.f90 modules. 
!!   The original implementation in CESM is due to Phil Rasch and Mathew Rothstein (2004-10) 
!!   based on http://www.fortran.com/kiss.f90 on 2004-03-16
!!   This code incorporates updates by Sean Santos (2013-08) to avoid compiler complaints 
!!   about integer overflows (required as part of the algorithm) in a platform-independent manner 
!!    
!!   The original KISS algorithm is described in posts to comp.lang.fortran on 
!!     2007-06-23 (e.g. https://groups.google.com/forum/#!msg/comp.lang.fortran/5Bi8cFoYwPE/pSFU7NaK224J) 
!!   The  KISS (Keep It Simple Stupid) random number generator. Combines:
!!    (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
!!    (2) A 3-shift shift-register generator, period 2^32-1,
!!    (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!!    Overall period>2^123; 
!!
!
MODULE mo_psrad_random
  USE mo_psrad_general, ONLY : dp, i8
  IMPLICIT NONE

  LOGICAL, PARAMETER:: big_endian = (transfer(1_i8, 1) == 0)

  INTEGER, PARAMETER :: state_size = 4
  INTEGER :: global_seed(state_size)  = &
             (/123456789,362436069,21288629,14921776/) 
  
  PRIVATE
  PUBLIC :: get_random, get_i_random, seed_size_random, set_seed_random
  INTERFACE get_random
    MODULE PROCEDURE kisssca, kiss_global, kissvec, kissvec_legacy, &
      kissvec_all, kissvec_global
  END INTERFACE get_random
  
CONTAINS
  ! -----------------------------------------------
  FUNCTION seed_size_random()
    INTEGER :: seed_size_random
    
    seed_size_random = state_size
  END FUNCTION seed_size_random   
  ! -----------------------------------------------
  SUBROUTINE set_seed_random(seed) 
    INTEGER,  INTENT(in) :: seed(:)  ! Dimension seed_size
    
    INTEGER :: mx
    mx = MIN(size(seed), state_size) 
    global_seed(1:mx) = seed(1:mx) 
  END SUBROUTINE set_seed_random
  ! -----------------------------------------------
  SUBROUTINE kissvec_all(kproma, kbdim, seed, harvest)
    INTEGER,  INTENT(in   ) :: kproma, kbdim
    INTEGER,  INTENT(inout) :: seed(:,:)  ! Dimension nproma, seed_size
    REAL(DP), INTENT(  out) :: harvest(:) ! Dimension nproma
  
    INTEGER :: mask(KBDIM) 
    
    mask(:) = 1
    CALL kissvec(kproma, KBDIM, seed, mask, harvest)
    
  END SUBROUTINE kissvec_all
  ! -----------------------------------------------
#define m(k,n) (ieor (k, ishft (k, n)))
  ! -----------------------------------------------
#ifdef BIG_ENDIAN
#define low_byte(i) transfer(ishft(i,bit_size(1)),1)
#else
#define low_byte(i) transfer(i,1)
#endif
  ! -----------------------------------------------

  SUBROUTINE kissvec(kproma, kbdim, seed, mask, harvest)
    INTEGER,  INTENT(in   ) :: kproma, kbdim
    INTEGER,  INTENT(inout) :: seed(:,:)      ! Dimension kbdim, seed_size or bigger
    INTEGER,  INTENT(in   ) :: mask(KBDIM)    
    REAL(DP), INTENT(  out) :: harvest(KBDIM) 
    
    INTEGER(i8) :: kiss(kproma) 
    INTEGER     :: jk 
    
    DO jk = 1, kproma
      IF(mask(jk) == 1) THEN  
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
  END SUBROUTINE kissvec

  SUBROUTINE kissvec_legacy(kproma, kbdim, seed, mask, harvest)
    INTEGER,  INTENT(in   ) :: kproma, kbdim
    INTEGER,  INTENT(inout) :: seed(:,:)      ! Dimension kbdim, seed_size or bigger
    LOGICAL,  INTENT(in   ) :: mask(KBDIM)    
    REAL(DP), INTENT(  out) :: harvest(KBDIM) 
    
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
  END SUBROUTINE kissvec_legacy
  ! -----------------------------------------------
  SUBROUTINE get_i_random(kproma, kbdim, seed, a, b, harvest)
    integer, intent(in   ) :: kproma, kbdim, a, b
    integer, intent(inout) :: seed(:,:)
    integer, intent(  out) :: harvest(KBDIM) 
    
    integer(i8) :: kiss(kproma) 
    
    seed(:,1) = low_byte(69069_i8 * seed(:,1) + 1327217885)
    seed(:,2) = m (seed(:,2), 13)
    seed(:,2) = m (seed(:,2), - 17)
    seed(:,2) = m (seed(:,2), 5)
    seed(:,3) = 18000 * iand (seed(:,3), 65535) + ishft (seed(:,3), - 16)
    seed(:,4) = 30903 * iand (seed(:,4), 65535) + ishft (seed(:,4), - 16)
    kiss(:) = int(seed(:,1), i8) + seed(:,2) + ishft (seed(:,3), 16) + seed(:,4)
    harvest(:) = int(mod(kiss(:), int(b-a+1,i8))) + a
  END SUBROUTINE get_i_random
  ! -----------------------------------------------
  SUBROUTINE kisssca(seed, harvest)
    INTEGER,       INTENT(inout) :: seed(:)
    REAL(KIND=DP), INTENT(  out) :: harvest
    
    INTEGER(i8) :: kiss

    seed(1) = low_byte(69069_i8 * seed(1) + 1327217885)
    seed(2) = m (m (seed(2), 13), - 17)
    seed(2) = m (seed(2), 5)
    seed(3) = 18000 * iand (seed(3), 65535) + ishft (seed(3), - 16)
    seed(4) = 30903 * iand (seed(4), 65535) + ishft (seed(4), - 16)
    kiss = int(seed(1), i8) + seed(2) + ishft (seed(3), 16) + seed(4)
    harvest = low_byte(kiss)*2.328306e-10_dp + 0.5_dp
  END SUBROUTINE kisssca
  ! -----------------------------------------------
  SUBROUTINE kiss_global(harvest)
    REAL(DP), INTENT(inout) :: harvest
    
    CALL kisssca(global_seed, harvest)
  END SUBROUTINE kiss_global
  ! -----------------------------------------------
  SUBROUTINE kissvec_global(harvest)
    REAL(DP), INTENT(inout) :: harvest(:)
    
    INTEGER :: i 
    DO i = 1, SIZE(harvest) 
      CALL kisssca(global_seed, harvest(i))
    END DO
  END SUBROUTINE kissvec_global
  
END MODULE mo_psrad_random
  ! -----------------------------------------------
