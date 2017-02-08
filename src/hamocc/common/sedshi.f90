!>
!! @file bgc.f90
!! @brief Shifting of solid components
!!
!! By this routine solid components are shifted (upward and) downward
!! to account for sedimant gain and loss. This includes a layer for
!! permanent burial which collects the partical matter (P, Si, C, clay)
!! over the full time of integration.
!!
!! Upward shift is currently disabled.
!!
#include "hamocc_omp_definitions.inc"
SUBROUTINE sedshi(start_idx,end_idx)

  USE mo_kind, ONLY       : wp
  USE mo_sedmnt, ONLY     : sedlay, seddw, burial, orgfa, oplfa, &
       &                    calfa, clafa, porsol, solfu, ks
  USE mo_memory_bgc, ONLY : bolay,  rcar
  USE mo_param1_bgc, ONLY : nsedtra, &
       &                    issso12, isssc12, issssil, issster
  USE mo_control_bgc, ONLY: bgc_nproma
  USE mo_hamocc_nml, ONLY: l_up_sedshi
 
  IMPLICIT NONE

  !! Arguments

  INTEGER, INTENT(in)  :: start_idx       !< 1st REAL of model grid.
  INTEGER, INTENT(in)  :: end_idx       !< 2nd REAL of model grid.

  !! Local variables

  INTEGER  :: j,k,iv
  REAL(wp) :: wsed(bgc_nproma)        !< shifting velocity for upward/downward shifts
  REAL(wp) :: fulsed(bgc_nproma)
  REAL(wp) :: sedlo,uebers
  REAL(wp) :: seddef                 !< sediment deficiency
  REAL(wp) :: spresent, buried
  REAL(wp) :: refill,frac
  !
  !----------------------------------------------------------------------
  !
  ! DOWNWARD SHIFTING
  ! shift solid sediment downwards, if layer is full, i.e., if
  ! the volume filled by the four constituents poc, opal, caco3, clay
  ! is more than porsol*seddw
  ! the outflow of layer i is given by sedlay(i)*porsol(i)*seddw(i), it is
  ! distributed in the layer below over a volume of porsol(i+1)*seddw(i+1)
  if (start_idx==0)RETURN

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(k,j,iv,uebers,sedlo,wsed) HAMOCC_OMP_DEFAULT_SCHEDULE
  DO k = 1, ks-1

     DO j = start_idx, end_idx
        
           IF (bolay(j) > 0._wp) THEN
              sedlo  = orgfa*rcar*sedlay(j,k,issso12)    &
                   & +      calfa*sedlay(j,k,isssc12)    &
                   & +      oplfa*sedlay(j,k,issssil)    &
                   & +      clafa*sedlay(j,k,issster)
              ! "full" sediment has sedlo=1. for sedlo>1., wsed is >0.
              wsed( j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp)) ! downward shifting velocity (?)
           ENDIF

     ENDDO !end j-loop

     ! filling downward  (accumulation)
     DO iv = 1, nsedtra
        DO j = start_idx, end_idx

              IF (bolay(j) > 0._wp) THEN
                 uebers = wsed(j)*sedlay(j,k,iv)                     ! 'uebersaettigung?'
                 sedlay(j,k  ,iv) = sedlay(j,k  ,iv) - uebers
                 sedlay(j,k+1,iv) = sedlay(j,k+1,iv) + uebers        &
                      &             *(seddw(k)*porsol(k))/(seddw(k+1)*porsol(k+1))
              ENDIF

        ENDDO !end j-loop
     ENDDO !end iv-loop

  ENDDO !end k-loop
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  ! store amount lost from bottom sediment layer - this is a kind of
  ! permanent burial in deep consolidated layer, and this stuff is
  ! effectively lost from the whole ocean+sediment(+atmosphere) system.
  ! Would have to be supplied by river runoff or simple addition e.g.
  ! to surface layers in the long range. Can be supplied again if a
  ! sediment column has a deficiency in volume.

  DO j = start_idx, end_idx

        IF (bolay(j) > 0._wp) THEN
           sedlo  = orgfa*rcar*sedlay(j,ks,issso12)    &
                & +      calfa*sedlay(j,ks,isssc12)    &
                & +      oplfa*sedlay(j,ks,issssil)    &
                & +      clafa*sedlay(j,ks,issster)
           wsed( j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))
        ENDIF

  ENDDO !end j-loop

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,sedlo,iv,uebers) HAMOCC_OMP_DEFAULT_SCHEDULE
  DO iv = 1, nsedtra
     DO j = start_idx, end_idx

           IF (bolay(j) > 0._wp) THEN
              uebers = wsed(j)*sedlay(j,ks,iv)
              sedlay(j,ks,iv) = sedlay(j,ks ,iv)-uebers
              burial(j,iv)    = burial(j,iv) + uebers*seddw(ks)*porsol(ks)
           ENDIF

     ENDDO !end j-loop
  ENDDO !end iv-loop
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

 IF(l_up_sedshi)THEN 

  ! UPWARD SHIFTING
  ! shift solid sediment upwards, if total sediment volume is less
  ! than required, i.e., if the volume filled by the four constituents
  ! poc, opal, caco3, clay (integrated over the total sediment column)
  ! is less than porsol*seddw (integrated over the total sediment column)
  ! first, the deepest box is filled from below with total required volume;
  ! then, successively, the following layers are filled upwards.
  ! if there is not enough solid matter to fill the column, add clay. 

  fulsed(:) = 0._wp

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(k,j,sedlo) HAMOCC_OMP_DEFAULT_SCHEDULE
  ! determine how the total sediment column is filled
  DO k = 1, ks
     DO j = start_idx, end_idx
           IF (bolay(j) > 0._wp) THEN
              sedlo  = orgfa*rcar*sedlay(j,k,issso12)        &
                   & +      calfa*sedlay(j,k,isssc12)        &
                   & +      oplfa*sedlay(j,k,issssil)        &
                   & +      clafa*sedlay(j,k,issster)
              fulsed(j) = fulsed(j) + porsol(k)*seddw(k)*sedlo
           ENDIF
     ENDDO !end j-loop
  ENDDO !end k-loop
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

  ! shift the sediment deficiency from the deepest (burial)
  ! layer into layer ks
!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,seddef,spresent,buried,refill,frac) HAMOCC_OMP_DEFAULT_SCHEDULE
  DO j = start_idx, end_idx

        IF (bolay(j) > 0._wp) THEN

           ! deficiency with respect to fully loaded sediment |packed in sedlay(i,j,ks) ??
           ! this is the volume of sediment shifted upwards from the burial layer

           ! 'sediment deficiency', solfu = total column inegrated solid fraction volume (bodensed)
           seddef = solfu-fulsed(j)

           ! total volume of solid constituents in buried layer
           spresent = orgfa*rcar*burial(j,issso12)             &
                &   +      calfa*burial(j,isssc12)             &
                &   +      oplfa*burial(j,issssil)             &
                &   +      clafa*burial(j,issster)

           ! determine whether an additional amount of clay is needed from the burial
           ! layer to fill the whole sediment; I assume that there is an infinite
           ! supply of clay from below
           burial(j,issster) = burial(j,issster)             &
                &              + MAX(0._wp, seddef - spresent) / clafa

           ! determine new volume of buried layer
           buried = orgfa*rcar*burial(j,issso12)               &
                & +      calfa*burial(j,isssc12)               &
                & +      oplfa*burial(j,issssil)               &
                & +      clafa*burial(j,issster)

           ! fill the deepest active sediment layer
           refill=seddef/buried
           frac = porsol(ks)*seddw(ks) !changed k to ks, ik

           sedlay(j,ks,issso12) = sedlay(j,ks,issso12)       &
                &                 + refill*burial(j,issso12)/frac
           sedlay(j,ks,isssc12) = sedlay(j,ks,isssc12)       &
                &                 + refill*burial(j,isssc12)/frac
           sedlay(j,ks,issssil) = sedlay(j,ks,issssil)       &
                &                 + refill*burial(j,issssil)/frac
           sedlay(j,ks,issster) = sedlay(j,ks,issster)       &
                &                 + refill*burial(j,issster)/frac

           ! account for losses in buried sediment
           burial(j,issso12) = burial(j,issso12)             &
                &              - refill*burial(j,issso12)
           burial(j,isssc12) = burial(j,isssc12)             &
                &              - refill*burial(j,isssc12)
           burial(j,issssil) = burial(j,issssil)             &
                &              - refill*burial(j,issssil)
           burial(j,issster) = burial(j,issster)             &
                &              - refill*burial(j,issster)
        ENDIF ! bolay >0

  ENDDO !end j-loop
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL

!HAMOCC_OMP_PARALLEL
!HAMOCC_OMP_DO PRIVATE(j,k,sedlo,iv,uebers,frac) HAMOCC_OMP_DEFAULT_SCHEDULE
  !     redistribute overload of deepest layer ks to layers 2 to ks
  DO  k = ks, 2, -1
     DO j = start_idx, end_idx

           IF (bolay(j) > 0._wp) THEN
              sedlo  = orgfa*rcar*sedlay(j,k,issso12)          &
                   & +      calfa*sedlay(j,k,isssc12)          &
                   & +      oplfa*sedlay(j,k,issssil)          &
                   & +      clafa*sedlay(j,k,issster)
              wsed(j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))
           ENDIF

     ENDDO !end j-loop

     DO iv = 1, 4
        DO j = start_idx, end_idx
              IF (bolay(j) > 0._wp) THEN
                 uebers = sedlay(j,k,iv)*wsed(j)
                 frac   = porsol(k)*seddw(k)/(porsol(k-1)*seddw(k-1))
                 sedlay(j,k,iv)   = sedlay(j,k,iv)   - uebers
                 sedlay(j,k-1,iv) = sedlay(j,k-1,iv) + uebers*frac ! note k-1 here = upward shift
              ENDIF
        ENDDO !end j-loop
     ENDDO !end iv-loop
  ENDDO !end k-loop
!HAMOCC_OMP_END_DO
!HAMOCC_OMP_END_PARALLEL
 ENDIF ! l_up_sedshi

END SUBROUTINE 
