!>
!! @file bgc.f90
!! @brief Shifting of solid components
!!
!! By this routine solid components are shifted (upward and) downward
!! to account for sedimant gain and loss. This includes a layer for
!! permanent local_sediment_mem%burial which collects the partical matter (P, Si, C, clay)
!! over the full time of integration.
!!
!! Upward shift is currently disabled.
!!
 

MODULE mo_sedshi

    USE mo_kind, ONLY       : wp
    USE mo_sedmnt, ONLY     : seddw, orgfa, oplfa, &
       &                    calfa, clafa, porsol, solfu
    USE mo_memory_bgc, ONLY : rcar
    USE mo_param1_bgc, ONLY : nsedtra, &
       &                    issso12, isssc12, issssil, issster
    USE mo_control_bgc, ONLY: bgc_nproma
    USE mo_hamocc_nml, ONLY: l_up_sedshi,ks
    
    USE mo_bgc_memory_types, ONLY: t_bgc_memory, t_sediment_memory
    

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: sedshi

CONTAINS

SUBROUTINE sedshi(local_bgc_mem, local_sediment_mem, start_idx,end_idx)

  
 
  IMPLICIT NONE

  !! Arguments
  TYPE(t_bgc_memory), POINTER :: local_bgc_mem
  TYPE(t_sediment_memory), POINTER :: local_sediment_mem

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
  ! the outflow of layer i is given by local_sediment_mem%sedlay(i)*porsol(i)*seddw(i), it is
  ! distributed in the layer below over a volume of porsol(i+1)*seddw(i+1)
  if (start_idx==0)RETURN

  DO k = 1, ks-1

     DO j = start_idx, end_idx
        
           IF (local_bgc_mem%bolay(j) > 0._wp) THEN
              sedlo  = orgfa*rcar*local_sediment_mem%sedlay(j,k,issso12)    &
                   & +      calfa*local_sediment_mem%sedlay(j,k,isssc12)    &
                   & +      oplfa*local_sediment_mem%sedlay(j,k,issssil)    &
                   & +      clafa*local_sediment_mem%sedlay(j,k,issster)
              ! "full" sediment has sedlo=1. for sedlo>1., wsed is >0.
              wsed( j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp)) ! downward shifting velocity (?)
           ENDIF

     ENDDO !end j-loop

     ! filling downward  (accumulation)
     DO iv = 1, nsedtra
        DO j = start_idx, end_idx

              IF (local_bgc_mem%bolay(j) > 0._wp) THEN
                 uebers = wsed(j)*local_sediment_mem%sedlay(j,k,iv)                     ! 'uebersaettigung?'
                 local_sediment_mem%sedlay(j,k  ,iv) = local_sediment_mem%sedlay(j,k  ,iv) - uebers
                 local_sediment_mem%sedlay(j,k+1,iv) = local_sediment_mem%sedlay(j,k+1,iv) + uebers        &
                      &             *(seddw(k)*porsol(k))/(seddw(k+1)*porsol(k+1))
              ENDIF

        ENDDO !end j-loop
     ENDDO !end iv-loop

  ENDDO !end k-loop
 
 

  ! store amount lost from bottom sediment layer - this is a kind of
  ! permanent local_sediment_mem%burial in deep consolidated layer, and this stuff is
  ! effectively lost from the whole ocean+sediment(+atmosphere) system.
  ! Would have to be supplied by river runoff or simple addition e.g.
  ! to surface layers in the long range. Can be supplied again if a
  ! sediment column has a deficiency in volume.

  DO j = start_idx, end_idx

        IF (local_bgc_mem%bolay(j) > 0._wp) THEN
           sedlo  = orgfa*rcar*local_sediment_mem%sedlay(j,ks,issso12)    &
                & +      calfa*local_sediment_mem%sedlay(j,ks,isssc12)    &
                & +      oplfa*local_sediment_mem%sedlay(j,ks,issssil)    &
                & +      clafa*local_sediment_mem%sedlay(j,ks,issster)
           wsed( j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))
        ENDIF

  ENDDO !end j-loop

  DO iv = 1, nsedtra
     DO j = start_idx, end_idx

           IF (local_bgc_mem%bolay(j) > 0._wp) THEN
              uebers = wsed(j)*local_sediment_mem%sedlay(j,ks,iv)
              local_sediment_mem%sedlay(j,ks,iv) = local_sediment_mem%sedlay(j,ks ,iv)-uebers
              local_sediment_mem%burial(j,iv)    = local_sediment_mem%burial(j,iv) + uebers*seddw(ks)*porsol(ks)
           ENDIF

     ENDDO !end j-loop
  ENDDO !end iv-loop
 
 

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

 
  ! determine how the total sediment column is filled
  DO k = 1, ks
     DO j = start_idx, end_idx
           IF (local_bgc_mem%bolay(j) > 0._wp) THEN
              sedlo  = orgfa*rcar*local_sediment_mem%sedlay(j,k,issso12)        &
                   & +      calfa*local_sediment_mem%sedlay(j,k,isssc12)        &
                   & +      oplfa*local_sediment_mem%sedlay(j,k,issssil)        &
                   & +      clafa*local_sediment_mem%sedlay(j,k,issster)
              fulsed(j) = fulsed(j) + porsol(k)*seddw(k)*sedlo
           ENDIF
     ENDDO !end j-loop
  ENDDO !end k-loop
 
 

  ! shift the sediment deficiency from the deepest (local_sediment_mem%burial)
  ! layer into layer ks
 
  DO j = start_idx, end_idx

        IF (local_bgc_mem%bolay(j) > 0._wp) THEN

           ! deficiency with respect to fully loaded sediment |packed in local_sediment_mem%sedlay(i,j,ks) ??
           ! this is the volume of sediment shifted upwards from the local_sediment_mem%burial layer

           ! 'sediment deficiency', solfu = total column inegrated solid fraction volume (bodensed)
           seddef = solfu-fulsed(j)

           ! total volume of solid constituents in buried layer
           spresent = orgfa*rcar*local_sediment_mem%burial(j,issso12)             &
                &   +      calfa*local_sediment_mem%burial(j,isssc12)             &
                &   +      oplfa*local_sediment_mem%burial(j,issssil)             &
                &   +      clafa*local_sediment_mem%burial(j,issster)

           ! determine whether an additional amount of clay is needed from the local_sediment_mem%burial
           ! layer to fill the whole sediment; I assume that there is an infinite
           ! supply of clay from below
           local_sediment_mem%burial(j,issster) = local_sediment_mem%burial(j,issster)             &
                &              + MAX(0._wp, seddef - spresent) / clafa

           ! determine new volume of buried layer
           buried = orgfa*rcar*local_sediment_mem%burial(j,issso12)               &
                & +      calfa*local_sediment_mem%burial(j,isssc12)               &
                & +      oplfa*local_sediment_mem%burial(j,issssil)               &
                & +      clafa*local_sediment_mem%burial(j,issster)

           ! fill the deepest active sediment layer
           refill=seddef/buried
           frac = porsol(ks)*seddw(ks) !changed k to ks, ik

           local_sediment_mem%sedlay(j,ks,issso12) = local_sediment_mem%sedlay(j,ks,issso12)       &
                &                 + refill*local_sediment_mem%burial(j,issso12)/frac
           local_sediment_mem%sedlay(j,ks,isssc12) = local_sediment_mem%sedlay(j,ks,isssc12)       &
                &                 + refill*local_sediment_mem%burial(j,isssc12)/frac
           local_sediment_mem%sedlay(j,ks,issssil) = local_sediment_mem%sedlay(j,ks,issssil)       &
                &                 + refill*local_sediment_mem%burial(j,issssil)/frac
           local_sediment_mem%sedlay(j,ks,issster) = local_sediment_mem%sedlay(j,ks,issster)       &
                &                 + refill*local_sediment_mem%burial(j,issster)/frac

           ! account for losses in buried sediment
           local_sediment_mem%burial(j,issso12) = local_sediment_mem%burial(j,issso12)             &
                &              - refill*local_sediment_mem%burial(j,issso12)
           local_sediment_mem%burial(j,isssc12) = local_sediment_mem%burial(j,isssc12)             &
                &              - refill*local_sediment_mem%burial(j,isssc12)
           local_sediment_mem%burial(j,issssil) = local_sediment_mem%burial(j,issssil)             &
                &              - refill*local_sediment_mem%burial(j,issssil)
           local_sediment_mem%burial(j,issster) = local_sediment_mem%burial(j,issster)             &
                &              - refill*local_sediment_mem%burial(j,issster)
        ENDIF ! local_bgc_mem%bolay >0

  ENDDO !end j-loop
 
 
  !     redistribute overload of deepest layer ks to layers 2 to ks
  DO  k = ks, 2, -1
     DO j = start_idx, end_idx

           IF (local_bgc_mem%bolay(j) > 0._wp) THEN
              sedlo  = orgfa*rcar*local_sediment_mem%sedlay(j,k,issso12)          &
                   & +      calfa*local_sediment_mem%sedlay(j,k,isssc12)          &
                   & +      oplfa*local_sediment_mem%sedlay(j,k,issssil)          &
                   & +      clafa*local_sediment_mem%sedlay(j,k,issster)
              wsed(j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))
           ENDIF

     ENDDO !end j-loop

     DO iv = 1, 4
        DO j = start_idx, end_idx
              IF (local_bgc_mem%bolay(j) > 0._wp) THEN
                 uebers = local_sediment_mem%sedlay(j,k,iv)*wsed(j)
                 frac   = porsol(k)*seddw(k)/(porsol(k-1)*seddw(k-1))
                 local_sediment_mem%sedlay(j,k,iv)   = local_sediment_mem%sedlay(j,k,iv)   - uebers
                 local_sediment_mem%sedlay(j,k-1,iv) = local_sediment_mem%sedlay(j,k-1,iv) + uebers*frac ! note k-1 here = upward shift
              ENDIF
        ENDDO !end j-loop
     ENDDO !end iv-loop
  ENDDO !end k-loop
 
 
 ENDIF ! l_up_sedshi

END SUBROUTINE 
END MODULE
