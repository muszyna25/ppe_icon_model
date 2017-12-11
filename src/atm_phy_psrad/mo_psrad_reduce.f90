MODULE mo_psrad_reduce

  USE mo_psrad_general, ONLY : wp, &
    nbndlw, nbndsw, ngptlw, ngptsw, ngpt_orig
  USE mo_psrad_lrtm_kgs, ONLY: ngpt_lrtm => ngpt
  USE mo_psrad_srtm_kgs, ONLY: ngpt_srtm => ngpt

  PRIVATE

  INTEGER, PARAMETER :: nbnd_max = MAX(nbndlw,nbndsw), &
    ngpt_max = MAX(ngptlw,ngptsw)
  INTEGER :: ngpt_current(nbnd_max), current
  REAL(wp) :: &
    rwgt(nbnd_max*ngpt_orig) ! Weights for reduced gpts

  ! RRTM weights for the original 16 g-intervals - same for LW and SW
  REAL(wp), PARAMETER :: wt(ngpt_orig) = (/& 
      0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
      0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
      0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
      0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
      0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
      0.0000750000_wp/)

  INTERFACE reduce
    MODULE PROCEDURE reduce0
    MODULE PROCEDURE reduce1
    MODULE PROCEDURE reduce2
    MODULE PROCEDURE reduce3
  END INTERFACE
  INTERFACE reduce_pack
    MODULE PROCEDURE reduce2_pack
    MODULE PROCEDURE reduce3_pack
  END INTERFACE

  PUBLIC :: init_rwgt, reduce, reduce_pack

  INTEGER, PARAMETER :: &
    ! # original gpt to combine in reduction
    ngn(ngpt_max,2) = RESHAPE((/& 
      !LW
      1,1,2,2,2,2,2,2,1,1, &                       ! band 1
      1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
      1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
      2,2,2,2,2,2,2,2, &                           ! band 6
      2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
      2,2,2,2,2,2,2,2, &                           ! band 8
      1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
      2,2,2,2,4,4, &                               ! band 10
      1,1,2,2,2,2,3,3, &                           ! band 11
      1,1,1,1,2,2,4,4, &                           ! band 12
      3,3,4,6, &                                   ! band 13
      8,8, &                                       ! band 14
      8,8, &                                       ! band 15
      4,12, &                                      ! band 16
      !SW
      2,2,2,2,4,4,             & ! band 16
      1,1,1,1,1,2,1,2,1,2,1,2, & ! band 17
      1,1,1,1,2,2,4,4,         & ! band 18
      1,1,1,1,2,2,4,4,         & ! band 19
      1,1,1,1,1,1,1,1,2,6,     & ! band 20
      1,1,1,1,1,1,1,1,2,6,     & ! band 21
      8,8,                     & ! band 22
      2,2,1,1,1,1,1,1,2,4,     & ! band 23
      2,2,2,2,2,2,2,2,         & ! band 24
      1,1,2,2,4,6,             & ! band 25
      1,1,2,2,4,6,             & ! band 26
      1,1,1,1,1,1,4,6,         & ! band 27
      1,1,2,2,4,6,             & ! band 28
      1,1,1,1,2,2,2,2,1,1,1,1, & ! band 29
      0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0/), SHAPE=(/ngpt_max,2/))

CONTAINS

  SUBROUTINE init_rwgt(which)
    INTEGER, INTENT(IN) :: which
    INTEGER :: nbnd, ntgt, i, offset, range(2), src, tgt
    REAL(wp) :: wtsum

    IF (which == 1) THEN
      nbnd = nbndlw
      ntgt = ngptlw
      ngpt_current(1:nbnd) = ngpt_lrtm
    ELSE
      nbnd = nbndsw
      ntgt = ngptsw
      ngpt_current(1:nbnd) = ngpt_srtm
    ENDIF
    current = which

    rwgt = 1.0_wp
    offset = 0
    range = 0
    src = 0
    DO tgt = 1,ntgt
      IF (range(2) == ngpt_orig) THEN
        offset = 0
        range = 0
      ENDIF
      offset = range(2)
      range = range(2) + (/1,ngn(tgt,current)/)
      wtsum = SUM(wt(range(1):range(2)))
      DO i = 1,ngn(tgt,current)
        src = src + 1
        rwgt(src) = wt(offset+i) / wtsum
      END DO
    END DO
  END SUBROUTINE init_rwgt

  SUBROUTINE reduce0(band, ref, weigh, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:)
    LOGICAL, INTENT(IN) :: weigh
    REAL(wp), INTENT(INOUT) :: k(:)

    INTEGER :: ipr, iprsm, iprsoff, igc, ngsoff
    REAL(wp) :: sumk

    k = 0
    iprsoff = (band-1) * ngpt_orig
    if (band == 1) then
      ngsoff = 0
    else
      ngsoff = SUM(ngpt_current(1:band-1))
    end if
    iprsm = 0
    DO igc = 1,ngpt_current(band)
      sumk = 0.0_wp
      DO ipr = 1, ngn(ngsoff + igc,current)
        iprsm = iprsm + 1
        IF (weigh) THEN
          sumk = sumk + ref(iprsm)*rwgt(iprsm+iprsoff)
        ELSE
          sumk = sumk + ref(iprsm)
        ENDIF
      ENDDO
      k(igc) = sumk
    ENDDO
  END SUBROUTINE reduce0

  SUBROUTINE reduce1(band, ref, weigh, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:,:)
    LOGICAL, INTENT(IN) :: weigh
    REAL(wp), INTENT(INOUT) :: k(:,:)

    INTEGER :: jt

    DO jt = 1,size(ref,1)
      CALL reduce0(band, ref(jt,:), weigh, k(jt,:))
    ENDDO
  END SUBROUTINE reduce1

  SUBROUTINE reduce2(band, ref, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:,:,:)
    REAL(wp), INTENT(INOUT) :: k(:,:,:)

    INTEGER :: jt, jp

    DO jt = 1,size(ref,1)
      DO jp = 1,size(ref,2)
        CALL reduce0(band, ref(jt,jp,:), .true., k(jt,jp,:))
      ENDDO
    ENDDO
  END SUBROUTINE reduce2

  SUBROUTINE reduce2_pack(band, ref, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:,:,:)
    REAL(wp), INTENT(INOUT) :: k(:,:)
    REAL(wp) :: tmp(size(ref,1),size(ref,2),size(k,2))

    CALL reduce2(band, ref, tmp)
    k = RESHAPE(tmp, SHAPE = (/size(k,1),size(k,2)/))
  END SUBROUTINE reduce2_pack

  SUBROUTINE reduce3(band, ref, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: k(:,:,:,:)

    INTEGER :: jn, jt, jp

    DO jn = 1,size(ref,1)
      DO jt = 1,size(ref,2)
        DO jp = 1,size(ref,3)
          CALL reduce0(band, ref(jn,jt,jp,:), .true., k(jn,jt,jp,:))
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE reduce3

  SUBROUTINE reduce3_pack(band, ref, k)
    INTEGER, INTENT(IN) :: band
    REAL(wp), INTENT(IN) :: ref(:,:,:,:)
    REAL(wp), INTENT(INOUT) :: k(:,:)
    REAL(wp) :: tmp(size(ref,1),size(ref,2),size(ref,3),size(k,2))

    CALL reduce3(band, ref, tmp)
    k = RESHAPE(tmp, SHAPE = (/size(k,1),size(k,2)/))
  END SUBROUTINE reduce3_pack
  
END MODULE mo_psrad_reduce
