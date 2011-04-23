PROGRAM sinterpw
  !
  USE mo_kind,   ONLY: dp

  INTEGER    :: nlon, nlat
  REAL (dp) ,ALLOCATABLE :: zglat(:)
  REAL (dp)  :: zglon, gilon, gilat
  REAL (dp) ,ALLOCATABLE :: sg1(:,:),sg2(:,:)
  REAL (dp)  :: pi
  CHARACTER(len=150) gfile, gifile
  REAL (dp ),ALLOCATABLE :: wt(:,:,:)

  INTEGER :: i, j

  NAMELIST /gaussgrid/ nlon, nlat, gfile
  NAMELIST /icongrid/ gifile

  !
  pi    = 2._dp*asin(1._dp)
  !
  !
  !nlon  = 640
  !nlat  = 320
  !gfile = "grid/hrefc6_0000_E000630.dat_42"
  !gifile= "grid/OPT_hex_grid.index.3"
  !gfile = "/uwork0/mripodas/stswm/mcase5/hmodc5_0360_E000550.dat"
  !gifile= "grid/h_test5T213ng8t360.gmt"
  !gifile= "grid/OPT_tri_grid.vert.3"
  !
  READ(5,gaussgrid)
  READ(5,icongrid)

  ALLOCATE (zglat(-1:nlat+2),sg1(nlon,nlat),sg2(nlon,nlat),wt(4, 2, 0:nlat))
  !
  OPEN (11,FILE=TRIM(gfile))
  !
  DO j = 1, nlat
    DO i = 1, nlon
      READ (11,'(4e30.16e3)') zglon, zglat(nlat-j+1), sg1(i,nlat-j+1),sg2(i,nlat-j+1)
      !WRITE(6,'(3e30.16e3)') zglon, zglat(nlat-j+1), sg(i,nlat-j+1)
      zglat(nlat-j+1)=zglat(nlat-j+1)*pi/180._dp
    ENDDO
  ENDDO
  !
  zglat(-1)       = -pi - zglat(2)
  zglat( 0)       = -pi - zglat(1)
  zglat(nlat+1)   =  pi - zglat(nlat)
  zglat(nlat+2)   =  pi - zglat(nlat-1)
  !
  !   Compute the weights for the horizontal interpolation from the
  !   Gaussian grid of the IFS to the triangular grid of the GME;
  !   weights for bicubic interpolation
  !
  DO j = 0,nlat
    CALL npr_lcbas (zglat(j-1), wt(1,1,j), wt(1,2,j))
  ENDDO
  !

  CALL npr_bicubicsw (nlon, nlat, zglat, wt, sg1, sg2, gifile)

  !
  CLOSE(11)
  !DO j = 1, nlat
    !DO i = 1, nlon
     ! WRITE (6,'(3e30.16e3)') zglon, zglat(j), sg(i,j)
    !ENDDO
  !ENDDO
  !
  !
  OPEN (11,FILE=TRIM(gifile))
  !
  reading_loop: DO
    READ (11,*,end=100) gilon, gilat
  !  WRITE(6,*) gilon, gilat
  ENDDO reading_loop
100 CLOSE(11)

  !


END PROGRAM sinterpw
