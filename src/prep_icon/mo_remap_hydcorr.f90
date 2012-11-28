!! This module contains subroutines for field adjustments
!! (e.g. hydrostatic correction) applied prior to horizontal
!! interpolation.
!!
MODULE mo_remap_hydcorr

  USE mo_kind,               ONLY: wp
  USE mo_parallel_config,    ONLY: nproma
  USE mo_physical_constants, ONLY: rd, grav, vtmpc1
  USE mo_nh_init_utils,      ONLY: virtual_temp
  USE mo_exception,          ONLY: message, finish
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_remap_config,       ONLY: dbg_level, MAX_NAME_LENGTH
  USE mo_mpi,                ONLY: get_my_mpi_work_id, p_comm_work,  &
    &                              p_bcast
  USE mo_remap_io,           ONLY: t_file_metadata, get_varID,       &
    &                              read_field2D, read_field3D
  USE mo_remap_sync,         ONLY: t_gather_c
  USE mo_util_binheap,       ONLY: t_heap_data
  USE mo_remap_intp,         ONLY: t_intp_data, interpolate_c
  USE mo_remap_shared,       ONLY: t_grid
  USE mo_ifs_coord,          ONLY: alloc_vct, auxhyb, geopot,        &
    &                              half_level_pressure, init_vct, vct

  IMPLICIT NONE
  INCLUDE 'cdi.inc'

  PUBLIC :: t_field_adjustment
  PUBLIC :: hydrostatic_correction

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_remap_hydcorr')  

  PRIVATE

  !> Meta-data for field adjustments (e.g. hydrostatic correction)
  !  applied prior to horizontal interpolation.
  TYPE t_field_adjustment
    LOGICAL                         :: lhydrostatic_correction
    CHARACTER (LEN=MAX_NAME_LENGTH) :: var_temp                  !< field name: "temperature"
    INTEGER                         :: code_temp
    CHARACTER (LEN=MAX_NAME_LENGTH) :: var_geosp                 !< field name: "surface geopotential"
    INTEGER                         :: code_geosp
    CHARACTER (LEN=MAX_NAME_LENGTH) :: var_qv                    !< field name: "specific humidity
    INTEGER                         :: code_qv
    REAL(wp)                        :: hpbl1                     !< height above ground of surface inversion top
    REAL(wp)                        :: hpbl2                     !< top of layer used to estimate the vertical 
                                                                 !  temperature gradient above the inversion
  END TYPE t_field_adjustment 


  INTEGER, PARAMETER :: NLEV_IFS = 91

  ! Pressure und height of 91 ECMWF model levels (main levels)
  ! for the US-Standard atmosphere (US 1976)(p_s=1013.25 hPa).
  !
  ! @author H. Frank, DWD
  !
  REAL(wp), PARAMETER :: IFS_LEVEL_DATA(4*NLEV_IFS) = &
    !    p_f [hPa]    z [m]   z-h [m]  Delta z
    &  (/     0.01 , 79303. ,  79303. , 4240.8 , &  !  1 
    &         0.03 , 72745. ,  72745. , 4240.8 , &  !  2 
    &         0.06 , 68689. ,  68689. , 3986.2 , &  !  3 
    &         0.10 , 64847. ,  64847. , 3776.6 , &  !  4 
    &         0.17 , 61204. ,  61204. , 3577.9 , &  !  5 
    &         0.28 , 57748. ,  57748. , 3389.6 , &  !  6 
    &         0.43 , 54470. ,  54470. , 3212.9 , &  !  7 
    &         0.64 , 51360. ,  51360. , 3041.5 , &  !  8 
    &         0.92 , 48442. ,  48442. , 2818.7 , &  !  9 
    &         1.30 , 45759. ,  45759. , 2571.0 , &  ! 10 
    &         1.78 , 43332. ,  43332. , 2319.6 , &  ! 11 
    &         2.38 , 41136. ,  41136. , 2101.7 , &  ! 12 
    &         3.12 , 39141. ,  39141. , 1912.2 , &  ! 13 
    &         4.02 , 37322. ,  37322. , 1746.3 , &  ! 14 
    &         5.09 , 35657. ,  35657. , 1600.4 , &  ! 15 
    &         6.34 , 34128. ,  34128. , 1471.4 , &  ! 16 
    &         7.80 , 32720. ,  32720. , 1356.9 , &  ! 17 
    &         9.47 , 31417. ,  31417. , 1260.2 , &  ! 18 
    &        11.37 , 30201. ,  30201. , 1179.7 , &  ! 19 
    &        13.50 , 29061. ,  29061. , 1106.4 , &  ! 20 
    &        15.88 , 27991. ,  27991. , 1039.4 , &  ! 21 
    &        18.52 , 26985. ,  26985. ,  978.1 , &  ! 22 
    &        21.41 , 26037. ,  26037. ,  921.9 , &  ! 23 
    &        24.57 , 25142. ,  25142. ,  870.2 , &  ! 24 
    &        27.99 , 24298. ,  24298. ,  822.5 , &  ! 25 
    &        31.67 , 23498. ,  23498. ,  778.4 , &  ! 26 
    &        35.63 , 22742. ,  22742. ,  737.6 , &  ! 27 
    &        39.85 , 22024. ,  22024. ,  699.8 , &  ! 28 
    &        44.33 , 21343. ,  21343. ,  664.7 , &  ! 29 
    &        49.07 , 20695. ,  20695. ,  632.0 , &  ! 30 
    &        54.07 , 20079. ,  20079. ,  601.6 , &  ! 31 
    &        59.31 , 19492. ,  19492. ,  574.2 , &  ! 32 
    &        64.80 , 18931. ,  18931. ,  548.3 , &  ! 33 
    &        70.51 , 18396. ,  18396. ,  523.6 , &  ! 34 
    &        76.43 , 17884. ,  17884. ,  500.5 , &  ! 35 
    &        82.57 , 17394. ,  17394. ,  480.9 , &  ! 36 
    &        88.96 , 16922. ,  16922. ,  464.6 , &  ! 37 
    &        95.62 , 16464. ,  16464. ,  451.3 , &  ! 38 
    &       102.58 , 16018. ,  16018. ,  440.7 , &  ! 39 
    &       109.89 , 15582. ,  15582. ,  432.6 , &  ! 40 
    &       117.59 , 15152. ,  15152. ,  426.8 , &  ! 41 
    &       125.75 , 14727. ,  14727. ,  423.3 , &  ! 42 
    &       134.40 , 14305. ,  14305. ,  420.8 , &  ! 43 
    &       143.59 , 13885. ,  13885. ,  418.4 , &  ! 44 
    &       153.35 , 13468. ,  13468. ,  416.0 , &  ! 45 
    &       163.72 , 13053. ,  13053. ,  413.6 , &  ! 46 
    &       174.72 , 12641. ,  12641. ,  411.2 , &  ! 47 
    &       186.38 , 12231. ,  12231. ,  408.8 , &  ! 48 
    &       198.76 , 11824. ,  11824. ,  406.4 , &  ! 49 
    &       211.87 , 11418. ,  11418. ,  404.1 , &  ! 50 
    &       225.77 , 11016. ,  11016. ,  402.2 , &  ! 51 
    &       240.48 , 10613. ,  10613. ,  404.0 , &  ! 52 
    &       256.07 , 10208. ,  10208. ,  406.5 , &  ! 53 
    &       272.56 ,  9800. ,   9800. ,  409.0 , &  ! 54 
    &       290.02 ,  9390. ,   9390. ,  411.4 , &  ! 55 
    &       308.48 ,  8977. ,   8977. ,  413.9 , &  ! 56 
    &       327.99 ,  8562. ,   8562. ,  416.3 , &  ! 57 
    &       348.62 ,  8144. ,   8144. ,  418.7 , &  ! 58 
    &       370.42 ,  7724. ,   7724. ,  421.1 , &  ! 59 
    &       393.44 ,  7302. ,   7302. ,  423.5 , &  ! 60 
    &       417.73 ,  6878. ,   6878. ,  425.6 , &  ! 61 
    &       443.34 ,  6451. ,   6451. ,  427.2 , &  ! 62 
    &       470.17 ,  6025. ,   6025. ,  424.3 , &  ! 63 
    &       497.96 ,  5605. ,   5605. ,  417.6 , &  ! 64 
    &       526.46 ,  5192. ,   5192. ,  407.3 , &  ! 65 
    &       555.40 ,  4792. ,   4792. ,  394.0 , &  ! 66 
    &       584.49 ,  4406. ,   4406. ,  378.1 , &  ! 67 
    &       613.50 ,  4036. ,   4036. ,  361.6 , &  ! 68 
    &       642.29 ,  3683. ,   3683. ,  345.0 , &  ! 69 
    &       670.73 ,  3347. ,   3347. ,  328.3 , &  ! 70 
    &       698.70 ,  3027. ,   3027. ,  311.7 , &  ! 71 
    &       726.07 ,  2724. ,   2724. ,  294.6 , &  ! 72 
    &       752.67 ,  2438. ,   2438. ,  277.4 , &  ! 73 
    &       778.40 ,  2169. ,   2169. ,  260.4 , &  ! 74 
    &       803.16 ,  1917. ,   1917. ,  243.6 , &  ! 75 
    &       826.81 ,  1682. ,   1682. ,  226.5 , &  ! 76 
    &       849.25 ,  1464. ,   1464. ,  209.4 , &  ! 77 
    &       870.38 ,  1264. ,   1264. ,  192.6 , &  ! 78 
    &       890.13 ,  1079. ,   1079. ,  176.1 , &  ! 79 
    &       908.44 ,   911. ,    911. ,  159.7 , &  ! 80 
    &       925.22 ,   760. ,    760. ,  143.4 , &  ! 81 
    &       940.44 ,   625. ,    625. ,  127.6 , &  ! 82 
    &       954.09 ,   505. ,    505. ,  112.4 , &  ! 83 
    &       966.17 ,   399. ,    399. ,   97.7 , &  ! 84 
    &       976.67 ,   309. ,    309. ,   83.3 , &  ! 85 
    &       985.63 ,   232. ,    232. ,   69.8 , &  ! 86 
    &       993.30 ,   167. ,    167. ,   60.4 , &  ! 87 
    &       999.84 ,   112. ,    112. ,   49.9 , &  ! 88 
    &      1005.12 ,    68. ,     68. ,   38.9 , &  ! 89 
    &      1009.15 ,    34. ,     34. ,   28.4 , &  ! 90 
    &      1012.05 ,    10. ,     10. ,   20.0 /)   ! 91 

  REAL(wp), PARAMETER :: IFS_LEVELS(NLEV_IFS,4) = &
    &  RESHAPE(IFS_LEVEL_DATA, (/ NLEV_IFS, 4 /) , order=(/ 2,1 /))

  INTEGER, PARAMETER :: TEMP   = 1
  INTEGER, PARAMETER :: IQV    = 2
  

CONTAINS

  SUBROUTINE compute_vtemp_and_height(file_metadata, fa, grid, gather_c, psfc, geosp, temp_v, z3d)
    TYPE (t_file_metadata),    INTENT(IN)    :: file_metadata  !< input file meta-data
    TYPE (t_field_adjustment), INTENT(IN)    :: fa             !< meta-data for hydrostatic correction
    TYPE (t_grid)  ,           INTENT(IN)    :: grid           !< source grid
    TYPE (t_gather_c),         INTENT(IN)    :: gather_c       !< communication pattern
    REAL(wp),                  INTENT(IN)    :: psfc(:,:)      !< surface orography height of input data [m] 
    REAL(wp),                  INTENT(IN)    :: geosp(:,:)     !< surface geopotential
    REAL(wp),                  INTENT(INOUT) :: temp_v(:,:,:)  !< virtual temperature (result)
    REAL(wp),                  INTENT(INOUT) :: z3d(:,:,:)     !< orography height (result)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &     TRIM(TRIM(modname)//'::compute_height')
    INTEGER :: nlen, jb,nlev, ierrstat, size_vct, &
      &        shape2D_glb(2), shape2D_loc(2), zaxisID
    REAL(wp), ALLOCATABLE :: temp(:,:,:), qv(:,:,:), &
      &                      rfield1D(:), rfield2D(:,:)
    REAL(wp), DIMENSION(nproma,NLEV_IFS+1) :: pres_ic, lnp_ic, geop_ic
    REAL(wp), DIMENSION(nproma,NLEV_IFS) :: delp, rdelp, rdlnpr, rdalpha, geop_mc

    ! allocate temporary fields
    nlev = NLEV_IFS
    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      shape2D_glb = (/ nproma,gather_c%nblks_c /)
    ELSE
      shape2D_glb = (/ 0,0 /)
    END IF
    shape2D_loc = (/ nproma, grid%p_patch%nblks_c /)
    ALLOCATE(rfield2D(    shape2D_glb(1), shape2D_glb(2)),  &
      &      rfield1D(gather_c%total_dim),                  &
      &      temp(nproma, nlev, grid%p_patch%nblks_c),      &
      &      qv(nproma, nlev, grid%p_patch%nblks_c),        &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! read 3D temperature field
    CALL read_field3D(file_metadata, fa%var_temp, fa%code_temp, nlev, gather_c, rfield1D, rfield2D, temp, zaxisID)

    ! read 3D specific humidity field
    CALL read_field3D(file_metadata, fa%var_qv, fa%code_qv, nlev, gather_c, rfield1D, rfield2D, qv)

    ! Compute virtual temperature of input data
    CALL virtual_temp(grid%p_patch, temp, qv, temp_v=temp_v )

    ! initialize vertical coordinate table
    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      size_vct = zaxisInqVctSize(zaxisID)
    END IF
    CALL p_bcast(size_vct, gather_c%rank0, p_comm_work)

    ALLOCATE(vct(size_vct), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      CALL zaxisInqVct(zaxisID, vct)
    END IF
    CALL p_bcast(vct, gather_c%rank0, p_comm_work)

    CALL alloc_vct(nlev)
    CALL init_vct(nlev)

    ! 1. Compute pressure and height of input data, using the IFS routines

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, pres_ic, lnp_ic, delp, rdelp, rdlnpr, &
!$OMP            rdalpha, geop_mc, geop_ic) 
    DO jb =1,grid%p_patch%nblks_c
      nlen = nproma
      IF (jb == grid%p_patch%nblks_c)  nlen = grid%p_patch%npromz_c
      
      CALL half_level_pressure(psfc(:,jb), nproma, nlen, nlev, pres_ic)
      CALL auxhyb(pres_ic, nproma, nlen, nlev,        & ! in
                  delp, rdelp, lnp_ic, rdlnpr, rdalpha) ! out

      CALL geopot(temp_v(:,:,jb), rdlnpr, rdalpha, geosp(:,jb), & ! in
                  nproma, 1, nlen, nlev, geop_mc, geop_ic    )    ! inout

      ! Compute 3D height coordinate field
      z3d(1:nlen,1:nlev,jb) = geop_mc(1:nlen,1:nlev)/grav
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! clean up
    DEALLOCATE(temp, qv, vct, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE compute_vtemp_and_height


  !> Compute virtual temperature.
  !
  !  Read temperature, qv, qc, qi, qr, qs from file and compute
  !  virtual temperature at lowest main level, and at hpbl1, hpbl2.
  !
  SUBROUTINE compute_vtemp_sfc(zsfc_in, nblks, npromz, temp_v, z3d, idx_hpbl1, idx_hpbl2, vtemp_sfc)
    REAL(wp),                  INTENT(IN)    :: zsfc_in(:,:)     !< surface orography height of input data [m] 
    INTEGER,                   INTENT(IN)    :: nblks, npromz
    REAL(wp),                  INTENT(IN)    :: temp_v(:,:,:), z3d(:,:,:)
    INTEGER,                   INTENT(IN)    :: idx_hpbl1, idx_hpbl2
    REAL(wp),                  INTENT(INOUT) :: vtemp_sfc(:,:)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &       TRIM(TRIM(modname)//'::hydrostatic_correction')

    INTEGER  :: i_startidx, i_endidx, jc,jb,          &
    &           idx_lowest_ml
    REAL(wp) :: dz, exdist, dtvdz

    IF (dbg_level >= 2)  CALL message(routine, "Start")
    ! linearly extrapolate surface virtual temperature from the
    ! temperature at the lowest main level using the vertical gradient
    ! between this level and the temperature at hpbl1
    idx_lowest_ml = UBOUND(z3d,2)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,dz,dtvdz,exdist) 
    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz
      DO jc=i_startidx, i_endidx
        ! extrapolation distance
        exdist = ( z3d(jc,idx_lowest_ml,jb) - zsfc_in(jc,jb) )
        ! temperature gradient:
        dz     = z3d(jc,idx_hpbl2,jb) - z3d(jc,idx_hpbl1,jb)
        dtvdz  = (temp_v(jc,idx_hpbl2,jb)  - temp_v(jc,idx_hpbl1,jb)) / dz
        vtemp_sfc(jc,jb) = temp_v(jc,idx_lowest_ml,jb) + dtvdz*exdist
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    IF (dbg_level >= 2)  CALL message(routine, "Done")
  END SUBROUTINE compute_vtemp_sfc


  !> Compute surface pressure from logarithm of surface pressure
  !
  SUBROUTINE compute_psfc(nblks, npromz, psfc)
    INTEGER,  INTENT(IN)    :: nblks, npromz  !< field blocks, length of last block
    REAL(wp), INTENT(INOUT) :: psfc(:,:)      !< logarithm of surface pressure (in) /surface pressure (result)
    ! local variables
    INTEGER :: jc, jb, i_startidx, i_endidx

!$OMP PARALLEL PRIVATE(i_startidx, i_endidx,jc)
!$OMP DO
    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx, i_endidx
        psfc(jc,jb) = EXP(psfc(jc,jb))
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE compute_psfc


  !> Compute orography height [m]
  !
  SUBROUTINE compute_zsfc(nblks, npromz, grav, zsfc)
    INTEGER,  INTENT(IN)    :: nblks, npromz  !< field blocks, length of last block
    REAL(wp), INTENT(IN)    :: grav           !< [m/s2] av. gravitational acceleration
    REAL(wp), INTENT(INOUT) :: zsfc(:,:)      !< surface geopotential (in) / surface height (result)
    ! local variables
    INTEGER :: jc, jb, i_startidx, i_endidx

!$OMP PARALLEL PRIVATE(i_startidx, i_endidx,jc)
!$OMP DO
    DO jb=1,nblks
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks) i_endidx = npromz

      DO jc=i_startidx, i_endidx
        zsfc(jc,jb) = zsfc(jc,jb)/grav
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE compute_zsfc


  !> Compute hydrostatic correction for surface pressure field, then perform
  !  horizontal interpolation operation.
  !
  !  @par Revision History
  !  Initial version:                    Guenther Zaengl, DWD (2012-05-20)
  !  Implementation for SCRIP algorithm: F. Prill,        DWD (2012-10-29)
  !
  SUBROUTINE hydrostatic_correction(intp_data, src_grid, dst_grid, file_metadata, &
    &                               gather_c, fa, psfc_in, psfc_out)
    TYPE(t_intp_data),         INTENT(IN)    :: intp_data       !< data structure with intp. weights
    TYPE (t_grid),             INTENT(IN)    :: src_grid        !< source grid
    TYPE (t_grid),             INTENT(IN)    :: dst_grid        !< destination grid
    TYPE (t_file_metadata),    INTENT(IN)    :: file_metadata   !< input file meta-data
    TYPE (t_gather_c),         INTENT(IN)    :: gather_c        !< communication pattern
    TYPE (t_field_adjustment), INTENT(IN)    :: fa              !< meta-data for hydrostatic correction
    REAL(wp),                  INTENT(INOUT) :: psfc_in(:,:)    !< logarithmic surface pressure (in)
    REAL(wp),                  INTENT(INOUT) :: psfc_out(:,:)   !< output surface pressure (result)

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = &
      &  TRIM(TRIM(modname)//':hydrostatic_correction')

    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    REAL(wp), PARAMETER :: dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

    INTEGER         :: i, jc,jb, ierrstat, nidx, i_startidx, i_endidx,   &
      &                ic, ib, js, shape2D_glb(2), shape2D_loc(2),       &
      &                nblks_in, npromz_in, nblks_out, npromz_out, nlev, &
      &                idx_hpbl1, idx_hpbl2
    REAL(wp)        :: val, rhpbl1, rhpbl2, dtvdz1, dtvdz2, exdist,    &
      &                tv_aux2, tv_aux3
    REAL(wp), ALLOCATABLE :: vtemp_sfc(:,:), zsfc_in(:,:),             &
      &                zsfc_out(:,:), rfield1D(:), rfield2D(:,:),      &
      &                p_aux1(:), p_aux2(:), p_aux3(:), z3d(:,:,:),    &
      &                temp_v(:,:,:)

    ! In general, the stencil size (and the stencil) depends on the current
    ! grid point:
    INTEGER                         :: nstpts 
    TYPE (t_heap_data), ALLOCATABLE :: stencil_data(:) 


    IF (dbg_level >= 11)  WRITE (0,*) "# perform hydrostatic correction."

    nblks_in   = src_grid%p_patch%nblks_c
    npromz_in  = src_grid%p_patch%npromz_c
    nblks_out  = dst_grid%p_patch%nblks_c
    npromz_out = dst_grid%p_patch%npromz_c

    ! allocate temporary fields:
    nlev = NLEV_IFS
    ALLOCATE(stencil_data(intp_data%max_totstencil),                              &
      &      p_aux1(intp_data%max_totstencil), p_aux2(intp_data%max_totstencil),  &
      &      p_aux3(intp_data%max_totstencil),                                    &
      &      vtemp_sfc(nproma, nblks_in),                                         &
      &      zsfc_in(nproma, nblks_in), zsfc_out(nproma, nblks_out),              &
      &      z3d(nproma,nlev,nblks_in),                                           &
      &      temp_v(nproma, nlev, nblks_in),                                      &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")
    ! allocate global input field only on PE rank0
    IF (get_my_mpi_work_id() == gather_c%rank0) THEN
      shape2D_glb = (/ nproma,gather_c%nblks_c /)
    ELSE
      shape2D_glb = (/ 0,0 /)
    END IF
    shape2D_loc = (/ nproma, nblks_in /)
    ALLOCATE(rfield2D(    shape2D_glb(1), shape2D_glb(2)), &
      &      rfield1D(gather_c%total_dim), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

    ! determine level indices nearest to hpbl1, hpbl2
    idx_hpbl1  = get_nearest_vlevel(fa%hpbl1, IFS_LEVELS(:,2))
    idx_hpbl2  = get_nearest_vlevel(fa%hpbl2, IFS_LEVELS(:,2))

    ! --- read logarithm of surface pressure and surface geopotential
    IF (dbg_level >= 5) &
      & WRITE (0,*) "read logarithm of surface pressure and surface geopotential"
    CALL read_field2D(file_metadata, fa%var_geosp, fa%code_geosp, gather_c, &
      &               rfield1D, rfield2D, zsfc_in)

    ! --- compute surface pressure from LNSP
    IF (dbg_level >= 5) &
      & WRITE (0,*) "compute surface pressure from LNSP"
    CALL compute_psfc(nblks_in, npromz_in, psfc_in)

    ! --- compute virtual temperature and height of input levels
    IF (dbg_level >= 5) &
      & WRITE (0,*) "compute virtual temperature and height of input levels"
    CALL compute_vtemp_and_height(file_metadata, fa, src_grid, gather_c, psfc_in, &
      &                           zsfc_in, temp_v, z3d)

    ! --- compute surface orography height [m]
    IF (dbg_level >= 5) &
      & WRITE (0,*) "compute surface orography height"
    CALL compute_zsfc(nblks_in, npromz_in, grav, zsfc_in)

    ! --- compute surface virtual temperature
    IF (dbg_level >= 5) &
      & WRITE (0,*) "compute surface virtual temperature"
    CALL compute_vtemp_sfc(zsfc_in, nblks_in, npromz_in, temp_v, z3d, &
      &                    idx_hpbl1, idx_hpbl2, vtemp_sfc)

    ! --- horizontal interpolation of surface orography height
    IF (dbg_level >= 5) &
      & WRITE (0,*) "horizontal interpolation of surface orography height"
    CALL interpolate_c(zsfc_in, zsfc_out, dst_grid, intp_data)

    ! --- loop over grid points in 2D field:

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,nidx,nstpts,stencil_data,js,ic,ib,rhpbl1,rhpbl2, &
!$OMP            dtvdz1,dtvdz2,p_aux1,exdist,tv_aux2,p_aux2,tv_aux3,p_aux3,i,val) 
    DO jb=1,nblks_out
      i_startidx = 1
      i_endidx   = nproma
      IF (jb == nblks_out) i_endidx = npromz_out

      DO jc=i_startidx, i_endidx
        ! --- loop over stencil points:
        nidx = intp_data%nidx(jc,jb)
        ! local stencil size
        nstpts = nidx  +                 &                            ! standard stencil
          &      intp_data%smax(jc,jb) - intp_data%smin(jc,jb) + 1    ! + sequential list
        ! build local stencil:
        stencil_data(1:nidx)%wgt      = intp_data%wgt (1:nidx, jc,jb) ! standard stencil part
        stencil_data(1:nidx)%sidx     = intp_data%iidx(1:nidx, jc,jb)
        stencil_data(1:nidx)%sblk     = intp_data%iblk(1:nidx, jc,jb)
        IF (intp_data%smax(jc,jb) > 0) THEN
          stencil_data((nidx+1):nstpts) = &                           ! sequential list part
            &    intp_data%sl(intp_data%smin(jc,jb):intp_data%smax(jc,jb))
        END IF

        ! --- estimate pressure at the height of the target point for each
        !     stencil point of the source data
        DO js=1,nstpts
          ic = stencil_data(js)%sidx
          ib = stencil_data(js)%sblk
          rhpbl1 = z3d(ic, idx_hpbl1, ib)
          rhpbl2 = z3d(ic, idx_hpbl2, ib)

          ! near-surface vertical temperature gradient
          dtvdz1 = (temp_v(ic,idx_hpbl1,ib) - vtemp_sfc(ic,ib)) / rhpbl1
          ! vertical temperature gradient above (possible) surface inversion
          dtvdz2 = (temp_v(ic,idx_hpbl2,ib) - temp_v(ic,idx_hpbl1,ib)) / (rhpbl2 - rhpbl1)

          ! 1a) integrate from sfc to hpbl1
          IF (ABS(dtvdz1) > dtvdz_thresh) THEN
            p_aux1(js) = psfc_in(ic,ib)*EXP(-grav/(rd*dtvdz1) * &
              &  LOG(temp_v(ic,idx_hpbl1,ib)/vtemp_sfc(ic,ib)))
          ELSE
            p_aux1(js) = psfc_in(ic,ib)*EXP(-grav*rhpbl1 / &
              &  (rd*0.5_wp*(temp_v(ic,idx_hpbl1,ib) + vtemp_sfc(ic,ib))) )
          ENDIF

          ! 1b) prolongate/shorten air column with the temperature gradient
          !     between hpbl1 and hpbl2 (i.e. dtvdz2)
          
          exdist  = zsfc_out(jc,jb) - zsfc_in(ic,ib)         ! extrapolation distance
          tv_aux2 = temp_v(ic,idx_hpbl1,ib) + dtvdz2*exdist ! estimated intermediate temperature
          IF (ABS(dtvdz2) > dtvdz_thresh) THEN
            p_aux2(js) = p_aux1(js)*EXP(-grav/(rd*dtvdz2)* &
                            LOG(tv_aux2/temp_v(ic,idx_hpbl1,ib)) )
          ELSE
            p_aux2(js) = p_aux1(js)*EXP(-grav*exdist / &
                           (rd*0.5_wp*(temp_v(ic,idx_hpbl2,ib)+tv_aux2)) )
          ENDIF

          ! 1c) integrate down to surface height of the target point

          tv_aux3 = tv_aux2 - dtvdz1*rhpbl1 ! estimated temperature at the height of the target point
          IF (ABS(dtvdz1) > dtvdz_thresh) THEN
            p_aux3(js) = LOG( p_aux2(js)*EXP(-grav/(rd*dtvdz1) * LOG(tv_aux3/tv_aux2) ) )
          ELSE
            p_aux3(js) = LOG( p_aux2(js)*EXP(grav*rhpbl1 / (rd*0.5_wp*(tv_aux2+tv_aux3)) ) )
          ENDIF

        END DO ! js

        ! --- execute horizontal interpolation
        val = 0._wp
        DO i=1,nstpts
          val = val + stencil_data(i)%wgt * p_aux3(i)
        END DO
        psfc_out(jc,jb) = val
      END DO ! jc
    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! clean up
    DEALLOCATE(stencil_data, p_aux1, p_aux2, p_aux3,                     &
      &        vtemp_sfc, zsfc_in, zsfc_out, rfield2D, rfield1D, z3d, temp_v, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE hydrostatic_correction


  !> @return Index of vertical level that is nearest to a given value.
  !
  FUNCTION get_nearest_vlevel(hvalue, levels) RESULT(idx)
    INTEGER :: idx
    REAL(wp), INTENT(IN) :: hvalue
    REAL(wp), INTENT(IN) :: levels(:)      !< vertical levels
    ! local variables:
    REAL(wp) :: mindist, dist
    INTEGER  :: nlev, i

    nlev = SIZE(levels)
    IF (nlev == 0) THEN
      idx = 0
      RETURN
    END IF
    idx     = 1
    mindist = ABS(levels(idx) - hvalue)
    DO i=2,nlev
      dist = ABS(levels(i) - hvalue)
      IF (dist < mindist) THEN
        idx     = i
        mindist = dist
      END IF
    END DO
  END FUNCTION get_nearest_vlevel

END MODULE mo_remap_hydcorr
