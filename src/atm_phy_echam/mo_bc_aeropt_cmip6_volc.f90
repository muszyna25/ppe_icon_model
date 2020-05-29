!>
!! @brief Read and apply optical properties of aerosol climatology 
!!        for volcanic stratospheric aerosols as provided for CMIP6
!!
!! @author J.S. Rast (MPI-M)
!!
!! @par Revision History
!!      2018-05-25: original version
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_bc_aeropt_cmip6_volc

  USE mo_kind,                   ONLY: wp, i8
  USE mo_mpi,                    ONLY: my_process_is_stdio, p_bcast, p_io
  USE mo_psrad_general,          ONLY: nbndlw, nbndsw
  USE mo_exception,              ONLY: finish
  USE mo_read_interface,         ONLY: openInputFile, closeFile, read_1D, &
    &                                  read_1D_extdim_time, &
    &                                  read_extdim_slice_extdim_extdim_extdim, &
    &                                  nf
  USE mo_latitude_interpolation, ONLY: latitude_weights_li
  USE mo_physical_constants,     ONLY: rgrav, rd
  USE mo_math_constants,         ONLY: deg2rad, pi_2
  USE mo_echam_phy_config,       ONLY: echam_phy_config
  USE mtime,                     ONLY: datetime 
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  PRIVATE
  PUBLIC                           :: read_bc_aeropt_cmip6_volc, &
                                      add_bc_aeropt_cmip6_volc

  REAL(wp), ALLOCATABLE, TARGET    :: aod_v_s(:,:,:),   & ! volcanic AOD solar
                                      ssa_v_s(:,:,:,:), & ! volcanic SSA solar
                                      asy_v_s(:,:,:,:), & ! volcanic ASY solar
                                      ext_v_s(:,:,:,:), & ! volcanic EXT solar
                                      aod_v_t(:,:,:),   & ! volcanic AOD therm
                                      ssa_v_t(:,:,:,:), & ! volcanic SSA therm
                                      ext_v_t(:,:,:,:), & ! volcanic EXT therm
                                      r_alt_clim(:),    & ! altitude (km)
                                      r_lat_clim(:)       ! latitudes
  REAL(wp)                         :: r_lat_shift, r_rdeltalat
  INTEGER, SAVE                    :: inm2_time_interpolation=-999999
  LOGICAL, SAVE                    :: l_mem_alloc=.FALSE.
  INTEGER, SAVE                    :: k_alt_clim, lat_clim
  INTEGER, PARAMETER               :: nmonths=2
!!$  REAL(wp), PARAMETER              :: rdog=rd*rgrav

CONTAINS
  !>
  !! SUBROUTINE su_bc_aeropt_cmip6_volc -- sets up the memory for fields in which
  !! the aerosol optical properties are stored when needed
SUBROUTINE su_bc_aeropt_cmip6_volc( p_patch_id,                          &
                        clat_dim,          calt_dim,          kyear      )
  INTEGER, INTENT(IN)          :: p_patch_id
  INTEGER(i8), INTENT(IN)      :: kyear
  CHARACTER(LEN=*), INTENT(IN) :: clat_dim, calt_dim

  CHARACTER(LEN=256)           :: cfname, csubprog_name
  CHARACTER(LEN=32)            :: ckyear

  INTEGER                      :: ifile_id, idim_id

  csubprog_name='su_bc_aeropt_cmip6_volc---of---mo_bc_aeropt_cmip6_volc.f90'
  IF (my_process_is_stdio()) THEN
    WRITE(ckyear,*) kyear
    ! for non-amip simulations, an arbitrary year may be stored in a file without
    ! the year in the filename. However, such a boundary condition file does not
    ! yet exist.
!!$    IF ( echam_phy_config(p_patch_id)%lamip ) THEN
      cfname='bc_aeropt_cmip6_volc_lw_b16_sw_b14_'//TRIM(ADJUSTL(ckyear))//'.nc'
!!$    ELSE
!!$      cfname='bc_aeropt_cmip6_volc_lw_b16_sw_b14.nc'
!!$    ENDIF
    ifile_id=openInputFile(cfname)
!!$write(0,*) 'cfname=',TRIM(ADJUSTL(cfname)),'ifile_id=',ifile_id
    CALL nf(nf_inq_dimid(ifile_id,TRIM(ADJUSTL(clat_dim)),idim_id), &
           TRIM(ADJUSTL(csubprog_name)))
    CALL nf(nf_inq_dimlen(ifile_id, idim_id, lat_clim), &
           TRIM(ADJUSTL(csubprog_name)))
    CALL nf(nf_inq_dimid(ifile_id,TRIM(ADJUSTL(calt_dim)),idim_id), &
           TRIM(ADJUSTL(csubprog_name)))
    CALL nf(nf_inq_dimlen(ifile_id, idim_id, k_alt_clim), &
           TRIM(ADJUSTL(csubprog_name)))
  END IF
  CALL p_bcast(lat_clim,p_io)
  CALL p_bcast(k_alt_clim,p_io)
! allocate memory for optical properties
  ALLOCATE(aod_v_s(nbndsw,0:lat_clim+1,nmonths))
  ALLOCATE(ext_v_s(nbndsw,k_alt_clim+1,0:lat_clim+1,nmonths))
  ALLOCATE(ssa_v_s(nbndsw,k_alt_clim+1,0:lat_clim+1,nmonths))
  ALLOCATE(asy_v_s(nbndsw,k_alt_clim+1,0:lat_clim+1,nmonths))
  ALLOCATE(aod_v_t(nbndlw,0:lat_clim+1,nmonths))
  ALLOCATE(ext_v_t(nbndlw,k_alt_clim+1,0:lat_clim+1,nmonths))
  ALLOCATE(ssa_v_t(nbndlw,k_alt_clim+1,0:lat_clim+1,nmonths))
  ALLOCATE(r_alt_clim(k_alt_clim))
  ALLOCATE(r_lat_clim(0:lat_clim+1))
  l_mem_alloc=.TRUE.
  
! initialize with zero
  aod_v_s(:,:,:)=0._wp
  ext_v_s(:,:,:,:)=0._wp
  ssa_v_s(:,:,:,:)=0._wp
  asy_v_s(:,:,:,:)=0._wp
  aod_v_t(:,:,:)=0._wp
  ext_v_t(:,:,:,:)=0._wp
  ssa_v_t(:,:,:,:)=0._wp
  r_alt_clim(:)=0._wp
  r_lat_clim(:)=0._wp
END SUBROUTINE su_bc_aeropt_cmip6_volc

  !> SUBROUTINE shift_months_bc_aeropt_cmip6_volc -- shifts months in order to read a new one.

SUBROUTINE shift_months_bc_aeropt_cmip6_volc

  aod_v_s(:,:,1)=aod_v_s(:,:,2)
  ext_v_s(:,:,:,1)=ext_v_s(:,:,:,2)
  ssa_v_s(:,:,:,1)=ssa_v_s(:,:,:,2)
  asy_v_s(:,:,:,1)=asy_v_s(:,:,:,2)
  aod_v_t(:,:,1)=aod_v_t(:,:,2)
  ext_v_t(:,:,:,1)=ext_v_t(:,:,:,2)
  ssa_v_t(:,:,:,1)=ssa_v_t(:,:,:,2)
  
END SUBROUTINE shift_months_bc_aeropt_cmip6_volc

  !> SUBROUTINE read_bc_aeropt_cmip6_volc -- read the aerosol optical properties 
  !! of the volcanic CMIP6 aerosols

SUBROUTINE read_bc_aeropt_cmip6_volc(current_date, p_patch_id)
  TYPE(datetime), POINTER, INTENT(in) :: current_date
  INTEGER, INTENT(in)                 :: p_patch_id

  !LOCAL VARIABLES
  INTEGER(i8) :: iyear(2)
  INTEGER :: imonth(2), nmonths, imonths

  TYPE(t_time_interpolation_weights) :: tiw

  tiw = calculate_time_interpolation_weights(current_date)  
  
  IF (tiw%month2_index == inm2_time_interpolation) RETURN

  IF (ALLOCATED(aod_v_s)) THEN

    ! skip shifting (and reading) of data if we have a turn of the year,
    ! but update control index. In general, reading takes place around the
    ! middle of a month, not when a month changes.
    ! This assumes that read_bc_aeropt_cmip6_volc
    ! is called with a sufficiently high period of less than half a month,
    ! half of the time interval at which data (here monthly) are provided.

    IF ( inm2_time_interpolation == 13 .AND. tiw%month2_index == 1 ) THEN
       inm2_time_interpolation=tiw%month2_index
       RETURN
    ENDIF

    CALL shift_months_bc_aeropt_cmip6_volc
    imonth(2)=tiw%month2_index
    iyear(2)=current_date%date%year
    IF (imonth(2) == 13 ) THEN
      imonth(2)=1
      iyear(2)=iyear(2)+1
    END IF
    CALL read_month_bc_aeropt_cmip6_volc (p_patch_id, &
         imonth(2),    iyear(2),    2   )
    
    inm2_time_interpolation=tiw%month2_index
    
  ELSE
    imonth(1)=tiw%month1_index
    imonth(2)=tiw%month2_index
    iyear(1)=current_date%date%year
    iyear(2)=current_date%date%year
    IF (imonth(1) == 0) THEN
      imonth(1)=12
      iyear(1)=iyear(1)-1
    END IF
    IF (imonth(2) == 13) THEN
      imonth(2)=1
      iyear(2)=iyear(2)+1
    END IF
    nmonths=2
    inm2_time_interpolation=tiw%month2_index
    CALL su_bc_aeropt_cmip6_volc(p_patch_id,                          &
                     'latitude',        'altitude',        iyear(1)   )
    DO imonths=1,nmonths
       CALL read_month_bc_aeropt_cmip6_volc (p_patch_id, &
         imonth(imonths),    iyear(imonths),    imonths   )
    END DO
  ENDIF
END SUBROUTINE read_bc_aeropt_cmip6_volc
!-------------------------------------------------------------------------
!> SUBROUTINE add_bc_aeropt_cmip6_volc

!! add aerosol optical properties for all wave length bands (solar and IR)
!! in the case of CMIP6 volcanic aerosols
!! The height profile is taken into account.
!!
!! !REVISION HISTORY:
!! original source by J.S. Rast (2010-02-19)
!! adapted to icon by J.S. Rast (2013-09-18)
SUBROUTINE add_bc_aeropt_cmip6_volc( &
     & current_date,          jg,              jcs,                  &
     & kproma,                kbdim,           klev,                 &
     & krow,                  nb_sw,           nb_lw,                &
     & zf,                    dz,                                    &
     & paer_tau_sw_vr,        paer_piz_sw_vr,  paer_cg_sw_vr,        &
     & paer_tau_lw_vr                                                )

  ! !INPUT PARAMETERS
  TYPE(datetime), POINTER, INTENT(in) :: current_date
  INTEGER,INTENT(in)  :: jg,     &! domain index
                         jcs,    &! start index in block
                         kproma, &! actual block length
                         kbdim,  &! maximum block length
                         krow,   &! block index
                         klev,   &! number of vertical levels
                         nb_lw,  &! number of wave length bands (far IR)
                         nb_sw    ! number of wave length bands (solar)
  REAL(wp),INTENT(in) :: dz(kbdim,klev),      & ! geometric height thickness [m]
                         zf(kbdim,klev)         ! geometric height 
! !OUTPUT PARAMETERS
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_lw):: &
   paer_tau_lw_vr      !aerosol optical depth (far IR)
  REAL(wp),INTENT(inout),DIMENSION(kbdim,klev,nb_sw):: &
   paer_tau_sw_vr,   & !aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   & !weighted sum of single scattering albedos, 
                       !sum_i(tau_i*omega_i)
   paer_cg_sw_vr       !weighted sum of asymmetry factors, 
                       !sum_i(tau_i*omega_i*g_i)

! !LOCAL VARIABLES
  
  INTEGER, PARAMETER                    :: norder=-1 ! latitudes in climatology order from N->S
  INTEGER, PARAMETER                    :: nm1=1, nm2=2 ! there are only two months of a year
                                                        ! stored in the field (saves memory)
  INTEGER                               :: jl,jk,jki,jwl
  INTEGER                               :: idx_lat_1, idx_lat_2, idx_lev
  REAL(wp)                              :: w1_lat, w2_lat
  INTEGER,  DIMENSION(kbdim,klev)       :: kindex ! index field for pressure interpolation
  LOGICAL,  DIMENSION(kbdim,klev)       :: l_kindex
  REAL(wp), DIMENSION(kbdim)            :: wgt1_lat,wgt2_lat
  INTEGER,  DIMENSION(kbdim)            :: inmw1_lat, inmw2_lat 
  REAL(wp), DIMENSION(kbdim,nb_sw)      :: zaod_s, zext_s_int, zfact_s
  REAL(wp), DIMENSION(kbdim,klev,nb_sw) :: zext_s, zomg_s, zasy_s
  REAL(wp), DIMENSION(kbdim,klev,nb_lw) :: zext_t, zomg_t
  REAL(wp), DIMENSION(kbdim,nb_lw)      :: zaod_t, zext_t_int, zfact_t 
  REAL(wp)                              :: p_lat_shift, p_rdeltalat
  INTEGER                               :: jc
  REAL(wp)                              :: dz_clim, z_max_lim_clim
  REAL(wp), DIMENSION(k_alt_clim+1)     :: r_alt_clim_q
  TYPE(t_time_interpolation_weights) :: tiw

  tiw = calculate_time_interpolation_weights(current_date)
  
! It is assumed that the pressure levels of the climatology do not change with time but
! are unequally spaced. Since the pressure of each icon level may change with time,
! each specific icon level may have its centre in a different level of the climatology at
! each time step.

! 1. calculate for each icon gridbox the index of the data set layer 
!     in which p_mid_icon is located and geometrical height of layers
  dz_clim=r_alt_clim(1)-r_alt_clim(2)
  z_max_lim_clim=r_alt_clim(1)+0.5_wp*dz_clim
  CALL altitude_index( &
       & jcs,          kproma,         kbdim,           klev,              &
       & zf,           dz_clim,        z_max_lim_clim,  k_alt_clim,        &
       & kindex,       l_kindex                                            )
!!$  IF (my_process_is_stdio()) THEN
!!$  WRITE(0,*) '======================================================================'
!!$  r_alt_clim_q(1)=r_alt_clim(1)+250._wp
!!$  DO jk=2,k_alt_clim
!!$    r_alt_clim_q(jk)=0.5_wp*(r_alt_clim(jk-1)+r_alt_clim(jk))
!!$  END DO
!!$  r_alt_clim_q(k_alt_clim+1)=r_alt_clim(k_alt_clim)-250._wp
!!$  WRITE(0,*) 'r_alt_clim_q=',r_alt_clim_q
!!$  DO jk=1,klev
!!$     IF (l_kindex(1,jk)) THEN
!!$        WRITE(0,*) 'kindex(1,jk)=',kindex(1,jk),'zf=',zf(1,jk),'r_alt_clim=',r_alt_clim_q(kindex(1,jk)),r_alt_clim_q(kindex(1,jk)+1)
!!$     ELSE
!!$        WRITE(0,*) 'kindex(1,jk)=',kindex(1,jk),'zf=',zf(1,jk)
!!$     END IF
!!$  END DO
!!$  WRITE(0,*) '======================================================================'
!!$  END IF
  p_lat_shift=r_lat_shift
  p_rdeltalat=r_rdeltalat
  CALL latitude_weights_li( &
       & jg,          jcs,            kproma,            kbdim,           &
       & krow,        wgt1_lat,       wgt2_lat,          inmw1_lat,       &
       & inmw2_lat,   p_lat_shift,    p_rdeltalat,       r_lat_clim,      &
       & lat_clim,    norder                                              )
! 2. Solar radiation
! 2.1 interpolate optical properties solar radiation
  DO jwl=1,nb_sw
     DO jk=1,klev
        DO jl=jcs,kproma
           idx_lat_1=inmw1_lat(jl)
           idx_lat_2=inmw2_lat(jl)
           w1_lat=wgt1_lat(jl)
           w2_lat=wgt2_lat(jl)
           idx_lev=kindex(jl,jk)
!!$           IF (my_process_is_stdio()) THEN
!!$              IF (jwl==9.and.jl==kproma) THEN
!!$                 write(0,*) '========================================'
!!$                 write(0,*) 'jk=',jk, 'kindex(jl,jk)=',idx_lev,'zf(jl,jk)=',zf(jl,jk)
!!$                 write(0,*) 'tiw%weight1=',tiw%weight1,'tiw%weight2=',tiw%weight2
!!$                 write(0,*) 'idx_lat_1=',idx_lat_1,'w1_lat=',w1_lat,'idx_lat_2=',idx_lat_2,'w2_lat=',w2_lat
!!$                 write(0,*) 'ext_v_s(jwl,idx_lev,idx_lat_1,nm1)=',ext_v_s(jwl,idx_lev,idx_lat_1,nm1)
!!$                 write(0,*) 'ext_v_s(jwl,idx_lev,idx_lat_2,nm1)=',ext_v_s(jwl,idx_lev,idx_lat_2,nm1)
!!$                 write(0,*) 'ext_v_s(jwl,idx_lev,idx_lat_1,nm2)=',ext_v_s(jwl,idx_lev,idx_lat_1,nm2)
!!$                 write(0,*) 'ext_v_s(jwl,idx_lev,idx_lat_2,nm2)=',ext_v_s(jwl,idx_lev,idx_lat_2,nm2)
!!$                 write(0,*) '========================================'
!!$              END IF
!!$           END IF
           zext_s(jl,jk,jwl)=tiw%weight1*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm2))
           zomg_s(jl,jk,jwl)=tiw%weight1*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm2))
           zasy_s(jl,jk,jwl)=tiw%weight1*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm2))
        END DO
     END DO
  END DO
  DO jwl=1,nb_sw
     DO jl=jcs,kproma
        idx_lat_1=inmw1_lat(jl)
        idx_lat_2=inmw2_lat(jl)
        w1_lat=wgt1_lat(jl)
        w2_lat=wgt2_lat(jl)
        zaod_s(jl,jwl)=tiw%weight1*(w1_lat*aod_v_s(jwl,idx_lat_1,nm1)+ &
                               w2_lat*aod_v_s(jwl,idx_lat_2,nm1))+ &
                       tiw%weight2*(w1_lat*aod_v_s(jwl,idx_lat_1,nm2)+ &
                               w2_lat*aod_v_s(jwl,idx_lat_2,nm2))
     END DO
  END DO
! 2.2 normalize zext to the correct total optical depth
!     the normalization factor generally depends on the wavelength if
!     the ratios of the extinction at different wavelengths are not 
!     independent of the height level. Generally, the aerosol composition
!     depends on height, this leads to different ratios of the extinction
!     between two given wavelengths at different heights.
  zext_s_int(jcs:kproma,1:nb_sw)=0._wp
  DO jwl=1,nb_sw
     DO jk=1,klev
        zext_s_int(jcs:kproma,jwl)=zext_s_int(jcs:kproma,jwl) + &
          zext_s(jcs:kproma,jk,jwl)*dz(jcs:kproma,jk)
     END DO
  END DO
  WHERE (zext_s_int(jcs:kproma,1:nb_sw) > 0._wp) 
     zfact_s(jcs:kproma,1:nb_sw)=zaod_s(jcs:kproma,1:nb_sw)/ &
                              zext_s_int(jcs:kproma,1:nb_sw)
  ELSEWHERE
     zfact_s(jcs:kproma,1:nb_sw)=1._wp
  END WHERE
  DO jwl=1,nb_sw
     DO jk=1,klev
        zext_s(jcs:kproma,jk,jwl)=zext_s(jcs:kproma,jk,jwl)* &
             dz(jcs:kproma,jk)*zfact_s(jcs:kproma,jwl)
!!$             zfact_s(jcs:kproma,jwl)*1000._wp
     END DO
  END DO
! 2.3 add optical parameters to the optical parameters of aerosols
!     inverse height profile
  DO jk=1,klev
     jki=klev-jk+1
     WHERE (zext_s(jcs:kproma,jki,1:nb_sw)>0._wp) 
     paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)*paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)+&
       zext_s(jcs:kproma,jki,1:nb_sw)*zomg_s(jcs:kproma,jki,1:nb_sw)*&
       zasy_s(jcs:kproma,jki,1:nb_sw)
     paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)+&
       zext_s(jcs:kproma,jki,1:nb_sw)*zomg_s(jcs:kproma,jki,1:nb_sw)
     paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)+&
       zext_s(jcs:kproma,jki,1:nb_sw)
     paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw)/&
          paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)
     paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)=paer_cg_sw_vr(jcs:kproma,jk,1:nb_sw)/&
          (paer_tau_sw_vr(jcs:kproma,jk,1:nb_sw)*paer_piz_sw_vr(jcs:kproma,jk,1:nb_sw))
     END WHERE
  END DO
! 3. far infrared
! 2.1 interpolate optical properties thermal radiation
  DO jwl=1,nb_lw
     DO jk=1,klev
        DO jl=jcs,kproma
           idx_lat_1=inmw1_lat(jl)
           idx_lat_2=inmw2_lat(jl)
           w1_lat=wgt1_lat(jl)
           w2_lat=wgt2_lat(jl)
           idx_lev=kindex(jl,jk)
           zext_t(jl,jk,jwl)=tiw%weight1*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm2))
           zomg_t(jl,jk,jwl)=tiw%weight1*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                             tiw%weight2*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm2))
        END DO
     END DO
  END DO
  DO jwl=1,nb_lw
     DO jl=jcs,kproma
        idx_lat_1=inmw1_lat(jl)
        idx_lat_2=inmw2_lat(jl)
        w1_lat=wgt1_lat(jl)
        w2_lat=wgt2_lat(jl)
        zaod_t(jl,jwl)=tiw%weight1*(w1_lat*aod_v_t(jwl,idx_lat_1,nm1)+ &
                               w2_lat*aod_v_t(jwl,idx_lat_2,nm1))+ &
                       tiw%weight2*(w1_lat*aod_v_t(jwl,idx_lat_1,nm2)+ &
                               w2_lat*aod_v_t(jwl,idx_lat_2,nm2))
     END DO
  END DO
! 2.2 normalize zext to the correct total optical depth
!     the normalization factor generally depends on the wavelength if
!     the ratios of the extinction at different wavelengths are not 
!     independent of the height level. Generally, the aerosol composition
!     depends on height, this leads to different ratios of the extinction
!     between two given wavelengths at different heights.
  zext_t_int(jcs:kproma,1:nb_lw)=0._wp
  DO jwl=1,nb_lw
     DO jk=1,klev
        zext_t_int(jcs:kproma,jwl)=zext_t_int(jcs:kproma,jwl) + &
          zext_t(jcs:kproma,jk,jwl)*dz(jcs:kproma,jk)
     END DO
  END DO
  WHERE (zext_t_int(jcs:kproma,1:nb_lw) > 0._wp) 
     zfact_t(jcs:kproma,1:nb_lw)=zaod_t(jcs:kproma,1:nb_lw)/ &
                              zext_t_int(jcs:kproma,1:nb_lw)
  ELSEWHERE
     zfact_t(jcs:kproma,1:nb_lw)=1._wp
  END WHERE
  DO jwl=1,nb_lw
     DO jk=1,klev
        zext_t(jcs:kproma,jk,jwl)=zext_t(jcs:kproma,jk,jwl)* &
             dz(jcs:kproma,jk)*zfact_t(jcs:kproma,jwl)
!!$             zfact_t(jcs:kproma,jwl)*1000._wp
     END DO
  END DO
! 2.3 add optical parameters to the optical parameters of aerosols
!     inverse height profile
  DO jk=1,klev
     jki=klev-jk+1
     ! use explicit DO loop to circumvent possible SXf90 compiler bug
     DO jc = jcs,kproma
        paer_tau_lw_vr(jc,jk,1:nb_lw)=paer_tau_lw_vr(jc,jk,1:nb_lw) + &
          zext_t(jc,jki,1:nb_lw)*(1._wp-zomg_t(jc,jki,1:nb_lw))
     ENDDO
  END DO  

END SUBROUTINE add_bc_aeropt_cmip6_volc

!------------------------------------------------------------------------
SUBROUTINE altitude_index( &
     & jcs,            kproma,         kbdim,             klev,              &
     & zf,             dz_clim,        z_max_lim_clim,    k_lev_clim,        &
     & kindex,         l_kindex                                              )
  INTEGER, INTENT(in)  :: jcs                  ! minimum block index
  INTEGER, INTENT(in)  :: kproma, kbdim, klev  ! dimensions of ICON fields
  REAL(wp), INTENT(in) :: zf(kbdim,klev),     &! mid-layer altitudes of ICON
                        & dz_clim,            &! layer thickness of climatology
                        & z_max_lim_clim       ! altitude of highest layer limit of clim.
  INTEGER, INTENT(in)  :: k_lev_clim           ! number of layers in climatology
  INTEGER, INTENT(out) :: kindex(kbdim,klev)   ! layer index of climatology in which
                                               ! icon layer is located
  LOGICAL, INTENT(out) :: l_kindex(kbdim,klev) ! .true. if index in range [1,k_lev_clim]
  REAL(wp)             :: dz_clim_inv
  INTEGER              :: jk
  dz_clim_inv=1._wp/dz_clim
  kindex(jcs:kproma,:)=FLOOR((z_max_lim_clim-zf(jcs:kproma,:))*dz_clim_inv)+1
  WHERE (kindex(jcs:kproma,:).LT.1.OR.kindex(jcs:kproma,:).GT.k_alt_clim)
    l_kindex(jcs:kproma,:)=.FALSE.
    kindex(jcs:kproma,:)=k_lev_clim+1
  END WHERE
!!$  IF (my_process_is_stdio()) THEN
!!$     do jk=1,SIZE(zf(1,:))
!!$        WRITE(0,*) 'jk=',jk,'zf(1,jk)=',zf(1,jk),'kindex(1,jk)=',kindex(1,jk),'z_max_lim_clim=',z_max_lim_clim
!!$     end do
!!$  END IF
END SUBROUTINE altitude_index
!-------------------------------------------------------------------------
! 
!> SUBROUTINE read_month_bc_aeropt_cmip6_volc -- reads optical aerosol parameters from 
!    original file containing the following quantities (attention: longitude means wavelengths!)
!!   tautl: monthly average of zonal mean total aod for all wavelengths tautl(time, latitude, longitude) 
!!   exts: extinction for each level corresp. to tautl in 1/m: exts(time, levels, latitude, longitude) 
!!   omega: single scattering albedo corresp. to tautl: omega(time, levels, latitude, longitude)
!!   asymm: asymmetry fractor corresp. to tautl: asymm(time, levels, latitude, longitude)
!!   levels: non-equidistant pressure levels corresp. to exts, omega, asymm: levels(levels).
!!     Only the mid-point pressures are given.
  SUBROUTINE read_month_bc_aeropt_cmip6_volc (p_patch_id,  &
    kmonth,           kyear,            ktime_step   )
  !
  INTEGER, INTENT(in)            :: p_patch_id   ! id number of the patch
  INTEGER, INTENT(in)            :: kmonth       ! number of month to be read
  INTEGER(i8), INTENT(in)        :: kyear        ! year of month
  INTEGER, INTENT(in)            :: ktime_step   ! month that has to be set by new data

  CHARACTER(len=256)             :: cfname   ! file name containing variables
  CHARACTER(len=32)              :: ckmonth,ckyear
  CHARACTER(len=256)             :: cdim_names(4)

  INTEGER                        :: ifile_id,lbnd,i
  REAL(wp), POINTER              :: zvar3d(:,:,:,:)
  REAL(wp), POINTER              :: zalt(:), zlat(:)
  REAL(wp)                       :: delta_alt

  IF (kmonth < 1 .OR. kmonth > 12 ) THEN
    WRITE(ckmonth,*) kmonth
    CALL finish ('read_month_bc_aeropt_cmip6_volc in mo_bc_aeropt_cmip6_volc', &
                 'month to be read outside valid range 1<=kmonth<=12, '// &
                 'kmonth='//TRIM(ADJUSTL(ckmonth))) 
  END IF
  WRITE(ckyear,*) kyear

  IF ( echam_phy_config(p_patch_id)%lamip ) THEN
    cfname='bc_aeropt_cmip6_volc_lw_b16_sw_b14_'//TRIM(ADJUSTL(ckyear))//'.nc'
  ELSE
    cfname='bc_aeropt_cmip6_volc_lw_b16_sw_b14.nc'
  ENDIF

  ifile_id=openInputFile(cfname)
  CALL read_1D(file_id=ifile_id, variable_name='altitude', return_pointer=zalt)
  ! order altitudes from top to bottom of atmosphere (should be outside reading routine)
  r_alt_clim=zalt
  DO lbnd=1,k_alt_clim
    zalt(k_alt_clim-lbnd+1)=r_alt_clim(lbnd)
  END DO
!convert from km to metres
  r_alt_clim(1:k_alt_clim)=1000._wp*zalt(1:k_alt_clim)
  delta_alt=r_alt_clim(1)-r_alt_clim(2)
  CALL read_1D(file_id=ifile_id, variable_name='latitude', return_pointer=zlat)
  r_lat_clim(1:lat_clim)=zlat(lat_clim:1:-1)*deg2rad
  r_lat_clim(0)=0.0_wp
  r_lat_clim(lat_clim+1)=-pi_2
  r_lat_shift=r_lat_clim(1)                      ! this is the value next to the N-pole (so +87.5 for example)
  r_rdeltalat=ABS(1.0_wp/(r_lat_clim(2)-r_lat_clim(1)))
  cdim_names(1)='month'
  cdim_names(2)='altitude'
  cdim_names(3)='latitude'
  cdim_names(4)='solar_bands'
  CALL read_extdim_slice_extdim_extdim_extdim( &
    &                             file_id=ifile_id, variable_name='ext_sun', &
    &                             return_pointer=zvar3d, dim_names=cdim_names, &
    &                             start_extdim1=kmonth, end_extdim1=kmonth)
  CALL rearrange_cmip6_volc_3d (zvar3d(1,:,:,:),'ext_v_s',ktime_step)
!convert units from 1/km to 1/m
  ext_v_s(:,:,:,ktime_step)=0.001_wp*ext_v_s(:,:,:,ktime_step)
  CALL read_extdim_slice_extdim_extdim_extdim( &
    &                             file_id=ifile_id, variable_name='omega_sun', &
    &                             return_pointer=zvar3d, dim_names=cdim_names, &
    &                             start_extdim1=kmonth, end_extdim1=kmonth)
  CALL rearrange_cmip6_volc_3d (zvar3d(1,:,:,:),'ssa_v_s',ktime_step)
  CALL read_extdim_slice_extdim_extdim_extdim( &
    &                             file_id=ifile_id, variable_name='g_sun', &
    &                             return_pointer=zvar3d, dim_names=cdim_names, &
    &                             start_extdim1=kmonth, end_extdim1=kmonth)
  CALL rearrange_cmip6_volc_3d (zvar3d(1,:,:,:),'asy_v_s',ktime_step)
  cdim_names(4)='terrestrial_bands'
  CALL read_extdim_slice_extdim_extdim_extdim( &
    &                             file_id=ifile_id, variable_name='ext_earth', &
    &                             return_pointer=zvar3d, dim_names=cdim_names, &
    &                             start_extdim1=kmonth, end_extdim1=kmonth)
  CALL rearrange_cmip6_volc_3d (zvar3d(1,:,:,:),'ext_v_t',ktime_step)
  !convert units from 1/km to 1/m
  ext_v_t(:,:,:,ktime_step)=0.001_wp*ext_v_t(:,:,:,ktime_step)
  CALL read_extdim_slice_extdim_extdim_extdim( &
    &                           file_id=ifile_id, variable_name='omega_earth', &
    &                           return_pointer=zvar3d, dim_names=cdim_names, &
    &                           start_extdim1=kmonth, end_extdim1=kmonth)
  CALL rearrange_cmip6_volc_3d (zvar3d(1,:,:,:),'ssa_v_t',ktime_step)
  ! calculate AOD of atmosphere
  DO lbnd=1,nbndsw
     CALL trapezoidal_rule(delta_alt,                       &
          & ext_v_s(lbnd,1:k_alt_clim,0:lat_clim+1,ktime_step), &
          & k_alt_clim, k_alt_clim, lat_clim+2,             &
          & aod_v_s(lbnd,0:lat_clim+1,ktime_step)               )
  END DO
  DO lbnd=1,nbndlw
     CALL trapezoidal_rule(delta_alt,                       &
          & ext_v_t(lbnd,1:k_alt_clim,0:lat_clim+1,ktime_step), &
          & k_alt_clim, k_alt_clim, lat_clim+2,             &
          & aod_v_t(lbnd,0:lat_clim+1,ktime_step)               )
  END DO
  DEALLOCATE(zalt,zlat,zvar3d)
  CALL closeFile(ifile_id)
  END SUBROUTINE read_month_bc_aeropt_cmip6_volc

!-------------------------------------------------------------------------
! 
!> SUBROUTINE rearrange_cmip6_volc -- rearrange dimensions and distribute on
!!   module variables for time step time_step

  SUBROUTINE rearrange_cmip6_volc_3d (array_in,cvar,time_step)

    REAL(wp), INTENT(in)          :: array_in(:,:,:)
    CHARACTER(LEN=*), INTENT(in)  :: cvar
    INTEGER, INTENT(in)           :: time_step

    REAL(wp), POINTER             :: var(:,:,:)
    INTEGER                       :: nbands, nlev, nlat, i
    REAL(wp), ALLOCATABLE         :: array_out(:,:,:)

    nlev=SIZE(array_in,1)
    nlat=SIZE(array_in,2)
    nbands=SIZE(array_in,3)
    ALLOCATE(array_out(nbands,nlev,nlat))
!!$    IF (my_process_is_stdio()) THEN
!!$       write(0,*) ' nbands=',nbands,' nlev=',nlev,' nlat=',nlat
!!$       write(0,*) 'rearrange_cmip6_volc_3d: zvar3d(:,16,5)'
!!$       DO i=1,nlev
!!$          write(0,*) i,array_in(i,16,5)
!!$       END DO
!!$    END IF 
    SELECT CASE (TRIM(cvar))
    CASE ('ext_v_s')
      var => ext_v_s(:,:,:,time_step)
    CASE ('ext_v_t')
      var => ext_v_t(:,:,:,time_step)
    CASE ('ssa_v_s')
      var => ssa_v_s(:,:,:,time_step)
    CASE ('ssa_v_t')
      var => ssa_v_t(:,:,:,time_step)
    CASE ('asy_v_s')
      var => asy_v_s(:,:,:,time_step)
    CASE DEFAULT 
      CALL finish('rearrange_cmip6_volc_2d of mo_bc_aeropt_cmip6_volc', &
                  'Variable '//TRIM(ADJUSTL(cvar))//' unknown')
    END SELECT
    ! attention: the indices of var are
    ! (1:bands[nbnd{sw,lw}],1:nlev[k_alt_clim],1:nlat+2[lat_clim+2])
    ! order latitudes from N->S and altitudes from bottom to top
    var(1:nbands,1:nlev,2:nlat+1)= &
         RESHAPE(array_in(nlev:1:-1,nlat:1:-1,1:nbands), &
         SHAPE = (/nbands,nlev,nlat/), ORDER=(/2,3,1/))
    var(1:nbands,1:nlev,1)=var(1:nbands,1:nlev,2)
    var(1:nbands,1:nlev,nlat+2)=var(1:nbands,1:nlev,nlat+1)
!!$    IF (my_process_is_stdio()) THEN
!!$       write(0,*) 'rearrange_cmip6_volc_3d, var(5,nlev-24+1,:)'
!!$       do i=2,nlat+1
!!$          write(0,*) i,var(5,nlev-24+1,i)
!!$       end do
!!$       write(0,*) 'rearrange_cmip6_volc_3d, var(5,:,nlat-16+1+1)'
!!$       do i=1,nlev
!!$          write(0,*) i, var(5,i,nlat-16+1+1)
!!$       end do
!!$    END IF
         NULLIFY(var)
  END SUBROUTINE rearrange_cmip6_volc_3d
!-------------------------------------------------------------------------
! 
SUBROUTINE trapezoidal_rule(dx,f,NMAX,n,m,f_ts)
  ! calculates trapezoidal sum for m functions f(NMAX,m) over the first n
  ! grid points.
  ! input: dx: interval length of equidistant grid of f(1:n,1:m)
  ! input: f(NMAX,m): m tabulated functions f(1:n,1:m) on n equidistant grid
  !                   points
  ! input: NMAX: first dimension of f as declared in calling (sub)program
  ! input: n: actual number of grid points tabulated for each function f(:,1:m)
  ! input: m: number of functions for which integral is calculated
  ! output: f_ts(m): the m integrals over f(:,1:m)
  INTEGER, INTENT(in)    :: NMAX,n,m
  REAL(wp), INTENT(in)   :: dx, f(NMAX,m)
  REAL(wp), INTENT(out)  :: f_ts(m)
  f_ts(:)=dx*(SUM(f(2:n-1,:),DIM=1)+0.5_wp*(f(1,:)+f(n,:)))
END SUBROUTINE trapezoidal_rule

END MODULE mo_bc_aeropt_cmip6_volc
