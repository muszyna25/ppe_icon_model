!>
!! @brief Read and apply optical properties of aerosol climatology 
!!        by G. Stenchikov (volcanic stratospheric aerosols)
!!        This is an adaption of mo_aero_volc of echam6 to icon
!!
!! @author J.S. Rast (MPI-M)
!!
!! @par Revision History
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!

MODULE mo_aero_stenchikov

  USE mo_kind,                 ONLY: wp
  USE mo_model_domain,         ONLY: t_patch, p_patch
  USE mo_parallel_config,      ONLY: nproma
  USE mo_lrtm_par,             ONLY: nbndlw
  USE mo_srtm_config,          ONLY: nbndsw=>jpsw
  USE mo_exception,            ONLY: finish
  USE mo_netcdf_read,          ONLY: netcdf_open_input, netcdf_close
  USE mo_netcdf_read,          ONLY: netcdf_read_2d_time, netcdf_read_3d_time, netcdf_read_1d
  USE mo_time_interpolation,   ONLY: wgt1_m=>wgt1_limm, wgt2_m=>wgt2_limm, &
                                     nmw1_m=>inm1_limm, nmw2_m=>inm2_limm
  USE mo_physical_constants,   ONLY: rgrav, rd
  USE mo_math_constants,       ONLY: deg2rad, pi_2
!!$  USE mo_echam_phy_memory,     ONLY: prm_field

  IMPLICIT NONE

  PRIVATE
  PUBLIC                           :: read_aero_stenchikov, add_aop_stenchikov

  INTERFACE reorder_stenchikov
    MODULE PROCEDURE reorder_stenchikov_2d
    MODULE PROCEDURE reorder_stenchikov_3d
  END INTERFACE reorder_stenchikov

  REAL(wp), ALLOCATABLE, TARGET    :: aod_v_s(:,:,:),   & ! volcanic AOD solar
                                      ssa_v_s(:,:,:,:), & ! volcanic SSA solar
                                      asy_v_s(:,:,:,:), & ! volcanic ASY solar
                                      ext_v_s(:,:,:,:), & ! volcanic EXT solar
                                      aod_v_t(:,:,:),   & ! volcanic AOD therm
                                      ssa_v_t(:,:,:,:), & ! volcanic SSA therm
                                      ext_v_t(:,:,:,:), & ! volcanic EXT therm
                                      p_lim_clim(:),    & ! limit lev press. data
                                      r_lat_clim(:)       ! latitudes
  REAL(wp)                         :: r_lat_shift, r_lat_max, r_rdeltalat
  INTEGER, SAVE                    :: inm2_time_interpolation=-999999
  INTEGER, PARAMETER               :: lev_clim=40, lat_clim=90, nmonths=2
  REAL(wp), PARAMETER              :: rdog=rd*rgrav
  LOGICAL                          :: laero_set=.false.

CONTAINS
  !>
  !! SUBROUTINE su_aero_kinne -- sets up the memory for fields in which
  !! the aerosol optical properties are stored when needed
SUBROUTINE su_aero_stenchikov

! allocate memory for optical properties

  ALLOCATE(aod_v_s(nbndsw,0:lat_clim+1,nmonths))
  ALLOCATE(ext_v_s(nbndsw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(ssa_v_s(nbndsw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(asy_v_s(nbndsw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(aod_v_t(nbndlw,0:lat_clim+1,nmonths))
  ALLOCATE(ext_v_t(nbndlw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(ssa_v_t(nbndlw,lev_clim,0:lat_clim+1,nmonths))
  ALLOCATE(p_lim_clim(lev_clim+1))
  ALLOCATE(r_lat_clim(0:lat_clim+1))

! initialize with zero
  aod_v_s(:,:,:)=0._wp
  ext_v_s(:,:,:,:)=0._wp
  ssa_v_s(:,:,:,:)=0._wp
  asy_v_s(:,:,:,:)=0._wp
  aod_v_t(:,:,:)=0._wp
  ext_v_t(:,:,:,:)=0._wp
  ssa_v_t(:,:,:,:)=0._wp
  p_lim_clim(lev_clim+1)=0._wp
END SUBROUTINE su_aero_stenchikov

  !> SUBROUTINE shift_months_aero_stenchikov -- shifts months in order to read a new one.

SUBROUTINE shift_months_aero_stenchikov

  aod_v_s(:,:,1)=aod_v_s(:,:,2)
  ext_v_s(:,:,:,1)=ext_v_s(:,:,:,2)
  ssa_v_s(:,:,:,1)=ssa_v_s(:,:,:,2)
  asy_v_s(:,:,:,1)=asy_v_s(:,:,:,2)
  aod_v_t(:,:,1)=aod_v_t(:,:,2)
  ext_v_t(:,:,:,1)=ext_v_t(:,:,:,2)
  ssa_v_t(:,:,:,1)=ssa_v_t(:,:,:,2)
  
END SUBROUTINE shift_months_aero_stenchikov

  !> SUBROUTINE read_aero_stenchikov -- read the aerosol optical properties 
  !! of the volcanic (Stratospheric) Stenchikov aerosols

SUBROUTINE read_aero_stenchikov(kyear, p_patch)
  
  INTEGER, INTENT(in)           :: kyear
  TYPE(t_patch), INTENT(in)     :: p_patch

  !LOCAL VARIABLES
  INTEGER                       :: imonth(2), iyear(2), nmonths, imonths
!!$  INTEGER                       :: zyrm1, zyr, zyrp1
!!$  LOGICAL                       :: lnewyear
  CHARACTER(len=20)             :: cfname_base,cyr
  CHARACTER(len=25)             :: cfname

  IF (nmw2_m == inm2_time_interpolation) RETURN
  IF (ALLOCATED(aod_v_s)) THEN
    CALL shift_months_aero_stenchikov
    imonth(1)=nmw2_m
    iyear(1)=kyear
    IF (imonth(1) == 13 ) THEN
      imonth(1)=1
      iyear(1)=iyear(1)+1
    END IF
    imonth(2)=imonth(1)
    iyear(2)=iyear(1)
    nmonths=1
  ELSE
    CALL su_aero_stenchikov
    imonth(1)=nmw1_m
    imonth(2)=nmw2_m
    iyear(1)=kyear
    iyear(2)=kyear
    IF (imonth(1) == 0) THEN
      imonth(1)=12
      iyear(1)=iyear(1)-1
    END IF
    IF (imonth(2) == 13) THEN
      imonth(2)=1
      iyear(2)=iyear(2)+1
    END IF
    nmonths=2
  ENDIF
  inm2_time_interpolation=nmw2_m
  DO imonths=1,nmonths
  CALL read_months_aero_stenchikov ( &
                     'tautl',            'exts',            'omega',              'asym',  &
                     'levels',           'longitude',       'latitude',           'levels',&
                     imonth(imonths),    iyear(imonths),    imonths                        )
  END DO
END SUBROUTINE read_aero_stenchikov
!-------------------------------------------------------------------------
!> SUBROUTINE add_aop_stenchikov

!! add aerosol optical properties for all wave length bands (solar and IR)
!! in the case of volcanic aerosols of Stenchikov
!! The height profile is taken into account.
!!
!! !REVISION HISTORY:
!! original source by J.S. Rast (2010-02-19)
!! adapted to icon by J.S. Rast (2013-09-18)
SUBROUTINE add_aop_stenchikov ( jg,                                       &
          & kproma,                 kbdim,              klev,             &
          & krow,                   nb_lw,              nb_sw,            &
          & paer_tau_lw_vr,         paer_tau_sw_vr,     paer_piz_sw_vr,   &
          & paer_cg_sw_vr,          ppd_hl,             pp_fl,            &
          & tk_fl                                                         )

! !INPUT PARAMETERS
  INTEGER,INTENT(in)  :: jg,     &! domain index
                         kproma, &! actual block length
                         kbdim,  &! maximum block length
                         krow,   &! block index
                         klev,   &! number of vertical levels
                         nb_lw,  &! number of wave length bands (far IR)
                         nb_sw    ! number of wave length bands (solar)
  REAL(wp),INTENT(in) :: ppd_hl(kbdim,klev)  ,& ! layer pressure thickness 
                         pp_fl(kbdim,klev)   ,& ! pressure at "full levels"
                         tk_fl(kbdim,klev)      ! temperature at "full lev."
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
  INTEGER, PARAMETER                    :: nm1=1, nm2=2 ! there are only two months stored
  INTEGER                               :: jl,jk,jki,jwl
  INTEGER                               :: idx_lat_1, idx_lat_2, idx_lev
  REAL(wp)                              :: w1_lat, w2_lat
  REAL(wp), DIMENSION(kbdim,klev)       :: zdeltag    ! layer thickness [m]
  INTEGER,  DIMENSION(kbdim,klev)       :: kindex ! index field for pressure interpolation
  REAL(wp), DIMENSION(kbdim)            :: wgt1_lat,wgt2_lat
  INTEGER,  DIMENSION(kbdim)            :: inmw1_lat, inmw2_lat 
  REAL(wp), DIMENSION(kbdim,nb_sw)      :: zaod_s, zext_s_int, zfact_s
  REAL(wp), DIMENSION(kbdim,klev,nb_sw) :: zext_s, zomg_s, zasy_s
  REAL(wp), DIMENSION(kbdim,klev,nb_lw) :: zext_t, zomg_t
  REAL(wp), DIMENSION(kbdim,nb_lw)      :: zaod_t, zext_t_int, zfact_t 
  REAL(wp)                              :: p_lat_shift, p_rdeltalat
  INTEGER                               :: n_lat

! It is assumed that the pressure levels of the climatology do not change with time but
! are unequally spaced. Since the pressure of each icon level may change with time,
! each specific icon level may have its centre in a different level of the climatology at
! each time step.

! 1. calculate for each icon gridbox the index of the data set layer 
!     in which p_mid_icon is located and geometrical height of layers
  CALL pressure_index(kproma,        kbdim,         klev,              &
                      pp_fl,         lev_clim,      p_lim_clim,        &
                      kindex)
  zdeltag(1:kproma,1:klev)= &
       & ppd_hl(1:kproma,1:klev)* &
       & tk_fl(1:kproma,1:klev)/pp_fl(1:kproma,1:klev)*rdog
  p_lat_shift=r_lat_shift
  p_rdeltalat=r_rdeltalat
  n_lat=lat_clim
  CALL latitude_weights_li(jg                                                   &
                        & ,kproma              ,kbdim            ,krow          &
                        & ,wgt1_lat            ,wgt2_lat         ,inmw1_lat     &
                        & ,inmw2_lat           ,p_lat_shift      ,p_rdeltalat   &
                        & ,norder                                               )
! 2. Solar radiation
! 2.1 interpolate optical properties solar radiation
  DO jwl=1,nb_sw
     DO jk=1,klev
        DO jl=1,kproma
           idx_lat_1=inmw1_lat(jl)
           idx_lat_2=inmw2_lat(jl)
           w1_lat=wgt1_lat(jl)
           w2_lat=wgt2_lat(jl)
           idx_lev=kindex(jl,jk)
           zext_s(jl,jk,jwl)=wgt1_m*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             wgt2_m*(w1_lat*ext_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ext_v_s(jwl,idx_lev,idx_lat_2,nm2))
           zomg_s(jl,jk,jwl)=wgt1_m*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             wgt2_m*(w1_lat*ssa_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ssa_v_s(jwl,idx_lev,idx_lat_2,nm2))
           zasy_s(jl,jk,jwl)=wgt1_m*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm1))+ &
                             wgt2_m*(w1_lat*asy_v_s(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*asy_v_s(jwl,idx_lev,idx_lat_2,nm2))
        END DO
     END DO
  END DO
  DO jwl=1,nb_sw
     DO jl=1,kproma
        idx_lat_1=inmw1_lat(jl)
        idx_lat_2=inmw2_lat(jl)
        w1_lat=wgt1_lat(jl)
        w2_lat=wgt2_lat(jl)
        zaod_s(jl,jwl)=wgt1_m*(w1_lat*aod_v_s(jwl,idx_lat_1,nm1)+ &
                               w2_lat*aod_v_s(jwl,idx_lat_2,nm1))+ &
                       wgt2_m*(w1_lat*aod_v_s(jwl,idx_lat_1,nm2)+ &
                               w2_lat*aod_v_s(jwl,idx_lat_2,nm2))
     END DO
  END DO
! 2.2 normalize zext to the correct total optical depth
!     the normalization factor generally depends on the wavelength if
!     the ratios of the extinction at different wavelengths are not 
!     independent of the height level. Generally, the aerosol composition
!     depends on height, this leads to different ratios of the extinction
!     between two given wavelengths at different heights.
  zext_s_int(1:kproma,1:nb_sw)=0._wp
  DO jwl=1,nb_sw
     DO jk=1,klev
        zext_s_int(1:kproma,jwl)=zext_s_int(1:kproma,jwl) + &
          zext_s(1:kproma,jk,jwl)*zdeltag(1:kproma,jk)
     END DO
  END DO
  WHERE (zext_s_int(1:kproma,1:nb_sw) > 0._wp) 
     zfact_s(1:kproma,1:nb_sw)=zaod_s(1:kproma,1:nb_sw)/ &
                              zext_s_int(1:kproma,1:nb_sw)
  ELSEWHERE
     zfact_s(1:kproma,1:nb_sw)=1._wp
  END WHERE
  DO jwl=1,nb_sw
     DO jk=1,klev
        zext_s(1:kproma,jk,jwl)=zext_s(1:kproma,jk,jwl)* &
             zdeltag(1:kproma,jk)*zfact_s(1:kproma,jwl)
     END DO
  END DO
! 2.3 add optical parameters to the optical parameters of aerosols
!     inverse height profile
  DO jk=1,klev
     jki=klev-jk+1
     WHERE (zext_s(1:kproma,jki,1:nb_sw)>0._wp) 
     paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw)*paer_cg_sw_vr(1:kproma,jk,1:nb_sw)+&
       zext_s(1:kproma,jki,1:nb_sw)*zomg_s(1:kproma,jki,1:nb_sw)*&
       zasy_s(1:kproma,jki,1:nb_sw)
     paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*&
       paer_piz_sw_vr(1:kproma,jk,1:nb_sw)+&
       zext_s(1:kproma,jki,1:nb_sw)*zomg_s(1:kproma,jki,1:nb_sw)
     paer_tau_sw_vr(1:kproma,jk,1:nb_sw)=paer_tau_sw_vr(1:kproma,jk,1:nb_sw)+&
       zext_s(1:kproma,jki,1:nb_sw)
     paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_piz_sw_vr(1:kproma,jk,1:nb_sw)/&
          paer_tau_sw_vr(1:kproma,jk,1:nb_sw)
     paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_cg_sw_vr(1:kproma,jk,1:nb_sw)/&
          (paer_tau_sw_vr(1:kproma,jk,1:nb_sw)*paer_piz_sw_vr(1:kproma,jk,1:nb_sw))
     END WHERE
  END DO
! 3. far infrared
! 2.1 interpolate optical properties solar radiation
  DO jwl=1,nb_lw
     DO jk=1,klev
        DO jl=1,kproma
           idx_lat_1=inmw1_lat(jl)
           idx_lat_2=inmw2_lat(jl)
           w1_lat=wgt1_lat(jl)
           w2_lat=wgt2_lat(jl)
           idx_lev=kindex(jl,jk)
           zext_t(jl,jk,jwl)=wgt1_m*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                             wgt2_m*(w1_lat*ext_v_t(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ext_v_t(jwl,idx_lev,idx_lat_2,nm2))
           zomg_t(jl,jk,jwl)=wgt1_m*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm1)+ &
                                     w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm1))+ &
                             wgt2_m*(w1_lat*ssa_v_t(jwl,idx_lev,idx_lat_1,nm2)+ &
                                     w2_lat*ssa_v_t(jwl,idx_lev,idx_lat_2,nm2))
        END DO
     END DO
  END DO
  DO jwl=1,nb_lw
     DO jl=1,kproma
        idx_lat_1=inmw1_lat(jl)
        idx_lat_2=inmw2_lat(jl)
        w1_lat=wgt1_lat(jl)
        w2_lat=wgt2_lat(jl)
        zaod_t(jl,jwl)=wgt1_m*(w1_lat*aod_v_t(jwl,idx_lat_1,nm1)+ &
                               w2_lat*aod_v_t(jwl,idx_lat_2,nm1))+ &
                       wgt2_m*(w1_lat*aod_v_t(jwl,idx_lat_1,nm2)+ &
                               w2_lat*aod_v_t(jwl,idx_lat_2,nm2))
     END DO
  END DO
! 2.2 normalize zext to the correct total optical depth
!     the normalization factor generally depends on the wavelength if
!     the ratios of the extinction at different wavelengths are not 
!     independent of the height level. Generally, the aerosol composition
!     depends on height, this leads to different ratios of the extinction
!     between two given wavelengths at different heights.
  zext_t_int(1:kproma,1:nb_lw)=0._wp
  DO jwl=1,nb_lw
     DO jk=1,klev
        zext_t_int(1:kproma,jwl)=zext_t_int(1:kproma,jwl) + &
          zext_t(1:kproma,jk,jwl)*zdeltag(1:kproma,jk)
     END DO
  END DO
  WHERE (zext_t_int(1:kproma,1:nb_lw) > 0._wp) 
     zfact_t(1:kproma,1:nb_lw)=zaod_t(1:kproma,1:nb_lw)/ &
                              zext_t_int(1:kproma,1:nb_lw)
  ELSEWHERE
     zfact_t(1:kproma,1:nb_lw)=1._wp
  END WHERE
  DO jwl=1,nb_lw
     DO jk=1,klev
        zext_t(1:kproma,jk,jwl)=zext_t(1:kproma,jk,jwl)* &
             zdeltag(1:kproma,jk)*zfact_t(1:kproma,jwl)
     END DO
  END DO
! 2.3 add optical parameters to the optical parameters of aerosols
!     inverse height profile
  DO jk=1,klev
     jki=klev-jk+1
     paer_tau_lw_vr(1:kproma,jk,1:nb_lw)=paer_tau_lw_vr(1:kproma,jk,1:nb_lw)+ &
          zext_t(1:kproma,jki,1:nb_lw)*(1._wp-zomg_t(1:kproma,jki,1:nb_lw))
  END DO  

END SUBROUTINE add_aop_stenchikov

!-------------------------------------------------------------------------
! 
!> SUBROUTINE pressure_index -- ! find index of pressure layer in which mid level 
!    pressures of icon are located
! !REVISION HISTORY:
! original source by J.S. Rast (2010-02-17)
! adapted to icon by J.S. Rast (2013-09-19)

SUBROUTINE pressure_index(kproma,        kbdim,         klev,              &
                          pp_mid,        klevels,       pp_bound,          &
                          kindex)

! !INPUT PARAMETERS
  INTEGER, INTENT(in)    :: kbdim, kproma, klev
  INTEGER, INTENT(in)    :: klevels !number of layers for indices are searched
  REAL(wp), INTENT(in)   :: pp_mid(kbdim,klev), & !echam midlevel pressures
                            pp_bound(klevels+1) !pressure at layer 
                                    !bounds of reference pressures
  INTEGER, INTENT(out)   :: kindex(kbdim,klev) !layer indices for echam press.

! !LOCAL VARIABLES

  LOGICAL                :: lp(kproma), lrepeat
  INTEGER                :: jk,il,kidx(kbdim)
  
  kidx(1:kproma)=2
  DO jk=1,klev
10   CONTINUE
     lrepeat=.FALSE.
     DO il=1,kproma
        lp(il)=pp_mid(il,jk).GT.pp_bound(kidx(il)).AND.kidx(il).LE.klevels
     END DO
     DO il=1,kproma
        IF (lp(il)) THEN
           kidx(il)=kidx(il)+1
           lrepeat=.TRUE.
        END IF
     END DO
     IF (lrepeat) THEN
        GOTO 10
     ELSE
        kindex(1:kproma,jk)=kidx(1:kproma)-1
     END IF
  END DO
END SUBROUTINE pressure_index

!-------------------------------------------------------------------------
! 
!> SUBROUTINE read_months_aero_stenchikov -- reads optical aerosol parameters from 
!    original file containing the following quantities (attention: longitude means wavelengths!)
!!   tautl: monthly average of zonal mean total aod for all wavelengths tautl(time, latitude, longitude) 
!!   exts: extinction for each level corresp. to tautl in 1/m: exts(time, levels, latitude, longitude) 
!!   omega: single scattering albedo corresp. to tautl: omega(time, levels, latitude, longitude)
!!   asymm: asymmetry fractor corresp. to tautl: asymm(time, levels, latitude, longitude)
!!   levels: non-equidistant pressure levels corresp. to exts, omega, asymm: levels(levels).
!!     Only the mid-point pressures are given.
  SUBROUTINE read_months_aero_stenchikov ( &
    caod,             cext,             cssa,               casy,            &
    clevels,          cwave_dim,        clat_dim,           clev_dim,        &
    kmonth,           kyear,            ktime_step                           )
!
  CHARACTER(len=*), INTENT(in)   :: caod,       &! name of variable containing optical depth of column
                                    cext,       &! name of variable containing level dependent extinction 
                                    cssa,       &! name of variable containing level dependent 
                                                 ! single scattering albedo
                                    casy,       &! name of variable containing level dependent 
                                                 ! asymmetry factor
                                    clevels,    &! mid point pressure of levels (not equidistant)
                                    cwave_dim,  &! name of wavelength dimension
                                    clat_dim,   &! name of latitude dimension
                                    clev_dim     ! name of level dimension
  INTEGER, INTENT(in)            :: kmonth       ! number of month to be read
  INTEGER, INTENT(in)            :: kyear        ! year of month
  INTEGER, INTENT(in)            :: ktime_step   ! month that has to be set by new data

  CHARACTER(len=256)             :: cfname   ! file name containing variables
  CHARACTER(len=32)              :: ckmonth,ckyear,ci_length,cj_length
  CHARACTER(len=256)             :: cdim_names(4)

!!$  TYPE(t_patch), INTENT(in)      :: p_patch
!!$  CHARACTER(len=*), INTENT(in), OPTIONAL     :: casl ! name of variable containing altitude of layer centres

  INTEGER                        :: ifile_id, kreturn
  REAL(wp), POINTER              :: zvar2d(:,:,:), zvar3d(:,:,:,:)
  REAL(wp), POINTER              :: zpmid(:), zlat(:)
!!$  REAL(wp), POINTER              :: zaod(:,:,:,:), zssa(:,:,:,:), zasy(:,:,:,:), zaer_ex(:,:,:,:)
!!$  CHARACTER(LEN=32)              :: cimnthb, cimnthe
!!$  CHARACTER(LEN=512)             :: cfnameyear,cyear

  IF (kmonth < 1 .OR. kmonth > 12 ) THEN
    WRITE(ckmonth,*) kmonth
    CALL finish ('read_months_aero_stenchikov in mo_aero_stenchikov', &
                 'month to be read outside valid range 1<=kmonth<=12, '// &
                 'kmonth='//TRIM(ADJUSTL(ckmonth))) 
  END IF
  WRITE(ckyear,*) kyear
  cfname='sato'//TRIM(ADJUSTL(ckyear))//'.nc'
  ifile_id=netcdf_open_input(cfname)
  cdim_names(1)=cwave_dim
  cdim_names(2)=clat_dim
  cdim_names(3)='time'
  zvar2d=>netcdf_read_2d_time(file_id=ifile_id, variable_name='tauttl', &
         dim_names=cdim_names(1:3), start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov (zvar2d(:,:,1),'aod',ktime_step)
  cdim_names(4)=cdim_names(3)
  cdim_names(3)=clev_dim
  zvar3d=>netcdf_read_3d_time(file_id=ifile_id, variable_name='exts', &
         dim_names=cdim_names, start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov (zvar3d(:,:,:,1),'ext',ktime_step)
  zvar3d=>netcdf_read_3d_time(file_id=ifile_id, variable_name='omega', &
         dim_names=cdim_names, start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov (zvar3d(:,:,:,1),'ssa',ktime_step)
  zvar3d=>netcdf_read_3d_time(file_id=ifile_id, variable_name='asymm', &
         dim_names=cdim_names, start_timestep=kmonth, end_timestep=kmonth)
  CALL reorder_stenchikov (zvar3d(:,:,:,1),'asy',ktime_step)
  zpmid=>netcdf_read_1d(file_id=ifile_id, variable_name=clev_dim)
  CALL p_lim_stenchikov(zpmid)
  zlat=>netcdf_read_1d(file_id=ifile_id, variable_name=clat_dim)
  IF (SIZE(zlat)/=lat_clim) THEN
    WRITE(ci_length,*) SIZE(zlat)
    WRITE(cj_length,*) lat_clim
    CALL finish ('read_months_aero_stenchikov of mo_aero_stenchikov','lat_clim= '// &
                 TRIM(ADJUSTL(cj_length))//' expected but found '//TRIM(ADJUSTL(ci_length))// &
                 ' elements.')
  END IF
  r_lat_clim(1:lat_clim)=zlat(lat_clim:1:-1)*deg2rad
  r_lat_clim(0)=0.0_wp
  r_lat_clim(lat_clim+1)=-pi_2
  r_lat_shift=r_lat_clim(1)                      ! this is the value at the N-pole (so +89)
  r_lat_max=r_lat_clim(lat_clim)               ! this is the value at the S-pole (so -89)
  r_rdeltalat=ABS(1.0_wp/(r_lat_clim(2)-r_lat_clim(1)))
  DEALLOCATE(zpmid,zlat,zvar2d,zvar3d)
  kreturn=netcdf_close(ifile_id)
  END SUBROUTINE read_months_aero_stenchikov
!-------------------------------------------------------------------------
! 
!> SUBROUTINE reorder_stenchikov -- reorder dimensions and distribute on
!!   module variables for time step time_step

  SUBROUTINE reorder_stenchikov_2d (pvar,cvar,time_step)

    REAL(wp), INTENT(in)          :: pvar(:,:)
    CHARACTER(LEN=*), INTENT(in)  :: cvar
    INTEGER, INTENT(in)           :: time_step

    REAL(wp), POINTER             :: var_solar(:,:), var_thermal(:,:)

    SELECT CASE (TRIM(cvar))
    CASE ('aod')
      var_solar => aod_v_s(:,:,time_step)
      var_thermal => aod_v_t(:,:,time_step)
    CASE DEFAULT 
      CALL finish('reorder_stenchikov_2d of mo_aero_stenchikov', &
                  'Variable '//TRIM(ADJUSTL(cvar))//' unknown')
    END SELECT
    ! attention: the indices of var_solar and var_thermal are (1:nbnd{sw,lw},1:lat_clim+2)
    var_solar(1:nbndsw-1,2:lat_clim+1)=pvar(nbndsw-1:1:-1,lat_clim:1:-1)
    var_solar(nbndsw,2:lat_clim+1)=pvar(nbndsw,lat_clim:1:-1)
    var_solar(1:nbndsw,1)=var_solar(1:nbndsw,1)
    var_solar(1:nbndsw,lat_clim+2)=var_solar(1:nbndsw,lat_clim)
    var_thermal(1:nbndlw-1,2:lat_clim+1)=pvar(nbndsw+nbndlw:nbndsw+2:-1,lat_clim:1:-1)
    var_thermal(nbndlw,2:lat_clim+1)=pvar(nbndsw-1,lat_clim:1:-1)
    var_thermal(1:nbndlw,1)=var_thermal(1:nbndlw,1)
    var_thermal(1:nbndlw,lat_clim+2)=var_thermal(1:nbndlw,lat_clim)
  END SUBROUTINE reorder_stenchikov_2d

!-------------------------------------------------------------------------
! 
!> SUBROUTINE reorder_stenchikov -- reorder dimensions and distribute on
!!   module variables for time step time_step

  SUBROUTINE reorder_stenchikov_3d (pvar,cvar,time_step)

    REAL(wp), INTENT(in)          :: pvar(:,:,:)
    CHARACTER(LEN=*), INTENT(in)  :: cvar
    INTEGER, INTENT(in)           :: time_step

    REAL(wp), POINTER             :: var_solar(:,:,:), var_thermal(:,:,:)

    ALLOCATE(var_solar(1:nbndsw,1:lev_clim,1:lat_clim+2), var_thermal(1:nbndlw,1:lev_clim,1:lat_clim+2))
    SELECT CASE (TRIM(cvar))
    CASE ('ext')
      var_solar => ext_v_s(:,:,:,time_step)
      var_thermal => ext_v_t(:,:,:,time_step)
    CASE ('ssa')
      var_solar => ssa_v_s(:,:,:,time_step)
      var_thermal => ssa_v_t(:,:,:,time_step)
    CASE ('asy')
      var_solar => asy_v_s(:,:,:,time_step)
      IF (ASSOCIATED(var_thermal)) NULLIFY(var_thermal)
    CASE DEFAULT 
      CALL finish('reorder_stenchikov_2d of mo_aero_stenchikov', &
                  'Variable '//TRIM(ADJUSTL(cvar))//' unknown')
    END SELECT
    ! attention: the indices of var_solar and var_thermal are (1:nbnd{sw,lw},1:lev_clim,1:lat_clim+2)
    var_solar(1:nbndsw-1,1:lev_clim,2:lat_clim+1)= &
      RESHAPE(pvar(nbndsw-1:1:-1,lat_clim:1:-1,lev_clim:1:-1), &
              (/nbndsw-1,lev_clim,lat_clim/),ORDER=(/1,3,2/))
    var_solar(nbndsw:nbndsw,1:lev_clim,2:lat_clim+1)= &
      RESHAPE(pvar(nbndsw,lat_clim:1:-1,lev_clim:1:-1),(/1,lev_clim,lat_clim/),ORDER=(/1,3,2/))
    var_solar(1:nbndsw,1:lev_clim,1)=var_solar(1:nbndsw,1:lev_clim,1)
    var_solar(1:nbndsw,1:lev_clim,lat_clim+2)=var_solar(1:nbndsw,1:lev_clim,lat_clim)
    IF (ASSOCIATED(var_thermal)) THEN
      var_thermal(1:nbndlw-1,1:lev_clim,2:lat_clim+1)= &
        RESHAPE(pvar(nbndsw+nbndlw:nbndsw+2:-1,lat_clim:1:-1,lev_clim:1:-1), &
                (/nbndlw-1,lev_clim,lat_clim/),ORDER=(/1,3,2/))
      var_thermal(nbndlw:nbndlw,1:lev_clim,2:lat_clim+1)= &
        RESHAPE(pvar(nbndsw-1,lat_clim:1:-1,lev_clim:1:-1),(/1,lev_clim,lat_clim/),ORDER=(/1,3,2/))
      var_thermal(1:nbndlw,1:lev_clim,1)=var_thermal(1:nbndlw,1:lev_clim,1)
      var_thermal(1:nbndlw,1:lev_clim,lat_clim+2)=var_thermal(1:nbndlw,1:lev_clim,lat_clim)
    END IF
    NULLIFY(var_solar,var_thermal)
!    DEALLOCATE(var_solar,var_thermal)
  END SUBROUTINE reorder_stenchikov_3d

!-------------------------------------------------------------------------
! 
!> SUBROUTINE p_lim_stenchikov -- calculate pressure at layer boundaries

  SUBROUTINE p_lim_stenchikov (p_pmid)

    REAL(wp), INTENT(inout)        :: p_pmid(lev_clim)

    INTEGER                        :: jk

    p_pmid(1:lev_clim)=p_pmid(lev_clim:1:-1)
    p_lim_clim(1)=0.0_wp
    DO jk=2,lev_clim+1
      p_lim_clim(jk)=2.0_wp*p_pmid(jk-1)-p_lim_clim(jk-1)
    END DO

  END SUBROUTINE p_lim_stenchikov
    
!-------------------------------------------------------------------------
! 
!> SUBROUTINE latitude_weights_li  -- calculate weights and indices for 
!             linear latitude interpolation.

  SUBROUTINE latitude_weights_li(jg                                                   &
                               & ,kproma              ,kbdim            ,krow         &
                               & ,wgt1_lat            ,wgt2_lat         ,inmw1_lat    &
                               & ,inmw2_lat           ,p_lat_shift      ,p_rdeltalat  &
                               & ,n_order                                             )

    ! n_order=1 if latitudes of climatology are in ascending (S->N), -1 if 
    ! latitudes are in descending (N->S) order.
    INTEGER, INTENT(in)               :: jg, kproma, kbdim, krow, n_order
    REAL(wp), INTENT(inout)           :: wgt1_lat(kbdim), wgt2_lat(kbdim)
    INTEGER, INTENT(inout)            :: inmw1_lat(kbdim), inmw2_lat(kbdim)
    REAL(wp), INTENT(in)              :: p_lat_shift, p_rdeltalat

    INTEGER                           :: jl
    REAL(wp)                          :: zlat(kbdim)
    
    zlat(1:kproma)=p_patch(jg)%cells%center(1:kproma,krow)%lat
    inmw1_lat(1:kproma)=MAX(INT(n_order*(zlat(1:kproma)-p_lat_shift)*p_rdeltalat+1),0)
    inmw2_lat(1:kproma)=inmw1_lat(1:kproma)+1
    wgt2_lat(1:kproma)=n_order*(zlat(1:kproma)-r_lat_clim(inmw1_lat(1:kproma)))*p_rdeltalat
    wgt1_lat(1:kproma)=1.0_wp-wgt2_lat(1:kproma)
!!$    write(0,*) '++++++++++++++++++++++++++++++'
!!$    write(0,*) 'latitudes=',MAXVAL(zlat(1:kproma))
!!$    write(0,*) (zlat(jl),inmw1_lat(jl),inmw2_lat(jl),wgt1_lat(jl),wgt2_lat(jl),jl=1,kproma)
!!$    write(0,*) '++++++++++++++++++++++++++++++'
  END SUBROUTINE latitude_weights_li

END MODULE mo_aero_stenchikov
