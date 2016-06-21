!>
!! @brief Read and apply monthly aerosol optical properties of S. Kinne
!! from yearly files.
!!
!! @author J.S. Rast (MPI-M)
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_bc_aeropt_kinne

  USE mo_kind,                 ONLY: wp, i8
  USE mo_model_domain,         ONLY: t_patch
  USE mo_parallel_config,      ONLY: nproma
  USE mo_lrtm_par,             ONLY: nbndlw
  USE mo_srtm_config,          ONLY: nbndsw=>jpsw
  USE mo_exception,            ONLY: finish, message
  USE mo_io_config,            ONLY: default_read_method
  USE mo_read_interface,       ONLY: openInputFile, closeFile, on_cells, &
    &                                t_stream_id, read_0D_real, read_3D_time
  USE mo_physical_constants,   ONLY: grav, rgrav, rd
  USE mo_echam_phy_memory,     ONLY: prm_field
  USE mo_echam_phy_config,     ONLY: echam_phy_config

  USE mo_bcs_time_interpolation, ONLY: tiw => current_time_interpolation_weights
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC                           :: read_bc_aeropt_kinne, set_bc_aeropt_kinne 

  REAL(wp), TARGET, ALLOCATABLE    :: aod_c_s(:,:,:,:), aod_f_s(:,:,:,:), &
                                      ssa_c_s(:,:,:,:), ssa_f_s(:,:,:,:), &
                                      asy_c_s(:,:,:,:), asy_f_s(:,:,:,:), &
                                      aod_c_f(:,:,:,:),                   &
                                      ssa_c_f(:,:,:,:),                   &
                                      asy_c_f(:,:,:,:),                   &
                                      z_km_aer_c_mo(:,:,:,:), &
                                      z_km_aer_f_mo(:,:,:,:)
  INTEGER(i8), SAVE                :: pre_year=-999999
  INTEGER, PARAMETER               :: lev_clim=40, nmonths=12
  REAL(wp)                         :: dz_clim
  REAL(wp)                         :: rdz_clim

CONTAINS
  !>
  !! SUBROUTINE su_bc_aeropt_kinne -- sets up the memory for fields in which
  !! the aerosol optical properties are stored when needed
SUBROUTINE su_bc_aeropt_kinne(p_patch)

  TYPE(t_patch), INTENT(in)       :: p_patch

  INTEGER                         :: nblks_len, nblks
  
  nblks=p_patch%nblks_c
  nblks_len=nproma
! allocate memory for optical properties
  ALLOCATE(aod_c_s(nblks_len,nbndsw,nblks,0:nmonths+1))
  ALLOCATE(aod_f_s(nblks_len,nbndsw,nblks,0:nmonths+1))
  ALLOCATE(ssa_c_s(nblks_len,nbndsw,nblks,0:nmonths+1))
  ALLOCATE(ssa_f_s(nblks_len,nbndsw,nblks,0:nmonths+1))
  ALLOCATE(asy_c_s(nblks_len,nbndsw,nblks,0:nmonths+1))
  ALLOCATE(asy_f_s(nblks_len,nbndsw,nblks,0:nmonths+1))
  ALLOCATE(aod_c_f(nblks_len,nbndlw,nblks,0:nmonths+1))
  ALLOCATE(ssa_c_f(nblks_len,nbndlw,nblks,0:nmonths+1))
  ALLOCATE(asy_c_f(nblks_len,nbndlw,nblks,0:nmonths+1))
  ALLOCATE(z_km_aer_c_mo(nblks_len,lev_clim,nblks,0:nmonths+1))
  ALLOCATE(z_km_aer_f_mo(nblks_len,lev_clim,nblks,0:nmonths+1))
! initialize with zero
  aod_c_s(:,:,:,:)=0._wp
  aod_f_s(:,:,:,:)=0._wp
  ssa_c_s(:,:,:,:)=0._wp
  ssa_f_s(:,:,:,:)=0._wp
  asy_c_s(:,:,:,:)=0._wp
  asy_f_s(:,:,:,:)=0._wp
  aod_c_f(:,:,:,:)=0._wp
  ssa_c_f(:,:,:,:)=0._wp
  asy_c_f(:,:,:,:)=0._wp
  z_km_aer_c_mo(:,:,:,:)=0._wp
  z_km_aer_f_mo(:,:,:,:)=0._wp
END SUBROUTINE su_bc_aeropt_kinne

  !> SUBROUTINE shift_months_bc_aeropt_kinne -- shifts December of current year into imonth=0 and 
  !! January of the following year into imonth=1 (these months do not need to be read again.

SUBROUTINE shift_months_bc_aeropt_kinne

  aod_c_s(:,:,:,0:1)=aod_c_s(:,:,:,12:13)
  aod_f_s(:,:,:,0:1)=aod_f_s(:,:,:,12:13)
  ssa_c_s(:,:,:,0:1)=ssa_c_s(:,:,:,12:13)
  ssa_f_s(:,:,:,0:1)=ssa_f_s(:,:,:,12:13)
  asy_c_s(:,:,:,0:1)=asy_c_s(:,:,:,12:13)
  asy_f_s(:,:,:,0:1)=asy_f_s(:,:,:,12:13)
  aod_c_f(:,:,:,0:1)=aod_c_f(:,:,:,12:13)
  ssa_c_f(:,:,:,0:1)=ssa_c_f(:,:,:,12:13)
  asy_c_f(:,:,:,0:1)=asy_c_f(:,:,:,12:13)
  
END SUBROUTINE shift_months_bc_aeropt_kinne

  !> SUBROUTINE read_bc_aeropt_kinne -- read the aerosol optical properties 
  !! of the Kinne aerosols

SUBROUTINE read_bc_aeropt_kinne(year, p_patch)
  
  INTEGER(i8), INTENT(in)       :: year
  TYPE(t_patch), INTENT(in)     :: p_patch

  !LOCAL VARIABLES
  INTEGER                       :: imonthb, imonthe

  IF (year > pre_year) THEN
    IF (ALLOCATED(aod_c_s)) THEN
      CALL shift_months_bc_aeropt_kinne
      imonthb=2
      imonthe=13
    ELSE
      CALL su_bc_aeropt_kinne(p_patch)
      imonthb=0
      imonthe=13
    ENDIF
    CALL read_months_bc_aeropt_kinne ( &
                     'aod',            'ssa',    'asy',                        'z_aer_coarse_mo',  &
                     'delta_z',        'lnwl',   'lev',                        imonthb,            &
                     imonthe,          year,     'bc_aeropt_kinne_sw_b14_coa', p_patch             )
    CALL read_months_bc_aeropt_kinne ( &
                     'aod',            'ssa',    'asy',                        'z_aer_coarse_mo',  &
                     'delta_z',        'lnwl',   'lev',                        imonthb,            &
                     imonthe,          year,     'bc_aeropt_kinne_lw_b16_coa', p_patch             )
    CALL read_months_bc_aeropt_kinne ( &
                     'aod',            'ssa',    'asy',                        'z_aer_fine_mo',    &
                     'delta_z',        'lnwl',   'lev',                        imonthb,            &
                     imonthe,          year,     'bc_aeropt_kinne_sw_b14_fin', p_patch             )
    rdz_clim=1._wp/dz_clim
    pre_year=year
  END IF    
END SUBROUTINE read_bc_aeropt_kinne
!-------------------------------------------------------------------------
!> SUBROUTINE set_bc_aeropt_kinne
!! set aerosol optical properties for all wave length bands (solar and IR)
!! in the case of the climatology of optical properties compiled by S.Kinne.
!! The height profile is taken into account.
!!
!! !REVISION HISTORY:
!! original source by J.S. Rast (2009-11-03) for echam6
!! adapted to icon by J.S. Rast (2013-08-28)
SUBROUTINE set_bc_aeropt_kinne (    jg,                                   &
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
  REAL(wp),INTENT(out),DIMENSION(kbdim,klev,nb_lw):: &
   paer_tau_lw_vr      !aerosol optical depth (far IR)
  REAL(wp),INTENT(out),DIMENSION(kbdim,klev,nb_sw):: &
   paer_tau_sw_vr,   & !aerosol optical depth (solar), sum_i(tau_i)
   paer_piz_sw_vr,   & !weighted sum of single scattering albedos, 
                       !sum_i(tau_i*omega_i)
   paer_cg_sw_vr       !weighted sum of asymmetry factors, 
                       !sum_i(tau_i*omega_i*g_i)

! !LOCAL VARIABLES
  
  INTEGER                     :: jl,jk,jwl
  REAL(wp), DIMENSION(kbdim,klev)   :: zh, &    ! altitude above surface
                                       zdeltag, & ! layer thickness
                                       zh_vr, &
                                       zdeltag_vr
  REAL(wp), DIMENSION(kbdim)        :: zq_int ! integral height profile
  REAL(wp), DIMENSION(kbdim,nb_lw)  :: zs_i
  REAL(wp), DIMENSION(kbdim,nb_sw)  :: zt_c, zt_f, &
                                       zs_c, zs_f, &
                                       zg_c, zg_f, & ! time interpolated
                                       ! aod, ssa ,ssa*asy 
                                       ! (coarse (c), fine natural (n), 
                                       !  fine anthropogenic (a))
                                       ztaua_c,ztaua_f ! optical depths
                                       ! at various altitudes
  REAL(wp), DIMENSION(kbdim,klev)   :: zq_aod_c, zq_aod_f ! altitude profile
                                       ! on echam grid (coarse and fine mode)
  INTEGER, DIMENSION(kbdim)         :: kindex ! index field
  REAL(wp), PARAMETER               :: rdog=rd/grav

! (i) calculate altitude above NN and layer thickness in 
!     echam for altitude profiles
     zdeltag(1:kproma,1:klev)= &
          & ppd_hl(1:kproma,1:klev)* &
          & tk_fl(1:kproma,1:klev)/pp_fl(1:kproma,1:klev)*rdog
     DO jk=1,klev
        zh(1:kproma,jk)=(prm_field(jg)%geom(1:kproma,jk,krow)+prm_field(jg)%geoi(1:kproma,klev+1,krow))*rgrav
     END DO
     DO jk=1,klev
        zdeltag_vr(1:kproma,jk)=zdeltag(1:kproma,klev-jk+1)
        zh_vr(1:kproma,jk)=zh(1:kproma,klev-jk+1)
     END DO
! (ii) calculate height profiles on echam grid for coarse and fine mode
     zq_aod_f(1:kproma,1:klev)=0._wp
     zq_aod_c(1:kproma,1:klev)=0._wp
     DO jk=1,klev
        kindex(1:kproma)=MAX(INT(zh_vr(1:kproma,jk)*rdz_clim+0.5_wp),1)
        DO jl=1,kproma
           IF (kindex(jl) > 0 .and. kindex(jl) <= lev_clim ) THEN
              zq_aod_c(jl,jk)= &
                & z_km_aer_c_mo(jl,kindex(jl),krow,tiw%month1)*tiw%weight1+ &
                & z_km_aer_c_mo(jl,kindex(jl),krow,tiw%month2)*tiw%weight2
              zq_aod_f(jl,jk)= &
               & z_km_aer_f_mo(jl,kindex(jl),krow,tiw%month1)*tiw%weight1+ &
               & z_km_aer_f_mo(jl,kindex(jl),krow,tiw%month2)*tiw%weight2
           END IF
        END DO
     END DO
! normalize height profile for coarse mode
     zq_int(1:kproma)=0._wp
     DO jk=1,klev
        zq_int(1:kproma)=zq_int(1:kproma)+ &
                       & zq_aod_c(1:kproma,jk)*zdeltag_vr(1:kproma,jk)
     ENDDO
     WHERE (zq_int(1:kproma) <= 0._wp)
        zq_int(1:kproma)=1._wp
     END WHERE
     DO jk=1,klev
        zq_aod_c(1:kproma,jk)=zdeltag_vr(1:kproma,jk)*zq_aod_c(1:kproma,jk)/ &
                            & zq_int(1:kproma)
     END DO
! normalize height profile for fine mode
     zq_int(1:kproma)=0._wp
     DO jk=1,klev
        zq_int(1:kproma)=zq_int(1:kproma)+ &
                       & zq_aod_f(1:kproma,jk)*zdeltag_vr(1:kproma,jk)
     ENDDO
     WHERE (zq_int(1:kproma) <= 0._wp)
        zq_int(1:kproma)=1._wp
     END WHERE
     DO jk=1,klev
        zq_aod_f(1:kproma,jk)=zdeltag_vr(1:kproma,jk)*zq_aod_f(1:kproma,jk)/ &
                            & zq_int(1:kproma)
     END DO

! (iii) far infrared
     zs_i(1:kproma,1:nb_lw)=1._wp-(tiw%weight1*ssa_c_f(1:kproma,1:nb_lw,krow,tiw%month1)+ &
                                   tiw%weight2*ssa_c_f(1:kproma,1:nb_lw,krow,tiw%month2))
     DO jk=1,klev
        DO jwl=1,nb_lw
           paer_tau_lw_vr(1:kproma,jk,jwl)=zq_aod_c(1:kproma,jk) * &
                zs_i(1:kproma,jwl) * &
                (tiw%weight1*aod_c_f(1:kproma,jwl,krow,tiw%month1) + &
                 tiw%weight2*aod_c_f(1:kproma,jwl,krow,tiw%month2)) 
        END DO
     END DO
! (iii) solar radiation
! time interpolated single scattering albedo (omega_f, omega_c)
     zs_c(1:kproma,1:nb_sw) = ssa_c_s(1:kproma,1:nb_sw,krow,tiw%month1)*tiw%weight1 + &
                              ssa_c_s(1:kproma,1:nb_sw,krow,tiw%month2)*tiw%weight2
     zs_f(1:kproma,1:nb_sw) = ssa_f_s(1:kproma,1:nb_sw,krow,tiw%month1)*tiw%weight1 + &
                              ssa_f_s(1:kproma,1:nb_sw,krow,tiw%month2)*tiw%weight2
! time interpolated asymmetry factor x ssa (omega_c*g_c, omega_{n,a}*g_{n,a})
     zg_c(1:kproma,1:nb_sw) = zs_c(1:kproma,1:nb_sw) * &
                              (asy_c_s(1:kproma,1:nb_sw,krow,tiw%month1)*tiw%weight1 + &
                               asy_c_s(1:kproma,1:nb_sw,krow,tiw%month2)*tiw%weight2)
     zg_f(1:kproma,1:nb_sw) = zs_f(1:kproma,1:nb_sw) * &
                              (asy_f_s(1:kproma,1:nb_sw,krow,tiw%month1)*tiw%weight1 + &
                               asy_f_s(1:kproma,1:nb_sw,krow,tiw%month2)*tiw%weight2)
! time interpolated aerosol optical depths
     zt_c(1:kproma,1:nb_sw)=tiw%weight1*aod_c_s(1:kproma,1:nb_sw,krow,tiw%month1) + &
                          & tiw%weight2*aod_c_s(1:kproma,1:nb_sw,krow,tiw%month2)
     zt_f(1:kproma,1:nb_sw)=tiw%weight1*aod_f_s(1:kproma,1:nb_sw,krow,tiw%month1) + &
                          & tiw%weight2*aod_f_s(1:kproma,1:nb_sw,krow,tiw%month2)
! height interpolation
! calculate optical properties
  DO jk=1,klev
! aerosol optical depth 
     DO jwl=1,nb_sw
        ztaua_c(1:kproma,jwl) = zt_c(1:kproma,jwl)*zq_aod_c(1:kproma,jk)
        ztaua_f(1:kproma,jwl) = zt_f(1:kproma,jwl)*zq_aod_f(1:kproma,jk)
     END DO
     paer_tau_sw_vr(1:kproma,jk,1:nb_sw) = ztaua_c(1:kproma,1:nb_sw) + &
                                         & ztaua_f(1:kproma,1:nb_sw) 
     paer_piz_sw_vr(1:kproma,jk,1:nb_sw) = &
                   & ztaua_c(1:kproma,1:nb_sw)*zs_c(1:kproma,1:nb_sw) + &
                   & ztaua_f(1:kproma,1:nb_sw)*zs_f(1:kproma,1:nb_sw)
     WHERE (paer_tau_sw_vr(1:kproma,jk,1:nb_sw) /= 0._wp) 
        paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=paer_piz_sw_vr(1:kproma,jk,1:nb_sw)&
                                           /paer_tau_sw_vr(1:kproma,jk,1:nb_sw)
     ELSEWHERE
        paer_piz_sw_vr(1:kproma,jk,1:nb_sw)=1._wp
     END WHERE
     paer_cg_sw_vr(1:kproma,jk,1:nb_sw)  = &
     &ztaua_c(1:kproma,1:nb_sw)*zs_c(1:kproma,1:nb_sw)*zg_c(1:kproma,1:nb_sw)+&
     &ztaua_f(1:kproma,1:nb_sw)*zs_f(1:kproma,1:nb_sw)*zg_f(1:kproma,1:nb_sw)
     WHERE (paer_tau_sw_vr(1:kproma,jk,1:nb_sw) /= 0._wp) 
        paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=paer_cg_sw_vr(1:kproma,jk,1:nb_sw)/&
                                          paer_piz_sw_vr(1:kproma,jk,1:nb_sw)/&
                                          paer_tau_sw_vr(1:kproma,jk,1:nb_sw)
     ELSEWHERE
        paer_cg_sw_vr(1:kproma,jk,1:nb_sw)=0._wp
     END WHERE
  ENDDO
!  WRITE(0,*) paer_tau_sw_vr(1:kproma,1,5)
END SUBROUTINE set_bc_aeropt_kinne
!-------------------------------------------------------------------------
! 
!> SUBROUTINE read_months_bc_aeropt_kinne -- reads optical aerosol parameters from file containing
!! aod, ssa, asy, aer_ex (altitude dependent extinction), dz_clim (layer 
!! thickness in meters), lev_clim (number of levels), and (optional) surface 
!! altitude in meters.
!!
SUBROUTINE read_months_bc_aeropt_kinne (                                   &
  caod,             cssa,             casy,               caer_ex,         &
  cdz_clim,         cwldim,           clevdim,            imnthb,          &
  imnthe,           iyear,            cfname,             p_patch          )
!
  CHARACTER(len=*), INTENT(in)   :: caod,    &! name of variable containing optical depth of column
                                    cssa,    &! name of variable containing single scattering albedo 
                                    casy,    &! name of variable containing asymmetry factor
                                              ! ssa and asy are assumed to be constant over column
                                    caer_ex, &! name of variable containing altitude dependent extinction
                                              ! aer_ex is normed to 1 (total over column is equal to 1)
                                    cdz_clim,&! layer thickness of climatology in meters
                                    cwldim,  &! name of wavelength dimension
                                    clevdim   ! name of level dimension in climatology
  INTEGER, INTENT(in)            :: imnthb, imnthe !begin and end month to be read
  INTEGER(i8), INTENT(in)        :: iyear ! base year. if month=0, month 12 of previous year is read, if month=13, month 1
                                          ! of subsequent year is read
  CHARACTER(len=*), INTENT(in)   :: cfname   ! file name containing variables
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch

  INTEGER                        :: ifile_id, kmonthb, kmonthe, ilen_cfname
  TYPE(t_stream_id)              :: stream_id
  REAL(wp), POINTER              :: zvar(:,:,:,:)
  REAL(wp), POINTER              :: zaod(:,:,:,:), zssa(:,:,:,:), zasy(:,:,:,:), zaer_ex(:,:,:,:)
  CHARACTER(LEN=32)              :: cimnthb, cimnthe
  CHARACTER(LEN=512)             :: cfnameyear,cyear

  IF (imnthb < 0 .OR. imnthe < imnthb .OR. imnthe > 13 ) THEN
    CALL finish ('read_months_bc_aeropt_kinne in mo_bc_aeropt_kinne', &
                 'months to be read outside valid range 0<=imnthb<=imnthe<=13, '// &
                 'imnthb='//TRIM(ADJUSTL(cimnthb))//', imnthe='//TRIM(ADJUSTL(cimnthe))) 
  END IF
  ilen_cfname=LEN_TRIM(cfname)
  IF (cfname(1:ilen_cfname) == 'bc_aeropt_kinne_sw_b14_coa') THEN
    zaod=>aod_c_s
    zssa=>ssa_c_s
    zasy=>asy_c_s
    zaer_ex=>z_km_aer_c_mo
  END IF
  IF (cfname(1:ilen_cfname) == 'bc_aeropt_kinne_lw_b16_coa') THEN
    zaod=>aod_c_f
    zssa=>ssa_c_f
    zasy=>asy_c_f
    zaer_ex=>z_km_aer_c_mo ! for the coarse mode, the altitude distribution is wavelength independent and
                           ! therefore for solar and long wave spectrum the same
  END IF
  IF (cfname(1:ilen_cfname) == 'bc_aeropt_kinne_sw_b14_fin') THEN
    zaod=>aod_f_s
    zssa=>ssa_f_s
    zasy=>asy_f_s
    zaer_ex=>z_km_aer_f_mo
  END IF
  IF (imnthb == 0) THEN
    WRITE(cyear,*) iyear-1

    IF ( echam_phy_config%lamip ) THEN
      cfnameyear=cfname//'_'//TRIM(ADJUSTL(cyear))//'.nc'
    ELSE
      cfnameyear=cfname//'.nc'
    ENDIF

    CALL message ('read_months_bc_aeropt_kinne of mo_bc_aeropt_kinne', &
   &              'reading from file '//TRIM(ADJUSTL(cfnameyear)))
    stream_id=openInputFile(cfnameyear, p_patch, default_read_method)
!    IF (ALLOCATED(zvar)) DEALLOCATE(zvar)
    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=caod, &
           &          return_pointer=zvar, start_timestep=12, end_timestep=12, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zaod(:,:,:,0:0)),SHAPE(zvar),cfnameyear,caod, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zaod(:,:,:,0)=zvar(:,:,:,1)
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=cssa, &
           &          return_pointer=zvar, start_timestep=12, end_timestep=12, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zssa(:,:,:,0:0)),SHAPE(zvar),cfnameyear,cssa, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zssa(:,:,:,0)=zvar(:,:,:,1)
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=casy, &
           &          return_pointer=zvar, start_timestep=12, end_timestep=12, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zasy(:,:,:,0:0)),SHAPE(zvar),cfnameyear,casy, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zasy(:,:,:,0)=zvar(:,:,:,1)
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=caer_ex, &
           &          return_pointer=zvar, start_timestep=12, end_timestep=12, &
           &          levelsDimName=clevdim)
    CALL shape_check_fields(SHAPE(zaer_ex(:,:,:,0:0)),SHAPE(zvar),cfnameyear,caer_ex, &
                                 'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zaer_ex(:,:,:,0)=zvar(:,:,:,1)
    DEALLOCATE(zvar)
    CALL closeFile(stream_id)
  END IF
  IF (imnthe > 0) THEN
    WRITE(cyear,*) iyear

    IF ( echam_phy_config%lamip ) THEN
      cfnameyear=cfname//'_'//TRIM(ADJUSTL(cyear))//'.nc'
    ELSE
      cfnameyear=cfname//'.nc'
    ENDIF

    stream_id=openInputFile(cfnameyear, p_patch, default_read_method)
    kmonthb=MAX(1,imnthb)
    kmonthe=MIN(12,imnthe)
!    IF (ALLOCATED(zvar)) DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=caod, &
           &          return_pointer=zvar, start_timestep=kmonthb, end_timestep=kmonthe, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zaod(:,:,:,kmonthb:kmonthe)),SHAPE(zvar),cfnameyear,caod, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zaod(:,:,:,kmonthb:kmonthe)=zvar
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=cssa, &
           &          return_pointer=zvar, start_timestep=kmonthb, end_timestep=kmonthe, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zssa(:,:,:,kmonthb:kmonthe)),SHAPE(zvar),cfnameyear,cssa, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zssa(:,:,:,kmonthb:kmonthe)=zvar
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=casy, &
           &          return_pointer=zvar, start_timestep=kmonthb, end_timestep=kmonthe, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zasy(:,:,:,kmonthb:kmonthe)),SHAPE(zvar),cfnameyear,casy, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zasy(:,:,:,kmonthb:kmonthe)=zvar
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=caer_ex, &
           &          return_pointer=zvar, start_timestep=kmonthb, end_timestep=kmonthe, &
           &          levelsDimName=clevdim)
    CALL shape_check_fields(SHAPE(zaer_ex(:,:,:,kmonthb:kmonthe)),SHAPE(zvar),cfnameyear,caer_ex, &
                                 'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zaer_ex(:,:,:,kmonthb:kmonthe)=zvar
    DEALLOCATE(zvar)
    CALL closeFile(stream_id)
  END IF
  IF (imnthe == 13) THEN
    WRITE(cyear,*) iyear+1

    IF ( echam_phy_config%lamip ) THEN
      cfnameyear=cfname//'_'//TRIM(ADJUSTL(cyear))//'.nc'
    ELSE
      cfnameyear=cfname//'.nc'
    ENDIF

    stream_id=openInputFile(cfnameyear, p_patch, default_read_method)
!    IF (ALLOCATED(zvar)) DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=caod, &
           &          return_pointer=zvar, start_timestep=1, end_timestep=1, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zaod(:,:,:,13:13)),SHAPE(zvar),cfnameyear,caod, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')    
    zaod(:,:,:,13)=zvar(:,:,:,1)
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=cssa, &
           &          return_pointer=zvar, start_timestep=1, end_timestep=1, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zssa(:,:,:,13:13)),SHAPE(zvar),cfnameyear,cssa, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zssa(:,:,:,13)=zvar(:,:,:,1)
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=casy, &
           &          return_pointer=zvar, start_timestep=1, end_timestep=1, &
           &          levelsDimName=cwldim)
    CALL shape_check_fields(SHAPE(zasy(:,:,:,13:13)),SHAPE(zvar),cfnameyear,casy, &
                                  'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zasy(:,:,:,13)=zvar(:,:,:,1)
    DEALLOCATE(zvar)

    CALL read_3D_time(stream_id=stream_id, location=on_cells, variable_name=caer_ex, &
           &          return_pointer=zvar, start_timestep=1, end_timestep=1, &
           &          levelsDimName=clevdim)
    CALL shape_check_fields(SHAPE(zaer_ex(:,:,:,13:13)),SHAPE(zvar),cfnameyear,caer_ex, &
                                 'read_months_bc_aeropt_kinne','mo_bc_aeropt_kinne')
    zaer_ex(:,:,:,13)=zvar(:,:,:,1)
    DEALLOCATE(zvar)
    CALL closeFile(stream_id)
    ifile_id = openInputFile(cfnameyear)
    dz_clim = read_0D_real (file_id=ifile_id, variable_name=cdz_clim)
    CALL closeFile(ifile_id)
  END IF
  END SUBROUTINE read_months_bc_aeropt_kinne
!-------------------------------------------------------------------------
! 
!> SUBROUTINE size_check_zerofields -- checks the shape of the shape of the 
!! fields read into the icon program

  SUBROUTINE shape_check_fields(kdim_icon,kdim_file,cfname,cvarname,croutine_name,cmodule_name)
    INTEGER,INTENT(in)            :: kdim_icon(:), kdim_file(:)
    CHARACTER(LEN=*), INTENT(in)  :: cfname, cvarname, croutine_name, cmodule_name
    INTEGER                       :: idim
    CHARACTER(LEN=2)              :: cidim 
    CHARACTER(LEN=32)             :: cidim_len_file, cidim_len_icon
    IF (SIZE(kdim_icon) /= SIZE(kdim_file)) THEN
      CALL finish(TRIM(ADJUSTL(croutine_name))//' of '//TRIM(ADJUSTL(cmodule_name )), &
                  'variable '//TRIM(ADJUSTL(cvarname))//' has wrong number of dimensions in file ' &
                  //TRIM(ADJUSTL(cfname)))
    END IF
    DO idim=1,SIZE(kdim_icon)
      IF (kdim_icon(idim) /= kdim_file(idim)) THEN
        WRITE(cidim,'(i2)') idim
        WRITE(cidim_len_icon,'(i32)') kdim_icon(idim)
        WRITE(cidim_len_file,'(i32)') kdim_file(idim)
        CALL finish(TRIM(ADJUSTL(croutine_name))//' of '//TRIM(ADJUSTL(cmodule_name )), &
                  'variable '//TRIM(ADJUSTL(cvarname))//' has wrong length in dimension ' &
                  //TRIM(ADJUSTL(cidim))//' length in icon model: '//TRIM(ADJUSTL(cidim_len_icon)) &
                  //' but length in file '//TRIM(ADJUSTL(cfname))//' is '//TRIM(ADJUSTL(cidim_len_file)))
      END IF
    END DO
  END SUBROUTINE shape_check_fields
END MODULE mo_bc_aeropt_kinne
