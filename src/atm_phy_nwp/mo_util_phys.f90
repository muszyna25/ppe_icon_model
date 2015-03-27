!>
!! Implementation of physics utility routines.
!!
!! @par Revision History
!!  Initial revision  :  F. Prill, DWD (2012-07-03)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_util_phys

  USE mo_kind,                  ONLY: wp
  USE mo_parallel_config,       ONLY: nproma
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 cpd, p0ref, rd,  &
    &                                 vtmpc1, t3,      &
    &                                 grav,            &
    &                                 tmelt
  USE mo_satad,                 ONLY: sat_pres_water, sat_pres_ice
  USE mo_fortran_tools,         ONLY: assign_if_present
  USE mo_impl_constants,        ONLY: min_rlcell_int
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs, nqtendphy, lart
  USE mo_ls_forcing_nml,        ONLY: is_ls_forcing
  USE mo_loopindices,           ONLY: get_indices_c
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_art_config,            ONLY: art_config
  USE mo_initicon_config,       ONLY: is_iau_active, iau_wgt_adv
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_satad,                 ONLY: qsat_rho

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_dyn_gust
  PUBLIC :: nwp_con_gust
  PUBLIC :: virtual_temp
  PUBLIC :: vap_pres
  PUBLIC :: calsnowlmt
  PUBLIC :: swdir_s
  PUBLIC :: rel_hum
  PUBLIC :: compute_field_rel_hum_wmo
  PUBLIC :: compute_field_rel_hum_ifs
  PUBLIC :: compute_field_omega
  PUBLIC :: nh_update_prog_phy

CONTAINS


  !-------------------------------------------------------------------------
  !!
  !! Calculate dynamic gusts as in near_surface of COSMO code
  !!     gust = ff10m + gust_factor * ustar
  !! where ff10m is the 10 m wind and the friction velocity ustar = SQRT(tcm)*ff1
  !!
  !! @par Revision History
  !! Developed by Helmut Frank, DWD (2013-03-13)
  !!
  ELEMENTAL FUNCTION nwp_dyn_gust( u_10m, v_10m, tcm, u1, v1, u_env, v_env) RESULT( vgust_dyn)

    REAL(wp), INTENT(IN) :: u_10m, &    ! zonal wind component at 10 m above ground [m/s]
      &                     v_10m, &    ! meridional wind component at 10 m above ground [m/s]
      &                     tcm  , &    ! transfer coefficient for momentum at surface
      &                     u1   , &    ! zonal wind at lowest model layer above ground [m/s]
      &                     v1   , &    ! meridional wind at lowest model layer above ground [m/s]
      &                     u_env, &    ! zonal wind at top of SSO envelope layer [m/s]
      &                     v_env       ! meridional wind at top of SSO envelope layer [m/s]

    REAL(wp) :: vgust_dyn               ! dynamic gust at 10 m above ground [m/s]

    REAL(wp) :: ff10m, ustar, uadd_sso
    REAL(wp), PARAMETER :: gust_factor = 8.0_wp ! previously 3.0_wp * 2.4_wp

    ff10m = SQRT( u_10m**2 + v_10m**2)
    uadd_sso = MAX(0._wp, SQRT(u_env**2 + v_env**2) - SQRT(u1**2 + v1**2))
    ustar = SQRT( MAX( tcm, 5.e-4_wp) * ( u1**2 + v1**2) )
    vgust_dyn = ff10m + uadd_sso + gust_factor*ustar

  END FUNCTION nwp_dyn_gust



  !-------------------------------------------------------------------------
  !!
  !! Calculate convective contribution to the wind gusts
  !!     gust_conv = \alpha MAX(0,U_850 - U_950)
  !! where \alpha=0.6 is a tunable constant and U_850-U_950 is the difference between
  !! the 850 hPa and 950 hPa wind speeds, which represents the low-level wind shear.
  !!
  !! @par Literature
  !! Bechthold, P. and J. Bidlot (2009): Parameterization of convective gusts. 
  !! ECMWF Newsletter No. 119
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-03-25)
  !!
  ELEMENTAL FUNCTION nwp_con_gust( u_850, u_950, v_850, v_950) RESULT(vgust_con)

    REAL(wp), INTENT(IN) :: u_850, &    ! zonal wind component at 850 hPa [m/s]
      &                     u_950, &    ! zonal wind component at 950 hPa [m/s]
      &                     v_850, &    ! meridional wind component at 850 hPa [m/s]
      &                     v_950       ! meridional wind component at 950 hPa [m/s]

    REAL(wp) :: vgust_con               ! convective contribution to the wind gusts [m/s]

    REAL(wp), PARAMETER :: alpha = 0.6_wp ! convective mixing parameter

    vgust_con = alpha * MAX(0._wp, SQRT((u_850**2 + v_850**2)) - SQRT((u_950**2 + v_950**2)))

  END FUNCTION nwp_con_gust



  !-------------
  !>
  !! SUBROUTINE virtual_temp
  !! Computes virtual temperature
  !!
  !! Required input fields: temperature, specific humidity, cloud and precipitation variables
  !! Output: virtual temperature
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  !!
  SUBROUTINE virtual_temp(p_patch, temp, qv, qc, qi, qr, qs, qg, qh, temp_v)


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN)                   :: temp(:,:,:) ! temperature (K)
    REAL(wp), INTENT(IN)                   :: qv  (:,:,:) ! specific humidity
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qc  (:,:,:) ! specific cloud water
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qi  (:,:,:) ! specific cloud ice
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qr  (:,:,:) ! specific rain water
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qs  (:,:,:) ! specific snow
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qg  (:,:,:) ! specific graupel
    REAL(wp), INTENT(IN), OPTIONAL, TARGET :: qh  (:,:,:) ! specific hail

    REAL(wp), INTENT(OUT) :: temp_v(:,:,:) ! virtual temperature (K)

    INTEGER :: jb, jk, jc, jt
    INTEGER :: nlen, nlev
    INTEGER :: num_qcpvars ! number of cloud or precipitation variables
    REAL(wp):: z_qsum(nproma,SIZE(temp,2))

    TYPE t_fieldptr
      REAL(wp), POINTER :: fld(:,:,:)
    END TYPE t_fieldptr
    TYPE(t_fieldptr) :: qptr(6)

    nlev = SIZE(temp,2) ! in order to be usable for input and output data

    num_qcpvars = 0
    IF (PRESENT(qc)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qc
    ENDIF
    IF (PRESENT(qr)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qr
    ENDIF
    IF (PRESENT(qi)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qi
    ENDIF
    IF (PRESENT(qs)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qs
    ENDIF
    IF (PRESENT(qg)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qg
    ENDIF
    IF (PRESENT(qh)) THEN
      num_qcpvars = num_qcpvars + 1
      qptr(num_qcpvars)%fld => qh
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jt,z_qsum) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF
      
      z_qsum(:,:) = 0._wp
      IF (num_qcpvars > 0) THEN
        DO jt = 1, num_qcpvars
          DO jk = 1, nlev
            DO jc = 1, nlen
              z_qsum(jc,jk) = z_qsum(jc,jk) + qptr(jt)%fld(jc,jk,jb)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          temp_v(jc,jk,jb) = temp(jc,jk,jb) * (1._wp + vtmpc1*qv(jc,jk,jb) - z_qsum(jc,jk))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE virtual_temp



  !> POINTWISE computation of shortwave direct downward flux as
  !! swdir_s = (1+ (\alpha/(1-\alpha))) * sobs - swdifd_s
  !! in W m**-2
  !!
  !! (domain independent and elemental)
  !!
  !! @par Revision History
  !! Initial revision by D. Reinert, DWD (2014-09-18) 
  ELEMENTAL FUNCTION swdir_s(albedo, swdifd_s, sobs)
    REAL(wp)             :: swdir_s
    REAL(wp), INTENT(IN) :: albedo      ! shortwave broadband albedo
    REAL(wp), INTENT(IN) :: swdifd_s    ! shortwave diffuse downward flux (sfc)
    REAL(wp), INTENT(IN) :: sobs        ! shortwave net flux (sfc)

    swdir_s = (1._wp + albedo/(1._wp - albedo)) * sobs - swdifd_s

  END FUNCTION swdir_s 


  !> POINTWISE computation of relative humidity as r=100. * e/e_sat,
  !! according to WMO standard
  !!
  !! (domain independent and elemental)
  !!
  !! @par Revision History
  !! Initial revision  :  F. Prill, DWD (2012-07-03) 
  ELEMENTAL FUNCTION rel_hum(temp, qv, p_ex)
    REAL(wp) :: rel_hum
    REAL(wp), INTENT(IN) :: temp, &  ! temperature
      &                     qv,   &  ! spec. water vapor content
      &                     p_ex     ! exner pressure
    ! local variables
    REAL(wp) :: pres, e_s, e

    ! compute dynamic pressure from Exner pressure:
    pres = p0ref * EXP((cpd/rd)*LOG(p_ex))
    ! approx. saturation vapor pressure:
    e_s = sat_pres_water(temp)
    ! compute vapor pressure from formula for specific humidity:
    e   = pres*qv / (rdv + o_m_rdv*qv)

    rel_hum = 100._wp * e/e_s

  END FUNCTION rel_hum


  !> POINTWISE computation of relative humidity as r=100. * e/e_sat, 
  !! according to IFS documentation
  !! I.e. For the temperature range 250.16<=T<=273.16, the saturation 
  !! vapour pressure is computed as a combination of the values over 
  !! water e_s_water and over ice e_s_ice.
  !!
  !! (domain independent and elemental)
  !!
  !! @par Revision History
  !! Initial revision  by Daniel Reinert, DWD (2013-07-15) 
  ELEMENTAL FUNCTION rel_hum_ifs(temp, qv, p_ex)
    REAL(wp) :: rel_hum_ifs
    REAL(wp), INTENT(IN) :: temp, &  ! temperature
      &                     qv,   &  ! spec. water vapor content
      &                     p_ex     ! exner pressure
    ! local variables
    REAL(wp) :: pres, e_s, e_s_water, e_s_ice, e
    REAL(wp), PARAMETER:: t_i = 250.16_wp  ! threshold value for mixed-phase clouds

    ! compute dynamic pressure from Exner pressure:
    pres = p0ref * EXP((cpd/rd)*LOG(p_ex))
    ! approx. saturation vapor pressure:
    IF (temp > t3) THEN
      e_s       = sat_pres_water(temp)
    ELSE IF (temp < (t3-23._wp)) THEN
      e_s       = sat_pres_ice(temp)
    ELSE
      e_s_water = sat_pres_water(temp)
      e_s_ice   = sat_pres_ice(temp)

      e_s       = e_s_ice + (e_s_water - e_s_ice) * ((temp - t_i)/(t3 - t_i))**2
    ENDIF

    ! compute vapor pressure from formula for specific humidity:
    e   = pres*qv / (rdv + o_m_rdv*qv)

    rel_hum_ifs = 100._wp * e/e_s

  END FUNCTION rel_hum_ifs



  !> computation of relative humidity as r=e/e_sat, according to WMO standard
  !!
  !! @par Revision History
  !! Initial revision  :  F. Prill, DWD (2012-07-04) 
  SUBROUTINE compute_field_rel_hum_wmo(ptr_patch, p_prog, p_diag, out_var, &
    &                              opt_slev, opt_elev, opt_rlstart, opt_rlend)

    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN) :: p_prog
    TYPE(t_nh_diag), INTENT(IN) :: p_diag
    ! output variable, dim: (nproma,nlev,nblks_c):
    REAL(wp),INTENT(INOUT) :: out_var(:,:,:)
    ! optional vertical start/end level:
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev, opt_elev
    ! start and end values of refin_ctrl flag:
    INTEGER, INTENT(in), OPTIONAL     :: opt_rlstart, opt_rlend
   
    ! local variables
    REAL(wp) :: temp, qv, p_ex
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, jk, jb

    ! default values
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(slev,     opt_slev)
    CALL assign_if_present(elev,     opt_elev)
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)
    ! values for the blocking
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,temp,qv,p_ex), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)
      
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
!CDIR NEXPAND
          out_var(jc,jk,jb) = rel_hum(temp, qv, p_ex)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_rel_hum_wmo



  !> computation of relative humidity as r=e/e_sat, according to IFS
  !!
  !! @par Revision History
  !! Initial revision  :  F. Prill, DWD (2012-07-04) 
  SUBROUTINE compute_field_rel_hum_ifs(ptr_patch, p_prog, p_diag, out_var, &
    &                              opt_lclip, opt_slev, opt_elev,          &
    &                              opt_rlstart, opt_rlend)

    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN) :: p_prog
    TYPE(t_nh_diag), INTENT(IN) :: p_diag
    ! output variable, dim: (nproma,nlev,nblks_c):
    REAL(wp),INTENT(INOUT) :: out_var(:,:,:)
    ! optional clipping to rh<=100%
    LOGICAL, INTENT(IN), OPTIONAL     :: opt_lclip
    ! optional vertical start/end level:
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev, opt_elev
    ! start and end values of refin_ctrl flag:
    INTEGER, INTENT(in), OPTIONAL     :: opt_rlstart, opt_rlend
   
    ! local variables
    REAL(wp) :: temp, qv, p_ex
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, jk, jb
    LOGICAL  :: lclip       ! clip rel. hum. to values <=100% 


    IF (PRESENT(opt_lclip)) THEN
      lclip = opt_lclip
    ELSE
      lclip = .FALSE.
    ENDIF


    ! default values
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(slev,     opt_slev)
    CALL assign_if_present(elev,     opt_elev)
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)
    ! values for the blocking
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,temp,qv,p_ex), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)
      
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=3
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
!CDIR NEXPAND
          out_var(jc,jk,jb) = rel_hum_ifs(temp, qv, p_ex)

          ! optional clipping, if lclip=.TRUE.
          out_var(jc,jk,jb) = MERGE(MIN(100._wp,out_var(jc,jk,jb)), out_var(jc,jk,jb), lclip)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_rel_hum_ifs



  !> computation of water vapour pressure
  !!
  !! water vapour pressure is computed as a function of specific humidity 
  !! qv and atmospheric pressure pres.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-07-25) 
  ELEMENTAL FUNCTION vap_pres(qv,pres)
  IMPLICIT NONE

    REAL(wp), INTENT(IN)  :: qv   ! specific humidity         [kg/kg]
    REAL(wp), INTENT(IN)  :: pres ! atmospheric pressure      [Pa]
    REAL(wp) :: vap_pres          ! water vapour pressure     [Pa]

    vap_pres = (qv * pres) / (rdv + O_m_rdv*qv)

  END FUNCTION vap_pres




  !------------------------------------------------------------------------------
  !>
  !! Description:
  !!   This subroutine calculates height of the snowfall limit (snowlmt).
  !!
  !! Method:
  !!   In a first step the wet bulb temperature is derived from pres, t and qv.
  !!   In a second step the snowfall limit is evaluated from 8000m down to the
  !!   the lowest model level (ke) and linearly interpolated to the height where
  !!   the wet bulb temperature is >= wbl (=+1.3C after P. Haechler, MeteoSwiss).
  !!   A flag (-999) is set to indicate that no snowlmt was found.
  !!
  !! @par Revision History
  !! Inherited from COSMO 5.0 by Daniel Reinert, DWD (2015-03-27)
  !! 
  !!
  SUBROUTINE calsnowlmt ( snowlmt, temp, pres, qv, hhl, hhlr, istart, iend, wbl)

    ! Parameter list:

    INTEGER, INTENT (IN)     ::  &
      istart, iend           ! loop start/end indices

    REAL(wp), INTENT (INOUT)   ::  &
      snowlmt(:)    ! height of the snowfall limit in m above sea level

    REAL(wp), INTENT (IN)    ::  &
      temp  (:,:), & ! temperature
      pres  (:,:), & ! pressure at full levels
      qv    (:,:), & ! specific humidity
      hhl   (:,:), & ! height of model half levels
      hhlr  (:)      ! height of model half levels resp. sea level

    REAL (wp), INTENT (IN)    ::  &
      wbl               ! (empirical) wet bulb temperature at snowfall limit (1.3C)
    !------------------------------------------------------------------------------
    ! Local variables

    INTEGER ::     i, k, ktopmin, nlev

    LOGICAL                  ::    &
      lfound(SIZE(temp,1))     ! Logical flag : =.TRUE when wet bulb temp corresponding to
                     !                  parameter "wbl" is found

    REAL (wp)       ::    &
      za = 0.78588481_wp,      & ! local storage
      zb = 7.567_wp,           &
      zc = 2066.92605_wp,      &
      zd = 33.45_wp,           &
      ze = 0.622_wp,           &
      zf = 0.378_wp,           &
      zg = 0.5_wp,             &
      zh = 0.6_wp,             &
      zi = 700._wp,            &
      zl = 0.1_wp,             &
      zm = 6400._wp,           &
      zn = 11.564_wp,          &
      zo = 1742._wp,           &
      td,tl,tp,ppp,            &
      zp,                      &  ! pressure in hPa
      ppp,                     &  ! pressure in dPa
      deltat,zt,               &
      ep,const,                &
      zh_bot, zh_top,          &
      zdt

    REAL(wp) :: wetblb(SIZE(temp,1),SIZE(temp,2))  ! wet-bulb temperature in Celsius

  !------------------------------------------------------------------------------

    ! Begin subroutine calsnowlmt

    ! number of vertical full levels
    nlev = SIZE(temp,2)

    ! Set the uppermost model level for the occurence of a wet bulb temperature (wbl)
    ! to about 8000m above surface
    ktopmin = 2
    DO k = nlev+1, 1, -1
      IF ( hhlr(k) < 8000.0_wp ) THEN
        ktopmin = k
      ENDIF
    ENDDO

    ! Initialize the definition mask and the output array snowlmt
    lfound (:) = .FALSE.
    snowlmt(:) = -999.0_wp

    DO k = ktopmin, nlev
      DO i = istart, iend
        zp     = (pres(i,k))/100._wp     ! in hPa
        ep     = MAX(1.0E-10_wp,qv(i,k))*zp /      &
                 (ze + zf*MAX(1.0E-10_wp,qv(i,k)))
        ep     = MAX(ep,1.0E-10_wp)
        CONST  = LOG10(ep) - za
        td     = (zd*CONST-zc) / (CONST-zb)              ! in Kelvin
        ! Wet bulb temperature after Egger/Joss
        tl     = (temp(i,k) - tmelt) *10._wp
        tp     = (td-tmelt) *10._wp
        ppp    = zp * 10._wp
        deltat = tl-tp
        zt     = tp + zg*deltat*(zh-tp/zi)
        wetblb(i,k) = zl * ( tp +                      & ! in Celsius
                      (deltat / (1._wp + zm*EXP(zn*zt/(zo+zt))/ppp)))

        IF ( wetblb(i,k) >= wbl ) THEN
          ! definition of snowlmt can be made in this column
          lfound (i) = .TRUE.
        ENDIF
      ENDDO
    ENDDO

    DO k = ktopmin+1, nlev
      DO i = istart, iend
        IF ( lfound(i) .AND. wetblb(i,k) >= wbl ) THEN
          ! definition of snowlmt is now made once
          lfound (i) = .FALSE.
          zh_bot     = 0.5_wp * ( hhl(i,k) + hhl(i,k+1) )
          zh_top     = 0.5_wp * ( hhl(i,k) + hhl(i,k-1) )
          zdt        = ( wbl - wetblb(i,k) ) /                 &
                       ( wetblb(i,k-1) - wetblb(i,k) )
          snowlmt(i) = zh_bot + (zh_top-zh_bot)*zdt
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE calsnowlmt






  !> computation of vertical velocity (dp/dt)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-03-28) 
  SUBROUTINE compute_field_omega(ptr_patch, p_prog, out_var, &
    &                            opt_slev, opt_elev, opt_rlstart, opt_rlend)

    TYPE(t_patch)        , INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    TYPE(t_nh_prog)      , INTENT(IN)    :: p_prog                 !< nonhydrostatic state
    REAL(wp)             , INTENT(INOUT) :: out_var(:,:,:)         !< output variable, dim: (nproma,nlev,nblks_c)
    INTEGER, INTENT(IN), OPTIONAL        :: opt_slev, opt_elev     !< optional vertical start/end level
    INTEGER, INTENT(IN), OPTIONAL        :: opt_rlstart, opt_rlend !< start and end values of refin_ctrl flag


    ! local
    REAL(wp):: w_avg(nproma)       ! vertical velocity averaged to full level
    INTEGER :: slev, elev          ! vertical start and end index
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_nchdom
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc, jk, jb  

    ! default values
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(slev,     opt_slev)
    CALL assign_if_present(elev,     opt_elev)
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)

    ! values for the blocking
    i_nchdom   = MAX(1,ptr_patch%n_childdom)
    i_startblk = ptr_patch%cells%start_blk(rl_start,1)
    i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx,w_avg), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, rl_start, rl_end)
      
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          ! half level to full level interpolation
          w_avg(jc) = 0.5_wp * (p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb))

          out_var(jc,jk,jb) = -p_prog%rho(jc,jk,jb)*grav*w_avg(jc)

        ENDDO
      ENDDO

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE compute_field_omega

  SUBROUTINE nh_update_prog_phy( pt_patch, pdtime, pt_diag, prm_nwp_tend, &
    &                            prm_diag, pt_prog_rcf, pt_prog )

    TYPE(t_patch),       INTENT(IN)   :: pt_patch     !!grid/patch info.
    TYPE(t_nh_diag)     ,INTENT(IN)   :: pt_diag      !<the diagnostic variables
    TYPE(t_nwp_phy_tend),TARGET, INTENT(IN):: prm_nwp_tend   !< atm tend vars
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag     !!the physics variables
    TYPE(t_nh_prog),     INTENT(INOUT):: pt_prog_rcf  !!the tracer field at
                                                   !!reduced calling frequency
    TYPE(t_nh_prog),     INTENT(IN)   :: pt_prog   !! NH prog state at dynamic time step
    REAL(wp),INTENT(in)            :: pdtime

    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !! blocks
    INTEGER :: i_startidx, i_endidx    !! slices
    INTEGER :: i_nchdom                !! domain index

    ! Local scalars:
    INTEGER  :: nlev        !< number of full levels
    INTEGER  :: jb          !block index
    INTEGER  :: jt          !tracers
    INTEGER  :: jk,jc,jg
    REAL(wp) :: zqc, zqi, zqcn, zqin

    REAL(wp) :: zrhw(nproma, pt_patch%nlev) ! relative humidity w.r.t. water

    jg = pt_patch%id

    ! number of vertical levels
    nlev = pt_patch%nlev

    ! local variables related to the blocking

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jt,i_startidx,i_endidx,zqc,zqi,zqcn,zqin,zrhw) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
           &                       i_startidx, i_endidx, rl_start, rl_end)

      ! add analysis increments from data assimilation to qv
      !
      IF (is_iau_active) THEN
        ! Compute relative humidity w.r.t. water
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            zrhw(jc,jk) = pt_prog_rcf%tracer(jc,jk,jb,iqv)/qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
          ENDDO
        ENDDO

        ! DA increments of humidity are limited to positive values if RH < 2% or QV < 2.5e-6
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            IF (zrhw(jc,jk) < 0.02_wp .OR. pt_prog_rcf%tracer(jc,jk,jb,iqv) < 2.5e-6_wp) THEN
              zqin = MAX(0._wp, pt_diag%qv_incr(jc,jk,jb))
            ELSE
              zqin = pt_diag%qv_incr(jc,jk,jb)
            ENDIF
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv * zqin
          ENDDO
        ENDDO
      ENDIF


! KF fix to positive values
      DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          zqc = pt_prog_rcf%tracer(jc,jk,jb,iqc)+pdtime*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc)
          zqi = pt_prog_rcf%tracer(jc,jk,jb,iqi)+pdtime*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)

          zqcn = MIN(0._wp,zqc)
          zqin = MIN(0._wp,zqi)

          pt_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp, zqc)
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp, zqi)

          ! Subtract moisture generated by artificial clipping of QC and QI from QV
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,iqv) + zqcn+zqin &
            &                       + pdtime*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqv))

        ENDDO
      ENDDO

      IF(lart .AND. art_config(jg)%lart_conv) THEN
! KL add convective tendency and fix to positive values
        DO jt=1,art_config(jg)%nconv_tracer  ! ASH
          DO jk = 1, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              pt_prog_rcf%conv_tracer(jb,jt)%ptr(jc,jk)=MAX(0._wp,pt_prog_rcf%conv_tracer(jb,jt)%ptr(jc,jk) &
                 +pdtime*prm_nwp_tend%conv_tracer_tend(jb,jt)%ptr(jc,jk))
            ENDDO
          ENDDO
        ENDDO
      ENDIF !lart

!DR additional clipping for qr, qs 
!DR (very small negative values may occur during the transport process (order 10E-15)) 
      DO jt=iqr, iqs  ! qr,qs
        DO jk = kstart_moist(jg), nlev
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO

      IF (atm_phy_nwp_config(jg)%lcalc_acc_avg) THEN

!DIR$ IVDEP
        prm_diag%rain_con(i_startidx:i_endidx,jb) =                                       &
          &                                  prm_diag%rain_con(i_startidx:i_endidx,jb)    &
          &                                  + pdtime                                     &
          &                                  * prm_diag%rain_con_rate(i_startidx:i_endidx,jb)
!DIR$ IVDEP
        prm_diag%snow_con(i_startidx:i_endidx,jb) =                                       &
          &                                  prm_diag%snow_con(i_startidx:i_endidx,jb)    &
          &                                  + pdtime                                     &
          &                                  * prm_diag%snow_con_rate(i_startidx:i_endidx,jb)

        !for grid scale part: see mo_nwp_gscp_interface/nwp_microphysics
!DIR$ IVDEP
        prm_diag%tot_prec(i_startidx:i_endidx,jb) =                                       &
          &                              prm_diag%tot_prec(i_startidx:i_endidx,jb)        &
          &                              +  pdtime                                        &
          &                              * (prm_diag%rain_con_rate(i_startidx:i_endidx,jb)&
          &                              +  prm_diag%snow_con_rate(i_startidx:i_endidx,jb))

      ENDIF


      IF(is_ls_forcing)THEN
        DO jt=1, nqtendphy  ! qv,qc,qi
          DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt)    &
                &                       + pdtime*prm_nwp_tend%ddt_tracer_ls(jk,jt))
            ENDDO
          ENDDO
        END DO
      ENDIF  ! is_ls_forcing


    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE nh_update_prog_phy


END MODULE mo_util_phys
