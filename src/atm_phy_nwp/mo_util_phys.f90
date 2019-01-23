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

  USE mo_kind,                  ONLY: vp, wp
  USE mo_parallel_config,       ONLY: nproma
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 cpd, p0ref, rd,  &
    &                                 vtmpc1, t3,      &
    &                                 grav,            &
    &                                 tmelt, earth_radius, &
    &                                 alvdcp, rd_o_cpd
  USE mo_exception,             ONLY: finish
  USE mo_satad,                 ONLY: sat_pres_water, sat_pres_ice
  USE mo_fortran_tools,         ONLY: assign_if_present
  USE mo_impl_constants,        ONLY: min_rlcell_int, min_rledge_int, &
    &                                 min_rlcell, dzsoil
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqni, ininact, &
       &                              iqm_max, nqtendphy, lart
  USE mo_nh_diagnose_pres_temp, ONLY: diag_pres, diag_temp
  USE mo_ls_forcing_nml,        ONLY: is_ls_forcing
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_nwp_tuning_config,     ONLY: tune_gust_factor
  USE mo_advection_config,      ONLY: advection_config
  USE mo_art_config,            ONLY: art_config
  USE mo_initicon_config,       ONLY: iau_wgt_adv, qcana_mode, qiana_mode
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_lnd_nwp_config,        ONLY: nlev_soil
  USE mo_nwp_lnd_types,         ONLY: t_lnd_diag
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_satad,                 ONLY: qsat_rho
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_nwp_sfc_interp,        ONLY: wsoil2smi
  USE mo_icon_interpolation_scalar,                     &
    &                           ONLY: edges2cells_scalar, cells2edges_scalar, &
    &                                 cells2verts_scalar, verts2edges_scalar, &
    &                                 cells2edges_scalar
  USE mo_math_gradients,        ONLY: grad_fd_norm, grad_fd_tang
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_edge
  USE mo_sync,                  ONLY: sync_patch_array, SYNC_C

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
  PUBLIC :: compute_field_pv
  PUBLIC :: compute_field_smi
  PUBLIC :: iau_update_tracer
  PUBLIC :: tracer_add_phytend
  PUBLIC :: cal_cape_cin

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
  ELEMENTAL FUNCTION nwp_dyn_gust( u_10m, v_10m, tcm, u1, v1, u_env, v_env, mtnmask) RESULT( vgust_dyn)

    REAL(wp), INTENT(IN) :: u_10m, &    ! zonal wind component at 10 m above ground [m/s]
      &                     v_10m, &    ! meridional wind component at 10 m above ground [m/s]
      &                     tcm  , &    ! transfer coefficient for momentum at surface
      &                     u1   , &    ! zonal wind at lowest model layer above ground [m/s]
      &                     v1   , &    ! meridional wind at lowest model layer above ground [m/s]
      &                     u_env, &    ! zonal wind at top of SSO envelope layer [m/s]
      &                     v_env, &    ! meridional wind at top of SSO envelope layer [m/s]
      &                     mtnmask     ! mask field for weighting SSO enhancement

    REAL(wp) :: vgust_dyn               ! dynamic gust at 10 m above ground [m/s]

    REAL(wp) :: ff10m, ustar, uadd_sso, gust_add

    ff10m = SQRT( u_10m**2 + v_10m**2)
    uadd_sso = MAX(0._wp, SQRT(u_env**2 + v_env**2) - SQRT(u1**2 + v1**2))
    ustar = SQRT( MAX( tcm, 5.e-4_wp) * ( u1**2 + v1**2) )
    gust_add = MAX(0._wp,MIN(2._wp,0.2_wp*(ff10m-10._wp)))*(1._wp+mtnmask)
    vgust_dyn = ff10m + mtnmask*uadd_sso + (tune_gust_factor+gust_add+2._wp*mtnmask)*ustar

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
      td,tl,tp,                &
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



  !> computation of soil mositure index (smi)
  !!
  !! Conversion of soil moisture into soil moisture index
  !! smi = (soil moisture - wilting point) / (field capacity - wilting point)
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2017-05-03) 
  !!
  SUBROUTINE compute_field_smi(ptr_patch, diag_lnd, ext_data, out_var, &
    &                            opt_rlstart, opt_rlend)

    TYPE(t_patch)        , INTENT(IN)    :: ptr_patch              !< patch on which computation is performed
    TYPE(t_lnd_diag)     , INTENT(IN)    :: diag_lnd               !< nwp diag land state
    TYPE(t_external_data), INTENT(IN)    :: ext_data               !< ext_data state
    REAL(wp)             , INTENT(INOUT) :: out_var(:,:,:)         !< output variable, dim: (nproma,nlev,nblks_c)
    INTEGER, INTENT(IN), OPTIONAL        :: opt_rlstart, opt_rlend !< start and end values of refin_ctrl flag


    ! local
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk
    INTEGER :: jc, jk, jb, ic  
    INTEGER :: i_count
    INTEGER :: ierr, ierr_wsoil2smi

    ! default values
    rl_start = 2
    rl_end   = min_rlcell_int-1
    ! check optional arguments
    CALL assign_if_present(rl_start, opt_rlstart)
    CALL assign_if_present(rl_end,   opt_rlend)

    ! values for the blocking
    i_startblk = ptr_patch%cells%start_block(rl_start)
    i_endblk   = ptr_patch%cells%end_block(rl_end)

!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,ic,i_count,ierr,ierr_wsoil2smi), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      ierr = 0

      ! loop over target (ICON) land points only
      i_count = ext_data%atm%lp_count(jb)

      DO ic = 1, i_count
        jc = ext_data%atm%idx_lst_lp(ic,jb)

        DO jk = 1, nlev_soil-1

          CALL wsoil2smi(wsoil   = diag_lnd%w_so(jc,jk,jb),     & !in
            &            dzsoil  = dzsoil(jk),                  & !in
            &            soiltyp = ext_data%atm%soiltyp(jc,jb), & !in
            &            smi     = out_var(jc,jk,jb),           & !out
            &            ierr    = ierr_wsoil2smi               ) !out
          !
          ierr = MIN(ierr, ierr_wsoil2smi)
        ENDDO
        ! assume no-gradient condition for soil moisture reservoir layer
        out_var(jc,nlev_soil,jb) = out_var(jc,nlev_soil-1,jb)
      ENDDO

      IF (ierr < 0) THEN
        CALL finish("compute_field_smi", "Landpoint has invalid soiltype (sea water or sea ice)")
      ENDIF

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_smi



  !> Computation of potential vorticity
  !! The full 3D-Ertel PV is calculated at the edges and interpolated to cells.
  !! The shallow atmosphere approximations are used.
  !!
  !! Implemented by Tobias Selz, LMU
  
  SUBROUTINE compute_field_pv(p_patch, p_int_state, p_metrics, p_prog, p_diag, out_var )

    TYPE(t_patch)        , INTENT(INOUT) :: p_patch              !< patch on which computation is performed
    TYPE(t_int_state)    , INTENT(IN)    :: p_int_state
    TYPE(t_nh_metrics)   , INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog)      , INTENT(IN)    :: p_prog                 !< nonhydrostatic state
    TYPE(t_nh_diag)      , INTENT(IN)    :: p_diag
    REAL(wp)             , INTENT(INOUT) :: out_var(:,:,:)         !< output variable, dim: (nproma,nlev,nblks_c)
  
    !Local variables
    !Indices
    INTEGER  :: slev, elev, rl_start, rl_end, i_nchdom,     &
      &         i_startblk, i_endblk, i_startidx, i_endidx, &
      &         jc, je, jk, jb, ivd1, ivd2
      
    REAL(wp) :: vdfac
        
    !temporary fields
    REAL(wp) :: pv_ef    (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                vt       (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                theta_cf (nproma,p_patch%nlev  ,p_patch%nblks_c),  &
                theta_vf (nproma,p_patch%nlev  ,p_patch%nblks_v),  &
                theta_ef (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                w_vh     (nproma,p_patch%nlev+1,p_patch%nblks_v),  & 
                w_eh     (nproma,p_patch%nlev+1,p_patch%nblks_e),  &
                ddtw_eh  (nproma,p_patch%nlev+1,p_patch%nblks_e),  &
                ddnw_eh  (nproma,p_patch%nlev+1,p_patch%nblks_e),  &
                ddtth_ef (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                ddnth_ef (nproma,p_patch%nlev  ,p_patch%nblks_e),  &
                vor_ef   (nproma,p_patch%nlev  ,p_patch%nblks_e)
                
    !Pointers to metric terms
    REAL(vp), POINTER :: ddnz(:,:,:), ddtz(:,:,:), gamma(:,:,:)
    
    ddnz  => p_metrics%ddxn_z_full
    ddtz  => p_metrics%ddxt_z_full
    gamma => p_metrics%ddqz_z_full_e


    ! Index bounds
    slev     = 1
    elev     = UBOUND(out_var,2)
    rl_start = 1
    rl_end   = min_rlcell

    ! values for the blocking
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk (rl_start,1)
    i_endblk   = p_patch%cells%end_blk   (rl_end,i_nchdom)
    
!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    !compute theta on cells
    DO jb = i_startblk, i_endblk
    
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          theta_cf(jc,jk,jb) = p_diag%temp(jc,jk,jb) / p_prog%exner(jc,jk,jb)
          
        ENDDO
      ENDDO
    
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! synchronize theta
    CALL sync_patch_array(SYNC_C, p_patch, theta_cf)

    !Get vt at edges (p_diag%vt is not up to date)
    CALL rbf_vec_interpol_edge( p_prog%vn, p_patch, p_int_state, vt)
    
    !Interpolate theta to vertices
    CALL cells2verts_scalar( theta_cf, p_patch, p_int_state%cells_aw_verts, theta_vf )
    
    !Interpolate theta to edges
    CALL cells2edges_scalar( theta_cf, p_patch, p_int_state%c_lin_e, theta_ef )
    
    !Interpolate w to vertices
    CALL cells2verts_scalar( p_prog%w, p_patch, p_int_state%cells_aw_verts, w_vh )
    
    !Interpolate w to edges
    CALL cells2edges_scalar( p_prog%w, p_patch, p_int_state%c_lin_e, w_eh )
    
    !Interpolate vorticity to edges
    CALL verts2edges_scalar( p_diag%omega_z, p_patch, p_int_state%v_1o2_e, vor_ef )
    
    !Calculate horizontal derivatives of w and theta
    CALL grad_fd_norm ( p_prog%w, p_patch, ddnw_eh  )
    CALL grad_fd_tang ( w_vh,     p_patch, ddtw_eh  )
    CALL grad_fd_norm ( theta_cf, p_patch, ddnth_ef )
    CALL grad_fd_tang ( theta_vf, p_patch, ddtth_ef )
    
    !Recompute loop indices for edges
    rl_start   = 3
    rl_end     = min_rledge_int-1
    i_startblk = p_patch%edges%start_blk (rl_start,1)
    i_endblk   = p_patch%edges%end_blk   (rl_end,i_nchdom)
    
!$OMP PARALLEL    
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,ivd1,ivd2,vdfac), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
     
      DO jk = slev, elev
      
        !Get indices for vertical derivatives of full level variables
        IF ( jk == slev ) THEN
          ivd1=slev
          ivd2=slev+1
          vdfac=1_wp
        ELSE IF ( jk == elev ) THEN
          ivd1=elev-1
          ivd2=elev
          vdfac=1_wp
        ELSE
          ivd1=jk-1
          ivd2=jk+1
          vdfac=2_wp
        END IF
          
        DO je = i_startidx, i_endidx

          !Ertel-PV calculation on edges
          pv_ef(je,jk,jb) =                                                                                     &
            &     (   0.5_wp*(ddnw_eh(je,jk,jb)+ddnw_eh(je,jk+1,jb))                                            &
            &       + ddnz(je,jk,jb)/gamma(je,jk,jb)*(w_eh(je,jk+1,jb)-w_eh(je,jk,jb))                          &
            &       + (p_prog%vn(je,ivd2,jb)-p_prog%vn(je,ivd1,jb))/vdfac/gamma(je,jk,jb)                       &
            &       - p_prog%vn(je,jk,jb)/earth_radius                                                          &
            &     )                                                                                             &
            &   * (   ddtth_ef(je,jk,jb)                                                                        &
            &       + ddtz(je,jk,jb)/gamma(je,jk,jb) * (theta_ef(je,ivd2,jb)-theta_ef(je,ivd1,jb))/vdfac        &
            &     )                                                                                             &
            &   + ( - (vt(je,ivd2,jb)-vt(je,ivd1,jb))/vdfac/gamma(je,jk,jb)                                     &
            &       - 0.5_wp*(ddtw_eh(je,jk,jb)+ddtw_eh(je,jk+1,jb))                                            &
            &       - ddtz(je,jk,jb)/gamma(je,jk,jb) * (w_eh(je,jk+1,jb)-w_eh(je,jk,jb))                        &
            &       + vt(je,jk,jb)/earth_radius                                                                 &
            &      )                                                                                            &
            &   * (   ddnth_ef(je,jk,jb)                                                                        &
            &       + ddnz(je,jk,jb)/gamma(je,jk,jb) * (theta_ef(je,ivd2,jb)-theta_ef(je,ivd1,jb))/vdfac        &
            &     )                                                                                             &
            &   + (   vor_ef(je,jk,jb)                                                                          &
            &       + ddtz(je,jk,jb)/gamma(je,jk,jb) * (p_prog%vn(je,ivd2,jb)-p_prog%vn(je,ivd1,jb))/vdfac      &
            &       - ddnz(je,jk,jb)/gamma(je,jk,jb) * (vt(je,ivd2,jb)-vt(je,ivd1,jb))/vdfac                    &
            &       + p_patch%edges%f_e(je,jb)                                                                  &
            &     )                                                                                             &
            &   * ( -(theta_ef(je,ivd2,jb)-theta_ef(je,ivd1,jb))/vdfac/gamma(je,jk,jb) )
                   
        ENDDO
      ENDDO
    ENDDO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
    !Interpolate to cells
    CALL edges2cells_scalar( pv_ef, p_patch, p_int_state%e_bln_c_s, out_var, opt_rlstart=2 )
    

    rl_start = 2
    rl_end   = min_rlcell_int

    ! values for the blocking
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%cells%start_blk (rl_start,1)
    i_endblk   = p_patch%cells%end_blk   (rl_end,i_nchdom)

    !Normalize with density
    !
!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
    
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
          out_var(jc,jk,jb) = out_var(jc,jk,jb) / p_prog%rho(jc,jk,jb)
        ENDDO
      ENDDO
    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL    

        
  END SUBROUTINE compute_field_pv




  !
  ! Add IAU increment to qv during IAU phase
  !
  ! Add analysis increments from data assimilation to qv
  !
  ! Initial revision by Daniel Reinert, DWD (2018-05-18)
  ! Previously, this code snippet was part of nh_update_tracer_phy
  ! 
  SUBROUTINE iau_update_tracer( pt_prog, p_metrics, pt_diag, pt_prog_rcf, &
    &                     jg, jb, i_startidx, i_endidx, kend )

    TYPE(t_nh_prog)    ,INTENT(IN)   :: pt_prog      !< NH prog state at dynamic time step
    TYPE(t_nh_metrics) ,INTENT(IN)   :: p_metrics    !< NH metrics variables
    TYPE(t_nh_diag)    ,INTENT(INOUT):: pt_diag      !< the diagnostic variables
    TYPE(t_nh_prog)    ,INTENT(INOUT):: pt_prog_rcf  !< the tracer field at
                                                      !< reduced calling frequency
    INTEGER            ,INTENT(IN)   :: jg           !< domain ID
    INTEGER            ,INTENT(IN)   :: jb           !< block index
    INTEGER            ,INTENT(IN)   :: i_startidx   !< hor. start idx
    INTEGER            ,INTENT(IN)   :: i_endidx     !< hor. end idx
    INTEGER            ,INTENT(IN)   :: kend         !< vert. end idx

    ! Local variables
    INTEGER  :: jk,jc
    REAL(wp) :: zqin
    REAL(wp) :: zrhw(nproma, kend) ! relative humidity w.r.t. water


    ! add analysis increments from data assimilation to qv
    !
    ! Diagnose pressure and temperature for subsequent calculations
    CALL diag_temp (pt_prog, pt_prog_rcf, advection_config(jg)%trHydroMass%list, pt_diag, &
                    jb, i_startidx, i_endidx, 1, kstart_moist(jg), kend)
    CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, kend)

    ! Compute relative humidity w.r.t. water
    DO jk = 1, kend
      DO jc = i_startidx, i_endidx
        zrhw(jc,jk) = pt_prog_rcf%tracer(jc,jk,jb,iqv)/qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
      ENDDO
    ENDDO

    DO jk = 1, kend
      DO jc = i_startidx, i_endidx
        IF (qcana_mode == 2 .AND. pt_prog_rcf%tracer(jc,jk,jb,iqc) > 0._wp) THEN
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + &
            iau_wgt_adv*pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
          pt_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqc) + &
            iau_wgt_adv*pt_diag%rhoc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ELSE 
          IF (qcana_mode >= 1) THEN
            zqin = (pt_diag%rhov_incr(jc,jk,jb)+pt_diag%rhoc_incr(jc,jk,jb))/pt_prog%rho(jc,jk,jb)
          ELSE
            zqin = pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
          ENDIF
          ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
          IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
            pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
        ENDIF

        IF (qiana_mode > 0) THEN
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqi) + &
            iau_wgt_adv*pt_diag%rhoi_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDIF
      ENDDO
    ENDDO


  END SUBROUTINE iau_update_tracer



  !
  ! Add slow-physics tendencies to tracer fields
  !
  ! Add slow-physics tendencies to tracer fields. Currently, 
  ! convection is the only slow-physics routine which provides tracer 
  ! tendencies.
  ! In addition, this routine
  ! - makes sure that tendencies from advection and/or convection 
  !   do not result in negative mass fractions. If negative values in qx  
  !   occur, these are clipped. The moisture which is spuriously created by this 
  !   clipping is substracted from qv.
  ! - Diagnoses amount of convective rain and snow (rain_con, snow_con), 
  !   as well as the total convective precipitation (prec_con).
  ! - applies large-scale-forcing tendencies, if ICON is run in single-column-mode.
  ! 
  SUBROUTINE tracer_add_phytend( p_rho_now, prm_nwp_tend, pdtime, prm_diag, &
    &                            pt_prog_rcf, jg, jb, i_startidx, i_endidx, kend)

    REAL(wp)             &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS         &
#endif
                        ,INTENT(IN)   :: p_rho_now(:,:)  !< total air density
    TYPE(t_nwp_phy_tend),INTENT(IN)   :: prm_nwp_tend    !< atm tend vars
    REAL(wp),            INTENT(IN)   :: pdtime          !< time step
    TYPE(t_nwp_phy_diag),INTENT(INOUT):: prm_diag        !< the physics variables
    TYPE(t_nh_prog),     INTENT(INOUT):: pt_prog_rcf     !< the tracer field at
                                                         !< reduced calling frequency
    INTEGER             ,INTENT(IN)   :: jg              !< domain ID
    INTEGER,             INTENT(IN)   :: jb              !< block index
    INTEGER,             INTENT(IN)   :: i_startidx, i_endidx
    INTEGER             ,INTENT(IN)   :: kend            !< vertical end index                             

    ! Local variables
    INTEGER  :: jt          ! tracer loop index
    INTEGER  :: idx         ! tracer position in container
    INTEGER  :: pos_qv      ! position of qv in local array zrhox
    INTEGER  :: jk,jc
    INTEGER  :: iq_start
    REAL(wp) :: zrhox(nproma,kend,5)
    REAL(wp) :: zrhox_clip(nproma,kend)
    !
    INTEGER, POINTER              :: ptr_conv_list(:)
    INTEGER, DIMENSION(3), TARGET :: conv_list_small
    INTEGER, DIMENSION(5), TARGET :: conv_list_large


    ! get list of water tracers which are affected by convection
    IF (atm_phy_nwp_config(jg)%ldetrain_conv_prec) THEN
      conv_list_large = (/iqv,iqc,iqi,iqr,iqs/)
      ptr_conv_list =>conv_list_large
    ELSE
      conv_list_small = (/iqv,iqc,iqi/)
      ptr_conv_list =>conv_list_small
    ENDIF

    zrhox_clip(:,:) = 0._wp

    ! add tendency due to convection
    DO jt=1,SIZE(ptr_conv_list)
      idx = ptr_conv_list(jt)
      DO jk = kstart_moist(jg), kend
        DO jc = i_startidx, i_endidx
          zrhox(jc,jk,jt) = p_rho_now(jc,jk)*pt_prog_rcf%tracer(jc,jk,jb,idx)  &
            &             + pdtime*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,idx)

          ! keep mass that is created due to artificial clipping
          zrhox_clip(jc,jk) = zrhox_clip(jc,jk) + MIN(0._wp,zrhox(jc,jk,jt))

          ! clip
          zrhox(jc,jk,jt) = MAX(0._wp, zrhox(jc,jk,jt))
        ENDDO
      ENDDO
      !
      ! Re-diagnose tracer mass fraction from partial mass
      IF (idx == iqv) THEN
        pos_qv = jt   ! store local qv-position for later use
        CYCLE         ! special treatment see below
      ENDIF
      !
      DO jk = kstart_moist(jg), kend
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,idx) = zrhox(jc,jk,jt)/p_rho_now(jc,jk)
        ENDDO
      ENDDO
    ENDDO ! jt
    !
    ! Special treatment for qv. 
    ! Rediagnose tracer mass fraction and substract mass created by artificial clipping.
    DO jk = kstart_moist(jg), kend
      DO jc = i_startidx, i_endidx
        pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp, &
          &                                    (zrhox(jc,jk,pos_qv) + zrhox_clip(jc,jk)) &
          &                                    /p_rho_now(jc,jk)                         &
          &                                   )
      ENDDO
    ENDDO



    IF(lart .AND. art_config(jg)%lart_conv) THEN
      ! add convective tendency and fix to positive values
      DO jt=1,art_config(jg)%nconv_tracer  ! ASH
        DO jk = 1, kend
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%conv_tracer(jb,jt)%ptr(jc,jk)=MAX(0._wp,pt_prog_rcf%conv_tracer(jb,jt)%ptr(jc,jk) &
               +pdtime*prm_nwp_tend%conv_tracer_tend(jb,jt)%ptr(jc,jk)/p_rho_now(jc,jk))
          ENDDO
        ENDDO
      ENDDO
    ENDIF !lart

    ! additional clipping for qr, qs, ... up to iqm_max
    ! (very small negative values may occur during the transport process (order 10E-15))
    IF (atm_phy_nwp_config(jg)%ldetrain_conv_prec) THEN
      iq_start = iqg  ! qr, qs already clipped above
    ELSE
      iq_start = iqr
    ENDIF 
    DO jt=iq_start, iqm_max  ! qr,qs,etc. 
      DO jk = kstart_moist(jg), kend
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,jt) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
        ENDDO
      ENDDO
    ENDDO
    
    ! clipping for number concentrations
    IF(ANY((/4,5,6/) == atm_phy_nwp_config(jg)%inwp_gscp))THEN
      DO jt=iqni, ininact  ! qni,qnr,qns,qng,qnh,qnc and ninact (but not yet ninpot)
        DO jk = kstart_moist(jg), kend
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) = MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt))
          ENDDO          
        ENDDO
      ENDDO
    END IF



    ! Diagnose convective precipitation amount
    IF (atm_phy_nwp_config(jg)%lcalc_acc_avg) THEN
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx

        prm_diag%rain_con(jc,jb) = prm_diag%rain_con(jc,jb)    &
          &                      + pdtime * prm_diag%rain_con_rate(jc,jb)

        prm_diag%snow_con(jc,jb) = prm_diag%snow_con(jc,jb)    &
          &                      + pdtime * prm_diag%snow_con_rate(jc,jb)

        prm_diag%prec_con(jc,jb) = prm_diag%rain_con(jc,jb) + prm_diag%snow_con(jc,jb)

      ENDDO
    ENDIF


    IF(is_ls_forcing)THEN
      DO jt=1, nqtendphy  ! qv,qc,qi
        DO jk = kstart_moist(jg), kend
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            pt_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, pt_prog_rcf%tracer(jc,jk,jb,jt)    &
              &                       + pdtime*prm_nwp_tend%ddt_tracer_ls(jk,jt))
          ENDDO
        ENDDO
      END DO
    ENDIF  ! is_ls_forcing

  END SUBROUTINE tracer_add_phytend




  SUBROUTINE cal_cape_cin ( i_startidx, i_endidx, kmoist, te, qve, prs, hhl,  &
                            cape_ml, cin_ml )

  !------------------------------------------------------------------------------
  !
  !>
  !! Description:
  !!  Computation of Convective Available Potential Energy CAPE,
  !!  Convective Inhibition CIN based on parcel theory.
  !!  This subroutine is based on COSMO code.
  !!        Helmut Frank
  !! 
  !! Input:  
  !!         - Temperature, specific humidity and pressure of environment
  !!
  !! Output: 
  !!         - cape_ml/cin_ml: CAPE/CIN based on a parcel with thermodynamical 
  !!                           properties of the lowest mean layer in the PBL (50hPa)
  !!      
  !! Motivation: 
  !!  Current parameter CAPE_CON is calculated in LM in the framework of the 
  !!  convective parametrisation scheme. Therefore this parameter is only available
  !!  at those gridpoints, where the scheme is called, but not continuously on the 
  !!  whole domain. This subroutine, on the other hand, provides continuous fields. 
  !!
  !! Method:
  !!  A dry/moist parcel ascent is performed following classic parcel theory.
  !!  Moist adiabatic ascent is calculated iteratively with an appropriate scheme.
  !!  Based on the temperature and moisture of the ascending parcel, CAPE and CIN
  !!  are computed, closely following the recommendations of Doswell and Rasmussen 
  !!  (1994), including a virtual temperature correction and searching for the 
  !!  most unstable parcel in the lower troposphere. Additionally, a mixed layer 
  !!  CAPE as well as the traditional Showalter Index and the surface lifted 
  !!  index are computed as further variables. 
  !!
  !!  References used during development: 
  !!  - C. A. Doswell and Rasmussen, E. N.: The Effect of Neglecting the 
  !!    Virtual Temperature Correction on CAPE Calculations. 
  !!    Weather and Forecasting, 9, 625-629.
  !!
  !!  - K. A. Emanuel (1994): Atmospheric Convection. Oxford University Press.
  !!
  !!  - H. Huntrieser et al. (1997): Comparison of Traditional and Newly Developed 
  !!    Thunderstorm Indices for Switzerland. Weather and Forecasting, 12, 
  !!    108-125.
  !!
  !!  - D. Bolton (1980): The Computation of Equivalent Potential Temperature. 
  !!    Monthly Weather Review, 108, 1046-1053
  !!
  !!  - Davies, J.M.,2002: On low-level thermodynamic parameters
  !!    associated with tornadic and nontornadic supercells.
  !!    Preprints, 21st Conf. On Severe Local Storms, San Antonio, Amer. Meteor. Soc.
  !!    http://members.cox.net/jondavies1/LLthermo.PDF
  !!
  !! @par Revision History
  !! Inherited from COSMO by Helmut Frank, DWD (2015-05-13)
  !! 
  !!

! Input data
!----------- 
  INTEGER, INTENT (IN) ::  &
    i_startidx, i_endidx,  &  ! start and end indices of loops in hoirozntal patch
    kmoist                    ! start index for moist processes

  REAL    (wp),    INTENT (IN) ::  &
    te  (:,:),   & ! environment temperature
    qve (:,:),   & ! environment specific humidity
    prs (:,:),   & ! full level pressure
    hhl (:,:)      ! height of half levels

! Output data
!------------ 
  REAL (wp), INTENT (OUT) :: &
    cape_ml  (:),   & ! mixed layer CAPE_ML
    cin_ml   (:)      ! mixed layer CIN_ML

! Local scalars and automatic arrays
!-----------------------------------
  INTEGER :: nlev

  REAL    (wp) ::           &      
    qvp_start(SIZE(qve,1)), & ! parcel initial specific humidity in mixed layer
    tp_start (SIZE(te,1))     ! parcel initial pot. temperature in mixed layer

  REAL (wp), PARAMETER :: p0 = 1.e5_wp   ! reference pressure for calculation of
  REAL (wp), PARAMETER :: missing_value  = -999.9_wp   ! Missing value for CIN (if no LFC/CAPE was found),

! Depth of mixed surface layer: 50hPa following Huntrieser, 1997.
! Other frequently used value is 100hPa.
  REAL (wp), PARAMETER :: ml_depth = 5000._wp

  INTEGER              :: &     
  i, k,                   & ! Indices of input/output fields
  k_ml(SIZE(te,1)),       & ! Index for calculation of mixed layer averages 
                            ! (potential temperature, moisture)
  kstart(SIZE(te,1)),     & ! Model level approx. corresponding to mixed layer mean pressure
  lcllev(SIZE(te,1)),     & ! Indices for Lifting Condensation Level LCL,
  lfclev(SIZE(te,1))        ! Level of Free Convection LFC
  
! The following parameters are help values for the iterative calculation 
! of the parcel temperature during the moist adiabatic ascent
  REAL    (wp)             :: esat,tguess1,tguess2,thetae1,thetae2
! REAL    (wp)             :: rp, r1,r2
  REAL    (wp)             :: q1, q2
  REAL    (wp), PARAMETER  :: eps=0.03

! this parameter helps to find the LFC above a capping inversion in cases, 
! where a LFC already was found in an unstable layer in the convective 
! boundary layer below. 
  REAL (wp), PARAMETER :: cc_comp    = 2.0_wp                      
      
  INTEGER ::    icount              ! counter for the iterative process
      
  REAL (wp) ::             &
    cin_help(SIZE(te,1)),  & ! help variable, the CIN above the LFC
    buo     (SIZE(te,1)),  & ! parcel buoyancy at level k
    tp      (SIZE(te,1)),  & ! temperature profile of ascending parcel
    qvp     (SIZE(te,1)),  & ! specific moisture profile of ascending parcel
    thp     (SIZE(te,1)),  & ! 1st guess theta_e of parcel for iterative 
    tvp,                   & ! virtual temperature of parcel at level k
    tve,                   & ! virtual temperature of environment at level k
    buo_belo,              & ! parcel buoyancy of level k+1 below
    esatp,                 & ! saturation vapour pressure at level k
    qvsp                     ! saturation specific humidity at level k
                             ! calculation of moist adiabatic ascent

  INTEGER :: lfcfound(SIZE(te,1))   ! flag indicating if a LFC has already been found
                                    ! below, in cases where several EL and LFC's occur
  
!------------------------------------------------------------------------------
! 
! A well mixed near surface layer is assumed (its depth is specified with 
! parameter ml_depth) Potential temperature and specific humidity are constant
! in this layer, they are calculated as arithmetical means of the corresponding
! variables of the environment (model) profile. The parcel starts from a level 
! approximately in the middle of this well mixed layer, with the average spec. 
! humidity and potential temperature as start values. 
!
!------------------------------------------------------------------------------
        
    nlev = SIZE( te,2)
    k_ml  (:)  = nlev  ! index used to step through the well mixed layer
    kstart(:)  = nlev  ! index of model level corresponding to average 
                       ! mixed layer pressure
    qvp_start(:) = 0.0_wp ! specific humidities in well mixed layer
    tp_start (:) = 0.0_wp ! potential temperatures in well mixed layer
            
    ! now calculate the mixed layer average potential temperature and 
    ! specific humidity
    DO k = nlev, kmoist, -1
      DO i = i_startidx, i_endidx

        IF ( prs(i,k) > (prs(i,nlev) - ml_depth)) THEN
          qvp_start(i) = qvp_start(i) + qve(i,k)
          tp_start (i) = tp_start (i) + te (i,k)*(p0/prs(i,k))**rd_o_cpd
             
          ! Find the level, where pressure approximately corresponds to the 
          ! average pressure of the well mixed layer. Simply assume a threshold
          ! of ml_depth/2 as average pressure in the layer, if this threshold 
          ! is surpassed the level with approximate mean pressure is found
          IF (prs(i,k) > prs(i,nlev) - ml_depth*0.5_wp) THEN
            kstart(i) = k
          ENDIF

          k_ml(i) = k - 1
        ELSE
          EXIT
        ENDIF

      ENDDO     
    ENDDO
        
    ! Calculate the start values for the parcel ascent, 
    DO i = i_startidx, i_endidx
      qvp_start(i) =  qvp_start(i) / (nlev-k_ml(i))
      tp_start (i) =  tp_start (i) / (nlev-k_ml(i))
    ENDDO     
  
  !------------------------------------------------------------------------------
  !
  ! Description:
  !   A single parcel ascent is performed, based on the given start 
  !   values kstart (level), tp_start (initial parcel temperature) and
  !   qvp_start (initial parcel specific humidity). 
  !
  !------------------------------------------------------------------------------
  
  ! Initialization
  
  cape_ml(:)  = 0.0_wp
  cin_ml(:)   = 0.0_wp
  
  lcllev  (:) = 0
  lfclev  (:) = 0
  lfcfound(:) = 0
  cin_help(:) = 0.0_wp
  tp (:)      = 0.0_wp
  qvp(:)      = 0.0_wp               
  buo(:)      = 0.0_wp
  
  ! Loop over all model levels above kstart
  kloop: DO k = nlev, kmoist, -1

    DO i = i_startidx, i_endidx
      IF ( k > kstart(i) ) CYCLE
         
      ! Dry ascent if below cloud base, assume first level is not saturated 
      ! (first approximation)
      IF (k > lcllev(i)) THEN
        tp (i)   = tp_start(i)*( prs(i,k)/p0)**rd_o_cpd   ! dry adiabatic process
        qvp(i)   = qvp_start(i)                           ! spec humidity conserved
            
        ! Calculate parcel saturation vapour pressure and saturation 
        ! specific humidity
        esatp = sat_pres_water( tp(i))
        qvsp  = fqvs( esatp, prs(i,k), qvp(i))
            
        ! Check whether parcel is saturated or not and 
        ! no LCL was already found below
        IF ( (qvp(i) >= qvsp) .AND. (lcllev(i) == 0) ) THEN  
          lcllev(i) = k                                    ! LCL is reached

          ! Moist ascent above LCL, first calculate an approximate thetae to hold 
          ! constant during the remaining ascent
!         rp      = qvp(i)/( 1._wp - qvp(i) )
!         thp(i)  = fthetae( tp(i),prs(i,k),rp )
          thp(i)  = fthetae( tp(i),prs(i,k), qvp(i) )

        ENDIF
      ENDIF
         
      ! Moist adiabatic process: the parcel temperature during this part of 
      ! the ascent is calculated iteratively using the iterative newton
      ! scheme, assuming the equivalent potential temperature of the parcel 
      ! at the LCL (thp) is held constant. The scheme converges usually within
      ! few (less than 10) iterations, its accuracy can be tuned with the 
      ! parameter "eps", a value of 0.03 is tested and recommended. 
            
      IF ( k <= lcllev(i) ) THEN                                
        ! The scheme uses a first guess temperature, which is the parcel 
        ! temperature at the level below. If it happens that the initial 
        ! parcel is already saturated, the environmental temperature 
        ! is taken as first guess instead
        IF (  k == kstart(i) ) THEN
          tguess1 = te(i,kstart(i))            
        ELSE
          tguess1 = tp(i)
        END IF
        icount = 0       ! iterations counter

        ! Calculate iteratively parcel temperature from 
        ! thp, prs and 1st guess tguess1
        DO
          esat     = sat_pres_water( tguess1)
!         r1       = rdv*esat/(prs(i,k)-esat)
!         thetae1  = fthetae( tguess1,prs(i,k),r1)
          q1       = fqvs( esat, prs(i,k), qvp(i) )
          thetae1  = fthetae( tguess1,prs(i,k),q1)

          tguess2  = tguess1 - 1.0_wp
          esat     = sat_pres_water( tguess2)
!         r2       = rdv*esat/(prs(i,k)-esat)
!         thetae2  = fthetae( tguess2,prs(i,k),r2)
          q2       = fqvs( esat, prs(i,k), qvp(i) )
          thetae2  = fthetae( tguess2,prs(i,k),q2)

          tguess1  = tguess1+(thetae1-thp(i))/(thetae2-thetae1)
          icount   = icount    + 1   

          IF ( ABS( thetae1-thp(i)) < eps .OR. icount > 20 ) THEN
            tp(i) = tguess1
            EXIT
          END IF
        END DO

        ! update specific humidity of the saturated parcel for new temperature
        esatp  = sat_pres_water( tp(i))
        qvp(i) = fqvs( esatp,prs(i,k),qvp(i))
      END IF
         
      ! Calculate virtual temperatures of parcel and environment
      tvp    = tp(i  ) * (1.0_wp + vtmpc1*qvp(i  )/(1.0_wp - qvp(i  )) )  
      tve    = te(i,k) * (1.0_wp + vtmpc1*qve(i,k)/(1.0_wp - qve(i,k)) ) 
         
      ! Calculate the buoyancy of the parcel at current level k, 
      ! save buoyancy from level k+1 below (buo_belo) to check if LFC or EL have been passed
      buo_belo = buo(i)
      buo(i)   = tvp - tve

      ! Check for level of free convection (LFC) and set flag accordingly. 
      ! Basic LFC condition is that parcel buoyancy changes from negative to 
      ! positive (comparison of buo with buo_belo). Tests showed that very 
      ! often the LFC is already found within the boundary layer below even if 
      ! significant capping inversions are present above (and since CIN is only
      ! defined below the LFC no CIN was accumulated in these cases.)
      ! To handle these situations in a meteorologically meaningful way an 
      ! additional flag "lfcfound" was introduced which is initially zero but 
      ! set to 1 if a second LFC was found, under the condition that the CIN 
      ! within the capping inversion is greater than the CAPE in the convective
      ! boundary layer below times the factor cc_comp (cc_comp = 1 - 2 
      ! recommended.)
      ! Help variable CIN_HELP saves all contributions to the total cin above 
      ! the LFC and has to be subtracted at the end from the final CIN in order
      ! to get the CIN only below the LFC (this is necessary since we do not 
      ! know yet where exactly we will find an LFC when performing the ascent
      ! from bottom to top in a stepwise manner.)

      ! Find the first LFC
      IF ( (buo(i) > 0.0_wp) .AND. (buo_belo <= 0.0_wp)            &
                             .AND. ( lfcfound(i)==0) ) THEN
            
        ! Check whether it is an LFC at one of the lowest model levels 
        ! (indicated by CAPE=0)
        IF ( (cape_ml(i) > 0.0_wp) .AND. ( lfcfound(i) == 0 ) ) THEN
          ! Check if there is a major capping inversion below, defined as 
          ! having CIN with an absolute value larger than the CAPE accumulated
          ! below times some arbitrary factor cc_comp - if this is the case the
          ! LFC index "lfclev" is updated to the current level k and 
          ! "lfcfound"-flag is now set to 1 assuming that we have found the 
          ! level of free convection finally. 
          IF ( cc_comp * ABS(cin_help(i)) > cape_ml(i) ) THEN
            lfclev(i)   = k
            cape_ml (i) = 0.0_wp
            cin_help(i) = 0.0_wp
            lfcfound(i) = 1
          ENDIF
        ELSE
          ! the LFC found is near the surface, set the LFC index to the current
          ! level k (lfclev) but do not set the flag "lfcfound" to zero to 
          ! indicate that a further LFC may be present above the boundary layer
          ! and an eventual capping inversion. Reset the CIN_HELP to zero to 
          ! store the contribution of CIN above this LFC.
          lfclev(i)   = k
          cin_help(i) = 0.0_wp
        ENDIF
      ENDIF
         
      ! Accumulation of CAPE and CIN according to definition given in Doswell 
      ! and Rasmussen (1994), 
      IF ( (buo(i) >= 0.0_wp) .AND. (k <= lfclev(i)) ) THEN   
        cape_ml(i)  = cape_ml(i)  + (buo(i)/tve)*grav*(hhl(i,k) - hhl(i,k+1))
      ELSEIF ( (buo(i) < 0.0) .AND. (k < kstart(i)) ) THEN  
        cin_ml(i)   = cin_ml(i)   + (buo(i)/tve)*grav*(hhl(i,k) - hhl(i,k+1))
        cin_help(i) = cin_help(i) + (buo(i)/tve)*grav*(hhl(i,k) - hhl(i,k+1))
      ENDIF

    ENDDO ! i = i_startidx, i_endidx
  ENDDO  kloop       ! End k-loop over levels
      
    ! Subtract the CIN above the LFC from the total accumulated CIN to 
    ! get only contriubtions from below the LFC as the definition demands.
  DO i = i_startidx, i_endidx

    ! make CIN positive
    cin_ml(i) = ABS (cin_ml(i) - cin_help(i))

    ! set the CIN to missing value if no LFC was found or no CAPE exists
    IF ( (lfclev(i) == 0) .OR. (cape_ml(i) == 0.0_wp)  ) cin_ml(i) = missing_value 
  ENDDO

CONTAINS

! Specific humidity at saturation as function of water vapor pressure zex,
! air pressure zpx, and specific humidity zqx.
  ELEMENTAL FUNCTION fqvs( zex, zpx, zqx)

    REAL(wp), INTENT(IN) :: zex   ! vapor pressure        [Pa]
    REAL(wp), INTENT(IN) :: zpx   ! atmospheric pressure  [Pa]
    REAL(wp), INTENT(IN) :: zqx   ! specific humidity     [kg/kg]
    REAL(wp)             :: fqvs  ! Equivalent potential temperature

    fqvs = zex/zpx *( rdv + o_m_rdv*zqx )        

  END FUNCTION fqvs

! Equivalent potential temperature to hold constant during ascent
! ELEMENTAL FUNCTION fthetae( ztx,zpx,zrx)
  ELEMENTAL FUNCTION fthetae( ztx,zpx,zqx)

    REAL(wp), INTENT(IN) :: ztx     ! air temperature       [K]
    REAL(wp), INTENT(IN) :: zpx     ! atmospheric pressure  [Pa]
!   REAL(wp), INTENT(IN) :: zrx     ! mixing ratio          [kg/kg]
    REAL(wp), INTENT(IN) :: zqx     ! specific humidity     [kg/kg]
    REAL(wp)             :: fthetae ! Equivalent potential temperature [K]

!   fthetae = (p0/zpx)**rd_o_cpd *ztx*exp( alvdcp*zrx/ztx)  
    fthetae = (p0/zpx)**rd_o_cpd *ztx*exp( alvdcp*zqx/(ztx*(1._wp-zqx)) )

  END FUNCTION fthetae

  END SUBROUTINE cal_cape_cin


END MODULE mo_util_phys
