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
    &                                 vtmpc1, t3
  USE mo_exception,             ONLY: finish, message
  USE mo_satad,                 ONLY: sat_pres_water, sat_pres_ice
  USE mo_fortran_tools,         ONLY: assign_if_present
  USE mo_impl_constants,        ONLY: min_rlcell_int, min_rledge_int
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,         ONLY: prm_diag
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqni, ininact, &
       &                              iqm_max, nqtendphy, lart, &
       &                              iqh, iqnc, iqnr, iqns, iqng, iqnh, msg_level
  USE mo_nh_diagnose_pres_temp, ONLY: diag_pres, diag_temp
  USE mo_ls_forcing_nml,        ONLY: is_ls_forcing
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_nwp_tuning_config,     ONLY: tune_gust_factor
  USE mo_advection_config,      ONLY: advection_config
  USE mo_art_config,            ONLY: art_config
  USE mo_initicon_config,       ONLY: iau_wgt_adv, qcana_mode, qiana_mode, qrsgana_mode, qnxana_2mom_mode
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_satad,                 ONLY: qsat_rho
  USE mo_upatmo_config,         ONLY: upatmo_config


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: nwp_dyn_gust
  PUBLIC :: nwp_con_gust
  PUBLIC :: virtual_temp
  PUBLIC :: vap_pres
  PUBLIC :: swdir_s
  PUBLIC :: rel_hum
  PUBLIC :: compute_field_rel_hum_wmo
  PUBLIC :: compute_field_rel_hum_ifs
  PUBLIC :: iau_update_tracer
  PUBLIC :: tracer_add_phytend

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_phys'

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
    !$acc routine seq

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
!$ACC ROUTINE SEQ

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
!$ACC ROUTINE SEQ
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
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

!!$ UB: do we need p_prog_rcf%tracer_ptr(iqv) instead of p_prog%tracer_ptr(iqv)?
          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
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
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif

!!$ UB: do we need p_prog_rcf%tracer_ptr(iqv) instead of p_prog%tracer_ptr(iqv)?
          ! get values for temperature, etc.:
          temp = p_diag%temp(jc,jk,jb)
          qv   = p_prog%tracer_ptr(iqv)%p_3d(jc,jk,jb)
          p_ex = p_prog%exner(jc,jk,jb)
          !-- compute relative humidity as r = e/e_s:
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
    CALL diag_pres (pt_prog, pt_diag, p_metrics, jb, i_startidx, i_endidx, 1, kend, &
      &             opt_lconstgrav=upatmo_config(jg)%nwp_phy%l_constgrav)

    ! Compute relative humidity w.r.t. water
    DO jk = 1, kend
      DO jc = i_startidx, i_endidx
        zrhw(jc,jk) = pt_prog_rcf%tracer(jc,jk,jb,iqv)/qsat_rho(pt_diag%temp(jc,jk,jb),pt_prog%rho(jc,jk,jb))
      ENDDO
    ENDDO

    ! GZ: This loop needs to be split for correct vectorization because rhoc_incr is allocated for qcana_mode >= 1 only;
    !     otherwise, the NEC runs into a segfault. Likewise, the remaining case selections need to be done outside the
    !     vectorized loops in order to avoid invalid memory accesses.
    DO jk = 1, kend
      IF (qcana_mode >= 1) THEN
        DO jc = i_startidx, i_endidx
          IF (qcana_mode == 2 .AND. pt_prog_rcf%tracer(jc,jk,jb,iqc) > 0._wp) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = pt_prog_rcf%tracer(jc,jk,jb,iqv) + &
              iau_wgt_adv*pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
            pt_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqc) + &
              iau_wgt_adv*pt_diag%rhoc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          ELSE 
            zqin = (pt_diag%rhov_incr(jc,jk,jb)+pt_diag%rhoc_incr(jc,jk,jb))/pt_prog%rho(jc,jk,jb)
            ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
            IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
              pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
            pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
          ENDIF
        ENDDO
      ELSE
        DO jc = i_startidx, i_endidx
          zqin = pt_diag%rhov_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb)
          ! DA increments of humidity are limited to positive values if p > 150 hPa and RH < 2% or QV < 5.e-7
          IF (pt_diag%pres(jc,jk,jb) > 15000._wp .AND. zrhw(jc,jk) < 0.02_wp .OR. &
            pt_prog_rcf%tracer(jc,jk,jb,iqv) < 5.e-7_wp) zqin = MAX(0._wp, zqin)
          pt_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqv) + iau_wgt_adv*zqin)
        ENDDO
      ENDIF

      IF (qiana_mode > 0) THEN
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqi) + &
            iau_wgt_adv*pt_diag%rhoi_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (qrsgana_mode > 0) THEN
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqr) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqr) + &
            iau_wgt_adv * pt_diag%rhor_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          pt_prog_rcf%tracer(jc,jk,jb,iqs) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqs) + &
            iau_wgt_adv * pt_diag%rhos_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (qrsgana_mode > 0 .AND. iqg <= iqm_max) THEN
        DO jc = i_startidx, i_endidx
          pt_prog_rcf%tracer(jc,jk,jb,iqg) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqg) + &
            iau_wgt_adv * pt_diag%rhog_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
        ENDDO
      ENDIF

      IF (atm_phy_nwp_config(jg)%l2moment) THEN
        DO jc = i_startidx, i_endidx
          IF (qcana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqnc) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnc) + &
                 iau_wgt_adv * pt_diag%rhonc_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
          IF (qiana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqni) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqni) + &
                 iau_wgt_adv * pt_diag%rhoni_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
          IF (qrsgana_mode > 0) THEN
            pt_prog_rcf%tracer(jc,jk,jb,iqh) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqh) + &
                 iau_wgt_adv * pt_diag%rhoh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqnr) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnr) + &
                 iau_wgt_adv * pt_diag%rhonr_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqns) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqns) + &
                 iau_wgt_adv * pt_diag%rhons_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqng) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqng) + &
                 iau_wgt_adv * pt_diag%rhong_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
            pt_prog_rcf%tracer(jc,jk,jb,iqnh) = MAX(0._wp,pt_prog_rcf%tracer(jc,jk,jb,iqnh) + &
                 iau_wgt_adv * pt_diag%rhonh_incr(jc,jk,jb)/pt_prog%rho(jc,jk,jb))
          END IF
        ENDDO
      ENDIF

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
    IF(atm_phy_nwp_config(jg)%l2moment)THEN
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


END MODULE mo_util_phys
