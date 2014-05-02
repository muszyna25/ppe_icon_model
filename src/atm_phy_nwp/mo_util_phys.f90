!>
!! Implementation of physics utility routines.
!!
!! @par Revision History
!!  Initial revision  :  F. Prill, DWD (2012-07-03)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
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
    &                                 grav
  USE mo_satad,                 ONLY: sat_pres_water, sat_pres_ice
  USE mo_fortran_tools,         ONLY: assign_if_present
  USE mo_impl_constants,        ONLY: min_rlcell_int
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag
  USE mo_run_config,            ONLY: iqv
  USE mo_loopindices,           ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = &
    &  '$Id$'


  PUBLIC :: nwp_dyn_gust
  PUBLIC :: nwp_con_gust
  PUBLIC :: virtual_temp
  PUBLIC :: vap_pres
  PUBLIC :: rel_hum
  PUBLIC :: compute_field_rel_hum_wmo
  PUBLIC :: compute_field_rel_hum_ifs
  PUBLIC :: compute_field_omega


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
  ELEMENTAL FUNCTION nwp_dyn_gust( u_10m, v_10m, tcm, u1, v1) RESULT( vgust_dyn)

    REAL(wp), INTENT(IN) :: u_10m, &    ! zonal wind component at 10 m above ground [m/s]
      &                     v_10m, &    ! meridional wind component at 10 m above ground [m/s]
      &                     tcm  , &    ! transfer coefficient for momentum at surface
      &                     u1   , &    ! zonal wind at lowest model layer above ground [m/s]
      &                     v1          ! meridional wind at lowest model layer above ground [m/s]
    REAL(wp) :: vgust_dyn               ! dynamic gust at 10 m above ground [m/s]

    REAL(wp) :: ff10m, ustar
    REAL(wp), PARAMETER :: gust_factor = 3.0_wp * 2.4_wp

    ff10m = SQRT( u_10m**2 + v_10m**2)
    ustar = SQRT( MAX( tcm, 5.e-4_wp) * ( u1**2 + v1**2) )
    vgust_dyn = ff10m + gust_factor*ustar

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
  SUBROUTINE virtual_temp(p_patch, temp, qv, qc, qi, qr, qs, temp_v)


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN)           :: temp(:,:,:) ! temperature (K)
    REAL(wp), INTENT(IN)           :: qv  (:,:,:) ! specific humidity
    REAL(wp), INTENT(IN), OPTIONAL :: qc  (:,:,:) ! specific cloud water
    REAL(wp), INTENT(IN), OPTIONAL :: qi  (:,:,:) ! specific cloud ice
    REAL(wp), INTENT(IN), OPTIONAL :: qr  (:,:,:) ! specific rain water
    REAL(wp), INTENT(IN), OPTIONAL :: qs  (:,:,:) ! specific snow

    REAL(wp), INTENT(OUT) :: temp_v(:,:,:) ! virtual temperature (K)

    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev
    LOGICAL :: l_cloud_precip

    nlev = SIZE(temp,2) ! in order to be usable for input and output data

    IF (PRESENT(qc) .AND. PRESENT(qi) .AND. PRESENT(qr) .AND. PRESENT(qs)) THEN
      l_cloud_precip = .TRUE.
    ELSE
      l_cloud_precip = .FALSE.
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      IF (l_cloud_precip) THEN
        DO jk = 1, nlev
          DO jc = 1, nlen
            temp_v(jc,jk,jb) = temp(jc,jk,jb) * (1._wp + vtmpc1*qv(jc,jk,jb) -      &
                              (qc(jc,jk,jb)+qi(jc,jk,jb)+qr(jc,jk,jb)+qs(jc,jk,jb)) )
          ENDDO
        ENDDO
      ELSE
        DO jk = 1, nlev
          DO jc = 1, nlen
            temp_v(jc,jk,jb) = temp(jc,jk,jb) * (1._wp + vtmpc1*qv(jc,jk,jb))
          ENDDO
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE virtual_temp



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


END MODULE mo_util_phys
