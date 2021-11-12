!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"
!>
!! Routines for optional diagnostic output variables in NWP
!! (formerly located in mo_util_phys)
!!
!! @par Revision History
!!  Initial revision  :  G. Zaengl, DWD (2020-02-17)
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

MODULE mo_opt_nwp_diagnostics

  USE mo_kind,                  ONLY: vp, wp
  USE mo_parallel_config,       ONLY: nproma
  USE mo_math_constants,        ONLY: pi
  USE mo_physical_constants,    ONLY: o_m_rdv        , & !! 1 - r_d/r_v &
    &                                 rdv,             & !! r_d / r_v
    &                                 vtmpc1,          &
    &                                 grav,            &
    &                                 tmelt, earth_radius, &
    &                                 alvdcp, rd_o_cpd, &
    &                                 rhoh2o, rhoice, K_w_0, K_i_0
  USE gscp_data,                ONLY: cloud_num, isnow_n0temp, zami, mu_rain, zams_ci, zams_gr, zbms, &
    &                                 znimax_Thom, zthn, mma, mmb, zcnue
  USE mo_2mom_mcrph_main,       ONLY: init_2mom_scheme,      &
    &                                 rain_coeffs  ! contains the parameters for the mue-Dm-relation
  USE mo_2mom_mcrph_types,      ONLY: particle, particle_frozen
  USE mo_2mom_mcrph_util,       ONLY: gfct
  USE mo_2mom_mcrph_processes,  ONLY: moment_gamma, rain_mue_dm_relation
  USE mo_exception,             ONLY: finish, message
  USE mo_satad,                 ONLY: sat_pres_water
  USE mo_fortran_tools,         ONLY: assign_if_present
  USE mo_impl_constants,        ONLY: min_rlcell_int, min_rledge_int, &
    &                                 min_rlcell, grf_bdywidth_c
  USE mo_impl_constants_grf,    ONLY: grf_bdyintp_start_c,  &
    &                                 grf_ovlparea_start_c, grf_fbk_start_c
  USE mo_model_domain,          ONLY: t_patch
  USE mo_nonhydro_types,        ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,         ONLY: t_nwp_phy_diag
  USE mo_run_config,            ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqni, &
       &                              iqh, iqnc, iqnr, iqns, iqng, iqnh, iqgl, iqhl, msg_level
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config
  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_io_config,             ONLY: echotop_meta
  USE mo_lnd_nwp_config,        ONLY: nlev_soil, dzsoil
  USE mo_nwp_lnd_types,         ONLY: t_lnd_diag
  USE mo_ext_data_types,        ONLY: t_external_data
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_nwp_sfc_interp,        ONLY: wsoil2smi
  USE mo_icon_interpolation_scalar,                     &
    &                           ONLY: edges2cells_scalar, cells2edges_scalar, &
    &                                 cells2verts_scalar, verts2edges_scalar, &
    &                                 cells2edges_scalar
  USE mo_math_gradients,        ONLY: grad_fd_norm, grad_fd_tang
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_edge
  USE mo_sync,                  ONLY: sync_patch_array, SYNC_C
  USE mo_grf_intp_data_strc,    ONLY: p_grf_state_local_parent
  USE mo_communication,         ONLY: exchange_data
  USE mo_grid_config,           ONLY: l_limited_area
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, get_my_mpi_work_id
#ifdef _OPENACC
  USE mo_mpi,                   ONLY: i_am_accel_node
#endif

  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calsnowlmt
  PUBLIC :: compute_field_omega
  PUBLIC :: compute_field_pv
  PUBLIC :: compute_field_sdi
  PUBLIC :: compute_field_lpi
  PUBLIC :: maximize_field_lpi
  PUBLIC :: compute_field_ceiling
  PUBLIC :: compute_field_hbas_sc
  PUBLIC :: compute_field_htop_sc
  PUBLIC :: compute_field_twater
  PUBLIC :: compute_field_q_sedim
  PUBLIC :: compute_field_dursun
  PUBLIC :: compute_field_tcond_max
  PUBLIC :: compute_field_uh_max
  PUBLIC :: compute_field_vorw_ctmax
  PUBLIC :: compute_field_w_ctmax
  PUBLIC :: compute_field_smi
  PUBLIC :: cal_cape_cin
  PUBLIC :: compute_field_dbz3d_lin
#ifdef HAVE_RADARFWO
  PUBLIC :: compute_field_dbz_1mom
  PUBLIC :: compute_field_dbz_2mom
#endif
  PUBLIC :: compute_field_dbzcmax
  PUBLIC :: compute_field_dbz850
  PUBLIC :: maximize_field_dbzctmax
  PUBLIC :: compute_field_echotop
  PUBLIC :: compute_field_echotopinm

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_opt_nwp_diagnostics'

CONTAINS



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
    REAL(wp):: w_avg               ! vertical velocity averaged to full level
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
      
!$ACC PARALLEL DEFAULT(PRESENT) PRIVATE(w_avg) IF( i_am_accel_node )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = slev, elev
#else
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx
#endif
          ! half level to full level interpolation
          w_avg = 0.5_wp * (p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb))

          out_var(jc,jk,jb) = -p_prog%rho(jc,jk,jb)*grav*w_avg

        ENDDO
      ENDDO
!$ACC END PARALLEL

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
      i_count = ext_data%atm%list_land%ncount(jb)

      DO ic = 1, i_count
        jc = ext_data%atm%list_land%idx(ic,jb)

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

!$ACC DATA CREATE( pv_ef, vt, theta_cf, theta_vf, theta_ef, w_vh, w_eh, ddtw_eh, ddnw_eh, &
!$ACC              ddtth_ef, ddnth_ef, vor_ef ) IF ( i_am_accel_node )
    
!$OMP PARALLEL    
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    !compute theta on cells
    DO jb = i_startblk, i_endblk
    
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      
!$ACC PARALLEL DEFAULT(PRESENT) IF( i_am_accel_node )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          theta_cf(jc,jk,jb) = p_diag%temp(jc,jk,jb) / p_prog%exner(jc,jk,jb)
          
        ENDDO
      ENDDO
!$ACC END PARALLEL
    
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

!$ACC PARALLEL DEFAULT(PRESENT) PRIVATE( ivd1, ivd2, vdfac )
!$ACC LOOP GANG VECTOR COLLAPSE(2)     
      DO jk = slev, elev
        DO je = i_startidx, i_endidx

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
!$ACC END PARALLEL

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

!$ACC END DATA
        
  END SUBROUTINE compute_field_pv


  !>
  !! compute_field_sdi
  !!
  !! Description:
  !!   calculation of the supercell detection indices (SDI1, SDI2)
  !!
  !! Method:
  !!   defined in:
  !!   Wicker L, J. Kain, S. Weiss and D. Bright, A Brief Description of the
  !!          Supercell Detection Index, (available from
  !!   http://www.spc.noaa.gov/exper/Spring_2005/SDI-docs.pdf)
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-05-13) 
  !!
  SUBROUTINE compute_field_sdi( ptr_patch, jg, ptr_patch_local_parent, p_int,    &
                                p_metrics, p_prog, p_diag,                 &
                                sdi_2 )

    IMPLICIT NONE

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch           !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg    ! domain ID of main grid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch_local_parent    !< parent grid for larger exchange halo
    TYPE(t_int_state),  INTENT(IN)    :: p_int   ! for reduced grid

    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag
 
    REAL(wp),           INTENT(OUT)   :: sdi_2(:,:)    !< output variable, dim: (nproma,nblks_c)

    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx

    REAL(wp) :: z_min, z_max

    REAL(wp) :: hsurf( nproma )
    REAL(wp) :: delta_z
    REAL(wp) :: vol( nproma )
    REAL(wp) :: vol_inv
    REAL(wp) :: area_norm

    TYPE(t_patch),  POINTER  :: p_pp

    REAL(wp) :: w_c

    ! vertical averages:
    REAL(wp) :: w_vmean        ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_vmean     ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_w_vmean      ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_zeta_vmean( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_zeta_vmean   ( nproma, ptr_patch%nblks_c)

    ! volume averages:
    REAL(wp) :: w_mean        ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_mean     ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_w_mean      ( nproma, ptr_patch%nblks_c)
    REAL(wp) :: zeta_zeta_mean( nproma, ptr_patch%nblks_c)
    REAL(wp) :: w_zeta_mean   ( nproma, ptr_patch%nblks_c)

    INTEGER, PARAMETER :: idx_w         = 1
    INTEGER, PARAMETER :: idx_zeta      = 2
    INTEGER, PARAMETER :: idx_w_w       = 3
    INTEGER, PARAMETER :: idx_zeta_zeta = 4
    INTEGER, PARAMETER :: idx_w_zeta    = 5

    ! the corresponding fields on the parent grid:
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: p_vmean  ! vertical means on the parent grid
    REAL(wp) :: p_mean( nproma, 5)                      ! means on the parent grid

    REAL(wp) :: w2_mean
    REAL(wp) :: zeta2_mean
    REAL(wp) :: helic_mean
    REAL(wp) :: helic_w_corr

    INTEGER :: jb, jc, jk
    INTEGER :: jb2, jc2
    INTEGER :: i, l, idx
    INTEGER :: ist
    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

    REAL(wp), POINTER :: p_fbkwgt(:,:,:)

    INTEGER :: nblks_c_lp

    REAL(wp) :: EPS = 1.0E-20_wp

    ! definition of the vertical integration limits:
    z_min = 1500.0  ! in m
    z_max = 5500.0  ! in m

    ! --- to prevent errors at the boundaries, set some field(s) to 0:

    i_rlstart = 1
    i_rlend   = grf_bdywidth_c

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        !sdi_1   ( jc, jb ) = 0.0_wp
        sdi_2 ( jc, jb ) = 0.0_wp

        w_vmean        (jc,jb) = 0.0_wp
        zeta_vmean     (jc,jb) = 0.0_wp
        w_w_vmean      (jc,jb) = 0.0_wp
        zeta_zeta_vmean(jc,jb) = 0.0_wp
        w_zeta_vmean   (jc,jb) = 0.0_wp

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ! --- calculate vertical averages ---

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,delta_z,vol,vol_inv,hsurf,w_c) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx

        w_vmean        (jc,jb) = 0.0_wp
        zeta_vmean     (jc,jb) = 0.0_wp
        w_w_vmean      (jc,jb) = 0.0_wp
        zeta_zeta_vmean(jc,jb) = 0.0_wp
        w_zeta_vmean   (jc,jb) = 0.0_wp

        vol  (jc) = 0.0_wp
        hsurf(jc) = p_metrics%z_ifc(jc,ptr_patch%nlev+1,jb)
      END DO

      DO jk = 1, ptr_patch%nlev
        DO jc = i_startidx, i_endidx

          IF ( ( p_metrics%z_ifc(jc,jk+1,jb) >= hsurf(jc) + z_min ) .AND.      &
            &  ( p_metrics%z_ifc(jc,jk  ,jb) <= hsurf(jc) + z_max ) ) THEN
            ! integrate only between z_min and z_max

            delta_z = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)

            w_c = 0.5_wp * ( p_prog%w(jc,jk  ,jb)    &
              &            + p_prog%w(jc,jk+1,jb) )

            w_vmean        (jc,jb) = w_vmean        (jc,jb) + delta_z * w_c
            zeta_vmean     (jc,jb) = zeta_vmean     (jc,jb) + delta_z * p_diag%vor(jc,jk,jb)
            w_w_vmean      (jc,jb) = w_w_vmean      (jc,jb) + delta_z * w_c * w_c
            zeta_zeta_vmean(jc,jb) = zeta_zeta_vmean(jc,jb) + delta_z * p_diag%vor(jc,jk,jb) * p_diag%vor(jc,jk,jb)
            w_zeta_vmean   (jc,jb) = w_zeta_vmean   (jc,jb) + delta_z * w_c * p_diag%vor(jc,jk,jb) 

            vol(jc) = vol(jc) + delta_z
          END IF

        END DO
      END DO

      DO jc = i_startidx, i_endidx
        vol_inv = 1.0_wp / vol(jc)
        w_vmean        (jc,jb) = w_vmean        (jc,jb) * vol_inv
        zeta_vmean     (jc,jb) = zeta_vmean     (jc,jb) * vol_inv
        w_w_vmean      (jc,jb) = w_w_vmean      (jc,jb) * vol_inv
        zeta_zeta_vmean(jc,jb) = zeta_zeta_vmean(jc,jb) * vol_inv
        w_zeta_vmean   (jc,jb) = w_zeta_vmean   (jc,jb) * vol_inv
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! --- average these values to the parent grid cells
    !
    !     motivation: enhance the horiz. averaging area, which is otherwise 
    !     too strongly limited by the limited halo of the domain decomposition.

    p_pp => ptr_patch_local_parent

    p_fbkwgt => p_grf_state_local_parent(jg)%fbk_wgt_bln

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx => p_pp%cells%child_idx
    iblk => p_pp%cells%child_blk

    nblks_c_lp = p_pp%cells%end_block( min_rlcell )

    ALLOCATE( p_vmean ( nproma, 5, nblks_c_lp), STAT=ist )
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_sdi', "allocate failed" )
    END IF

    ! first nullify all lateral grid points

    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_bdyintp_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = grf_fbk_start_c + 1

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,idx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO idx=1, 5
        DO jc = i_startidx, i_endidx
          p_vmean ( jc, idx, jb ) = 0.0_wp
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    
    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_fbk_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = min_rlcell_int

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        p_vmean(jc, idx_w, jb) =                                                 &
            w_vmean        ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + w_vmean        ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + w_vmean        ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + w_vmean        ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_zeta, jb) =                                              &
            zeta_vmean     ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + zeta_vmean     ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + zeta_vmean     ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + zeta_vmean     ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_w_w, jb) =                                               &
            w_w_vmean      ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + w_w_vmean      ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + w_w_vmean      ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + w_w_vmean      ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_zeta_zeta, jb) =                                         &
            zeta_zeta_vmean( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + zeta_zeta_vmean( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + zeta_zeta_vmean( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + zeta_zeta_vmean( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

        p_vmean(jc, idx_w_zeta, jb) =                                            &
            w_zeta_vmean   ( iidx(jc,jb,1), iblk(jc,jb,1) ) * p_fbkwgt(jc,jb,1)  &
          + w_zeta_vmean   ( iidx(jc,jb,2), iblk(jc,jb,2) ) * p_fbkwgt(jc,jb,2)  &
          + w_zeta_vmean   ( iidx(jc,jb,3), iblk(jc,jb,3) ) * p_fbkwgt(jc,jb,3)  &
          + w_zeta_vmean   ( iidx(jc,jb,4), iblk(jc,jb,4) ) * p_fbkwgt(jc,jb,4)

       END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! --- Exchange of these fields on the parent grid

    CALL exchange_data( p_pp%comm_pat_c, recv=p_vmean )

    ! --- Average over the neighbouring parent grid cells

    ! consistency check
    IF ( p_int%cell_environ%max_nmbr_iter /= 1 ) THEN
      CALL finish( modname//':compute_field_sdi', "cell_environ is not built with the right number of iterations" )
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,p_mean,  &
!$OMP            l,jc2,jb2,area_norm,i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,             &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        p_mean( jc, 1:5 ) = 0.0_wp
      END DO

      DO l=1, p_int%cell_environ%max_nmbr_nghbr_cells
        DO jc = i_startidx, i_endidx

          jc2       = p_int%cell_environ%idx      ( jc, jb, l)
          jb2       = p_int%cell_environ%blk      ( jc, jb, l)
          area_norm = p_int%cell_environ%area_norm( jc, jb, l)

          p_mean( jc, idx_w        ) = p_mean( jc, idx_w        ) + area_norm * p_vmean( jc2, idx_w,         jb2)
          p_mean( jc, idx_zeta     ) = p_mean( jc, idx_zeta     ) + area_norm * p_vmean( jc2, idx_zeta,      jb2)
          p_mean( jc, idx_w_w      ) = p_mean( jc, idx_w_w      ) + area_norm * p_vmean( jc2, idx_w_w,       jb2)
          p_mean( jc, idx_zeta_zeta) = p_mean( jc, idx_zeta_zeta) + area_norm * p_vmean( jc2, idx_zeta_zeta, jb2)
          p_mean( jc, idx_w_zeta   ) = p_mean( jc, idx_w_zeta   ) + area_norm * p_vmean( jc2, idx_w_zeta,    jb2)

        END DO
      END DO

      ! write back to the 4 child grid cells:
      DO i = 1, 4
        DO jc = i_startidx, i_endidx
          w_mean        ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_w )
          zeta_mean     ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_zeta )
          w_w_mean      ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_w_w )
          zeta_zeta_mean( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_zeta_zeta )
          w_zeta_mean   ( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_mean( jc, idx_w_zeta )
        END DO
      END DO

    END DO
!$OMP END PARALLEL

    DEALLOCATE( p_vmean, STAT=ist ) 
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_sdi', "deallocate failed" )
    END IF

    ! --- calculate SDI_1, SDI_2 (now again on the child-grid)

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,helic_mean,w2_mean,zeta2_mean,helic_w_corr) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,      &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx

        helic_mean  = w_zeta_mean(jc,jb)    - w_mean(jc,jb) * zeta_mean(jc,jb)
        zeta2_mean  = zeta_zeta_mean(jc,jb) - zeta_mean(jc,jb) * zeta_mean(jc,jb)
        w2_mean     = w_w_mean(jc,jb)       - w_mean(jc,jb) * w_mean(jc,jb)

        IF ( ( w2_mean > EPS ) .AND. ( zeta2_mean > EPS ) ) THEN

          helic_w_corr = helic_mean / SQRT( w2_mean * zeta2_mean )

          !sdi_1(jc,jb) = helic_w_corr * zeta_vmean(jc,jb)

          IF ( w_vmean(jc,jb) > 0.0_wp ) THEN
            sdi_2(jc,jb) = helic_w_corr * ABS( zeta_vmean(jc,jb) )
          ELSE
            sdi_2(jc,jb) = 0.0_wp
          END IF

        ELSE
         !sdi_1(jc,jb) = 0.0_wp
          sdi_2(jc,jb) = 0.0_wp
        END IF

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_sdi


  !>
  !! Calculate the Lightning Potential Index (LPI)
  !!
  !! Literature
  !!       - B. Lynn, Y. Yair, 2010: Prediction of lightning flash density with the WRF model,
  !!           Adv. Geosci., 23, 11-16
  !! adapted from the COSMO-implementation by Uli Blahak.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-05-27) 
  !!
  SUBROUTINE compute_field_LPI( ptr_patch, jg, ptr_patch_local_parent, p_int,   &
                                p_metrics, p_prog, p_prog_rcf, p_diag,          &
                                lpi )

    IMPLICIT NONE

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg                ! domain ID of main grid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch_local_parent  !< parent grid for larger exchange halo
    TYPE(t_int_state),  INTENT(IN)    :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    REAL(wp),           INTENT(OUT)   :: lpi(:,:)          !< output variable, dim: (nproma,nblks_c)

    TYPE(t_patch),  POINTER  :: p_pp

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx

    REAL(wp) :: delta_z
    REAL(wp) :: q_liqu, q_solid, epsw
    REAL(wp) :: lpi_incr
    REAL(wp) :: w_c
    REAL(wp) :: Tmelt_m_20K

    REAL(wp) :: w_updraft_crit
    REAL(wp) :: qg_low_limit, qg_upp_limit, qg_scale_inv, q_i, q_s, q_g
    REAL(wp) :: w_thresh, w_frac_tresh

    REAL(wp) :: vol( nproma )

    INTEGER :: jc, jk, jb
    INTEGER :: jc2, jb2
    INTEGER :: l, i
    INTEGER :: ist

    INTEGER :: nmbr_w  ( nproma, ptr_patch%nblks_c )
    !INTEGER :: nmbr_all( nproma, ptr_patch%nblks_c )  ! only necessary for 'variant 1 of updraft in environment crit.'

    INTEGER, ALLOCATABLE :: p_nmbr_w(:,:)
    INTEGER  :: p_nmbr_w_sum  (nproma)
    INTEGER  :: p_nmbr_all_sum(nproma)
    REAL(wp) :: p_frac_w      (nproma)

    REAL(wp) :: frac_w( nproma, ptr_patch%nblks_c )

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
    INTEGER :: nblks_c_lp


    IF (.not.atm_phy_nwp_config(jg)%lhave_graupel) THEN
      CALL finish( modname//'compute_field_LPI',  &
        &     "no graupel available! Either switch off LPI output or change the microphysics scheme" )
    END IF

    Tmelt_m_20K = Tmelt - 20.0_wp

    ! --- Thresholds, limits, ... for several criteria: ---

    ! for the 'updraft criterion' in the LPI integral:
    w_updraft_crit = 0.5_wp  ! in m/s; threshold for w

    ! for the 'Graupel-criterion' in the LPI integral:
    qg_low_limit = 0.0002_wp  ! in kg/kg; lower limit for the 'Graupel-criterion'
    qg_upp_limit = 0.001_wp   ! in kg/kg; upper limit for the 'Graupel-criterion'

    qg_scale_inv = 1.0_wp / ( qg_upp_limit - qg_low_limit )


!!$ UB: need to check this threshold for ICON-D2!
    ! for the 'updraft in environment'-criterion:
    w_thresh = 0.5_wp  ! in m/s; threshold for w
                       ! see Lynn, Yair (2010)
                       ! however, Uli Blahak determined 1.1 m/s for COSMO

    w_frac_tresh = 0.5_wp   ! in 100%; threshold for the fraction of points with 'w>w_thresh'
                            ! to identify growing thunderstorms
                            ! Lynn, Yair (2010) propose 0.5 (however, the whole criterion is unclear)


    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ! nullify every grid point (lateral boundary, too)
    lpi(:,:) = 0.0_wp

    ! --- calculation of the LPI integral ---

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,vol,delta_z,w_c,q_liqu,q_i,q_s,q_g,q_solid,epsw,lpi_incr), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        vol     (jc)     = 0.0_wp
        nmbr_w  (jc, jb) = 0
       !nmbr_all(jc, jb) = 0  ! only necessary for 'variant 1 of updraft in environment crit.'
      END DO

      DO jk = kstart_moist(jg), ptr_patch%nlev
        DO jc = i_startidx, i_endidx

          IF (      (p_diag%temp(jc,jk,jb) <= Tmelt)     &
              .AND. (p_diag%temp(jc,jk,jb) >= Tmelt_m_20K) ) THEN

            delta_z = p_metrics%ddqz_z_full(jc,jk,jb)
            w_c     = 0.5_wp * ( p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb) )

            q_liqu  = p_prog_rcf%tracer(jc,jk,jb,iqc)   &
                    + p_prog_rcf%tracer(jc,jk,jb,iqr)

            IF (atm_phy_nwp_config(jg)%l2moment) THEN
              ! sometimes during IAU slightly negative values of the hydrometeors were encountered, which
              ! lead to a crash in the sqrt() below. Therefore we clip all hydrometeors at 0.0 for security:
              q_i = MAX(p_prog_rcf%tracer(jc,jk,jb,iqi), 0.0_wp)
              q_s = MAX(p_prog_rcf%tracer(jc,jk,jb,iqs), 0.0_wp)
              q_g = MAX(p_prog_rcf%tracer(jc,jk,jb,iqg) + p_prog_rcf%tracer(jc,jk,jb,iqh), 0.0_wp)
            ELSE
              q_i = MAX(p_prog_rcf%tracer(jc,jk,jb,iqi), 0.0_wp)
              q_s = MAX(p_prog_rcf%tracer(jc,jk,jb,iqs), 0.0_wp)
              q_g = MAX(p_prog_rcf%tracer(jc,jk,jb,iqg), 0.0_wp)
            END IF

            IF (atm_phy_nwp_config(jg)%inwp_gscp == 7) THEN
              q_g = q_g + p_prog_rcf%tracer(jc,jk,jb,iqgl) + p_prog_rcf%tracer(jc,jk,jb,iqhl)
            END IF
            
            q_solid = q_g *                                                 &
                 &    ( SQRT( q_i * q_g  ) / MAX( q_i + q_g, 1.0e-20_wp) +  &
                 &      SQRT( q_s * q_g  ) / MAX( q_s + q_g, 1.0e-20_wp) )

            q_solid = q_g *                                                 &
                 &    ( SQRT( q_i * q_g  ) / MAX( q_i + q_g, 1.0e-20_wp) +  &
                 &      SQRT( q_s * q_g  ) / MAX( q_s + q_g, 1.0e-20_wp) )

            q_solid = q_g *                                                 &
                 &    ( SQRT( q_i * q_g  ) / MAX( q_i + q_g, 1.0e-20_wp) +  &
                 &      SQRT( q_s * q_g  ) / MAX( q_s + q_g, 1.0e-20_wp) )

            epsw = 2.0_wp * SQRT( q_liqu * q_solid ) / MAX( q_liqu + q_solid, 1.0e-20_wp)

            ! 'updraft-criterion' in the LPI integral
            IF ( w_c >= w_updraft_crit ) THEN
              ! only (strong enough) updrafts should be counted, no downdrafts
              lpi_incr = w_c * w_c * epsw * delta_z

              ! additional 'Graupel-criterion' in the LPI integral
              lpi_incr = lpi_incr * MAX( MIN( (q_g - qg_low_limit) * qg_scale_inv, 1.0_wp ), 0.0_wp )

            ELSE
              lpi_incr = 0.0_wp
            END IF

            lpi(jc,jb) = lpi(jc,jb) + lpi_incr
            vol(jc)    = vol(jc)    + delta_z

            ! sum up for the later test on growth of thunderstorm ('updraft in environment crit.'):

            ! 'variant 1 of updraft in environment crit.' (3D-criterion, MB)
            !IF (w_c >= w_thresh) THEN
            !  nmbr_w(jc,jb) = nmbr_w(jc,jb) + 1
            !END IF
            !nmbr_all(jc,jb) = nmbr_all(jc,jb) + 1

            ! 'variant 2 of updraft in environment crit.' (2D-criterion, UB)
            IF (w_c >= w_thresh) THEN
              nmbr_w(jc,jb) = 1
            END IF

          END IF

        END DO
      END DO

      ! normalization
      DO jc = i_startidx, i_endidx
        IF ( vol(jc) > 1.0e-30_wp ) THEN
          lpi(jc,jb) = lpi(jc,jb) / vol(jc)
        ELSE
          lpi(jc,jb) = 0.0_wp
        END IF
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! --- 'updraft in environment criterion' ---


    ! spatial filtering, i.e. 
    ! test for growth phase of the thunderstorm
    ! (remark: currently only 'variant 2' (2D-criterion) is implemented here)

    ! --- average nmbr_w to the parent grid cells

    p_pp => ptr_patch_local_parent

    ! Set pointers to index and coefficient fields for cell-based variables
    iidx => p_pp%cells%child_idx
    iblk => p_pp%cells%child_blk

    nblks_c_lp = p_pp%cells%end_block( min_rlcell )

    ALLOCATE( p_nmbr_w ( nproma, nblks_c_lp), STAT=ist )
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_lpi', "allocate failed" )
    END IF

    ! first nullify all lateral grid points

    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_bdyintp_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = grf_fbk_start_c + 1

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        p_nmbr_w( jc, jb) = 0.0_wp
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    ! Start/End block in the parent domain
    IF (jg == 1 .AND. l_limited_area) THEN
      i_rlstart = grf_fbk_start_c
    ELSE
      i_rlstart = grf_ovlparea_start_c
    ENDIF
    i_rlend   = min_rlcell_int

    i_startblk = p_pp%cells%start_block( i_rlstart )
    i_endblk   = p_pp%cells%end_block  ( i_rlend   )


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        p_nmbr_w(jc, jb) =                            &
              nmbr_w( iidx(jc,jb,1), iblk(jc,jb,1) )  &
            + nmbr_w( iidx(jc,jb,2), iblk(jc,jb,2) )  &
            + nmbr_w( iidx(jc,jb,3), iblk(jc,jb,3) )  &
            + nmbr_w( iidx(jc,jb,4), iblk(jc,jb,4) ) 
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! --- Exchange of these fields on the parent grid
    CALL exchange_data( p_pp%comm_pat_c, recv=p_nmbr_w )

    ! --- Average over the neighbouring parent grid cells

    ! consistency check
    IF ( p_int%cell_environ%max_nmbr_iter /= 1 ) THEN
      CALL finish( modname//':compute_field_lpi', "cell_environ is not built with the right number of iterations" )
    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,  &
!$OMP            p_nmbr_w_sum, p_nmbr_all_sum,  &
!$OMP            l,jc2,jb2,p_frac_w,i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_pp, jb, i_startblk, i_endblk,           &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        p_nmbr_w_sum  (jc)  = 0
        p_nmbr_all_sum(jc) = 0
      END DO

      DO l=1, p_int%cell_environ%max_nmbr_nghbr_cells
        DO jc = i_startidx, i_endidx

          jc2 = p_int%cell_environ%idx( jc, jb, l)
          jb2 = p_int%cell_environ%blk( jc, jb, l)

          IF ( p_int%cell_environ%area_norm( jc, jb, l) > 1.0e-7_wp ) THEN
            p_nmbr_w_sum  (jc) = p_nmbr_w_sum  (jc) + p_nmbr_w(jc2, jb2)
            p_nmbr_all_sum(jc) = p_nmbr_all_sum(jc) + 4
          END IF
        END DO
      END DO

      DO jc = i_startidx, i_endidx
        p_frac_w(jc) = DBLE( p_nmbr_w_sum(jc) ) / DBLE( p_nmbr_all_sum(jc) )
      END DO

      ! write back to the 4 child grid cells:
      DO i = 1, 4
        DO jc = i_startidx, i_endidx
          frac_w( iidx(jc,jb,i), iblk(jc,jb,i) ) = p_frac_w(jc)
        END DO
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE( p_nmbr_w, STAT=ist )
    IF ( ist /= 0 ) THEN
      CALL finish( modname//':compute_field_lpi', "deallocate failed" )
    END IF

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc = i_startidx, i_endidx
        ! finally, this is the 'updraft in environment criterion':
        IF ( frac_w(jc,jb) < w_frac_tresh ) THEN
          lpi(jc, jb) = 0.0_wp
        END IF

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_LPI


  !>
  !! Do a maximization step for the calculation of the LPI_MAX.
  !!
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-09-17) 
  !!
  SUBROUTINE maximize_field_LPI( ptr_patch, jg, ptr_patch_local_parent, p_int,   &
                                p_metrics, p_prog, p_prog_rcf, p_diag,           &
                                lpi_max )

    IMPLICIT NONE

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg                ! domain ID of main grid
    TYPE(t_patch), TARGET, INTENT(IN) :: ptr_patch_local_parent  !< parent grid for larger exchange halo
    TYPE(t_int_state),  INTENT(IN)    :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    REAL(wp), INTENT(INOUT) :: lpi_max( nproma, ptr_patch%nblks_c )

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc

    REAL(wp) :: lpi( nproma, ptr_patch%nblks_c )

    CALL compute_field_LPI( ptr_patch, jg, ptr_patch_local_parent, p_int,   &
                            p_metrics, p_prog, p_prog_rcf, p_diag,                &
                            lpi )

    ! Maximize

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int
    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        lpi_max(jc,jb) = MAX( lpi_max(jc,jb), lpi(jc,jb) )
      END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE maximize_field_LPI


  !>
  !! Calculate the ceiling height
  !! = height above MSL, for which cloud coverage > 4/8
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-21) 
  !!
  SUBROUTINE compute_field_ceiling( ptr_patch, jg,    &
                                p_metrics, prm_diag,  &
                                ceiling_height )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg            ! domain ID of main grid
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag

    REAL(wp),             INTENT(OUT) :: ceiling_height(:,:)    !< output variable, dim: (nproma,nblks_c)

    LOGICAL ::  cld_base_found( nproma ) 

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,cld_base_found), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        cld_base_found(jc) = .FALSE.
        ceiling_height(jc,jb) = p_metrics%z_mc(jc,1,jb)  ! arbitrary default value
      END DO

      DO jk = ptr_patch%nlev, kstart_moist(jg), -1

        DO jc = i_startidx, i_endidx

          IF ( .NOT.(cld_base_found(jc)) .AND. (prm_diag%clc(jc,jk,jb) > 0.5_wp) ) THEN
            ceiling_height(jc,jb) = p_metrics%z_mc(jc,jk,jb)
            cld_base_found(jc) = .TRUE.
          ENDIF

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_ceiling


  !>
  !! Calculate the height of base over MSL from the shallow convection parameterization.
  !!
  !! This subroutine is quite similar to compute_field_htop_sc.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-22) 
  !!
  SUBROUTINE compute_field_hbas_sc( ptr_patch,        &
                                p_metrics, prm_diag,  &
                                hbas_sc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag

    REAL(wp),             INTENT(OUT) :: hbas_sc(:,:)    !< output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc, idx

    REAL(wp), PARAMETER :: undefValue = 0.0_wp

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,idx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF ( prm_diag%ktype(jc,jb) == 2 ) THEN
          ! the column is identified as shallow convection
          idx = prm_diag%mbas_con( jc, jb)
          hbas_sc(jc,jb) = p_metrics%z_mc( jc, idx, jb)
       ELSE
          hbas_sc(jc,jb) = undefValue
        END IF

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_hbas_sc


  !>
  !! Calculate the height of top over MSL from the shallow convection parameterization.
  !!
  !! This subroutine is quite similar to compute_field_hbas_sc.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-22) 
  !!
  SUBROUTINE compute_field_htop_sc( ptr_patch,        &
                                p_metrics, prm_diag,  &
                                htop_sc )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag

    REAL(wp),             INTENT(OUT) :: htop_sc(:,:)    !< output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jc, idx

    REAL(wp), PARAMETER :: undefValue = 0.0_wp

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,idx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx

        IF ( prm_diag%ktype(jc,jb) == 2 ) THEN
          ! the column is identified as shallow convection
          idx = prm_diag%mtop_con( jc, jb)
          htop_sc(jc,jb) = p_metrics%z_mc( jc, idx, jb)
       ELSE
          htop_sc(jc,jb) = undefValue
        END IF

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_htop_sc


  !>
  !! Calculate total column integrated water (twater)
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-23) 
  !!
  SUBROUTINE compute_field_twater( ptr_patch, jg,      &
                                   p_metrics, p_prog, p_prog_rcf,  &
                                   twater )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg            ! domain ID of main grid
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics
    TYPE(t_nh_prog),      INTENT(IN)  :: p_prog, p_prog_rcf

    REAL(wp),             INTENT(OUT) :: twater(:,:)    !< output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: q_water( nproma, ptr_patch%nlev )

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    twater( :, 1:i_startblk-1 ) = 0.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,q_water), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, ptr_patch%nlev
        DO jc = i_startidx, i_endidx
          q_water(jc,jk) = p_prog_rcf%tracer(jc,jk,jb,iqv)   &
            &            + p_prog_rcf%tracer(jc,jk,jb,iqc)
        END DO
      END DO

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqi)%p_3d ) ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_water(jc,jk) = q_water(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqi)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqr)%p_3d ) ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_water(jc,jk) = q_water(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqr)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqs)%p_3d ) ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_water(jc,jk) = q_water(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqs)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_water(jc,jk) = q_water(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqg)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_water(jc,jk) = q_water(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqh)
          END DO
        END DO
      END IF

      ! calculate vertically integrated mass

      twater(:, jb) = 0.0_wp
      DO jk = 1, ptr_patch%nlev
        DO jc = i_startidx, i_endidx
          twater(jc,jb) = twater(jc,jb)       &
            &            + p_prog%rho(jc,jk,jb) * q_water(jc,jk) * p_metrics%ddqz_z_full(jc,jk,jb)
        END DO
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_twater


  !>
  !! Calculate specific content of precipitation particles
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-23) 
  !!
  SUBROUTINE compute_field_q_sedim( ptr_patch, jg, p_prog_rcf, q_sedim )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch     !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg            ! domain ID of main grid
    TYPE(t_nh_prog),      INTENT(IN)  :: p_prog_rcf

    REAL(wp),             INTENT(OUT) :: q_sedim(:,:,:)  !< output variable, dim: (nproma,nlev,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    q_sedim( :, :, 1:i_startblk-1 ) = 0.0_wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! it is assumed that at least a warm rain Kessler scheme is used, i.e. qr is available too
      DO jk = kstart_moist(jg), ptr_patch%nlev
        DO jc = i_startidx, i_endidx
          q_sedim(jc,jk,jb) = p_prog_rcf%tracer(jc,jk,jb,iqr)
        END DO
      END DO

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqs)%p_3d ) ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_sedim(jc,jk,jb) = q_sedim(jc,jk,jb) + p_prog_rcf%tracer(jc,jk,jb,iqs)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_sedim(jc,jk,jb) = q_sedim(jc,jk,jb) + p_prog_rcf%tracer(jc,jk,jb,iqg)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_sedim(jc,jk,jb) = q_sedim(jc,jk,jb) + p_prog_rcf%tracer(jc,jk,jb,iqh)
          END DO
        END DO
      END IF

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_q_sedim


  !>
  !! Calculate 
  !!     TCOND_MAX   (total column-integrated condensate, max. during the last hour)
  !! and TCOND10_MAX (total column-integrated condensate above z(T=-10 degC), max. during the last hour)
  !! Here, compute columnwise amximum of these input fields and the newly computed fields.
  !! 
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-23) 
  !!
  SUBROUTINE compute_field_tcond_max( ptr_patch, jg,                      &
                                   p_metrics, p_prog, p_prog_rcf, p_diag, &
                                   flag_tcond_max, flag_tcond10_max,      &
                                   tcond_max,      tcond10_max )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    INTEGER,            INTENT(IN)    :: jg                ! domain ID of main grid
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    LOGICAL,            INTENT(IN)    :: flag_tcond_max    ! if true, then calculate tcond_max
    LOGICAL,            INTENT(IN)    :: flag_tcond10_max  ! if true, then calculate tcond10_max

    REAL(wp),           INTENT(INOUT) :: tcond_max  (:,:)    !< output variable, dim: (nproma,nblks_c)
    REAL(wp),           INTENT(INOUT) :: tcond10_max(:,:)    !< output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: q_cond( nproma, ptr_patch%nlev )
    REAL(wp) :: tcond( nproma )

    REAL(wp), PARAMETER :: Tzero_m_10K = 263.15    ! in K

    ! Consistency check:
    IF ( (.NOT. flag_tcond_max) .AND. (.NOT. flag_tcond10_max) ) THEN
      CALL finish( "compute_field_tcond_max", "at least one of the two flags must be set to .TRUE." )
    END IF

    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,q_cond,tcond), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      ! it is assumed that at least a warm rain Kessler scheme is used, i.e. qr is available too
      DO jk = kstart_moist(jg), ptr_patch%nlev
        DO jc = i_startidx, i_endidx
          q_cond(jc,jk) = p_prog_rcf%tracer(jc,jk,jb,iqc)  &
            &           + p_prog_rcf%tracer(jc,jk,jb,iqr)
        END DO
      END DO

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqi)%p_3d ) ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqi)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED( p_prog_rcf%tracer_ptr(iqs)%p_3d ) ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqs)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqg)
          END DO
        END DO
      END IF

      IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            q_cond(jc,jk) = q_cond(jc,jk) + p_prog_rcf%tracer(jc,jk,jb,iqh)
          END DO
        END DO
      END IF

      ! calculate vertically integrated mass

      IF ( flag_tcond_max ) THEN
        tcond( i_startidx: i_endidx) = 0.0_wp
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            tcond(jc) = tcond(jc)       &
              &       + p_prog%rho(jc,jk,jb) * q_cond(jc,jk) * p_metrics%ddqz_z_full(jc,jk,jb)
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          tcond_max( jc,jb ) = MAX( tcond_max( jc,jb ), tcond(jc) )
        END DO
      END IF

      IF ( flag_tcond10_max ) THEN
        tcond( i_startidx: i_endidx) = 0.0_wp
        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            IF ( p_diag%temp( jc,jk,jb) <= Tzero_m_10K ) THEN
              tcond(jc) = tcond(jc)       &
                &       + p_prog%rho(jc,jk,jb) * q_cond(jc,jk) * p_metrics%ddqz_z_full(jc,jk,jb)
            END IF
          END DO
        END DO

        DO jc = i_startidx, i_endidx
          tcond10_max( jc,jb ) = MAX( tcond10_max( jc,jb ), tcond(jc) )
        END DO
      END IF

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE compute_field_tcond_max


  !>
  !! Calculate UH_MAX (updraft helicity, max.  during the last hour)
  !! For the definition see: Kain et al. (2008) Wea. Forecasting
  !!
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-23) 
  !! Inserted variable boundaries for vertical integration, Vera Maurer, DWD (2021-03-10)
  !! Vertical integration changed by Uli Blahak, DWD (2021-03-20)
  !!
  SUBROUTINE compute_field_uh_max( ptr_patch,                 &
                                   p_metrics, p_prog, p_diag, &
                                   zmin_in, zmax_in,          &
                                   uh_max )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag

    REAL(wp),           INTENT(IN)    :: zmin_in        !< lower boundary for vertical integration
    REAL(wp),           INTENT(IN)    :: zmax_in        !< upper boundary for vertical integration

    REAL(wp),           INTENT(INOUT) :: uh_max(:,:)    !< input/output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: zmin( nproma )
    REAL(wp) :: zmax( nproma )
    REAL(wp) :: uhel( nproma )
    REAL(wp) :: w_c

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zmin,zmax,uhel,w_c), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        IF (zmin_in < 500._wp) THEN
          zmin(jc) = MAX( p_metrics%z_ifc( jc, ptr_patch%nlev+1, jb), zmin_in )
        ELSE
          zmin(jc) = MAX( p_metrics%z_ifc( jc, ptr_patch%nlev+1, jb) + 500.0_wp, zmin_in )
        END IF
        zmax(jc) = zmin(jc) + zmax_in - zmin_in
      END DO

      uhel( i_startidx:i_endidx ) = 0.0_wp
      DO jk = 1, ptr_patch%nlev
        DO jc = i_startidx, i_endidx

          ! Parts of the grid boy are within the bounds, integrate over the exact bounds [zmin,zmax]:
          !  (It also works if the integration layer is so narrow that the bounds are in the same grid box)
          IF ( ( p_metrics%z_ifc( jc, jk+1, jb) <= zmax(jc) ) .AND.     &
            &  ( p_metrics%z_ifc( jc, jk, jb)   >= zmin(jc) ) ) THEN

            w_c = 0.5_wp * ( p_prog%w(jc,jk,jb) + p_prog%w(jc,jk+1,jb) )
            
            ! a simple box-integration in the vertical, but honouring the exact integration bounds zmin, zmax;
            ! only updrafts are counted:
            uhel(jc) = uhel(jc) + MAX( w_c, 0.0_wp) * p_diag%vor(jc,jk,jb) * &
                 ( MIN(p_metrics%z_ifc(jc,jk,jb), zmax(jc)) - MAX(p_metrics%z_ifc( jc, jk+1, jb), zmin(jc)) )
            
          END IF

        END DO
      END DO

      DO jc = i_startidx, i_endidx
        uh_max(jc,jb) = MERGE(uhel(jc), uh_max(jc,jb), ABS(uhel(jc)) > ABS(uh_max(jc,jb)) )
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_uh_max


  !>
  !! Calculate VORW_CTMAX (Maximum rotation amplitude during the last hour)
  !!
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-23)
  !! Vertical integration changed by Uli Blahak, DWD (2021-03-20)
  !!
  SUBROUTINE compute_field_vorw_ctmax( ptr_patch,          &
                                       p_metrics, p_diag,  &
                                       vorw_ctmax )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),    INTENT(IN)    :: p_diag

    REAL(wp),           INTENT(INOUT) :: vorw_ctmax(:,:)  !< input/output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    REAL(wp) :: zmin( nproma )
    REAL(wp) :: zmax( nproma )
    REAL(wp) :: vort( nproma )

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zmin,zmax,vort), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        zmin(jc) = p_metrics%z_ifc( jc, ptr_patch%nlev+1, jb)
        zmax(jc) = MAX( 3000.0_wp,  zmin(jc) + 1500.0_wp )
      END DO

      vort( i_startidx:i_endidx ) = 0.0_wp
      DO jk = 1, ptr_patch%nlev
        DO jc = i_startidx, i_endidx

          ! Parts of the grid box are within the bounds, integrate over the exact bounds [zmin,zmax]:
          !  (It also works if the integration layer is so narrow that the bounds are in the same grid box)
          IF ( ( p_metrics%z_ifc( jc, jk+1, jb) <= zmax(jc) ) .AND.     &
            &  ( p_metrics%z_ifc( jc, jk, jb)   >= zmin(jc) ) ) THEN

            ! simple box-integration in the vertical, but honouring the exact integration bounds zmin, zmax:
            vort(jc) = vort(jc) + p_diag%vor(jc,jk,jb) * &
                 ( MIN(p_metrics%z_ifc(jc,jk,jb), zmax(jc)) - MAX(p_metrics%z_ifc( jc, jk+1, jb), zmin(jc)) )
            
          END IF

        END DO
      END DO

      DO jc = i_startidx, i_endidx
        vort(jc) = vort(jc) / (zmax(jc) - zmin(jc))
        vorw_ctmax(jc,jb) = MERGE(vort(jc), vorw_ctmax(jc,jb), ABS(vort(jc)) > ABS(vorw_ctmax(jc,jb)) )
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_vorw_ctmax


  !>
  !! Calculate W_CTMAX (Maximum updraft track during the last hour)
  !!
  !! Implementation analogous to those of Uli Blahak in COSMO.
  !!
  !! @par Revision History
  !! Initial revision by Michael Baldauf, DWD (2019-10-23) 
  !!
  SUBROUTINE compute_field_w_ctmax( ptr_patch,             &
                                    p_metrics, p_prog,     &
                                    w_ctmax )

    TYPE(t_patch),      INTENT(IN)    :: ptr_patch         !< patch on which computation is performed
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog

    REAL(wp),           INTENT(INOUT) :: w_ctmax(:,:)    !< input/output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, ptr_patch%nlev
        DO jc = i_startidx, i_endidx

          IF ( p_metrics%z_mc( jc, jk, jb) <= 10000.0_wp ) THEN
            w_ctmax(jc,jb) = MAX( w_ctmax(jc,jb), p_prog%w(jc,jk,jb) )
          END IF

        END DO
      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_w_ctmax


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
    i_startidx, i_endidx,  &  ! start and end indices of loops in horizontal patch
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
#ifdef __SX__
  REAL    (wp)             :: tguess1v(SIZE(te,1))
  LOGICAL                  :: lcalc(SIZE(te,1))
#endif
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
  LOGICAL :: lexit(SIZE(te,1))
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
    lexit(:)     = .FALSE.

    ! now calculate the mixed layer average potential temperature and 
    ! specific humidity
    DO k = nlev, kmoist, -1
      IF (ALL(lexit(i_startidx:i_endidx))) EXIT
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
          lexit(i) = .TRUE.
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

#ifdef __SX__
    ENDDO ! i = i_startidx, i_endidx

    ! Vectorized version
    DO i = i_startidx, i_endidx
      lcalc(i) = .FALSE.
      IF ( k > kstart(i) ) CYCLE
      IF ( k <= lcllev(i) ) THEN
        ! The scheme uses a first guess temperature, which is the parcel
        ! temperature at the level below. If it happens that the initial
        ! parcel is already saturated, the environmental temperature
        ! is taken as first guess instead
        IF (  k == kstart(i) ) THEN
          tguess1v(i) = te(i,kstart(i))
        ELSE
          tguess1v(i) = tp(i)
        END IF
        lcalc(i) = .TRUE.
      ENDIF
    ENDDO ! i = i_startidx, i_endidx


    ! Calculate iteratively parcel temperature from thp, prs and 1st guess tguess1
    DO icount = 1, 21
      IF (COUNT(lcalc(i_startidx:i_endidx)) > 0) THEN
        DO i = i_startidx, i_endidx
          IF ( lcalc(i) ) THEN
            esat     = sat_pres_water( tguess1v(i))
            q1       = fqvs( esat, prs(i,k), qvp(i) )
            thetae1  = fthetae( tguess1v(i),prs(i,k),q1)

            tguess2  = tguess1v(i) - 1.0_wp
            esat     = sat_pres_water( tguess2)
            q2       = fqvs( esat, prs(i,k), qvp(i) )
            thetae2  = fthetae( tguess2,prs(i,k),q2)

            tguess1v(i)  = tguess1v(i)+(thetae1-thp(i))/(thetae2-thetae1)

            IF ( ABS( thetae1-thp(i)) < eps .OR. icount > 20) THEN
              tp(i) = tguess1v(i)
              lcalc(i) = .false.
            END IF
          END IF
        ENDDO ! i = i_startidx, i_endidx
      ELSE
        EXIT
      ENDIF
    END DO

    ! update specific humidity of the saturated parcel for new temperature
    DO i = i_startidx, i_endidx
      IF ( k > kstart(i) ) CYCLE
      IF ( k <= lcllev(i) ) THEN
        esatp  = sat_pres_water( tp(i))
        qvp(i) = fqvs( esatp,prs(i,k),qvp(i))
      END IF
    ENDDO ! i = i_startidx, i_endidx

    DO i = i_startidx, i_endidx
      IF ( k > kstart(i) ) CYCLE

#else

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
#endif       
  
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
    fthetae = (p0/zpx)**rd_o_cpd *ztx*EXP( alvdcp*zqx/(ztx*(1._wp-zqx)) )

  END FUNCTION fthetae

  END SUBROUTINE cal_cape_cin

  !> Wrapper routine to get the 3D radar reflectivity field depending on the microphysics scheme
  !  and store it in p_diag%dbz3d(:,:,:)
  !!
  !! @par Revision History
  !! Initial revision  :  U. Blahak, DWD (2020-01-20) 
  SUBROUTINE compute_field_dbz3d_lin(jg, ptr_patch, p_prog,  p_prog_rcf, p_diag, prm_diag, dbz3d_lin)

    ! Domain index for later use with EMVORADO calc_dbz_vec():
    INTEGER, INTENT(in)  :: jg
    ! patch on which computation is performed:
    TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch
    ! nonhydrostatic state
    TYPE(t_nh_prog), INTENT(IN)       :: p_prog            !< at timelevel nnow(jg)
    TYPE(t_nh_prog), INTENT(IN)       :: p_prog_rcf        !< at timelevel nnow_rcf(jg)
    TYPE(t_nh_diag), INTENT(INOUT)    :: p_diag
    TYPE(t_nwp_phy_diag), INTENT(IN)  :: prm_diag
    REAL(wp),        INTENT(OUT)      :: dbz3d_lin(:,:,:)  !< reflectivity in mm^6/m^3
    
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_dbz3d_lin'
    REAL(wp) :: rho, qnc_s(nproma,ptr_patch%nblks_c)
    INTEGER  :: i_rlstart, i_rlend, i_startblk, i_endblk, i_startidx, i_endidx, i_startidx_1, i_endidx_2, &
      &         jc, jk, jb

    ! without halo or boundary points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    ! Determine first and last index for first and last block, respectively
    CALL get_indices_c( ptr_patch, i_startblk, i_startblk, i_endblk, i_startidx_1, i_endidx, i_rlstart, i_rlend)
    CALL get_indices_c( ptr_patch, i_endblk, i_startblk, i_endblk, i_startidx, i_endidx_2, i_rlstart, i_rlend)

    SELECT CASE ( atm_phy_nwp_config(jg)%inwp_gscp )
    CASE ( 1 )

      IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 2) THEN
        ! Not yet implemented in microphysics! We give a dummy value here.
        qnc_s(:,:) = cloud_num               ! 1/kg
      ELSE IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 1) THEN
        qnc_s(:,:) = prm_diag%cloud_num(:,:) ! neglect difference of 1/m^3 and 1/kg for this near-surface value
      ELSE
        qnc_s(:,:) = cloud_num               ! 1/kg
      END IF
      
      CALL compute_field_dbz_1mom( npr       = nproma,                           &
                                   nlev      = ptr_patch%nlev,                   &
                                   nblks     = ptr_patch%nblks_c,                &
                                   startblk  = i_startblk,                       &
                                   endblk    = i_endblk,                         &
                                   jk_start  = kstart_moist(jg),                 &
                                   startidx1 = i_startidx_1,                     &
                                   endidx2   = i_endidx_2,                       &
                                   lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
                                   lmessage_full  = (msg_level > 15),            &
                                   my_id_for_message = get_my_mpi_work_id(),     &
                                   rho_w     = rhoh2o,                           &
                                   rho_ice   = rhoice,                           &
                                   K_w       = K_w_0,                            &
                                   K_ice     = K_i_0,                            &
                                   T_melt    = Tmelt,                            &
                                   igscp     = atm_phy_nwp_config(jg)%inwp_gscp, &
                                   q_crit_radar = 1e-8_wp,                       &
                                   T         = p_diag%temp(:,:,:),               &
                                   rho       = p_prog%rho(:,:,:),                &
                                   q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
                                   q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
                                   q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
                                   q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
                                   n_cloud_s = qnc_s(:,:),                       &  ! 1/kg
                                   z_radar   = dbz3d_lin(:,:,:)                  )

    CASE ( 2 )

      IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 2) THEN
        ! Not yet implemented in microphysics! We give a dummy value here.
        qnc_s(:,:) = cloud_num               ! 1/kg
      ELSE IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 1) THEN
        qnc_s(:,:) = prm_diag%cloud_num(:,:) ! neglect difference of 1/m^3 and 1/kg for this near-surface value
      ELSE
        qnc_s(:,:) = cloud_num               ! 1/kg
      END IF
      
      CALL compute_field_dbz_1mom( npr       = nproma,                           &
                                   nlev      = ptr_patch%nlev,                   &
                                   nblks     = ptr_patch%nblks_c,                &
                                   startblk  = i_startblk,                       &
                                   endblk    = i_endblk,                         &
                                   jk_start  = kstart_moist(jg),                 &
                                   startidx1 = i_startidx_1,                     &
                                   endidx2   = i_endidx_2,                       &
                                   lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
                                   lmessage_full  = (msg_level > 15),            &
                                   my_id_for_message = get_my_mpi_work_id(),     &
                                   rho_w     = rhoh2o,                           &
                                   rho_ice   = rhoice,                           &
                                   K_w       = K_w_0,                            &
                                   K_ice     = K_i_0,                            &
                                   T_melt    = Tmelt,                            &
                                   igscp     = atm_phy_nwp_config(jg)%inwp_gscp, &
                                   q_crit_radar = 1e-8_wp,                       &
                                   T         = p_diag%temp(:,:,:),               &
                                   rho       = p_prog%rho(:,:,:),                &
                                   q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
                                   q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
                                   q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
                                   q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
                                   q_graupel = p_prog_rcf%tracer(:,:,:,iqg),     &
                                   n_cloud_s = qnc_s(:,:),                       &  ! 1/kg
                                   z_radar   = dbz3d_lin(:,:,:)                  )
      
    CASE ( 4, 5, 6 )

      CALL compute_field_dbz_2mom( npr       = nproma,                           &
                                   nlev      = ptr_patch%nlev,                   &
                                   nblks     = ptr_patch%nblks_c,                &
                                   startblk  = i_startblk,                       &
                                   endblk    = i_endblk,                         &
                                   jk_start  = kstart_moist(jg),                 &
                                   startidx1 = i_startidx_1,                     &
                                   endidx2   = i_endidx_2,                       &
                                   lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
                                   lmessage_full  = (msg_level > 15),            &
                                   my_id_for_message = get_my_mpi_work_id(),     &
                                   rho_w     = rhoh2o,                           &
                                   rho_ice   = rhoice,                           &
                                   K_w       = K_w_0,                            &
                                   K_ice     = K_i_0,                            &
                                   T_melt    = Tmelt,                            &
                                   q_crit_radar = 1e-8_wp,                       &
                                   T         = p_diag%temp(:,:,:),               &
                                   rho       = p_prog%rho(:,:,:),                &
                                   q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
                                   q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
                                   q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
                                   q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
                                   q_graupel = p_prog_rcf%tracer(:,:,:,iqg),     &
                                   q_hail    = p_prog_rcf%tracer(:,:,:,iqh),     &
                                   n_cloud   = p_prog_rcf%tracer(:,:,:,iqnc),    &
                                   n_ice     = p_prog_rcf%tracer(:,:,:,iqni),    &
                                   n_rain    = p_prog_rcf%tracer(:,:,:,iqnr),    &
                                   n_snow    = p_prog_rcf%tracer(:,:,:,iqns),    &
                                   n_graupel = p_prog_rcf%tracer(:,:,:,iqng),    &
                                   n_hail    = p_prog_rcf%tracer(:,:,:,iqnh),    &
                                   z_radar   = dbz3d_lin(:,:,:)                  )

      
    CASE ( 7 )

      CALL compute_field_dbz_2mom( npr       = nproma,                           &
                                   nlev      = ptr_patch%nlev,                   &
                                   nblks     = ptr_patch%nblks_c,                &
                                   startblk  = i_startblk,                       &
                                   endblk    = i_endblk,                         &
                                   jk_start  = kstart_moist(jg),                 &
                                   startidx1 = i_startidx_1,                     &
                                   endidx2   = i_endidx_2,                       &
                                   lmessage_light = (msg_level > 12 .AND. my_process_is_mpi_workroot()), &
                                   lmessage_full  = (msg_level > 15),            &
                                   my_id_for_message = get_my_mpi_work_id(),     &
                                   rho_w     = rhoh2o,                           &
                                   rho_ice   = rhoice,                           &
                                   K_w       = K_w_0,                            &
                                   K_ice     = K_i_0,                            &
                                   T_melt    = Tmelt,                            &
                                   q_crit_radar = 1e-8_wp,                       &
                                   T         = p_diag%temp(:,:,:),               &
                                   rho       = p_prog%rho(:,:,:),                &
                                   q_cloud   = p_prog_rcf%tracer(:,:,:,iqc),     &
                                   q_ice     = p_prog_rcf%tracer(:,:,:,iqi),     &
                                   q_rain    = p_prog_rcf%tracer(:,:,:,iqr),     &
                                   q_snow    = p_prog_rcf%tracer(:,:,:,iqs),     &
                                   q_graupel = p_prog_rcf%tracer(:,:,:,iqg),     &
                                   q_hail    = p_prog_rcf%tracer(:,:,:,iqh),     &
                                   n_cloud   = p_prog_rcf%tracer(:,:,:,iqnc),    &
                                   n_ice     = p_prog_rcf%tracer(:,:,:,iqni),    &
                                   n_rain    = p_prog_rcf%tracer(:,:,:,iqnr),    &
                                   n_snow    = p_prog_rcf%tracer(:,:,:,iqns),    &
                                   n_graupel = p_prog_rcf%tracer(:,:,:,iqng),    &
                                   n_hail    = p_prog_rcf%tracer(:,:,:,iqnh),    &
                                   ql_graupel= p_prog_rcf%tracer(:,:,:,iqgl),    &
                                   ql_hail   = p_prog_rcf%tracer(:,:,:,iqhl),    &
                                   z_radar   = dbz3d_lin(:,:,:)                  )

      
    CASE DEFAULT
      
      CALL finish( modname//': get_field_dbz3d_lin',  &
        &     "dbz3d-computation not available for this microphysics scheme! Available for inwp_gscp=1,2,4,5,6 or 7" )

    END SELECT

!#endif
    
  END SUBROUTINE compute_field_dbz3d_lin

  !>
  !! Calculate radar reflectivity for the 1-moment microphysics scheme in linear units mm^6/m^3
  !!
  !! @par Revision History
  !! Initial revision by U. Blahak, DWD (2020-01-20) 
  !!
  SUBROUTINE compute_field_dbz_1mom( npr, nlev, nblks, startblk, endblk, jk_start,       &
                                     startidx1, endidx2,                                 &
                                     lmessage_light, lmessage_full, my_id_for_message,   &
                                     rho_w, rho_ice,                                     &
                                     K_w, K_ice, T_melt, igscp, q_crit_radar,            &
                                     T, rho, q_cloud, q_rain, q_ice, q_snow, z_radar,    &
                                     q_graupel, n_cloud_s )
 
   !------------------------------------------------------------------------------
    !
    ! Description:  Calculation of grid point values for effective radar
    !               reflectivity factor Z in dBZ.
    !
    ! Method:       Rayleigh-Approximation for the Back-Scattering
    !               (no attenuation!), Debye-Approximation for the
    !               effective refractive index of two-component ice-air-
    !               mixture particles (dry ice, snow, and graupel).
    !               Melting particles by substitution of the ice substance
    !               with water when T_air > 273.16 K. For the applied
    !               Rayleigh scattering, this is equal to assuming instantaneous
    !               melting, because only the square of the total water mass of the particle
    !               counts ( = the total square of its dipole moment).
    !               
    ! Inputs:       npr, nlev, nblks : field dimensions
    !               startblk, endblk, jk_start, startidx1, endidx2   : loop start and end indices
    !               rho_w        : bulk density of pure water [kg/m**3]
    !               rho_ice      : bulk density of pure ice   [kg/m**3]
    !               K_w          : dielectric constant of water
    !               K_ice        : dielectric constant of ice
    !               T_melt       : melting temperature of ice
    !               igscp        : itype of grid scale precip (1,2)
    !               q_crit_radar : threshold for the q's to compute reflectivity  [kg/m**3]
    !               T            : temperature field          [K]
    !               rho          : air density                [kg/m**3]
    !               q_cloud      : cloud water mixing ratio   [kg/kg] 
    !               q_rain       : rain water mixing ratio    [kg/kg] 
    !               q_ice        : cloud ice mixing ratio     [kg/kg]
    !               q_snow       : snow mixing ratio          [kg/kg]  
    !   OPTIONAL:   q_graupel    : graupel mixing ratio       [kg/kg] 
    !   OPTIONAL:   n_cloud_s    : surface value of cloud droplet number concentration [1/kg] 
    !
    ! Output:       z_radar      : 3D field of Z              [mm^6/m^3]
    ! 
    !
    !------------------------------------------------------------------------------
    
    ! Input/Output parameters:
    !-------------------------

    INTEGER , INTENT(IN)       :: npr, nlev, nblks
    INTEGER,  INTENT(IN)       :: startblk, endblk, jk_start, startidx1, endidx2, my_id_for_message
    LOGICAL,  INTENT(in)       :: lmessage_light, lmessage_full
    REAL(wp), INTENT(IN)       :: rho_w, rho_ice, K_w, K_ice, T_melt, q_crit_radar
    INTEGER , INTENT(IN)       :: igscp
    REAL(wp), INTENT(IN)       :: T(:,:,:),          &
                                  rho(:,:,:),        &
                                  q_cloud(:,:,:),    &
                                  q_rain(:,:,:),     &
                                  q_ice(:,:,:),      &
                                  q_snow(:,:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: q_graupel(:,:,:), n_cloud_s(:,:)
    REAL(wp), INTENT(OUT)      :: z_radar(:,:,:)

    ! Local Variables
    !----------------

    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_dbz_1mom'
    
    REAL(wp) ::  rho_c, rho_i, rho_r, rho_s, rho_g, z_cloud, z_rain, z_snow, z_ice, z_grau

    INTEGER        :: jc, jk, jb, i_startidx, i_endidx
    LOGICAL, SAVE  :: firstcall = .TRUE.
    logical        :: lqnc_input

    REAL(wp)       :: z_fac_ice_dry, z_fac_ice_wet
    REAL(wp)       :: ztc, m2s, m3s, alf, bet, hlp, zn0s, x_c, x_i_mono, n_i
    INTEGER        :: nn, pe_center


    REAL(wp), SAVE :: z_r, z_g, z_s, p_r, p_s, p_g, mom_fac, z_fac_c
    REAL(wp), SAVE :: nor, nog, amg, bmg, nos, ams, bms, &
                      ami, bmi, D_c_fix, x_c_fix, mue_rain_c

    REAL(wp), PARAMETER :: zn0s1       = 13.5_wp * 5.65E5_wp, & ! parameter in N0S(T)
                           zn0s2       = -0.107_wp,           & ! parameter in N0S(T), Field et al
                           eps         = 1.00E-15_wp,         &
                           n_cloud_min = 1.00E7_wp,           & ! min. allowed cloud number conc. in externally provided n_cloud_s
                           convfac     = 1.E18_wp,            &
                           z10olog10   = 10.0_wp / LOG(10.0_wp)

    TYPE(particle) :: cloud_tmp

    ! Variables for cloud ice parameterization 
    REAL(wp)                             ::  fxna               ! statement function for ice crystal number, Cooper(1986)
    REAL(wp)                             ::  fxna_cooper        ! statement function for ice crystal number, Cooper(1986)
    REAL(wp)                             ::  ztx, ztmelt        ! dummy arguments for statement functions
    REAL(wp), SAVE                       ::  znimax             ! Maximum ice concentration
    ! This is constant in both cloudice and graupel
    LOGICAL, PARAMETER                   ::  lsuper_coolw = .TRUE.

    REAL(wp) :: zdebug

    CHARACTER(len=250) :: message_text
    
!------------------------------------------------------------------------------

    ! Number of activate ice crystals;  ztx is temperature. Statement functions
    fxna       (ztx,ztmelt) = 1.0E2_wp  * EXP(0.2_wp   * (ztmelt - ztx))  ! 1/m^3
    fxna_cooper(ztx,ztmelt) = 5.0E+0_wp * EXP(0.304_wp * (ztmelt - ztx))  ! 1/m^3

!------------------------------------------------------------------------------

    IF (lmessage_light .OR. lmessage_full) THEN
      message_text(:) = ' '
      WRITE(message_text,'(a,i0)') 'Computing dbz3d_lin for inw_gscp = ', igscp
      CALL message(TRIM(routine), TRIM(message_text), all_print=.TRUE.)
    ENDIF
      
    lqnc_input = PRESENT(n_cloud_s)

    z_fac_ice_dry = (rho_w/rho_ice)**2 * K_ice/K_w
    z_fac_ice_wet = 1.0_wp
    
    IF (firstcall) THEN

      mom_fac = (6.0_wp / (pi * rho_w))**2
      
      ! Parameters for cloud droplets:
      !   PSD consistent to autoconversion parameterization of SB2001:
      !   gamma distribution w.r.t. mass with gamma shape parameter zcnue from gscp_data.f90
      D_c_fix = 2.0E-5_wp                         ! constant mean mass diameter of cloud droplets (if n_cloud_s(:,:)
      x_c_fix = pi * rho_w / 6.0_wp * D_c_fix**3  !   is not present in argument list)
      cloud_tmp%nu = zcnue     ! gamma shape parameter mu for cloud droplets, in Axel's notation this is called nu instead!
      cloud_tmp%mu = 1.0_wp    ! 2nd shape parameter of generalized gamma distribution 
      z_fac_c = moment_gamma(cloud_tmp,2) * mom_fac

      ! Parameters for cloud ice:
      ami = zami              ! mass-size-relation prefactor
      bmi = 3.0_wp            ! mass-size-relation exponent
      IF( lsuper_coolw) THEN
        znimax = znimax_Thom        ! znimax_Thom = 250.E+3_wp 1/m^3
      ELSE
        znimax = fxna(zthn,T_melt)  ! Maximum number of cloud ice crystals in 1/m^3
      END IF

      ! Parameters for rain:
      mue_rain_c = mu_rain
      nor = 8.0e6_wp * EXP(3.2_wp*mue_rain_c) * (0.01_wp)**(-mue_rain_c)
      p_r = (7.0_wp+mue_rain_c) / (4.0_wp+mue_rain_c)
      z_r = nor*gfct(7.0_wp+mue_rain_c) * (pi*rho_w*nor*gfct(4.0_wp+mue_rain_c)/6.0_wp)**(-p_r)
      
      ! Parameters for snow and graupel:
      IF (igscp == 1) THEN
        
        ams = zams_ci
        bms = zbms
        p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
        z_s = mom_fac*ams**2 * gfct(2.0_wp*bms+1.0_wp) * (ams*gfct(bms+1.0_wp))**(-p_s)

        IF (lmessage_light) THEN
          WRITE (*, *) TRIM(routine)//": cloud ice scheme (using rain and snow)"
          WRITE (*,'(A,F10.3)') '     p_r = ', p_r
          WRITE (*,'(A,F10.3)') '     z_r = ', z_r
          WRITE (*,'(A,F10.3)') '     p_s = ', p_s
          WRITE (*,'(A,F10.3)') '     z_s = ', z_s
        ENDIF
      
        firstcall = .FALSE.

      ELSEIF (igscp == 2) THEN

        IF ( .NOT. PRESENT(q_graupel) ) THEN
          message_text(:) = ' '
          WRITE(message_text, '(a,i0,a)') 'optional argument q_graupel needed for inwp_gscp = ',igscp,' but is not present!'
          CALL finish(TRIM(routine), TRIM(message_text))
        END IF

        ams = zams_gr         
        bms = zbms        
        amg = 169.6_wp
        bmg = 3.1_wp
        nog = 4.E6_wp
        p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
        p_g = (2.0_wp*bmg+1.0_wp)/(bmg+1.0_wp)
        z_s = mom_fac*ams**2 * gfct(2.0_wp*bms+1.0_wp) *       &
                         (ams*gfct(bms+1.0_wp))**(-p_s)
        z_g = mom_fac*amg**2 * nog*gfct(2.0_wp*bmg+1.0_wp) *       &
                         (amg*nog*gfct(bmg+1.0_wp))**(-p_g)
        
        IF (lmessage_light) THEN
          WRITE (*, *) TRIM(routine)//": graupel scheme (using rain, snow, graupel)"
          WRITE (*,'(A,F10.3)') '     p_r = ', p_r
          WRITE (*,'(A,F10.3)') '     z_r = ', z_r
          WRITE (*,'(A,F10.3)') '     p_s = ', p_s
          WRITE (*,'(A,F10.3)') '     z_s = ', z_s
          WRITE (*,'(A,F10.3)') '     p_g = ', p_g
          WRITE (*,'(A,F10.3)') '     z_g = ', z_g
        ENDIF
        firstcall = .FALSE.
        
      ELSE
        
        message_text(:) = ' '
        WRITE(message_text, '(a,i0,a)') 'inwp_gscp = ',igscp,' not implemented for DBZ calculation!'
        CALL finish(TRIM(routine), trim(message_text))

      ENDIF
      
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ztc,zn0s,nn,hlp,alf,bet,m2s,m3s, &
!$OMP            rho_c,rho_i,rho_r,rho_s,rho_g,z_cloud,z_rain,z_snow,z_ice,z_grau,x_c,n_i,x_i_mono), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = startblk, endblk

      IF (jb == startblk) THEN
        i_startidx = startidx1
      ELSE
        i_startidx = 1
      ENDIF

      IF (jb == endblk) THEN
        i_endidx = endidx2
      ELSE
        i_endidx = npr
      ENDIF

      DO jk = jk_start, nlev
        DO jc = i_startidx, i_endidx

          ! convert mixing ratios into densities
          rho_c = q_cloud(jc,jk,jb) * rho(jc,jk,jb)
          rho_i = q_ice(jc,jk,jb)   * rho(jc,jk,jb)
          rho_r = q_rain(jc,jk,jb)  * rho(jc,jk,jb)
          rho_s = q_snow(jc,jk,jb)  * rho(jc,jk,jb)

          ! .. cloud droplets (gamma distribution w.r.t. mass x and mu_x=zcnue from gscp_data.f90):
          IF (rho_c >= q_crit_radar) THEN
            IF (lqnc_input) THEN
              x_c = q_cloud(jc,jk,jb) / MAX(n_cloud_s(jc,jb), n_cloud_min)
            ELSE
              x_c = x_c_fix
            END IF
            z_cloud = z_fac_c * rho_c * x_c
            z_radar(jc,jk,jb) = z_cloud * convfac
          ELSE
            z_radar(jc,jk,jb) = 0._wp
          END IF
          
          ! .. rain:
          IF (rho_r >= q_crit_radar) THEN
            z_rain  = z_r * EXP(p_r * LOG(rho_r))
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_rain * convfac
          ENDIF

          ! .. cloud ice (monodisperse size distribution):
          IF (rho_i >= q_crit_radar) THEN
            IF( lsuper_coolw) THEN
              n_i   = MIN( fxna_cooper(T(jc,jk,jb),T_melt), znimax )  ! 1/m^3
            ELSE
              n_i   = MIN( fxna(T(jc,jk,jb),T_melt), znimax )         ! 1/m^3
            END IF
            x_i_mono = rho_i / n_i
            z_ice = MERGE(z_fac_ice_dry, z_fac_ice_wet, T(jc,jk,jb) < T_melt) * mom_fac * rho_i * x_i_mono
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_ice * convfac
          END IF

          ! .. snow:
          IF (rho_s >= q_crit_radar) THEN

            IF (isnow_n0temp == 1) THEN
              ! Calculate n0s using the temperature-dependent
              ! formula of Field et al. (2005)
              ztc = T(jc,jk,jb) - T_melt
              ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
              zn0s = zn0s1*EXP(zn0s2*ztc)
              zn0s = MIN(zn0s,1e9_wp)
              zn0s = MAX(zn0s,1e6_wp)
            ELSEIF (isnow_n0temp == 2) THEN
              ! Calculate n0s using the temperature-dependent moment
              ! relations of Field et al. (2005)
              ztc = T(jc,jk,jb) - T_melt
              ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
              nn  = 3
              hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
                  + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
              alf = EXP(hlp * LOG(10.0_wp))
              bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
                  + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
!!$ UB: caution! Here is the exponent bms=2.0 hardwired! not ideal!
              m2s = rho_s / ams
              m3s = alf*EXP(bet * LOG(m2s))
              hlp  = zn0s1*EXP(zn0s2*ztc)
!!$ UB: the 13.5 is actually 3^3 / gamma(3) ...
              zn0s = 13.50_wp * m2s*(m2s/m3s)**3
              zn0s = MAX(zn0s,0.5_wp*hlp)
              zn0s = MIN(zn0s,1e2_wp*hlp)
              zn0s = MIN(zn0s,1e9_wp)
              zn0s = MAX(zn0s,8e5_wp)
            ELSE
              ! Old constant n0s
              zn0s = 8.0e5_wp
            ENDIF
            hlp = z_s * EXP((1.0_wp-p_s) * LOG(zn0s))

            z_snow = hlp * EXP(p_s*LOG(rho_s)) * MERGE(z_fac_ice_dry, z_fac_ice_wet, T(jc,jk,jb) < T_melt)
            
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_snow * convfac

          END IF

        ENDDO
      ENDDO

      IF (igscp == 2) THEN

        DO jk = jk_start, nlev
          DO jc = i_startidx, i_endidx
            rho_g = q_graupel(jc,jk,jb) * rho(jc,jk,jb)
            IF (rho_g >= q_crit_radar) THEN
              z_grau = z_g * EXP(p_g * LOG(rho_g)) * MERGE(z_fac_ice_dry, z_fac_ice_wet, T(jc,jk,jb) < T_melt)
              z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_grau * convfac
            ENDIF
          ENDDO
        ENDDO

      END IF

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF ( lmessage_light .OR. lmessage_full ) THEN
      zdebug = MAXVAL(z10olog10 * LOG(z_radar + eps))
      message_text(:) = ' '
      WRITE (message_text, '(a,i4,a,1x,f6.1)') 'reflectivity statistics on proc ',my_id_for_message, &
           ': MAX dBZ tot : ', zdebug
      CALL message (TRIM(routine), TRIM(message_text), all_print=.TRUE.)
    ENDIF

  END SUBROUTINE compute_field_dbz_1mom


  !>
  !! Calculate radar reflectivity for the 2-moment microphysics scheme in linear units mm^6/m^3
  !!
  !! @par Revision History
  !! Initial revision by U. Blahak, DWD (2020-01-20) 
  !!
  SUBROUTINE compute_field_dbz_2mom( npr, nlev, nblks, startblk, endblk, jk_start,          &
                                     startidx1, endidx2,                                    &
                                     lmessage_light, lmessage_full, my_id_for_message,      &
                                     rho_w, rho_ice,                                        &
                                     K_w, K_ice, T_melt, q_crit_radar, T, rho,              &
                                     q_cloud, q_rain, q_ice, q_snow, q_graupel, q_hail,     &
                                     n_cloud, n_rain, n_ice, n_snow, n_graupel, n_hail,     &
                                     ql_graupel, ql_hail, z_radar )

    !------------------------------------------------------------------------------
    !
    ! Description:  Calculation of grid point values for effective radar
    !               reflectivity factor Z in dBZ.
    !
    ! Method:       Rayleigh-Approximation for the Back-Scattering
    !               (no attenuation!), Debye-Approximation for the
    !               effective refractive index of two-component ice-air-
    !               mixture particles (dry ice, snow, and graupel).
    !               Melting particles by substitution of the ice substance
    !               with water when T_air > 273.16 K. For the applied
    !               Rayleigh scattering, this is equal to assuming instantaneous
    !               melting, because only the square of the total water mass of the particle
    !               counts ( = the total square of its dipole moment).
    !               
    ! Inputs:       npr, nlev, nblks : field dimensions
    !               startblk, endblk, jk_start, startidx1, endidx2   : loop start and end indices
    !               rho_w        : bulk density of pure water [kg/m**3]
    !               rho_ice      : bulk density of pure ice   [kg/m**3]
    !               K_w          : dielectric constant of water
    !               K_ice        : dielectric constant of ice
    !               T_melt       : melting temperature of ice
    !               q_crit_radar : threshold for the q's to compute reflectivity  [kg/m**3]
    !               T            : temperature field          [K]
    !               rho          : air density                [kg/m**3]
    !               q_cloud      : cloud water mixing ratio   [kg/kg] 
    !               q_rain       : rain water mixing ratio    [kg/kg] 
    !               q_ice        : cloud ice mixing ratio     [kg/kg]
    !               q_snow       : snow mixing ratio          [kg/kg]  
    !               q_graupel    : graupel mixing ratio       [kg/kg]
    !               q_hail       : hail mixing ratio          [kg/kg]
    !               n_cloud      : cloud water number density [1/kg] 
    !               n_rain       : rain water number density  [1/kg] 
    !               n_ice        : cloud ice number density   [1/kg] 
    !               n_snow       : snow number density        [1/kg] 
    !               n_graupel    : graupel number density     [1/kg]
    !               n_hail       : hail number density        [1/kg]
    !   OPTIONAL:   ql_graupel   : liquid water on graupel mixing ratio  [kg/kg]
    !   OPTIONAL:   ql_hail      : liquid water on hail mixing ratio     [kg/kg]
    !
    ! Output:       z_radar      : 3D field of Z              [mm^6/m^3]
    ! 
    !
    !------------------------------------------------------------------------------

    ! Input/Output parameters:
    !-------------------------

    INTEGER,  INTENT(IN) :: npr, nlev, nblks
    INTEGER,  INTENT(IN) :: startblk, endblk, jk_start, startidx1, endidx2, my_id_for_message
    LOGICAL,  INTENT(in) :: lmessage_light, lmessage_full
    REAL(wp), INTENT(in) :: K_w, K_ice, T_melt, rho_w, rho_ice, q_crit_radar

    REAL(wp), INTENT(IN) :: T(:,:,:),          &
                            rho(:,:,:),        &
                            q_cloud(:,:,:),    &
                            q_rain(:,:,:),     &
                            q_ice(:,:,:),      &
                            q_snow(:,:,:),     &
                            q_graupel(:,:,:),  &
                            q_hail(:,:,:),     &
                            n_cloud(:,:,:),    &
                            n_rain(:,:,:),     &
                            n_ice(:,:,:),      &
                            n_snow(:,:,:),     &
                            n_graupel(:,:,:),  &
                            n_hail(:,:,:)

    REAL(wp), INTENT(IN), OPTIONAL ::          &
                            ql_graupel(:,:,:), &
                            ql_hail(:,:,:)
    
    REAL(wp), INTENT(OUT) :: z_radar(:,:,:)

    ! Local Variables
    !----------------

    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_dbz_2mom'

    CHARACTER(len=250) :: message_text

    REAL(wp), PARAMETER :: eps  = 1.00E-15_wp
    REAL(wp), PARAMETER :: eps2 = 1.00E-20_wp
    REAL(wp), PARAMETER :: convfac = 1.E18_wp

    INTEGER       :: jc,jk, jb, pe_center, i_startidx, i_endidx
    LOGICAL, SAVE :: firstcall = .TRUE.
    LOGICAL       :: llwf_scheme

    REAL(wp)      :: T_a,                            &
                     q_c, n_c, x_c,                  &
                     q_r, n_r, x_r,                  &
                     q_g, n_g, x_g,                  &
                     q_h, n_h, x_h,                  &
                     q_s, n_s, x_s,                  &
                     q_i, n_i, x_i,                  &
                     z_fac_ice_dry, z_fac_ice_wet,   &
                     d_r, muD, z_fac_r_muD

    REAL(wp), SAVE :: z_fac_c, z_fac_r, z_fac_i, z_fac_s, z_fac_g, z_fac_h, mom_fac

    TYPE(particle)        :: cloud, rain
    TYPE(particle_frozen) :: ice, snow, graupel, hail
!!$ LWF scheme not yet implemented:    CLASS(particle_lwf)    :: graupel_lwf, hail_lwf

    IF (lmessage_light .OR. lmessage_full) THEN
      message_text(:) = ' '
      WRITE(message_text,'(a,i0)') 'Computing dbz3d_lin for inw_gscp = 4, 5, 6 or 7'
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

    IF (PRESENT(ql_graupel) .AND. PRESENT(ql_hail)) THEN
      llwf_scheme = .TRUE.
    ELSE
      llwf_scheme = .FALSE.
    END IF

    CALL init_2mom_scheme(cloud, rain, ice, snow, graupel, hail)

    ! .. K_ice and K_w  might change from call to call, so not in firstcall only:
    z_fac_ice_dry = (rho_w/rho_ice)**2 * K_ice/K_w
    z_fac_ice_wet = 1.0_wp

    IF (firstcall) THEN
      mom_fac = (6.0_wp / (pi * rho_w))**2
      z_fac_c = moment_gamma(cloud,2)   * mom_fac
      z_fac_r = moment_gamma(rain,2)    * mom_fac
      z_fac_i = moment_gamma(ice,2)     * mom_fac
      z_fac_s = moment_gamma(snow,2)    * mom_fac
      z_fac_g = moment_gamma(graupel,2) * mom_fac
      z_fac_h = moment_gamma(hail,2)    * mom_fac
      IF (lmessage_light) THEN
        WRITE (*, *) TRIM(routine)//":"
        WRITE (*,'(A,D10.3)') "     z_fac_c = ",z_fac_c
        WRITE (*,'(A,D10.3)') "     z_fac_r = ",z_fac_r
        WRITE (*,'(A,D10.3)') "     z_fac_i = ",z_fac_i
        WRITE (*,'(A,D10.3)') "     z_fac_s = ",z_fac_s
        WRITE (*,'(A,D10.3)') "     z_fac_g = ",z_fac_g
        WRITE (*,'(A,D10.3)') "     z_fac_h = ",z_fac_h
        WRITE (*,'(A,D10.3)') "     fak_ice = ",z_fac_ice_dry
      ENDIF
      firstcall = .FALSE.
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,T_a,q_c,n_c,q_r,n_r,q_i,n_i,q_s,n_s,q_g,n_g,q_h,n_h, &
!$OMP&           x_c,x_r,x_i,x_s,x_g,x_h,d_r,muD,z_fac_r_muD), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = startblk, endblk

      IF (jb == startblk) THEN
        i_startidx = startidx1
      ELSE
        i_startidx = 1
      ENDIF

      IF (jb == endblk) THEN
        i_endidx = endidx2
      ELSE
        i_endidx = npr
      ENDIF

      DO jk = jk_start, nlev
        DO jc = i_startidx, i_endidx

          z_radar(jc,jk,jb) = 0.0_wp

          T_a = T(jc,jk,jb)
          q_c = MAX(q_cloud(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
          n_c = MAX(n_cloud(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
          q_r = MAX(q_rain(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          n_r = MAX(n_rain(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          q_i = MAX(q_ice(jc,jk,jb)  , 0.0_wp) * rho(jc,jk,jb)
          n_i = MAX(n_ice(jc,jk,jb)  , 0.0_wp) * rho(jc,jk,jb)
          q_s = MAX(q_snow(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          n_s = MAX(n_snow(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          IF (llwf_scheme) THEN
            q_g = MAX(q_graupel(jc,jk,jb) + ql_graupel(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
            q_h = MAX(q_hail(jc,jk,jb)    + ql_hail(jc,jk,jb)   , 0.0_wp) * rho(jc,jk,jb)
          ELSE
            q_g = MAX(q_graupel(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
            q_h = MAX(q_hail(jc,jk,jb)   , 0.0_wp) * rho(jc,jk,jb)
          END IF
          n_g = MAX(n_graupel(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
          n_h = MAX(n_hail(jc,jk,jb)   , 0.0_wp) * rho(jc,jk,jb)

          x_c = MIN( MAX(q_c/(n_c+eps2),cloud%x_min),cloud%x_max )
          x_r = MIN( MAX(q_r/(n_r+eps2),rain%x_min),rain%x_max )
          x_i = MIN( MAX(q_i/(n_i+eps2),ice%x_min),ice%x_max )
          x_s = MIN( MAX(q_s/(n_s+eps2),snow%x_min),snow%x_max )
          x_g = MIN( MAX(q_g/(n_g+eps2),graupel%x_min),graupel%x_max )
          x_h = MIN( MAX(q_h/(n_h+eps2),hail%x_min),hail%x_max )

          ! .. Cloud water reflectivity:
          IF (q_c < q_crit_radar) x_c = 0.0_wp
          z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_fac_c * q_c * x_c

          ! .. Rain water reflectivity:
          IF (q_r >= q_crit_radar) THEN
            IF (q_c > q_crit_radar) THEN
              ! Inside of cloud cores assume generalized gamma DSD:
              z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_fac_r * q_r * x_r
            ELSE
              ! Outside of cloud cores assume mu-D-relation Seifert (2008):
              d_r = rain%a_geo * EXP(rain%b_geo * LOG(x_r))
              muD = rain_mue_dm_relation(rain_coeffs,d_r)
              z_fac_r_muD = mom_fac*(muD+6.0)*(muD+5.0)*(muD+4.0)/((muD+3.0)*(muD+2.0)*(muD+1.0))
              z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_fac_r_muD * q_r * x_r
            END IF
          END IF

          ! .. Ice species reflectivity:
          IF (q_i < q_crit_radar) x_i = 0.0_wp
          IF (q_s < q_crit_radar) x_s = 0.0_wp
          IF (q_g < q_crit_radar) x_g = 0.0_wp
          IF (q_h < q_crit_radar) x_h = 0.0_wp

          IF (T_a < T_melt) THEN
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_i * q_i * x_i * z_fac_ice_dry
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_s * q_s * x_s * z_fac_ice_dry
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_g * q_g * x_g * z_fac_ice_dry
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_h * q_h * x_h * z_fac_ice_dry
          ELSE
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_i * q_i * x_i * z_fac_ice_wet
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_s * q_s * x_s * z_fac_ice_wet
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_g * q_g * x_g * z_fac_ice_wet
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_h * q_h * x_h * z_fac_ice_wet
          ENDIF

          ! conversion of output unit
          z_radar(jc,jk,jb) = z_radar(jc,jk,jb) * convfac

        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (lmessage_light .OR. lmessage_full ) THEN
      message_text(:) = ' '
      WRITE (message_text, '(A,i4,2(A,F10.1))') 'on proc ',my_id_for_message,': '// &
           'MAX dBZ = ', &
           MAXVAL(10.0_wp / LOG(10.0_wp) * LOG(z_radar + eps)), &
           '  MIN dBZ = ', &
           MINVAL(10.0_wp / LOG(10.0_wp) * LOG(z_radar + eps))
      CALL message(TRIM(routine), TRIM(message_text), all_print=.TRUE.)
    END IF
  
  END SUBROUTINE compute_field_dbz_2mom

  !>
  !! Compute column maximum reflectivity from dbz3d_lin
  !!
  !! @par Revision History
  !! Initial revision by Ulrich Blahak, DWD (2020-01-23) 
  !!
  SUBROUTINE compute_field_dbzcmax( ptr_patch, jg, dbz3d_lin, dbz_cmax )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               ! domain ID of main grid
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(OUT) :: dbz_cmax(:,:)  !< output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

        dbz_cmax(:,jb) = -HUGE(1.0_wp)

        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx

            dbz_cmax(jc,jb) = MAX (dbz_cmax(jc,jb), dbz3d_lin(jc,jk,jb))

          END DO
        END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_dbzcmax

  !>
  !! Compute column maximum radar reflectivity from dbz3d_lin and maximize over time
  !!
  !! @par Revision History
  !! Initial revision by Ulrich Blahak, DWD (2020-01-23) 
  !!
  SUBROUTINE maximize_field_dbzctmax( ptr_patch, jg, dbz3d_lin, dbz_ctmax )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               ! domain ID of main grid
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(INOUT) :: dbz_ctmax(:,:)  !< input/output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jk = kstart_moist(jg), ptr_patch%nlev
          DO jc = i_startidx, i_endidx

            dbz_ctmax(jc,jb) = MAX (dbz_ctmax(jc,jb), dbz3d_lin(jc,jk,jb))

          END DO
        END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE maximize_field_dbzctmax

  !>
  !! Compute radar reflectivity in approx. 850 hPa from dbz3d_lin
  !!
  !! @par Revision History
  !! Initial revision by Ulrich Blahak, DWD (2020-01-23) 
  !!
  SUBROUTINE compute_field_dbz850( ptr_patch, k850, dbz3d_lin, dbz_850 )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: k850(:,:)        !< level index field indicating 850 hPa
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(OUT) :: dbz_850(:,:)  !< output variable, dim: (nproma,nblks_c)

    INTEGER :: i_rlstart,  i_rlend
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jb, jk, jc

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      dbz_850( :,jb) = 0.0_wp

      DO jc = i_startidx, i_endidx

        jk = k850(jc,jb)

        ! Just overtake the values from dbz3d_lin(:,:,:) in linear space. The conversion to dBZ
        ! will be done by a post_op ("post operation") right before output to file. See the
        ! corresponding add_var(..., post_op=(...) ) in mo_nwp_phy_state.f90!
        dbz_850(jc,jb) = dbz3d_lin(jc,jk,jb)

      END DO
      
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_dbz850

  !>
  !! Compute ECHOTOPs in Pa from linear dbz3d_lin
  !!
  !! @par Revision History
  !! Initial revision by Ulrich Blahak, DWD (2020-01-23) 
  !!
  SUBROUTINE compute_field_echotop( ptr_patch, jg, p_diag, dbz3d_lin, echotop_p )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               !< domain ID of main grid
    TYPE(t_nh_diag),      INTENT(IN)  :: p_diag           !< type which contains the dbz3d_lin(:,:,:) field
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(INOUT) :: echotop_p(:,:,:)  !< input/output variable, dim: (nproma,nechotop,nblks_c)

    INTEGER               :: i_rlstart,  i_rlend
    INTEGER               :: i_startblk, i_endblk
    INTEGER               :: i_startidx, i_endidx
    INTEGER               :: jb, jc, jk, lev_etop
    INTEGER               :: jk_echotop(1:nproma)
    REAL(wp)              :: pechotop
    REAL(wp)              :: zthresh, zzthresh, zpA, zpB, zzdbzA,  zzdbzB

    REAL(wp), PARAMETER   :: repsilon = 1.0E8_wp*TINY(1.0_wp)  ! To prevent numerical division by 0 below
    
    ! NOTE: pressure does not have to be recomputed/diagnosed here because this was already done when computing
    !       dbz3d_lin in the call to compute_field_dbz3d_lin() in mo_nh_stepping().
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    DO lev_etop = 1, echotop_meta(jg)%nechotop

      zthresh  = 10.0_wp ** ( 0.1_wp*echotop_meta(jg)%dbzthresh(lev_etop))
      zzthresh = (zthresh+repsilon) ** 0.66666_wp  ! convert to something ~rain rate for interpolation,
                                                   ! take into account repsilon in the same way than in the EXP(LOG()) below

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zzdbzA,zzdbzB,zpA,zpB,pechotop,jk_echotop), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        ! Find the model level just below the echotop:
        jk_echotop(:)  = -999
        DO jk = MAX(kstart_moist(jg),2), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            IF ( jk_echotop(jc) < -900 .AND. dbz3d_lin(jc,jk,jb) >= zthresh ) THEN
              jk_echotop(jc) = jk
            END IF
          END DO
        END DO

        ! Interpolate the exact echotop pressure log-linearily and take the min to the pre-existing "old" value:
        DO jc = i_startidx, i_endidx
          IF (jk_echotop(jc) >= -900) THEN
            jk = jk_echotop(jc)
            ! Interpolation to the pressure where the reflectivity-equivalent normalized rain rate Z^(2/3)
            !  equals the equivalently transformed threshold:
            ! The reflectivity in dbz3d_lin(:,:,:) is linear:
            zzdbzA = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk  ,jb)+repsilon) )
            zzdbzB = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk-1,jb)+repsilon) ) ! Should be < zzdbzA, according to the logic above
            zzdbzB = MIN(zzdbzB - zzdbzA, -repsilon) ! To prevent numerical division by "almost" 0 in the next line
            ! This is the logarithmically interpolated pressure:
            zpA  = LOG( p_diag%pres(jc,jk  ,jb) )
            zpB  = LOG( p_diag%pres(jc,jk-1,jb) )
            pechotop = EXP(zpA + (zpB-zpA) / zzdbzB * (zzthresh-zzdbzA))
            IF ( echotop_p(jc,lev_etop,jb) < 0.0_wp ) THEN
              echotop_p (jc,lev_etop,jb) = pechotop
            ELSE
              echotop_p (jc,lev_etop,jb) = MIN(echotop_p(jc,lev_etop,jb), pechotop)
            END IF
          END IF
        END DO
        
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO

  END SUBROUTINE compute_field_echotop

  !>
  !! Compute ECHOTOPs in m MSL from linear dbz3d_lin
  !!
  !! @par Revision History
  !! Initial revision by Ulrich Blahak, DWD (2020-01-23) 
  !!
  SUBROUTINE compute_field_echotopinm( ptr_patch, jg, p_metrics, dbz3d_lin, echotop_z )

    IMPLICIT NONE

    TYPE(t_patch),        INTENT(IN)  :: ptr_patch        !< patch on which computation is performed
    INTEGER,              INTENT(IN)  :: jg               !< domain ID of grid
    TYPE(t_nh_metrics),   INTENT(IN)  :: p_metrics 
    REAL(wp),             INTENT(IN)  :: dbz3d_lin(:,:,:) !< reflectivity in mm^6/m^3

    REAL(wp),             INTENT(INOUT) :: echotop_z(:,:,:)  !< input/output variable, dim: (nproma,nechotop,nblks_c)

    INTEGER               :: i_rlstart,  i_rlend
    INTEGER               :: i_startblk, i_endblk
    INTEGER               :: i_startidx, i_endidx
    INTEGER               :: jb, jc, jk, lev_etop
    INTEGER               :: jk_echotop(1:nproma)
    REAL(wp)              :: zechotop
    REAL(wp)              :: zthresh, zzthresh, zA, zB, zzdbzA,  zzdbzB

    REAL(wp), PARAMETER   :: repsilon = 1.0E8_wp*TINY(1.0_wp)  ! To prevent numerical division by 0 below

    
    
    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = ptr_patch%cells%start_block( i_rlstart )
    i_endblk   = ptr_patch%cells%end_block  ( i_rlend   )

    DO lev_etop = 1, echotop_meta(jg)%nechotop

      zthresh  = 10.0_wp ** ( 0.1_wp*echotop_meta(jg)%dbzthresh(lev_etop))
      zzthresh = (zthresh+repsilon) ** 0.66666_wp  ! convert to something ~rain rate for interpolation,
                                                   ! take into account repsilon in the same way than in the EXP(LOG()) below

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,zzdbzA,zzdbzB,zA,zB,zechotop,jk_echotop), ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( ptr_patch, jb, i_startblk, i_endblk,     &
                            i_startidx, i_endidx, i_rlstart, i_rlend)

        ! Find the model level just below the echotop:
        jk_echotop(:)  = -999
        DO jk = MAX(kstart_moist(jg),2), ptr_patch%nlev
          DO jc = i_startidx, i_endidx
            IF ( jk_echotop(jc) < -900 .AND. dbz3d_lin(jc,jk,jb) >= zthresh ) THEN
              jk_echotop(jc) = jk
            END IF
          END DO
        END DO

        ! Interpolate the exact echotop height linearily and take the max to the pre-existing "old" value:
        DO jc = i_startidx, i_endidx
          IF (jk_echotop(jc) >= -900) THEN
            jk = jk_echotop(jc)
            ! Interpolation to the height where the reflectivity-equivalent normalized rain rate Z^(2/3)
            !  equals the equivalently transformed threshold:
            ! The reflectivity in dbz3d_lin(:,:,:) is linear:
            zzdbzA = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk  ,jb)+repsilon) )
            zzdbzB = EXP( 0.66666_wp * LOG(dbz3d_lin(jc,jk-1,jb)+repsilon) ) ! Should be < zzdbzA, according to the logic above
            zzdbzB = MIN(zzdbzB - zzdbzA, -repsilon) ! To prevent numerical division by "almost" 0 in the next line
            ! This is the interpolated height:
            zA  = p_metrics%z_mc( jc, jk  , jb)  ! lower bound for linear interpolation
            zB  = p_metrics%z_mc( jc, jk-1, jb)  ! upper bound
            zechotop = zA + (zB-zA) / zzdbzB * (zzthresh-zzdbzA)
            echotop_z (jc,lev_etop,jb) = MAX(echotop_z(jc,lev_etop,jb), zechotop)
          END IF
        END DO
        
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END DO

  END SUBROUTINE compute_field_echotopinm

  !>
  !! Calculate sunshine duration
  !!
  !!  According to WMO (2003),2 sunshine duration during a given period is defined as
  !!  the sum of that sub-period for which the perpendicular direct solar irradiance exceeds 120 W m–2
  !!  WMO-No. 8 Guide to Meteorological Instruments and Methods of Observation
  !!
  !! The direct solar irradiance at the surface is calculated from the shortwave net flux at surface, 
  !! the shortwave upward flux and the shortwave diffuse downward radiative flux. It is divided
  !! by the cosine of solar zenith angle to get the perpendicular solar irradiance.
  !! If the direct solar irradiance exeeds 120 Wm-2 the sunshine duration is extended by the fast physics timestep
  !!
  !! settings for sunshine duration
  !!
  !! dursun_thresh is the threshold for solar direct irradiance in W/m2
  !!       above which the sunshine duration is increased (default 120 W/m2)
  !!
  !! dursun_thresh_width is the smoothness / width of the threshold
  !!       function (e.g. if equal to 60 W/m2 the sunshine duration will
  !!       increase from zero at 170 W/m2 to dt at 230 W/m2)
  !!
  !! Note: MeteoSwiss uses 200 W/m2 instead of WMO value of 120 W/m2 and a width of 60 W/m2
  !!
  !! Note: use dursun_thresh=120.0_wp and dursun_thresh_width=0.01_wp to
  !!       reproduce original behaviour, according to WMO
  !!
  !! @par Revision History
  !! Initial revision by Burkhardt Rockel, Hereon (2021-06-17) 
  !! More sunshine duration fields by Guy de Morsier, MeteoSwiss (2021-11-02)
  !!
  SUBROUTINE compute_field_dursun( pt_patch, dt_phy, dursun,                    &
    &                              swflxsfc, swflx_up_sfc, swflx_dn_sfc_diff,   &
    &                              cosmu0, dursun_thresh, dursun_thresh_width,  &
    &                              dursun_m, dursun_r, pi0, pres, twater)

    TYPE(t_patch),      INTENT(IN)    :: pt_patch              !< patch on which computation is performed
    REAL(wp),           INTENT(IN)    :: dt_phy                !< time interval for fast physics
    REAL(wp), INTENT(INOUT)           :: dursun(:,:)           !< sunshine duration (s)
    REAL(wp),           INTENT(IN)    :: swflxsfc(:,:)         !< shortwave net flux at surface [W/m2]
    REAL(wp),           INTENT(IN)    :: swflx_up_sfc(:,:)     !< shortwave upward flux at the surface [W/m2]
    REAL(wp),           INTENT(IN)    :: swflx_dn_sfc_diff(:,:)!< shortwave diffuse downward radiative flux at the surface [W/m2]
    REAL(wp),           INTENT(IN)    :: cosmu0(:,:)           !< cosine of solar zenith angle
    REAL(wp),           INTENT(IN)    :: dursun_thresh         !< threshold for solar direct irradiance in W/m2
    REAL(wp),           INTENT(IN)    :: dursun_thresh_width   !< smoothness / width of the threshold
    REAL(wp), INTENT(INOUT), OPTIONAL :: dursun_m(:,:)         !< maximum sunshine duration (s)
    REAL(wp), INTENT(INOUT), OPTIONAL :: dursun_r(:,:)         !< relative sunshine duration (s)
    REAL(wp), INTENT(IN), OPTIONAL    :: pi0(:,:)              !< local solar incoming flux at TOA [W/m2]
    REAL(wp), INTENT(IN), OPTIONAL    :: pres(:,:)             !< pressure
    REAL(wp), INTENT(IN), OPTIONAL    :: twater(:,:)           !< total column water

    ! Use a minimum value to avoid div0: The exact value does not make a difference
    ! as the dursun_thresh [W/m2] will not be hit for such small values anyway.
    REAL(wp) :: cosmu0_dark = 1.e-9_wp

    ! local variables
    LOGICAL  :: l_present_dursun_m, l_present_dursun_r
    REAL(wp) :: sun_el, swrad_dir, theta_sun, xval, &
                zsct ! solar constant (at time of year)
    INTEGER  :: i_rlstart,  i_rlend
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jb, jc

    l_present_dursun_m = .FALSE.
    l_present_dursun_r = .FALSE.
    IF (PRESENT(dursun_m)) l_present_dursun_m=.TRUE.
    IF (PRESENT(dursun_r)) l_present_dursun_r=.TRUE.

    ! without halo or boundary  points:
    i_rlstart = grf_bdywidth_c + 1
    i_rlend   = min_rlcell_int

    i_startblk = pt_patch%cells%start_block( i_rlstart )
    i_endblk   = pt_patch%cells%end_block  ( i_rlend   )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,sun_el,swrad_dir,theta_sun,xval,zsct), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( pt_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        IF(cosmu0(jc,jb)>cosmu0_dark) THEN

          ! compute sunshine duration
          xval = (swflxsfc(jc,jb) + swflx_up_sfc(jc,jb) - swflx_dn_sfc_diff(jc,jb))/cosmu0(jc,jb)
          xval = (xval - dursun_thresh)/(dursun_thresh_width/pi)
          IF (xval > 0.5_wp*pi) THEN
            dursun(jc,jb) = dursun(jc,jb) + dt_phy
          ELSEIF (xval > -0.5_wp*pi) THEN
            dursun(jc,jb) = dursun(jc,jb) + dt_phy* 0.5_wp*(SIN(xval) + 1.0_wp)
          ENDIF
        ENDIF

        IF (l_present_dursun_m .AND. l_present_dursun_r) THEN
          ! estimate direct solar radiation for cloud free conditions 
          ! (after R. G. Allen et al. 2006, Agricultural and Forest Meteorology 
          !  doi:10.1016/j.agrformet.2006.05.012                               )
          ! The prefactor in eq. (17) is 0.94 instead of 0.98 in order to better
          ! fit the COSMO clear sky radiation. The turbidity factor $K_t$ is set
          ! to 0.8 (between 0.5 for extremely turbid air and 1.0 for clean air)

          ! sun elevation angle
          sun_el = ASIN(cosmu0(jc,jb))
          ! sun elevation angle in radians
          theta_sun = MAX(0.0_wp, sun_el)
          ! from "calculate solar incoming flux at TOA" in mo_nh_interface_nwp.f90 line 1308
          ! get solar constant 
          zsct = pi0(jc,jb)/cosmu0(jc,jb)
          IF ( swflxsfc(jc,jb) > 0.0001_wp ) THEN
            swrad_dir = zsct * 0.94_wp * EXP(                                    &
                 - 0.00146_wp * pres(jc,jb) / 1.0E3_wp / 0.8_wp / SIN(theta_sun) &
                 - 0.075_wp * (twater(jc,jb)/SIN(theta_sun))**0.4_wp             )
          ELSE
            swrad_dir = 0.0_wp
          ENDIF
  
          ! maximum possible sunshine duration (same formula as for SSD above)
          xval = (swrad_dir-dursun_thresh)/(dursun_thresh_width/pi)
          IF (xval > 0.5_wp*pi) THEN
            dursun_m(jc,jb) = dursun_m(jc,jb) + dt_phy
          ELSEIF (xval > -0.5_wp*pi) THEN
            dursun_m(jc,jb) = dursun_m(jc,jb) + dt_phy* 0.5_wp*(SIN(xval) + 1.0_wp)
          ENDIF
  
          ! relative sunshine duration (%)
          IF (dursun_m(jc,jb) > 0.0_wp) THEN
            dursun_r(jc,jb) = 100.0_wp*dursun(jc,jb)/dursun_m(jc,jb)
          ENDIF
          dursun_r(jc,jb) = MIN(100.0_wp, MAX(0.0_wp, dursun_r(jc,jb)))
        ENDIF

      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_field_dursun

END MODULE mo_opt_nwp_diagnostics

