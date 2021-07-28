!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_rte_rrtmgp_interface
  USE mo_kind,                       ONLY: wp
  USE mo_math_constants,             ONLY: pi
  USE mo_physical_constants,         ONLY: rhoh2o
  USE mo_exception,                  ONLY: finish, warning
  USE mo_radiation_random,           ONLY: seed_size
  USE mo_cloud_optics_parameters,    ONLY: rad_perm
  USE mo_radiation_cloud_optics,     ONLY: cloud_optics
  USE mo_radiation_cld_sampling,     ONLY: sample_cld_state
  USE mo_bc_aeropt_kinne,            ONLY: set_bc_aeropt_kinne
  USE mo_bc_aeropt_stenchikov,       ONLY: add_bc_aeropt_stenchikov
  USE mo_bc_aeropt_splumes,          ONLY: add_bc_aeropt_splumes

  USE mo_optical_props,              ONLY: ty_optical_props_1scl, &
                                           ty_optical_props_2str
  USE mo_gas_concentrations,         ONLY: ty_gas_concs
  USE mo_source_functions,           ONLY: ty_source_func_lw
  USE mo_rte_lw,                     ONLY: rte_lw
  USE mo_rte_sw,                     ONLY: rte_sw
  USE mo_fluxes,                     ONLY: ty_fluxes_broadband
  USE mo_icon_fluxes_sw,             ONLY: ty_icon_fluxes_sw, set_fractions
  USE mo_rte_rrtmgp_setup,           ONLY: k_dist_lw, k_dist_sw, &
                                           cloud_optics_lw, cloud_optics_sw, &
                                           stop_on_err

  USE mo_rad_diag,                   ONLY: rad_aero_diag
  USE mo_timer,                      ONLY: ltimer, timer_start, timer_stop, &
   &                                       timer_lrtm, timer_srtm
  USE mo_radiation_config,           ONLY: lrad_aero_diag
  USE mo_echam_rad_config,           ONLY: echam_rad_config
  USE mtime,                         ONLY: datetime


#ifdef RRTMGP_MERGE_DEBUG
  USE mo_rte_rrtmgp_merge_debug, ONLY: write_record_interface_echam
#endif

! These need to be sent once in the init phae from the atmo to the ps_rad
!          zf(kbdim,klev),               & !< geometric height at full level in m
!          zh(kbdim,klev+1),             & !< geometric height at half level in m
!          dz(kbdim,klev),               & !< geometric height thickness in m

  IMPLICIT NONE

  PRIVATE
#ifdef _OPENACC
  LOGICAL, PARAMETER :: use_acc  = .TRUE.
#else
  LOGICAL, PARAMETER :: use_acc  = .FALSE.
#endif

  LOGICAL, PARAMETER :: top_at_1 = .true.
  LOGICAL            :: lneed_aerosols

  REAL(wp), PARAMETER :: pressure_scale = 100._wp,         &
                         inverse_pressure_scale = 1._wp/pressure_scale, &
                         droplet_scale = 1.0e2, &
                         nir_vis_boundary   = 14500._wp

  PUBLIC :: rte_rrtmgp_interface, pressure_scale, inverse_pressure_scale, &
            droplet_scale

CONTAINS

  !-------------------------------------------------------------------

  !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !!
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Some cloud physical properties are prescribed, which are
  !!   required to derive cloud optical properties
  !!
  !! @par The gases are passed into RRTM via two multi-constituent arrays:


  ! TODO/BUG?
  !-------------------------------------------------------------------
  SUBROUTINE rte_rrtmgp_interface(                                               &
      jg, jb, jcs, jce, nproma, klev,                                       &
      & irad_aero     , ktype                                              ,&
      & psctm, ssi_factor,                                                  &
      & loland          ,loglac          ,this_datetime                    ,&
      & pcos_mu0        ,daylght_frc                                       ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & emissivity                                                         ,&
      & zf              ,zh              ,dz                               ,&
      & pp_sfc          ,pp_fl           ,pp_hl                            ,&
      & tk_sfc          ,tk_fl           ,tk_hl                            ,&
      & xm_dry          ,xvmr_vap        ,xm_liq          ,xm_ice          ,&
      & cdnc            ,xc_frc                                            ,&
      & xvmr_co2        ,xvmr_ch4        ,xvmr_n2o        ,xvmr_cfc        ,&
      & xvmr_o3         ,xvmr_o2                                           ,&
      & lw_upw          ,lw_upw_clr      ,lw_dnw          ,lw_dnw_clr      ,&
      & sw_upw          ,sw_upw_clr      ,sw_dnw          ,sw_dnw_clr      ,&
      & vis_dn_dir_sfc  ,par_dn_dir_sfc  ,nir_dn_dir_sfc                   ,&
      & vis_dn_dff_sfc  ,par_dn_dff_sfc  ,nir_dn_dff_sfc                   ,&
      & vis_up_sfc      ,par_up_sfc      ,nir_up_sfc                       ,&
      & aer_aod_533     ,aer_ssa_533     ,aer_asy_533                      ,&
      & aer_aod_2325    ,aer_ssa_2325    ,aer_asy_2325                     ,&
      & aer_aod_9731                                                       )
#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:1
#endif
     !-------------------------------------------------------------------

    INTEGER,INTENT(IN) :: &
         jg,           & !< domain index
         jb,           & !< block index
         jcs, jce,     & !< starting and ending columns
         nproma, klev, & !< array dimensions(?)
         irad_aero,    & !< aerosol control
         ktype(:)        !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(:)                 !< fraction of TSI in the 14 RRTM SW bands

    LOGICAL,INTENT(IN) ::              &
         loland(:),                & !< land sea mask, land=.true.
         loglac(:)                   !< glacier mask, glacier=.true.

    TYPE(datetime), POINTER ::  this_datetime !< actual time step

    REAL(WP),INTENT(IN)  :: &
         pcos_mu0(:),     & !< mu0 for solar zenith angle
         daylght_frc(:),  & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:),  & !< surface albedo for vis range and dir light
         alb_nir_dir(:),  & !< surface albedo for NIR range and dir light
         alb_vis_dif(:),  & !< surface albedo for vis range and dif light
         alb_nir_dif(:),  & !< surface albedo for NIR range and dif light
         emissivity(:),   & !< sufrace emissivity
         zf(:,:),         & !< geometric height at full level in m
         zh(:,:),         & !< geometric height at half level in m
         dz(:,:),         & !< geometric height thickness in m
         pp_sfc(:),       & !< surface pressure in Pa
         pp_fl(:,:),      & !< full level pressure in Pa
         pp_hl(:,:),      & !< half level pressure in Pa
         tk_sfc(:),       & !< surface temperature in K
         tk_fl(:,:),      & !< full level temperature in K
         tk_hl(:,:),      & !< half level temperature in K
         xm_dry(:,:),     & !< dry air     mass in kg/m2
         xvmr_vap(:,:),   & !< water vapor volume mixing ratio 
         xm_liq(:,:),     & !< cloud water mass in kg/m2
         xm_ice(:,:),     & !< cloud ice   mass in kg/m2
         cdnc(:,:),       & !< cloud nuclei concentration
         xc_frc(:,:),     & !< fractional cloud cover
         xvmr_co2(:,:),   & !< co2 volume mixing ratio
         xvmr_ch4(:,:),   & !< ch4 volume mixing ratio
         xvmr_n2o(:,:),   & !< n2o volume mixing ratio
         xvmr_cfc(:,:,:), & !< cfc volume mixing ratio
         xvmr_o3(:,:),    & !< o3  volume mixing ratio
         xvmr_o2(:,:)       !< o2  volume mixing ratio

    REAL(wp), INTENT(OUT)   :: &
      & lw_dnw_clr(:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:),& !< Clear-sky upward   shortwave at all levels
      & lw_dnw(:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:)       !< All-sky   upward   shortwave at all levels

    REAL(wp), INTENT(OUT) ::  &
         vis_dn_dir_sfc(:),   & !< Diffuse downward flux surface visible radiation
         par_dn_dir_sfc(:),   & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:),   & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:),   & !< Direct  downward flux surface visible radiation
         par_dn_dff_sfc(:),   & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:),   & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:),   & !< Upward  flux surface visible radiation
         par_up_sfc    (:),   & !< Upward  flux surface PAR
         nir_up_sfc    (:),   & !< Upward  flux surface near-infrared radiation
         aer_aod_533   (:,:), & !< Aerosol optical density at 533 nm
         aer_ssa_533   (:,:), & !< Single scattering albedo at 533 nm
         aer_asy_533   (:,:), & !< Asymmetry factor at 533 nm
         aer_aod_2325  (:,:), & !< Aerosol optical density at 2325 nm
         aer_ssa_2325  (:,:), & !< Single scattering albedo at 2325 nm
         aer_asy_2325  (:,:), & !< Asymmetry factor at 2325 nm
         aer_aod_9731  (:,:)    !< Aerosol optical density at 9731 nm

    LOGICAL :: lclearsky

    ! --------------------------------------------------------------------------
    INTEGER :: ncol_supplied, ncol_needed, ncol_chunk, jchunk_start, jchunk_end
    INTEGER :: nbndsw, nbndlw, i, ilev, iband,  jl, jk, jband
    ! --- Aerosol optical properites - vertically reversed fields
    REAL(wp) ::                &
         x_cdnc       (nproma)  !< Scale factor for Cloud Droplet Number Concentration
                                !<  baustelle - x_cndc should be used in cloud optics but isn't
    REAL(wp), ALLOCATABLE :: &
         aer_tau_lw(:,:,:),  & !< LW optical thickness of aerosols
         aer_tau_sw(:,:,:),  & !< aerosol optical thickness
         aer_ssa_sw(:,:,:),  & !< aerosol single scattering albedo
         aer_asy_sw(:,:,:)     !< aerosol asymmetry factor
    ! --------------------------------------------------------------------------
    !
    ! Aerosol optical properties are computed at this level because they require
    !   geographic and temporal information; the geographic information is lost
    !   when the data provided to RTE+RRTMGP are extracted from larger arrarys
    !
    ! IF (aero == ...) THEN
    ! iaero=0: No aerosol
    ! iaero=13: only tropospheric Kinne aerosols
    ! iaero=14: only Stenchikov's volcanic aerosols
    ! iaero=15: tropospheric Kinne aerosols + volcanic Stenchikov's aerosols
    ! set all aerosols to zero first

    lneed_aerosols = (irad_aero /= 0)
    IF(lneed_aerosols) THEN
      nbndlw = k_dist_lw%get_nband()
      nbndsw = k_dist_sw%get_nband()
      
      ALLOCATE( aer_tau_lw(nproma,klev,nbndlw), &
                aer_tau_sw(nproma,klev,nbndsw), &
                aer_ssa_sw(nproma,klev,nbndsw), &
                aer_asy_sw(nproma,klev,nbndsw)  )
         
      !$ACC enter data create(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw) 

      !$ACC kernels default(present) async(1)
      aer_tau_lw(:,:,:) = 0.0_wp
      aer_tau_sw(:,:,:) = 0.0_wp
      aer_ssa_sw(:,:,:) = 1.0_wp
      aer_asy_sw(:,:,:) = 0.0_wp
      !$ACC end kernels

      IF (irad_aero==13 .OR. irad_aero==15 .OR. irad_aero==18) THEN
      ! iaero=13: only Kinne aerosols are used
      ! iaero=15: Kinne aerosols plus Stenchikov's volcanic aerosols are used
      ! iaero=18: Kinne background aerosols (of natural origin, 1850) are set
        write (0,*) "Kinne"
        CALL set_bc_aeropt_kinne(this_datetime,                       &
              & jg,                                                    &
              & jcs, nproma,    nproma,                klev,           &
              & jb,             nbndsw,                nbndlw,         &
              & zf,             dz,                                    &
              & aer_tau_sw,     aer_ssa_sw,            aer_asy_sw,     &
              & aer_tau_lw                                              )
      END IF
      IF (irad_aero==14 .OR. irad_aero==15 .OR. irad_aero==18) THEN
      ! iaero=14: only Stechnikov's volcanic aerosols are used (added to zero)
      ! iaero=15: Stenchikov's volcanic aerosols are added to Kinne aerosols
      ! iaero=18: Stenchikov's volcanic aerosols are added to Kinne background
      !           aerosols (of natural origin, 1850)
        write (0,*) "Stenchikov"
#ifdef _OPENACC
        CALL warning('mo_rte_rrtmgp_interface/rte_rrtmgp_interface','Stenchikov aerosols ACC not implemented')
#endif
        !$acc update host(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw, dz, pp_fl)
        CALL add_bc_aeropt_stenchikov(this_datetime,    jg,               &
              & jcs, nproma,      nproma,                 klev,       &
              & jb,               nbndsw,                nbndlw,           &
              & dz,               pp_fl,                                   &
              & aer_tau_sw,    aer_ssa_sw,         aer_asy_sw,     &
              & aer_tau_lw                                              )
        !$acc update device(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw)
      END IF
      !!$    IF (irad_aero==16) THEN
      !!$      CALL add_aop_volc_ham( &
      !!$           & nproma,           kbdim,                 klev,             &
      !!$           & jb,             nbndlw,                nbndsw,           &
      !!$           & aer_tau_lw,    aer_tau_sw,         aer_ssa_sw,    &
      !!$           & aer_asy_sw                                               )
      !!$    END IF
      !!$    IF (irad_aero==17) THEN
      !!$      CALL add_aop_volc_crow( &
      !!$           & nproma,           kbdim,                 klev,             &
      !!$           & jb,             nbndlw,                nbndsw,           &
      !!$           & aer_tau_lw,    aer_tau_sw,         aer_ssa_sw,    &
      !!$           & aer_asy_sw                                               )
      !!$    END IF
      IF (irad_aero==18) THEN
      ! iaero=18: Simple plumes are added to Stenchikov's volcanic aerosols
      !           and Kinne background aerosols (of natural origin, 1850)
        write (0,*) "Plumes"
#ifdef _OPENACC
        CALL warning('mo_rte_rrtmgp_interface/rte_rrtmgp_interface','Plumes ACC not implemented')
#endif
        !$acc update host(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw, zf, dz, zh(:,klev+1))
        CALL add_bc_aeropt_splumes(jg,                                     &
              & jcs, nproma,           nproma,                 klev,             &
              & jb,             nbndsw,                this_datetime,    &
              & zf,               dz,                    zh(:,klev+1),     &
              & aer_tau_sw,    aer_ssa_sw,         aer_asy_sw,     &
              & x_cdnc                                                     )
        !$acc update device(aer_tau_sw, aer_ssa_sw, aer_asy_sw)
      END IF

      ! this should be decativated in the concurrent version and make the aer_* global variables for output
      IF (lrad_aero_diag) THEN
        CALL rad_aero_diag (                                  &
          & 1,               nproma,          nproma,         &
          & klev,            nbndlw,          nbndsw,         &
          & aer_tau_lw,      aer_tau_sw,      aer_ssa_sw,     &
          & aer_asy_sw,                                       &
          & aer_aod_533,     aer_ssa_533,     aer_asy_533,    &
          & aer_aod_2325,    aer_ssa_2325,    aer_asy_2325,   &
          & aer_aod_9731,                                     &
          & opt_use_acc = use_acc                             )
      ENDIF
      !
      ! Map solar bands from RRTMG to GP order
      !
      CALL rearrange_bands2rrtmgp(nproma,klev,nbndsw, aer_tau_sw)
      CALL rearrange_bands2rrtmgp(nproma,klev,nbndsw, aer_ssa_sw)
      CALL rearrange_bands2rrtmgp(nproma,klev,nbndsw, aer_asy_sw)

      ! Aerosol optical properties have reverse vertical orientation - reorient
      !DA TODO
      CALL reorient_3d_wrt2(aer_tau_lw)
      CALL reorient_3d_wrt2(aer_tau_sw)
      CALL reorient_3d_wrt2(aer_ssa_sw)
      CALL reorient_3d_wrt2(aer_asy_sw)
    ELSE
      ! allocate dummy zero-size arrays
      ALLOCATE( aer_tau_lw(1,1,0), &
                aer_tau_sw(1,1,0), &
                aer_ssa_sw(1,1,0), &
                aer_asy_sw(1,1,0)  )
    END IF
    !
    ! Turn geography and number concentrations into effective radii
    !
    ! --------------------------------------------------------------------------
    ! Set flag for the optional computation of clear-sky fluxes
    lclearsky     = echam_rad_config(jg)%lclearsky
    ! --------------------------------------------------------------------------
    !
    !
    ! Assumption: all arrays share the "column" dimension
    !
    ncol_supplied = size(pcos_mu0) ! baustelle - this should be = nproma?
    ncol_needed   = jce-jcs+1
    ncol_chunk    = echam_rad_config(jg)%rrtmgp_columns_chunk
    !
    ! RTE+RRTMGP process all columns supplied and assume a starting index of 1.
    !   If these conditions are satisfied we can call the interface directly...
    !
    IF (jcs==1 .and. ncol_needed == ncol_supplied .and. ncol_chunk == ncol_needed) THEN

       CALL rte_rrtmgp_interface_onBlock(                              &
          & lclearsky,                                                 &
          & ncol_needed,       klev,                                   &
          & ktype(:),                                                  &
          & psctm,             ssi_factor,                             &
          & loland(:),         loglac(:),                              &
          & pcos_mu0(:),       daylght_frc(:),                         &
          & alb_vis_dir(:),    alb_nir_dir(:),                         &
          & alb_vis_dif(:),    alb_nir_dif(:),                         &
          & emissivity(:),                                             &
          & zf(:,:),           zh(:,:),           dz(:,:),             &
          & pp_sfc(:),         pp_fl(:,:),        pp_hl(:,:),          &
          & tk_sfc(:),         tk_fl(:,:),        tk_hl(:,:),          &
          & xm_dry(:,:),       xvmr_vap(:,:),     xm_liq(:,:),         &
          & xm_ice(:,:),       cdnc(:,:),         xc_frc(:,:),         &
          & xvmr_co2(:,:),     xvmr_ch4(:,:),     xvmr_n2o (:,:),      &
          & xvmr_cfc(:,:,:),   xvmr_o3(:,:),      xvmr_o2(:,:),        &
          & aer_tau_lw(:,:,:),                                         &
          & aer_tau_sw(:,:,:), aer_ssa_sw(:,:,:), aer_asy_sw(:,:,:),   &
          !
          & lw_upw(:,:),       lw_upw_clr (:,:),                       &
          & lw_dnw(:,:),       lw_dnw_clr (:,:),                       &
          & sw_upw(:,:),       sw_upw_clr(:,:),                        &
          & sw_dnw(:,:),       sw_dnw_clr(:,:),                        &
          & vis_dn_dir_sfc(:), par_dn_dir_sfc(:), nir_dn_dir_sfc(:),   &
          & vis_dn_dff_sfc(:), par_dn_dff_sfc(:), nir_dn_dff_sfc(:),   &
          & vis_up_sfc(:),     par_up_sfc(:),     nir_up_sfc(:)        )
       !
    ELSE
       !
       ! ... but if the starting column is > 1 and/or there are trailing columns
       !   which shouldn't be processed we make copies to supply contiguous memory
       !   to RTE+RRMTPG
       !

       DO jchunk_start = jcs,jce, ncol_chunk
        jchunk_end = MIN(jchunk_start + ncol_chunk - 1, jce)
        CALL shift_and_call_rte_rrtmgp_interface_onBlock(                &
            & lclearsky,                                                 &
            & jchunk_start,      jchunk_end,                             &
            & klev,                                                      &
            & ktype(:),                                                  &
            & psctm,             ssi_factor,                             &
            & loland(:),         loglac(:),                              &
            & pcos_mu0(:),       daylght_frc(:),                         &
            & alb_vis_dir(:),    alb_nir_dir(:),                         &
            & alb_vis_dif(:),    alb_nir_dif(:),                         &
            & emissivity(:),                                             &
            & zf(:,:),           zh(:,:),           dz(:,:),             &
            & pp_sfc(:),         pp_fl(:,:),        pp_hl(:,:),          &
            & tk_sfc(:),         tk_fl(:,:),        tk_hl(:,:),          &
            & xm_dry(:,:),       xvmr_vap(:,:),     xm_liq(:,:),         &
            & xm_ice(:,:),       cdnc(:,:),         xc_frc(:,:),         &
            & xvmr_co2(:,:),     xvmr_ch4(:,:),     xvmr_n2o (:,:),      &
            & xvmr_cfc(:,:,:),   xvmr_o3(:,:),      xvmr_o2(:,:),        &
            & aer_tau_lw(:,:,:),                                         &
            & aer_tau_sw(:,:,:), aer_ssa_sw(:,:,:), aer_asy_sw(:,:,:),   &
            !
            & lw_upw(:,:),       lw_upw_clr (:,:),                       &
            & lw_dnw(:,:),       lw_dnw_clr (:,:),                       &
            & sw_upw(:,:),       sw_upw_clr(:,:),                        &
            & sw_dnw(:,:),       sw_dnw_clr(:,:),                        &
            & vis_dn_dir_sfc(:), par_dn_dir_sfc(:), nir_dn_dir_sfc(:),   &
            & vis_dn_dff_sfc(:), par_dn_dff_sfc(:), nir_dn_dff_sfc(:),   &
            & vis_up_sfc(:),     par_up_sfc(:),     nir_up_sfc(:)        )
       END DO
       !
    END IF

  !$ACC wait
  !$ACC exit data delete(aer_tau_lw, aer_tau_sw, aer_ssa_sw, aer_asy_sw) if( lneed_aerosols )

  END SUBROUTINE rte_rrtmgp_interface
 ! -------------------------------------------------------------------------------------
  SUBROUTINE clamp_pressure(src, tgt, low, high)
    REAL(wp), INTENT(IN) :: src(:,:)
    REAL(wp), INTENT(OUT) :: tgt(:,:)
    REAL(wp), INTENT(IN) :: low, high

    !tgt(:,:) = min(high, max(low, src))
    INTEGER :: i, j, m, n

    ! min and max are level-dependent
    REAL(wp), DIMENSION(SIZE(src,2)) :: tgt_min, tgt_max 
    
    !$ACC data create(tgt_min, tgt_max) present(src,tgt)
    
    m = SIZE(src,1)
    n = SIZE(src,2)
    !$ACC parallel loop default(none) gang vector async(1)
    DO j = 1, n
      tgt_min(j)=low  + (j-1)*epsilon(tgt_min)
      tgt_max(j)=high - (n-j)*epsilon(tgt_max)
    END DO

    !$ACC parallel loop default(none) gang vector collapse(2) async(1)
    DO j = 1, n
      DO i = 1 , m
         tgt(i,j) = min(tgt_max(j), max(tgt_min(j), src(i,j)))
      ENDDO
    ENDDO

    !$ACC end data
  END SUBROUTINE clamp_pressure

  SUBROUTINE clamp_temperature(src, tgt, low, high)
    REAL(wp), INTENT(IN) :: src(:,:)
    REAL(wp), INTENT(OUT) :: tgt(:,:)
    REAL(wp), INTENT(IN) :: low, high

    !$ACC kernels default(none) present(tgt,src) async(1)
    tgt(:,:) = min(high, max(low, src(:,:)))
    !$ACC end kernels
  END SUBROUTINE clamp_temperature
  !----------------------------------------------- !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !!
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of thise routine is to reorder the input in the vertical.  In
  !!   addition some cloud physical properties are prescribed, which are
  !!   required to derive cloud optical properties
  !!

  SUBROUTINE rte_rrtmgp_interface_onBlock(                   &
       & lclearsky,                                          &
       & ncol,           klev,                               &
       & ktype,                                              &
       & psctm,          ssi_factor,                         &
       & laland,         laglac,                             &
       & pcos_mu0,       daylght_frc,                        &
       & alb_vis_dir,    alb_nir_dir,                        &
       & alb_vis_dif,    alb_nir_dif,                        &
       & emissivity,                                         &
       & zf,             zh,             dz,                 &
       & pp_sfc,         pp_fl,          pp_hl,              &
       & tk_sfc,         tk_fl,          tk_hl,              &
       & xm_dry,         xvmr_vap,       xm_liq,             &
       & xm_ice,         cdnc,           cld_frc,            &
       & xvmr_co2,       xvmr_ch4,       xvmr_n2o ,          &
       & xvmr_cfc ,      xvmr_o3,        xvmr_o2,            &
       & aer_tau_lw,                                         &
       & aer_tau_sw,     aer_ssa_sw  ,   aer_asy_sw,         &
       & flx_uplw,       flx_uplw_clr,                       &
       & flx_dnlw,       flx_dnlw_clr,                       &
       & flx_upsw,       flx_upsw_clr,                       &
       & flx_dnsw,       flx_dnsw_clr,                       &
       & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc,     &
       & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc,     &
       & vis_up_sfc,     par_up_sfc,     nir_up_sfc          )

#ifdef __INTEL_COMPILER
!DIR$ OPTIMIZE:1
#endif

    LOGICAL,INTENT(IN)  :: lclearsky                     !< flag for clear-sky computations

    INTEGER,INTENT(IN)  :: &
         ncol,             & !< number of columns
         klev,             & !< number of levels
         ktype(:)            !< type of convection

    REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp),INTENT(IN) :: ssi_factor(:)                 !< fraction of TSI in the 14 RRTM SW bands

    LOGICAL,INTENT(IN) :: &
         laland(:),   & !< land sea mask, land=.true.
         laglac(:)      !< glacier mask, glacier=.true.

    REAL(WP),INTENT(IN)  ::    &
         pcos_mu0(:),      & !< mu0 for solar zenith angle
         daylght_frc(:),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:),   & !< surface albedo for vis range and dir light
         alb_nir_dir(:),   & !< surface albedo for NIR range and dir light
         alb_vis_dif(:),   & !< surface albedo for vis range and dif light
         alb_nir_dif(:),   & !< surface albedo for NIR range and dif light
         emissivity(:),    & !< surface longwave emissivity
         zf(:,:),          & !< geometric height at full level in m
         zh(:,:),          & !< geometric height at half level in m
         dz(:,:),          & !< geometric height thickness in m
         
         pp_sfc(:),        & !< surface pressure in Pa
         pp_fl(:,:),       & !< full level pressure in Pa
         pp_hl(:,:),       & !< full level pressure in Pa
         tk_sfc(:),        & !< surface temperature in K
         tk_fl(:,:),       & !< full level temperature in K
         tk_hl(:,:),       & !< half level temperature in K
         xm_dry(:,:),      & !< dry air     mass in kg/m2
         xvmr_vap(:,:),    & !< water vapor volume mixing ratio
         xm_liq(:,:),      & !< cloud water mass in kg/m2
         xm_ice(:,:),      & !< cloud ice   mass in kg/m2
         aer_tau_lw(:,:,:),& !< aerosol optical depth, longwave (ncol, nlay, nbndlw)
         aer_tau_sw(:,:,:),& !< aerosol optical depth,            shortwave (ncol, nlay, nbndlw)
         aer_ssa_sw(:,:,:),& !< aerosol single-scattering albedo, shortwave (ncol, nlay, nbndlw)
         aer_asy_sw(:,:,:),& !< aerosol asymetry parameter,       shortwave (ncol, nlay, nbndlw)
         cdnc(:,:),        & !< cloud nuclei concentration
         cld_frc(:,:),     & !< fractional cloud cover
         xvmr_co2(:,:),    & !< co2 volume mixing ratio
         xvmr_ch4(:,:),    & !< ch4 volume mixing ratio
         xvmr_n2o(:,:),    & !< n2o volume mixing ratio
         xvmr_cfc(:,:,:),  & !< cfc volume mixing ratio (kbdim,klev,2)
         xvmr_o3(:,:),     & !< o3  volume mixing ratio
         xvmr_o2(:,:)        !< o2  volume mixing ratio

    REAL (wp), TARGET, INTENT (OUT) ::       &
         flx_uplw    (:,:), & !<   upward LW flux profile, all sky
         flx_uplw_clr(:,:), & !<   upward LW flux profile, clear sky
         flx_dnlw    (:,:), & !< downward LW flux profile, all sky
         flx_dnlw_clr(:,:), & !< downward LW flux profile, clear sky
         flx_upsw    (:,:), & !<   upward SW flux profile, all sky
         flx_upsw_clr(:,:), & !<   upward SW flux profile, clear sky
         flx_dnsw    (:,:), & !< downward SW flux profile, all sky
         flx_dnsw_clr(:,:)    !< downward SW flux profile, clear sky

    REAL (wp), TARGET, INTENT (OUT) :: &
         vis_dn_dir_sfc(:) , & !< Diffuse downward flux surface visible radiation
         par_dn_dir_sfc(:) , & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc(:) , & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc(:) , & !< Direct  downward flux surface visible radiation
         par_dn_dff_sfc(:) , & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc(:) , & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc    (:) , & !< Upward  flux surface visible radiation
         par_up_sfc    (:) , & !< Upward  flux surface PAR
         nir_up_sfc    (:)     !< Upward  flux surface near-infrared radiation

    ! -----------------------------------------------------------------------
    ! At this stage, columns are stored in the arrays in a contiguous form, i.e.
    ! all loops run over (1:ncol, 1:klev) where ncol=jce-jcs+1 and arrays are copied
    ! into the arrays before accordingly.

    INTEGER  :: jk, jl !< loop indices
    INTEGER  :: nbndlw, nbndsw, ngptsw
    REAL(wp) ::                      &
         zsemiss(k_dist_lw%get_nband(),ncol) !< LW surface emissivity by band
    ! --- local scaled variables
    REAL(wp) ::                   &
         cld_frc_loc, & !< secure cloud fraction
         ziwp       (ncol,klev), & !< in cloud ice water content       [g/m2]
         zlwp       (ncol,klev), & !< in cloud liquid water content    [g/m2]
         ziwc, zlwc                !< in cloud water concentration     [g/m3]

    REAL(wp) ::                &
         re_drop (ncol,klev), & !< effective radius of liquid
         re_cryst(ncol,klev)
    !
    ! Random seeds for sampling. Needs to get somewhere upstream
    !
    INTEGER :: rnseeds(ncol,seed_size), gpt_lims(2), gpt_start, gpt_end
    INTEGER :: band, gpt, i, j
    REAL(wp) :: low, high

    TYPE(ty_source_func_lw)     :: source_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: atmos_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: aerosol_lw !check types regarding acc later
    TYPE(ty_optical_props_1scl) :: clouds_lw, clouds_bnd_lw !check types regarding acc later
    TYPE(ty_optical_props_2str) :: atmos_sw !check types regarding acc later
    TYPE(ty_optical_props_2str) :: aerosol_sw !check types regarding acc later
    TYPE(ty_optical_props_2str) :: clouds_sw, clouds_bnd_sw !check types regarding acc later

    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_mask
    TYPE(ty_fluxes_broadband) :: fluxes_lw !check acc
    TYPE(ty_icon_fluxes_sw  ) :: fluxes_sw !check acc
    TYPE(ty_fluxes_broadband) :: fluxes_lwcs, fluxes_swcs !check acc

    TYPE(ty_gas_concs) :: gas_concs !check acc
    REAL(wp), DIMENSION(ncol       ) :: tsi_norm_factor, mu0
    REAL(wp), DIMENSION(ncol,klev  ) :: tlay, play, dummy2d
    REAL(wp), DIMENSION(ncol,klev+1) :: tlev, plev
    REAL(wp), DIMENSION(k_dist_sw%get_nband(),ncol) :: albdir, albdif
    REAL(wp) :: band_lims(2,k_dist_sw%get_nband()), delwave, frc_vis
    REAL(wp) :: toa_flux(ncol,k_dist_sw%get_ngpt())
    ! Threshold within which a cloud fraction is considered = 0 or 1.
    REAL(wp) :: cld_frc_thresh
    LOGICAL  :: do_frac_cloudiness
    !--------------------------------
    INTEGER, PARAMETER :: n_gas_names = 8
    CHARACTER(len=5), PARAMETER :: gas_names(n_gas_names) = (/ &
       'h2o  ', 'co2  ', 'ch4  ', 'o2   ', 'o3   ', 'n2o  ','cfc11', 'cfc12'/)
    !--------------------------------
    ! Variables for effective radii computations
    REAL (wp), PARAMETER :: &
       ccwmin = 1.e-7_wp, &    ! min condensate for lw cloud opacity
       zkap_cont = 1.143_wp, & ! continental (Martin et al. ) breadth param
       zkap_mrtm = 1.077_wp    ! maritime (Martin et al.) breadth parameter
    REAL (wp) :: effective_radius
    REAL (wp) :: reimin, reimax, relmin, relmax, re_cryst_scal, re_drop_scal, zkap
    LOGICAL   :: lcldlyr
    !
    !DA TODO: rearrange the data section to reduce memory consumption
    !
    !$ACC data present(cld_frc, xm_ice, xm_liq, dz, pcos_mu0, emissivity,            &
    !$ACC              alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,           &
    !$ACC              daylght_frc, laland, laglac, dz, cdnc )                       &
    !$ACC      create (ziwp, zlwp, mu0, zsemiss, albdif, re_cryst, re_drop,          &
    !$ACC              albdir, rnseeds, tsi_norm_factor, toa_flux,                   &
    !$ACC              plev, play, tlev, tlay)

    nbndlw = k_dist_lw%get_nband()
    nbndsw = k_dist_sw%get_nband()
    ngptsw = k_dist_sw%get_ngpt()
    cld_frc_thresh = 4._wp*spacing(1._wp)

    ! 1.0 Constituent properties
    !--------------------------------
    !
    ! Is there fractional cloudiness i.e. differences from 0 or 1? If not we can skip McICA sampling
    !
    do_frac_cloudiness = .false.
    !$ACC parallel loop default(none) copy(do_frac_cloudiness) &
    !$ACC               gang vector collapse(2) async(1)
    DO jk = 1, klev
      DO jl = 1, ncol
         IF (min(abs(cld_frc(jl,jk) - 1._wp), abs(cld_frc(jl,jk))) > cld_frc_thresh) THEN
            do_frac_cloudiness = .true.
#ifndef _OPENACC
            EXIT
#endif
         END IF
      END DO
    END DO
    !
    ! Keep the fractional cloudiness turned off for now
    do_frac_cloudiness = .false. 
    !
    ! 1.1 Cloud condensate and cloud fraction
    !
    effective_radius = &
      1.0e6_wp * droplet_scale * (3.0e-9_wp / (4.0_wp * pi * rhoh2o))**(1.0_wp/3.0_wp) 

    reimin = MAX(cloud_optics_lw%get_min_radius_ice(), cloud_optics_sw%get_min_radius_ice()) ! 10.0_wp  !
    reimax = MIN(cloud_optics_lw%get_max_radius_ice(), cloud_optics_sw%get_max_radius_ice()) ! 124.0_wp !
    
    relmin = MAX(cloud_optics_lw%get_min_radius_liq(), cloud_optics_sw%get_min_radius_liq()) ! 2.5_wp  ! 
    relmax = MIN(cloud_optics_lw%get_max_radius_liq(), cloud_optics_sw%get_max_radius_liq()) ! 21.5_wp ! 

    IF (relmax <= relmin .OR. reimax <= reimin) THEN
      CALL finish('rte_rrtmgp_interface_onBlock (mo_rte_rrtmgp_interface.f90)', &
                  'Droplet minimun size required is bigger than maximum')
    END IF

    !$ACC parallel loop default(none) gang vector collapse(2) async(1)
    DO jk = 1, klev
      DO jl = 1, ncol
        !
        ! --- Cloud liquid and ice mass: [kg/m2 in cell] --> [g/m2 in cloud]
        !
        cld_frc_loc = MAX(EPSILON(1.0_wp),cld_frc(jl,jk))
        ziwp(jl,jk) = xm_ice(jl,jk)*1000.0_wp/cld_frc_loc
        zlwp(jl,jk) = xm_liq(jl,jk)*1000.0_wp/cld_frc_loc

        ! Mask which tells cloud optics that this cell is clear
        lcldlyr = cld_frc(jl,jk) > cld_frc_thresh !!!
        IF (.NOT. lcldlyr) THEN
          ziwp(jl,jk) = 0.0_wp
          zlwp(jl,jk) = 0.0_wp
        END IF
        !
        ! --- cloud water and ice concentrations [g/m3]
        !
        ziwc = ziwp(jl,jk)/dz(jl,jk) !!!
        zlwc = zlwp(jl,jk)/dz(jl,jk)
        !
        IF (lcldlyr .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN

          zkap = zkap_mrtm
          IF ( laland(jl) .AND. .NOT.laglac(jl) ) zkap = zkap_cont
          
          re_cryst_scal = MAX(reimin, MIN(reimax,83.8_wp*ziwc**0.216_wp))
          re_drop_scal  = MAX(relmin, MIN(relmax, &
            effective_radius * zkap * (zlwc / cdnc(jl,jk))**(1.0_wp/3.0_wp) ))

          re_cryst(jl,jk) = re_cryst_scal
          re_drop (jl,jk) = re_drop_scal
        ELSE
          re_cryst(jl,jk) = reimin
          re_drop (jl,jk) = relmin
        END IF
      END DO
    END DO
    !
    ! RRTMGP cloud optics: here compute effective radius of liquid and ice from formulae in
    !   mo_radiation_cloud_optics.
    !
    !--------------------------------
    !
    ! 1.2 Gas concentrations
    !

    ! At this stage, columns are stored in the arrays in a contiguous form, i.e.
    ! all loops run over (1:ncol, 1:klev) where ncol=jce-jcs+1 and arrays are copied
    ! into the arrays before accordingly.
    ! The gas profile routine provides all gas concentrations in volume mixing ratios
    !
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC wait
    CALL stop_on_err(gas_concs%init(gas_names))
    CALL stop_on_err(gas_concs%set_vmr('h2o',   xvmr_vap))
    CALL stop_on_err(gas_concs%set_vmr('co2',   xvmr_co2))
    CALL stop_on_err(gas_concs%set_vmr('ch4',   xvmr_ch4))
    CALL stop_on_err(gas_concs%set_vmr('o2',    xvmr_o2))
    CALL stop_on_err(gas_concs%set_vmr('o3',    xvmr_o3))
    CALL stop_on_err(gas_concs%set_vmr('n2o',   xvmr_n2o))
    CALL stop_on_err(gas_concs%set_vmr('cfc11', xvmr_cfc(:,:,1)))
    CALL stop_on_err(gas_concs%set_vmr('cfc12', xvmr_cfc(:,:,2)))

    !--------------------------------
    !
    ! Restrict out-of-bounds temperatures and pressures
    !
    ! The air pressure on levels plev, on the upper and lower boundaries of a layer,
    ! is used to determine the air mass transferred by radiation. For safety
    ! reasons plev is limited to values >= 0 Pa (and <=10**6 Pa so that high is defined).
    low = 0._wp
    high = 1000000._wp
    CALL clamp_pressure(pp_hl, plev, low, high)
    !
    ! The pressure in the layer play is used for optical properties and
    ! is limited here to the range for which tables are defined in the
    ! RRTMGP data file and stored in k_dist_lw.
    ! Thus radiation can be computed for pressure values lower or higher
    ! than the defined range, albeit at lower precision.

    low =  k_dist_lw%get_press_min()
    high = k_dist_lw%get_press_max()
    CALL clamp_pressure(pp_fl, play, low, high)

    low =  k_dist_lw%get_temp_min()
    high = k_dist_lw%get_temp_max()
    CALL clamp_temperature(tk_hl, tlev, low, high)
    CALL clamp_temperature(tk_fl, tlay, low, high)

    !--------------------------------
    !
    ! Boundary conditions
    !
    ! baustelle - shouldn't the min solar zenith cosine be parameterized?
!!debug++
    !$ACC kernels default(none) async(1)
    mu0(:) = MAX(1.e-10_wp,MIN(1.0_wp,pcos_mu0(:)))
    !$ACC end kernels
!!debug--


    ! 2.0 Surface Properties
    ! --------------------------------
   !$ACC parallel loop gang vector default(none) collapse(2) async(1)
   DO j=1,ncol
      DO i=1,nbndlw
        zsemiss(i,j)=emissivity(j)
      END DO
    END DO

    !
    ! Surface albedo interpolation
    !
    !DA TODO: the next function has to run on GPU
    band_lims = k_dist_sw%get_band_lims_wavenumber()

    !$ACC parallel default(none) copyin(band_lims) async(1)
    !$ACC loop collapse(2)
    DO j=1,ncol
      DO band=1,nbndsw
        delwave = band_lims(2,band) - band_lims(1,band)
        frc_vis = MAX(0.0_wp, MIN(1.0_wp, &
          (band_lims(2,band) - nir_vis_boundary) / delwave))

        albdif(band,j) = alb_vis_dif(j) * frc_vis + &
                         alb_nir_dif(j) * (1.0_wp - frc_vis)
        albdir(band,j) = alb_vis_dir(j) * frc_vis + &
                         alb_nir_dir(j) * (1.0_wp - frc_vis)
      END DO
    END DO
    !$ACC end parallel

    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
    !
    ! 3.1 Aerosols - moved to a layer above
    !

    !
!!    ! 3.2 Clouds optical properties
!!    call stop_on_err(clouds_bnd_lw%alloc_1scl(ncol, klev, &
!!                     k_dist_lw%get_band_lims_wavenumber()))
!!    call stop_on_err(clouds_bnd_sw%alloc_2str(ncol, klev, &
!!                     k_dist_sw%get_band_lims_wavenumber()))
!!    !$ACC enter data create (clouds_bnd_lw, clouds_bnd_sw,   clouds_bnd_lw%tau,    &
!!    !$ACC                    clouds_bnd_sw%tau, clouds_bnd_sw%ssa, clouds_bnd_sw%g)
!!
!!    !$ACC update host(ktype, cld_frc, laglac, laland, icldlyr, zlwp, ziwp, zlwc, ziwc, cdnc)
!!    !DA TODO: this -> to GPU 
!!    CALL cloud_optics(                                                     &
!!         & laglac,         laland,         ncol,         ncol,             &
!!         & klev,           ktype,                                          &
!!         & icldlyr,    zlwp,       ziwp,       zlwc,                       &
!!         & ziwc,       cdnc,   clouds_bnd_lw%tau, clouds_bnd_sw%tau,       &
!!         & clouds_bnd_sw%ssa, clouds_bnd_sw%g,  re_drop,        re_cryst        )

!!    ! write(0,*) "Clouds checksum: ", sum(ktype), sum(cld_frc), sum(merge(1,0,laglac)), sum(merge(1,0,laland)),&
!!    ! sum(icldlyr), sum(zlwp), sum(ziwp), &
!!    ! sum(ziwc), sum(cdnc), sum(clouds_bnd_lw%tau), sum(clouds_bnd_sw%tau), sum(clouds_bnd_sw%ssa), sum(clouds_bnd_sw%g)

!!    !$ACC update device(clouds_bnd_lw%tau, clouds_bnd_sw%tau, clouds_bnd_sw%ssa, &
!!    !$ACC               clouds_bnd_sw%g)
!!    ! baustelle : the above statement has to go once cloud_optics on GPU 
!!    ! remember to set a "end data" further down...
!!    !
!!    ! RRTMG band order differs from GP
!!    !
!!    CALL rearrange_bands2rrtmgp(ncol,klev,nbndsw, clouds_bnd_sw%tau)
!!    CALL rearrange_bands2rrtmgp(ncol,klev,nbndsw, clouds_bnd_sw%ssa)
!!    CALL rearrange_bands2rrtmgp(ncol,klev,nbndsw, clouds_bnd_sw%g)
!!    !
!!    ! new cloud optics: delete section 3.2, since LW and SW cloud optics can be computed
!!    !    separately within each section, reducing memory footprint
!!    !

    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    IF (ltimer) CALL timer_start(timer_lrtm)
    !
    ! 4.1 Longwave radiative Transfer

    !
    ! 4.1.2 Gas optics
    !
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC wait
    CALL stop_on_err(source_lw%alloc    (ncol, klev, k_dist_lw))
    CALL stop_on_err(atmos_lw%alloc_1scl(ncol, klev, k_dist_lw))
    ! baustelle - surface and lowest layer temperature are considered the same
    !$ACC data create(source_lw, source_lw%lay_source,     source_lw%lev_source_inc, &
    !$ACC                        source_lw%lev_source_dec, source_lw%sfc_source,     &
    !$ACC                        source_lw%sfc_source_Jac,                           &
    !$ACC             atmos_lw,  atmos_lw%tau)

    CALL stop_on_err( &
           k_dist_lw%gas_optics(play, plev, tlay, tk_sfc, &
                                gas_concs, atmos_lw, source_lw, &
                                tlev = tlev))
    !
    ! 4.1.2 Aerosol optical depth: add to clear-sky
    !  If irad_aero == 0, aer_tau_lw will not be allocated here
    !  and we need to skip this step
    !
    IF ( lneed_aerosols ) THEN
      CALL stop_on_err(aerosol_lw%alloc_1scl(ncol, klev, &
                                             k_dist_lw%get_band_lims_wavenumber()))
      !$ACC data present(aer_tau_lw) create(aerosol_lw, aerosol_lw%tau)
      !
      !DA TODO: this can be just a pointer assignment
      !$ACC parallel loop default(none) gang vector collapse(3) async(1)
      DO band = 1, nbndlw
        DO j = 1, klev
          DO i = 1, ncol
            aerosol_lw%tau(i,j,band) = aer_tau_lw(i,j,band)
          END DO
        END DO
      END DO
      !
      ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
      !$ACC wait
      CALL stop_on_err(aerosol_lw%increment(atmos_lw))
      ! aerosols
      !$ACC end data
      DEALLOCATE(aerosol_lw%tau)
      CALL aerosol_lw%finalize()
    END IF
    !
    !
    ! 4.1.3 Longwave clear-sky fluxes
    !
    IF (lclearsky) THEN
       !
       fluxes_lwcs%flux_up => flx_uplw_clr
       fluxes_lwcs%flux_dn => flx_dnlw_clr
       !
       ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
       !$ACC wait
       CALL stop_on_err(rte_lw(atmos_lw, top_at_1, source_lw, zsemiss, fluxes_lwcs))
       !
    END IF

    ! new cloud optics: allocate memory for cloud optical properties:
    CALL stop_on_err(clouds_bnd_lw%alloc_1scl(ncol, klev, &
                     k_dist_lw%get_band_lims_wavenumber()))
    !$ACC data create (clouds_bnd_lw, clouds_bnd_lw%tau)    
    ! then compute cloud optics

    ! !$ACC update host(zlwp,     ziwp,    re_drop,    re_cryst)
    ! write (0,*) "newcloudsss", sum(zlwp),     sum(ziwp),    sum(re_drop),    sum(re_cryst)

    CALL stop_on_err(cloud_optics_lw%cloud_optics( &
                     zlwp,     ziwp,    re_drop,    re_cryst,   clouds_bnd_lw ))
    ! This will require computing logical masks for ice and liquid clouds
    !   nrghice (ice roughness) is 1, 2, or 3; probably any values is fine
    !
    ! 4.1.4 McICA sampling of cloud optics for LW
    !

    IF (do_frac_cloudiness) THEN
      !
      ! Seeds for random numbers come from least significant digits of
      ! pressure field
      !
      !DA TODO: remove plev_vr alltogether
      !$ACC parallel loop default(none) gang vector collapse(2) async(1)
      DO jk=1,seed_size
        DO jl=1,ncol
          rnseeds(jl,jk) = &
            int((inverse_pressure_scale * plev(jl,klev+2-jk) -  &
             int(inverse_pressure_scale * plev(jl,klev+2-jk)))* 1E9 + rad_perm)
        END DO
      END DO

      ALLOCATE(cloud_mask(ncol,klev,k_dist_lw%get_ngpt()))
      CALL stop_on_err(clouds_lw%alloc_1scl(ncol, klev, k_dist_lw))
      !
      ! Do not change the order of these data sections!!
      !
      !$ACC data create(clouds_lw, clouds_lw%tau)
      !$ACC data create(cloud_mask)

      CALL sample_cld_state(ncol, ncol, klev,     &
                            k_dist_lw%get_ngpt(), &
                            rnseeds, cld_frc, cloud_mask)

      !DA todo: restructure loop such that not so many gpu kernels are created (so in one construct)
      DO band = 1, nbndlw
        gpt_lims(:) = clouds_lw%convert_band2gpt(band)
        gpt_start = gpt_lims(1) ! avoid copying array gpt_lims
        gpt_end   = gpt_lims(2)

        !$ACC parallel loop default(none) present(clouds_bnd_lw) &
        !$ACC               gang vector collapse(3) async(1)
        DO gpt = gpt_start, gpt_end
          DO j = 1, klev
            DO i = 1, ncol
              IF (cloud_mask(i,j,gpt)) THEN
                clouds_lw%tau(i,j,gpt) = clouds_bnd_lw%tau(i,j,band)
              ELSE
                clouds_lw%tau(i,j,gpt) = 0.0_wp
              END IF
            ENDDO
          ENDDO
        ENDDO

      ENDDO
      !$ACC end data
      DEALLOCATE(cloud_mask)

      ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
      !$ACC wait
      CALL stop_on_err(clouds_lw%increment(atmos_lw))

      !$ACC end data
      DEALLOCATE(clouds_lw%tau)
      CALL clouds_lw%finalize()
    ELSE
      CALL stop_on_err(clouds_bnd_lw%increment(atmos_lw))
    END IF

    !$ACC end data
    DEALLOCATE(clouds_bnd_lw%tau)
    CALL clouds_bnd_lw%finalize()

    !
    ! 4.1.5 Longwave all-sky fluxes
    !
    fluxes_lw%flux_up => flx_uplw
    fluxes_lw%flux_dn => flx_dnlw
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC wait
    CALL stop_on_err(rte_lw(atmos_lw, top_at_1, source_lw, &
                            zsemiss, fluxes_lw))
    !
    ! 4.1.6 End of longwave calculations - free memory
    !
    !$ACC end data
    DEALLOCATE(atmos_lw%tau)
    CALL source_lw%finalize()
    CALL atmos_lw%finalize()

    IF (ltimer) CALL timer_stop(timer_lrtm)
    !
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !
    ! 4.2 Shortwave calculations
    !
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !
    IF (ltimer) CALL timer_start(timer_srtm)
    !
    ! 4.2.1 Array and type allocation for shortwave
    !--------------------------------
    !
    ! Shortwave gas optical properties and source functions
    !
    CALL stop_on_err(atmos_sw%alloc_2str(ncol, klev, k_dist_sw))
    !$ACC data create(atmos_sw, atmos_sw%tau, atmos_sw%ssa, atmos_sw%g, &
    !$ACC             toa_flux)

    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC wait
    CALL stop_on_err(&
       k_dist_sw%gas_optics(play, plev, tlay, &
                            gas_concs, atmos_sw, &
                            toa_flux))
    !
    ! Normalize incident radiation
    !
!    tsi_norm_factor(:) = 1._wp/SUM(toa_flux, dim=2)
    !$ACC kernels default(none) async(1)
    tsi_norm_factor(:) = 0._wp
    !$ACC end kernels

    !$ACC parallel default(none) async(1)
    !$ACC loop seq
    DO gpt = 1, ngptsw
      !$ACC loop gang vector
      DO i = 1, ncol
        tsi_norm_factor(i) = tsi_norm_factor(i) + toa_flux(i,gpt)
      END DO
    END DO
    !$ACC end parallel

    !$ACC parallel loop default(none) gang vector collapse(2) async(1)
    DO gpt = 1, ngptsw 
      DO i = 1, ncol
        toa_flux(i,gpt) = toa_flux(i,gpt) * daylght_frc(i) * psctm / tsi_norm_factor(i)
        !
        ! baustelle - should be more general
        !
        toa_flux(i,gpt) = toa_flux(i,gpt) * (1360.9_wp/1368.22_wp)
      ENDDO
    ENDDO
    !
    ! 4.2.2 Aerosol optical depth: add to clear-sky, reorder bands
    !
    IF ( lneed_aerosols ) THEN
      CALL stop_on_err(aerosol_sw%alloc_2str(ncol, klev, &
                                            k_dist_sw%get_band_lims_wavenumber()))
      !$ACC data create (aerosol_sw, aerosol_sw%tau, aerosol_sw%ssa, aerosol_sw%g) &
      !$ACC      present(aer_tau_sw, aer_ssa_sw, aer_asy_sw)
      !DA TODO: this could be just a pointer assignment
      !$ACC parallel loop default(none) gang vector collapse(3) async(1)
      DO band = 1, nbndsw
        DO j = 1, klev
          DO i = 1, ncol
            aerosol_sw%tau(i,j,band) = aer_tau_sw(i,j,band)
            aerosol_sw%ssa(i,j,band) = aer_ssa_sw(i,j,band)
            aerosol_sw%g  (i,j,band) = aer_asy_sw(i,j,band)
          END DO
        END DO
      END DO
      ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
      !$ACC wait
      CALL stop_on_err(aerosol_sw%increment(atmos_sw))
      ! aerosol_sw
      !$ACC end data
      CALL aerosol_sw%finalize()
    END IF
    !
    ! 4.2.3 Shortwave clear-sky fluxes
    !
    IF (lclearsky) THEN
       !
       fluxes_swcs%flux_up => flx_upsw_clr
       fluxes_swcs%flux_dn => flx_dnsw_clr
       !
       ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
       !$ACC wait
       CALL stop_on_err(rte_sw(atmos_sw, top_at_1, mu0, toa_flux, albdir, albdif, fluxes_swcs))
       !
    END IF
    
    ! new cloud optics: allocate memory for cloud optical properties:
    CALL stop_on_err(clouds_bnd_sw%alloc_2str(ncol, klev, &
                     k_dist_sw%get_band_lims_wavenumber()))
    !$ACC data create (clouds_bnd_sw, clouds_bnd_sw%tau, clouds_bnd_sw%ssa, clouds_bnd_sw%g)    
    ! then compute cloud optics
    CALL stop_on_err(cloud_optics_sw%cloud_optics( &
                     zlwp,     ziwp,    re_drop,    re_cryst,   clouds_bnd_sw ))
    ! 4.2.4 McICA sampling of cloud optics for size(source, dim=1)W
    !    Reset random seeds so SW doesn't depend on what's happened in LW but is also independent
    !
    IF (do_frac_cloudiness) THEN

      !$ACC parallel loop default(none) gang vector collapse(2) async(1)
      DO jk=1,seed_size
        DO jl=1,ncol
          rnseeds(jl,jk) = &
            int((inverse_pressure_scale * plev(jl,klev+1-seed_size+jk) -  &
             int(inverse_pressure_scale * plev(jl,klev+1-seed_size+jk)))* 1E9 + rad_perm)
        END DO
      END DO

      ALLOCATE(cloud_mask(ncol,klev,k_dist_sw%get_ngpt()))
      CALL stop_on_err(clouds_sw%alloc_2str(ncol, klev, k_dist_sw))
      !
      ! Do not change the order of these data sections!!
      !
      !$ACC data create(clouds_sw, clouds_sw%tau, clouds_sw%ssa, clouds_sw%g)
      !$ACC data create(cloud_mask)

      CALL sample_cld_state(ncol, ncol, klev, &
                            k_dist_sw%get_ngpt(), &
                            rnseeds, cld_frc, cloud_mask)

      !## !$ACC update wait host(cloud_mask, rnseeds, cld_frc)
      !## write(0,*) "cld_state: ", sum(merge(1,0, cloud_mask)), sum(rnseeds), sum(cld_frc)

      DO band = 1, nbndsw
        gpt_lims(:) = clouds_sw%convert_band2gpt(band)
        gpt_start = gpt_lims(1) ! avoid copying array gpt_lims
        gpt_end   = gpt_lims(2)

        !$ACC parallel loop default(none) present(clouds_bnd_sw) &
        !$ACC               gang vector collapse(3) async(1)
        DO gpt = gpt_start, gpt_end
          DO j = 1, klev
            DO i = 1, ncol
              IF (cloud_mask(i,j,gpt)) THEN
                clouds_sw%g  (i,j,gpt) = clouds_bnd_sw%g  (i,j,band)
                clouds_sw%ssa(i,j,gpt) = clouds_bnd_sw%ssa(i,j,band)
                clouds_sw%tau(i,j,gpt) = clouds_bnd_sw%tau(i,j,band)
              ELSE
                clouds_sw%g  (i,j,gpt) = 0.0_wp
                clouds_sw%ssa(i,j,gpt) = 1.0_wp
                clouds_sw%tau(i,j,gpt) = 0.0_wp
              ENDIF
            ENDDO
          ENDDO
        ENDDO

      ENDDO
      !$ACC end data
      DEALLOCATE(cloud_mask)
      
      ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
      !$ACC wait
!!$      Delta scaling should be activated here when fractional clouds are used but is 
!!$      not yet tested in this case
!!$      CALL stop_on_err(clouds_bnd_sw%delta_scale()) ! necessary for cases w=g near 1
      CALL finish ('rte_rrtmgp_interface_onBlock, mo_rte_rrtmgp_interface.f90', &
     &             'Delta scaling for fractional clouds must be activated and tested, '// &
     &             'see comment lines in code')
      CALL stop_on_err(clouds_sw%increment(atmos_sw))

      ! clouds_sw
      !$ACC end data
      CALL clouds_sw%finalize()
    ELSE
      CALL stop_on_err(clouds_bnd_sw%delta_scale()) ! necessary for cases w=g near 1
      CALL stop_on_err(clouds_bnd_sw%increment(atmos_sw))
    END IF

    !$ACC end data
    CALL clouds_bnd_sw%finalize()
    !
    ! 4.2.5 Shortwave all-sky fluxes
    !
    fluxes_sw%flux_up => flx_upsw
    fluxes_sw%flux_dn => flx_dnsw
    fluxes_sw%vis_dn_dir_sfc => vis_dn_dir_sfc
    fluxes_sw%par_dn_dir_sfc => par_dn_dir_sfc
    fluxes_sw%nir_dn_dir_sfc => nir_dn_dir_sfc
    fluxes_sw%vis_dn_dff_sfc => vis_dn_dff_sfc
    fluxes_sw%par_dn_dff_sfc => par_dn_dff_sfc
    fluxes_sw%nir_dn_dff_sfc => nir_dn_dff_sfc
    fluxes_sw%vis_up_sfc => vis_up_sfc
    fluxes_sw%par_up_sfc => par_up_sfc
    fluxes_sw%nir_up_sfc => nir_up_sfc

    CALL set_fractions(fluxes_sw, atmos_sw, psctm, ssi_factor)
    ! RTE-RRTMGP ACC code is synchronous, so need to wait before calling it
    !$ACC wait
    CALL stop_on_err(rte_sw(atmos_sw, top_at_1, &
                            mu0, toa_flux, albdir, albdif, &
                            fluxes_sw))

    !
    ! 4.2.6 End of shortwave calculations - free memory
    !
    !$ACC end data
    CALL atmos_sw%finalize()
        
    IF (ltimer) CALL timer_stop(timer_srtm)

#ifdef RRTMGP_MERGE_DEBUG
!$OMP CRITICAL (write_record)
    CALL write_record_interface_echam(nproma, pcos_mu0, daylght_frc, &
      rnseeds1, rnseeds2, alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
      tk_sfc, zf, zh, dz, pp_fl, pp_hl, tk_fl, tk_hl, &
      play, plev, tlay, tlev, &
      xm_dry, xvmr_vap, xvmr_co2, xvmr_ch4, xvmr_o2, xvmr_o3, xvmr_n2o, cdnc, &
      cld_frc, &
      flx_dnlw_clr, flx_uplw_clr, flx_dnsw_clr, flx_upsw_clr, &
      flx_dnlw, flx_uplw, flx_dnsw, flx_upsw, &
      vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
      vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
      vis_up_sfc,     par_up_sfc,     nir_up_sfc      )
!$OMP END CRITICAL (write_record)
#endif

  !$ACC end data
  END SUBROUTINE rte_rrtmgp_interface_onBlock
  ! ----------------------------------------------------------------------------
  SUBROUTINE shift_and_call_rte_rrtmgp_interface_onBlock(    &
    & lclearsky,                                      &
    & jcs,            jce,                            &
    &                 klev,                           &
    !
    & ktype,                                          &
    & psctm,          ssi_factor,                     &
    & laland,         laglac,                         &
    & pcos_mu0,       daylght_frc,                    &
    & alb_vis_dir,    alb_nir_dir,                    &
    & alb_vis_dif,    alb_nir_dif,                    &
    & emissivity,                                     &
    & zf,             zh,             dz,             &
    & pp_sfc,         pp_fl,          pp_hl,          &
    & tk_sfc,         tk_fl,          tk_hl,          &
    & xm_dry,         xvmr_vap,       xm_liq,         &
    & xm_ice,         cdnc,           xc_frc,         &
    & xvmr_co2,       xvmr_ch4,       xvmr_n2o,       &
    & xvmr_cfc,       xvmr_o3,        xvmr_o2,        &
    & aer_tau_lw,                                     &
    & aer_tau_sw,     aer_ssa_sw,     aer_asy_sw,     &
    !
    & lw_upw,         lw_upw_clr,                     &
    & lw_dnw,         lw_dnw_clr,                     &
    & sw_upw,         sw_upw_clr,                     &
    & sw_dnw,         sw_dnw_clr,                     &
    & vis_dn_dir_sfc, par_dn_dir_sfc, nir_dn_dir_sfc, &
    & vis_dn_dff_sfc, par_dn_dff_sfc, nir_dn_dff_sfc, &
    & vis_up_sfc,     par_up_sfc,     nir_up_sfc      )

 LOGICAL,INTENT(IN)  :: lclearsky                     !< flag for clear-sky computations

 INTEGER,INTENT(IN)  :: &
      & jcs,            & !< cell/column index, start
      & jce,            & !< cell/column index, end
      & klev,           & !< number of full levels
      & ktype(:)            !< type of convection

 REAL(wp),INTENT(IN) :: psctm                         !< orbit and time dependent solar constant for radiation time step
 REAL(wp),INTENT(IN) :: ssi_factor(:)                 !< fraction of TSI in the 14 RRTM SW bands

 LOGICAL,INTENT(IN) :: &
      & laland(:),   & !< land sea mask, land=.true.
      & laglac(:)      !< glacier mask, glacier=.true.

 REAL(WP),INTENT(IN)  ::    &
      & pcos_mu0(:),      & !< mu0 for solar zenith angle
      & daylght_frc(:),   & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
      & alb_vis_dir(:),   & !< surface albedo for vis range and dir light
      & alb_nir_dir(:),   & !< surface albedo for NIR range and dir light
      & alb_vis_dif(:),   & !< surface albedo for vis range and dif light
      & alb_nir_dif(:),   & !< surface albedo for NIR range and dif light
      & emissivity(:),    & !< surface longwave emissivity
      & zf(:,:),          & !< geometric height at full level in m
      & zh(:,:),          & !< geometric height at half level in m
      & dz(:,:),          & !< geometric height thickness in m
      & pp_sfc(:),        & !< surface pressure in Pa
      & pp_fl(:,:),       & !< full level pressure in Pa
      & pp_hl(:,:),       & !< full level pressure in Pa
      & tk_sfc(:),        & !< surface temperature in K
      & tk_fl(:,:),       & !< full level temperature in K
      & tk_hl(:,:),       & !< half level temperature in K
      & xm_dry(:,:),      & !< dry air     mass in kg/m2
      & xvmr_vap(:,:),    & !< water vapor volume mixing ratio
      & xm_liq(:,:),      & !< cloud water mass in kg/m2
      & xm_ice(:,:),      & !< cloud ice   mass in kg/m2
      & cdnc(:,:),        & !< cloud nuclei concentration
      & xc_frc(:,:),      & !< fractional cloud cover
      & xvmr_co2(:,:),    & !< co2 volume mixing ratio
      & xvmr_ch4(:,:),    & !< ch4 volume mixing ratio
      & xvmr_n2o(:,:),    & !< n2o volume mixing ratio
      & xvmr_cfc(:,:,:),  & !< cfc volume mixing ratio (kbdim,klev,2)
      & xvmr_o3(:,:),     & !< o3  volume mixing ratio
      & xvmr_o2(:,:),     & !< o2  volume mixing ratio
      & aer_tau_lw(:,:,:),& !< aerosol optical depth, longwave (ncol, nlay, nbndlw)
      & aer_tau_sw(:,:,:),& !< aerosol optical depth,            shortwave (ncol, nlay, nbndlw)
      & aer_ssa_sw(:,:,:),& !< aerosol single-scattering albedo, shortwave (ncol, nlay, nbndlw)
      & aer_asy_sw(:,:,:)   !< aerosol asymetry parameter,       shortwave (ncol, nlay, nbndlw)

 REAL (wp), TARGET, INTENT (OUT) ::       &
      & lw_upw    (:,:), & !<   upward LW flux profile, all sky
      & lw_upw_clr(:,:), & !<   upward LW flux profile, clear sky
      & lw_dnw    (:,:), & !< downward LW flux profile, all sky
      & lw_dnw_clr(:,:), & !< downward LW flux profile, clear sky
      & sw_upw    (:,:), & !<   upward SW flux profile, all sky
      & sw_upw_clr(:,:), & !<   upward SW flux profile, clear sky
      & sw_dnw    (:,:), & !< downward SW flux profile, all sky
      & sw_dnw_clr(:,:)    !< downward SW flux profile, clear sky

 REAL (wp), TARGET, INTENT (OUT) :: &
      & vis_dn_dir_sfc(:) , & !< Diffuse downward flux surface visible radiation
      & par_dn_dir_sfc(:) , & !< Diffuse downward flux surface PAR
      & nir_dn_dir_sfc(:) , & !< Diffuse downward flux surface near-infrared radiation
      & vis_dn_dff_sfc(:) , & !< Direct  downward flux surface visible radiation
      & par_dn_dff_sfc(:) , & !< Direct  downward flux surface PAR
      & nir_dn_dff_sfc(:) , & !< Direct  downward flux surface near-infrared radiation
      & vis_up_sfc    (:) , & !< Upward  flux surface visible radiation
      & par_up_sfc    (:) , & !< Upward  flux surface PAR
      & nir_up_sfc    (:)     !< Upward  flux surface near-infrared radiation

 INTEGER :: ncol !< number of columns needed

 ! Shifted input arguments
 !
 REAL(wp)  ::                                    &
      & s_zf             (jce-jcs+1,klev),       & !< geometric height at full level in m
      & s_zh             (jce-jcs+1,klev+1),     & !< geometric height at half level in m
      & s_dz             (jce-jcs+1,klev),       & !< geometric height thickness in m
      & s_pp_fl          (jce-jcs+1,klev),       & !< full level pressure in Pa
      & s_pp_hl          (jce-jcs+1,klev+1),     & !< full level pressure in Pa
      & s_tk_fl          (jce-jcs+1,klev),       & !< full level temperature in K
      & s_tk_hl          (jce-jcs+1,klev+1),     & !< half level temperature in K
      & s_xm_dry         (jce-jcs+1,klev),       & !< dry air     mass in kg/m2
      & s_xvmr_vap       (jce-jcs+1,klev),       & !< water vapor volume mixing ratio
      & s_xm_liq         (jce-jcs+1,klev),       & !< cloud water mass in kg/m2
      & s_xm_ice         (jce-jcs+1,klev),       & !< cloud ice   mass in kg/m2
      & s_cdnc           (jce-jcs+1,klev),       & !< cloud nuclei concentration
      & s_xc_frc         (jce-jcs+1,klev),       & !< fractional cloud cover
      & s_xvmr_co2       (jce-jcs+1,klev),       & !< co2 volume mixing ratio
      & s_xvmr_ch4       (jce-jcs+1,klev),       & !< ch4 volume mixing ratio
      & s_xvmr_n2o       (jce-jcs+1,klev),       & !< n2o volume mixing ratio
      & s_xvmr_cfc       (jce-jcs+1,klev,2),     & !< cfc volume mixing ratio
      & s_xvmr_o3        (jce-jcs+1,klev),       & !< o3  volume mixing ratio
      & s_xvmr_o2        (jce-jcs+1,klev)          !< o2  volume mixing ratio

 REAL(wp), ALLOCATABLE ::    & 
      & s_aer_tau_lw(:,:,:), &
      & s_aer_tau_sw(:,:,:), &
      & s_aer_ssa_sw(:,:,:), &
      & s_aer_asy_sw(:,:,:)
 ! Shifted output arguments
 !
 REAL(wp)  ::                              &
      & s_lw_upw       (jce-jcs+1,klev+1), & !<   upward LW flux profile, all sky
      & s_lw_upw_clr   (jce-jcs+1,klev+1), & !<   upward LW flux profile, clear sky
      & s_lw_dnw       (jce-jcs+1,klev+1), & !< downward LW flux profile, all sky
      & s_lw_dnw_clr   (jce-jcs+1,klev+1), & !< downward LW flux profile, clear sky
      & s_sw_upw       (jce-jcs+1,klev+1), & !<   upward SW flux profile, all sky
      & s_sw_upw_clr   (jce-jcs+1,klev+1), & !<   upward SW flux profile, clear sky
      & s_sw_dnw       (jce-jcs+1,klev+1), & !< downward SW flux profile, all sky
      & s_sw_dnw_clr   (jce-jcs+1,klev+1)    !< downward SW flux profile, clear sky

  ! Shift input arguments that would be non-contiguous when sliced
  !
  ncol = jce-jcs+1

  !$ACC data create( &
  !$ACC              s_zf,             s_zh,             s_dz,         &
  !$ACC              s_pp_fl,          s_pp_hl,                        &
  !$ACC              s_tk_fl,          s_tk_hl,                        &
  !$ACC              s_xm_dry,         s_xvmr_vap,       s_xm_liq,     &
  !$ACC              s_xm_ice,         s_cdnc,           s_xc_frc,     &
  !$ACC              s_xvmr_co2,       s_xvmr_ch4,       s_xvmr_n2o,   &
  !$ACC              s_xvmr_cfc,       s_xvmr_o3,        s_xvmr_o2,    &
  !$ACC              s_lw_upw,         s_lw_upw_clr,                   &
  !$ACC              s_lw_dnw,         s_lw_dnw_clr,                   &
  !$ACC              s_sw_upw,         s_sw_upw_clr,                   &
  !$ACC              s_sw_dnw,         s_sw_dnw_clr                    )

  ! (ncol, klev)
  !$ACC kernels default(present) async(1)
  s_zf          (1:ncol,:)   = zf          (jcs:jce,:)
  s_dz          (1:ncol,:)   = dz          (jcs:jce,:)
  s_pp_fl       (1:ncol,:)   = pp_fl       (jcs:jce,:)
  s_tk_fl       (1:ncol,:)   = tk_fl       (jcs:jce,:)
  s_xm_dry      (1:ncol,:)   = xm_dry      (jcs:jce,:)
  s_xvmr_vap    (1:ncol,:)   = xvmr_vap    (jcs:jce,:)
  s_xm_liq      (1:ncol,:)   = xm_liq      (jcs:jce,:)
  s_xm_ice      (1:ncol,:)   = xm_ice      (jcs:jce,:)
  s_cdnc        (1:ncol,:)   = cdnc        (jcs:jce,:)
  s_xc_frc      (1:ncol,:)   = xc_frc      (jcs:jce,:)
  s_xvmr_co2    (1:ncol,:)   = xvmr_co2    (jcs:jce,:)
  s_xvmr_ch4    (1:ncol,:)   = xvmr_ch4    (jcs:jce,:)
  s_xvmr_n2o    (1:ncol,:)   = xvmr_n2o    (jcs:jce,:)
  s_xvmr_o3     (1:ncol,:)   = xvmr_o3     (jcs:jce,:)
  s_xvmr_o2     (1:ncol,:)   = xvmr_o2     (jcs:jce,:)
  !$ACC end kernels

  ! (ncol, klev+1)
  !$ACC kernels default(present) async(1)
  s_zh          (1:ncol,:)   = zh          (jcs:jce,:)
  s_pp_hl       (1:ncol,:)   = pp_hl       (jcs:jce,:)
  s_tk_hl       (1:ncol,:)   = tk_hl       (jcs:jce,:)
  !$ACC end kernels

  ! (ncol, klev, 2)
  !$ACC kernels default(present) async(1)
  s_xvmr_cfc    (1:ncol,:,:) = xvmr_cfc    (jcs:jce,:,:)
  !$ACC end kernels

  IF ( lneed_aerosols ) THEN
    ! Aerosols are present, irad_aero /= 0
    !      
    ALLOCATE( s_aer_tau_lw(jce-jcs+1,klev,k_dist_lw%get_nband()), &
              s_aer_tau_sw(jce-jcs+1,klev,k_dist_sw%get_nband()), &
              s_aer_ssa_sw(jce-jcs+1,klev,k_dist_sw%get_nband()), &
              s_aer_asy_sw(jce-jcs+1,klev,k_dist_sw%get_nband())  )
    !
    !$ACC enter data create(s_aer_tau_lw, s_aer_tau_sw, s_aer_ssa_sw, s_aer_asy_sw)
    !
    ! (ncol, klev, nbndlw)
    !$ACC kernels default(present) async(1)
    s_aer_tau_lw  (1:ncol,:,:) = aer_tau_lw(jcs:jce,:,:)
    !$ACC end kernels

    ! (ncol, klev, nbndsw)
    !$ACC kernels default(present) async(1)
    s_aer_tau_sw  (1:ncol,:,:) = aer_tau_sw(jcs:jce,:,:)
    s_aer_ssa_sw  (1:ncol,:,:) = aer_ssa_sw(jcs:jce,:,:)
    s_aer_asy_sw  (1:ncol,:,:) = aer_asy_sw(jcs:jce,:,:)
    !$ACC end kernels
  ELSE
    ! allocate dummy zero-size arrays
    ALLOCATE( s_aer_tau_lw(1,1,0), &
              s_aer_tau_sw(1,1,0), &
              s_aer_ssa_sw(1,1,0), &
              s_aer_asy_sw(1,1,0)  )
  END IF

  ! Call radiation with shifted input arguments and receive shifted output arguments
  !
  CALL rte_rrtmgp_interface_onBlock(                                                 &
      & lclearsky,                                                                   &
      &   ncol,                klev,                                                 &
      !
      &   ktype    (jcs:jce),                                                        &
      &   psctm,                  ssi_factor,                                        &
      & laland     (jcs:jce),     laglac     (jcs:jce),                              &
      & pcos_mu0   (jcs:jce),     daylght_frc(jcs:jce),                              &
      & alb_vis_dir(jcs:jce),     alb_nir_dir(jcs:jce),                              &
      & alb_vis_dif(jcs:jce),     alb_nir_dif(jcs:jce),                              &
      & emissivity (jcs:jce),                                                        &
      & s_zf(:,:),                s_zh(:,:),                s_dz(:,:),               &
      & pp_sfc     (jcs:jce),     s_pp_fl(:,:),             s_pp_hl(:,:),            &
      & tk_sfc     (jcs:jce),     s_tk_fl(:,:),             s_tk_hl(:,:),            &
      & s_xm_dry(:,:),            s_xvmr_vap(:,:),          s_xm_liq(:,:),           &
      & s_xm_ice(:,:),            s_cdnc(:,:),              s_xc_frc(:,:),           &
      & s_xvmr_co2(:,:),          s_xvmr_ch4(:,:),          s_xvmr_n2o(:,:),         &
      & s_xvmr_cfc(:,:,:),        s_xvmr_o3(:,:),           s_xvmr_o2(:,:),          &
      & s_aer_tau_lw(:,:,:),                                                         &
      & s_aer_tau_sw(:,:,:),      s_aer_ssa_sw(:,:,:),      s_aer_asy_sw(:,:,:),     &
      !     
      & s_lw_upw(:,:),            s_lw_upw_clr(:,:),                                 &
      & s_lw_dnw(:,:),            s_lw_dnw_clr(:,:),                                 &
      & s_sw_upw(:,:),            s_sw_upw_clr(:,:),                                 &
      & s_sw_dnw(:,:),            s_sw_dnw_clr(:,:),                                 &
      & vis_dn_dir_sfc(jcs:jce),  par_dn_dir_sfc(jcs:jce),  nir_dn_dir_sfc(jcs:jce), &
      & vis_dn_dff_sfc(jcs:jce),  par_dn_dff_sfc(jcs:jce),  nir_dn_dff_sfc(jcs:jce), &
      & vis_up_sfc    (jcs:jce),  par_up_sfc    (jcs:jce),  nir_up_sfc    (jcs:jce)  )

  ! Shift output arguments
  !
  ! (ncol, klev+1)
  !$ACC kernels default(present) async(1)
  lw_upw         (jcs:jce,:) = s_lw_upw         (1:ncol,:)
  lw_upw_clr     (jcs:jce,:) = s_lw_upw_clr     (1:ncol,:)
  lw_dnw         (jcs:jce,:) = s_lw_dnw         (1:ncol,:)
  lw_dnw_clr     (jcs:jce,:) = s_lw_dnw_clr     (1:ncol,:)
  sw_upw         (jcs:jce,:) = s_sw_upw         (1:ncol,:)
  sw_upw_clr     (jcs:jce,:) = s_sw_upw_clr     (1:ncol,:)
  sw_dnw         (jcs:jce,:) = s_sw_dnw         (1:ncol,:)
  sw_dnw_clr     (jcs:jce,:) = s_sw_dnw_clr     (1:ncol,:)
  !$ACC end kernels

  !$ACC exit data delete(s_aer_tau_lw, s_aer_tau_sw, s_aer_ssa_sw, s_aer_asy_sw) if (lneed_aerosols)
  !$ACC end data
END SUBROUTINE shift_and_call_rte_rrtmgp_interface_onBlock

SUBROUTINE reorient_3d_wrt2 (field)
  REAL(wp), INTENT(INOUT)  :: field(:,:,:)
  INTEGER                  :: nl1, nl2, nl3, il1, il2, il3, idx1, idx2
  REAL(WP)                 :: tmp

  nl1=SIZE(field,1)
  nl2=SIZE(field,2)
  nl3=SIZE(field,3)

  !$ACC parallel loop present(field) gang vector collapse(3) async(1)
  DO il3 = 1, nl3
    DO il2 = 1, nl2/2
      DO il1 = 1, nl1
        idx1 = il2
        idx2 = nl2 - il2 + 1
        tmp = field(il1,idx1,il3)
        field(il1,idx1,il3) = field(il1,idx2,il3)
        field(il1,idx2,il3) = tmp
      END DO
    END DO
  END DO
END SUBROUTINE reorient_3d_wrt2

SUBROUTINE rearrange_bands2rrtmgp(nproma, klev, nbnd, field)
  INTEGER,  INTENT(IN)    :: nproma, klev, nbnd
  REAL(wp), INTENT(INOUT) :: field(nproma,klev,nbnd) 

  REAL(wp) :: last
  INTEGER  :: jk, jl, jband, i

#ifndef _OPENACC
  field(:,:,:) = field(:,:,[nbnd, (i, i = 1, nbnd-1)])
#else
  !$ACC parallel default(present) async(1)
  !$ACC loop gang vector collapse(2)
  DO jk=1,klev
    DO jl=1,nproma
      last = field(jl,jk,nbnd)

      !$ACC loop seq
      DO jband=nbnd, 2, -1
        field(jl,jk,jband) = field(jl,jk,jband-1)
      END DO

      field(jl,jk,1) = last
    END DO
  END DO
  !$ACC end parallel
#endif
END SUBROUTINE rearrange_bands2rrtmgp

END MODULE mo_rte_rrtmgp_interface
