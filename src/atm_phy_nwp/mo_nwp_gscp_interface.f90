!OPTION! -cont
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! microphysical parameterisations:
!!
!! inwp_gscp == 1 : one_moment bulk microphysics by Doms and Schaettler(2004) 
!!                                                and Seifert and Beheng(2006)
!! inwp_gscp == 2 : one-moment graupel scheme
!!
!! inwp_gscp == 3 : two-moment cloud ice scheme of Koehler (2013)
!!
!! inwp_gscp == 4 : two-moment bulk microphysics by Seifert and Beheng (2006)
!!                  with prognostic cloud droplet number
!!
!! inwp_gscp == 5 : two-moment bulk microphysics by Seifert and Beheng (2006)
!!                  with prognostic cloud droplet number and some aerosol,
!!                  CCN and IN tracers
!!
!! inwp_gscp == 6 : two-moment bulk microphysics by Seifert and Beheng (2006)
!!                  incorporating prognostic aerosol as CCN and IN from the
!!                  ART extension
!!
!! inwp_gscp == 9 : a simple Kessler-type warm rain scheme
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
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
#include "icon_contiguous_defines.inc"
!----------------------------

MODULE mo_nwp_gscp_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message
  USE mo_parallel_config,      ONLY: nproma

  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, iss, iorg, iso4, idu, iedmf
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, iqr, iqs,       &
                                     iqni, iqni_nuc, iqg, iqh, iqnr, iqns,     &
                                     iqng, iqnh, iqnc, inccn, ininpot, ininact,&
                                     iqtvar, iqgl, iqhl
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, iprog_aero
  USE gscp_kessler,            ONLY: kessler
  USE gscp_cloudice,           ONLY: cloudice
  USE gscp_graupel,            ONLY: graupel
  USE gscp_hydci_pp_ice,       ONLY: hydci_pp_ice
  USE mo_exception,            ONLY: finish
  USE mo_2mom_mcrph_driver,    ONLY: two_moment_mcrph, set_qnc, &
       &                             set_qnr,set_qni,set_qns,set_qng,set_qnh
#ifdef __ICON_ART
  USE mo_art_clouds_interface, ONLY: art_clouds_interface_2mom
#endif
  USE mo_nwp_diagnosis,        ONLY: nwp_diag_output_minmax_micro
  USE gscp_data,               ONLY: cloud_num
  USE mo_cpl_aerosol_microphys,ONLY: specccn_segalkhain, ncn_from_tau_aerosol_speccnconst, &
                                     specccn_segalkhain_simple
  USE mo_grid_config,          ONLY: l_limited_area
  USE mo_satad,                ONLY: satad_v_3D, satad_v_3D_gpu

  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop,    &
      &                              timer_phys_micro_specific,                &
      &                              timer_phys_micro_satad                               

  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_microphysics

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_microphysics(  tcall_gscp_jg,                & !>input
                            &   lsatad,                       & !>input
                            &   p_patch,p_metrics,            & !>input
                            &   p_prog,                       & !>inout
                            &   ptr_tracer,                   & !>inout
                            &   ptr_tke,                      & !>in
                            &   p_diag ,                      & !>inout
                            &   prm_diag,prm_nwp_tend,        & !>inout
                            &   ext_data,                     & !>in
                            &   lcompute_tt_lheat             ) !>in 



    TYPE(t_patch)          , INTENT(in)   :: p_patch        !!<grid/patch info.
    TYPE(t_nh_metrics)     , INTENT(in)   :: p_metrics
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog          !<the dyn prog vars
    REAL(wp), CONTIGUOUS_ARGUMENT(inout) :: ptr_tracer(:,:,:,:)
    REAL(wp), CONTIGUOUS_ARGUMENT(in) :: ptr_tke(:,:,:)
    TYPE(t_nh_diag)        , INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag        !<the atm phys vars
    TYPE(t_nwp_phy_tend)   , TARGET, INTENT(inout):: prm_nwp_tend    !< atm tend vars
    TYPE(t_external_data)  , INTENT(in)   :: ext_data

    REAL(wp)               , INTENT(in)   :: tcall_gscp_jg   !< time interval for 
                                                             !< microphysics
    LOGICAL                , INTENT(in)   :: lsatad          !< satad on/off

    LOGICAL                , INTENT(in)   :: lcompute_tt_lheat !< TRUE: store temperature tendency
                                                               ! due to microphysics for latent heat nudging

    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels !CK<
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_rlstart, i_rlend

    ! Variables for tendencies
    REAL(wp), DIMENSION(nproma,p_patch%nlev) :: ddt_tend_t , ddt_tend_qv, ddt_tend_qc, &
                                                ddt_tend_qi, ddt_tend_qr, ddt_tend_qs

    ! Local scalars:

    INTEGER :: jc,jb,jg,jk               !<block indices

    REAL(wp) :: zncn(nproma,p_patch%nlev),qnc(nproma,p_patch%nlev),qnc_s(nproma),rholoc,rhoinv
    LOGICAL  :: l_nest_other_micro
    LOGICAL  :: ldiag_ttend, ldiag_qtend



    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! domain ID
    jg = p_patch%id


    IF ( ASSOCIATED(prm_nwp_tend%ddt_temp_gscp)   ) THEN
      ldiag_ttend = .TRUE.
    ELSE
      ldiag_ttend = .FALSE.
    ENDIF
    IF ( ASSOCIATED(prm_nwp_tend%ddt_tracer_gscp) ) THEN
      ldiag_qtend = .TRUE.
    ELSE
      ldiag_qtend = .FALSE.
    ENDIF

    ! boundary conditions for number densities
    IF (jg > 1) THEN
       IF (atm_phy_nwp_config(jg)%inwp_gscp .ne. atm_phy_nwp_config(jg-1)%inwp_gscp) THEN
          l_nest_other_micro = .true.
       ELSE
          l_nest_other_micro = .false.
       END IF
    ELSE
      l_nest_other_micro = .false. 
    END IF

    !$acc data create(ddt_tend_t, ddt_tend_qv, ddt_tend_qc, ddt_tend_qi, ddt_tend_qr, ddt_tend_qs, &
    !$acc             zncn, qnc, qnc_s)

    SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
    CASE(4,5,6,7)

#ifdef _OPENACC
        CALL finish('mo_nwp_gscp_interface:','only graupel microphysics (inwp_gscp=2) is supported on GPU!')
#endif
       ! Update lateral boundaries of nested domains
       IF ( (l_limited_area.AND.jg==1) .OR. l_nest_other_micro) THEN

          IF (msg_level > 10) &
               & CALL message('mo_nwp_gscp_interface: ',"lateral boundaries for number densities")

          i_rlstart  = 1
          i_rlend    = grf_bdywidth_c
          i_startblk = p_patch%cells%start_blk(i_rlstart,1)
          i_endblk   = p_patch%cells%end_blk(i_rlend,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc,rholoc,rhoinv) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk, i_endblk

             CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                  i_startidx, i_endidx, i_rlstart, i_rlend)
             
             DO jk = 1, nlev
                DO jc = i_startidx, i_endidx
                   rholoc = p_prog%rho(jc,jk,jb)
                   rhoinv = 1.0_wp / rholoc
                   ptr_tracer(jc,jk,jb,iqnc) = set_qnc(ptr_tracer(jc,jk,jb,iqc)*rholoc)*rhoinv
                   ptr_tracer(jc,jk,jb,iqnr) = set_qnr(ptr_tracer(jc,jk,jb,iqr)*rholoc)*rhoinv
                   ptr_tracer(jc,jk,jb,iqni) = set_qni(ptr_tracer(jc,jk,jb,iqi)*rholoc)*rhoinv
                   ptr_tracer(jc,jk,jb,iqns) = set_qns(ptr_tracer(jc,jk,jb,iqs)*rholoc)*rhoinv
                   ptr_tracer(jc,jk,jb,iqng) = set_qng(ptr_tracer(jc,jk,jb,iqg)*rholoc)*rhoinv
                   ptr_tracer(jc,jk,jb,iqnh) = set_qnh(ptr_tracer(jc,jk,jb,iqh)*rholoc)*rhoinv
                ENDDO
             ENDDO
          ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
       ENDIF
    CASE (3)
       ! do something for QNI in case of gscp=3
    CASE DEFAULT
       ! Nothing to do for other schemes
    END SELECT

   
    ! exclude boundary interpolation zone of nested domains
    i_rlstart = grf_bdywidth_c+1
    i_rlend   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    ! Some run time diagnostics (can also be used for other schemes)
    IF (msg_level>14 .AND. atm_phy_nwp_config(jg)%l2moment) THEN
       CALL nwp_diag_output_minmax_micro(p_patch, p_prog, p_diag, ptr_tracer)
    END IF
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,zncn,qnc,qnc_s,ddt_tend_t,ddt_tend_qv,        &
!$OMP            ddt_tend_qc,ddt_tend_qi,ddt_tend_qr,ddt_tend_qs) ICON_OMP_GUIDED_SCHEDULE

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)


        IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 2) THEN

          ! Preparation for coupling of more advanced microphysics schemes (inwp_gscp>=2) with aerosol climatology
          ! Not yet implemented
          CALL ncn_from_tau_aerosol_speccnconst (nproma, nlev, i_startidx, i_endidx, kstart_moist(jg), nlev, &
            p_metrics%z_ifc(:,:,jb), prm_diag%aerosol(:,iss,jb), prm_diag%aerosol(:,iso4,jb),                &
            prm_diag%aerosol(:,iorg,jb), prm_diag%aerosol(:,idu,jb), zncn)

          CALL specccn_segalkhain (nproma, nlev, i_startidx, i_endidx, kstart_moist(jg), nlev, zncn,         &
            p_prog%w(:,:,jb), ptr_tracer(:,:,jb,iqc), p_prog%rho(:,:,jb), p_metrics%z_ifc(:,:,jb), qnc)

        ELSE IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 1) THEN

          IF (iprog_aero == 0) THEN
            !$acc parallel default(present)
            !$acc loop gang vector
            DO jc=i_startidx,i_endidx
              qnc_s(jc) = prm_diag%cloud_num(jc,jb)
            END DO
            !$acc end parallel
          ELSE

            CALL ncn_from_tau_aerosol_speccnconst (nproma, nlev, i_startidx, i_endidx, nlev, nlev, &
              p_metrics%z_ifc(:,:,jb), prm_diag%aerosol(:,iss,jb), prm_diag%aerosol(:,iso4,jb),    &
              prm_diag%aerosol(:,iorg,jb), prm_diag%aerosol(:,idu,jb), zncn)

            CALL specccn_segalkhain_simple (nproma, i_startidx, i_endidx, zncn(:,nlev), prm_diag%cloud_num(:,jb))

            ! Impose lower limit on cloud_num over land
            DO jc = i_startidx, i_endidx
              IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. ext_data%atm%llake_c(jc,jb)) &
                prm_diag%cloud_num(jc,jb) = MAX(175.e6_wp,prm_diag%cloud_num(jc,jb))
            ENDDO
!!$ UB: formally qnc_s is in the wrong unit (1/m^3) for the 1-moment schemes. Should be 1/kg.
!!$   However: since only the near-surface value of level nlev is used and the vertical profile is disregarded
!!$            anyways, we neglect this small near-surface difference and assume rho approx. 1.0. 
            qnc_s(i_startidx:i_endidx) = prm_diag%cloud_num(i_startidx:i_endidx,jb)

          ENDIF

        ELSE

          !$acc parallel default(present)
          !$acc loop gang vector
          DO jc=i_startidx,i_endidx
            qnc_s(jc) = cloud_num
          END DO
          !$acc end parallel

        ENDIF

        ! tt_lheat to be in used LHN
        ! lateron the updated p_diag%temp is added again
        IF (lcompute_tt_lheat) THEN
          !$acc parallel default(present)
          !$acc loop gang vector collapse(2)
          DO jk=1,nlev
            DO jc=i_startidx,i_endidx
              prm_diag%tt_lheat(jc,jk,jb) = prm_diag%tt_lheat(jc,jk,jb) - p_diag%temp(jc,jk,jb)
            ENDDO
          ENDDO
          !$acc end parallel
        ENDIF

        IF (timers_level > 10) CALL timer_start(timer_phys_micro_specific) 
        SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)

          
        CASE(0)  ! no microphysics scheme - in this case, this interface should not be called anyway
          
          WRITE(0,*) "                           "

        CASE(1)  ! COSMO-EU scheme (2-cat ice: cloud ice, snow)
                 ! version unified with COSMO scheme
                 ! unification version: COSMO_V4_23

          CALL cloudice (                                   &
            & nvec   =nproma                           ,    & !> in:  actual array size
            & ke     =nlev                             ,    & !< in:  actual array size
            & ivstart=i_startidx                       ,    & !< in:  start index of calculation
            & ivend  =i_endidx                         ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                 ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                    ,    & !< in:  timestep
            & qi0    =atm_phy_nwp_config(jg)%qi0       ,    & 
            & qc0    =atm_phy_nwp_config(jg)%qc0       ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)    ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)           ,    & !< inout:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)           ,    & !< in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )         ,    & !< in:  density
            & qv     =ptr_tracer (:,:,jb,iqv)   ,    & !< inout:  spec. humidity
            & qc     =ptr_tracer (:,:,jb,iqc)   ,    & !< inout:  cloud water
            & qi     =ptr_tracer (:,:,jb,iqi)   ,    & !< inout:  cloud ice
            & qr     =ptr_tracer (:,:,jb,iqr)   ,    & !< inout:  rain water
            & qs     =ptr_tracer (:,:,jb,iqs)   ,    & !< inout:  snow
            & qnc    = qnc_s                           ,    & !< cloud number concentration
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)    ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)    ,    & !< out: precipitation rate of snow
            & pri_gsp=prm_diag%ice_gsp_rate (:,jb)     ,    & !< out: precipitation rate of cloud ice
            & qrsflux= prm_diag%qrs_flux (:,:,jb)      ,    & !< out: precipitation flux
            & ldiag_ttend = ldiag_ttend                ,    & !< in:  if temp. tendency shall be diagnosed
            & ldiag_qtend = ldiag_qtend                ,    & !< in:  if moisture tendencies shall be diagnosed
            & ddt_tend_t  = ddt_tend_t                 ,    & !< out: tendency temperature
            & ddt_tend_qv = ddt_tend_qv                ,    & !< out: tendency QV
            & ddt_tend_qc = ddt_tend_qc                ,    & !< out: tendency QC
            & ddt_tend_qi = ddt_tend_qi                ,    & !< out: tendency QI
            & ddt_tend_qr = ddt_tend_qr                ,    & !< out: tendency QR
            & ddt_tend_qs = ddt_tend_qs                ,    & !< out: tendency QS
            & idbg=msg_level/2                         ,    &
            & l_cv=.TRUE.                              ,    &
            & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) !< in: latent heat choice
          
        CASE(2)  ! COSMO-DE (3-cat ice: snow, cloud ice, graupel)

          CALL graupel (                                     &
            & nvec   =nproma                            ,    & !> in:  actual array size
            & ke     =nlev                              ,    & !< in:  actual array size
            & ivstart=i_startidx                        ,    & !< in:  start index of calculation
            & ivend  =i_endidx                          ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                  ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                     ,    & !< in:  timestep
            & qi0    =atm_phy_nwp_config(jg)%qi0        ,    & 
            & qc0    =atm_phy_nwp_config(jg)%qc0        ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)     ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)            ,    & !< in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)            ,    & !< in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )          ,    & !< in:  density
            & qv     =ptr_tracer (:,:,jb,iqv)    ,    & !< in:  spec. humidity
            & qc     =ptr_tracer (:,:,jb,iqc)    ,    & !< in:  cloud water
            & qi     =ptr_tracer (:,:,jb,iqi)    ,    & !< in:  cloud ice
            & qr     =ptr_tracer (:,:,jb,iqr)    ,    & !< in:  rain water
            & qs     =ptr_tracer (:,:,jb,iqs)    ,    & !< in:  snow
            & qg     =ptr_tracer (:,:,jb,iqg)    ,    & !< in:  graupel
            & qnc    = qnc_s                            ,    & !< cloud number concentration
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)     ,    & !< out: precipitation rate of snow
            & pri_gsp=prm_diag%ice_gsp_rate (:,jb)      ,    & !< out: precipitation rate of cloud ice
            & prg_gsp=prm_diag%graupel_gsp_rate (:,jb)  ,    & !< out: precipitation rate of graupel
            & qrsflux= prm_diag%qrs_flux (:,:,jb)       ,    & !< out: precipitation flux
            & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
            & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
            & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
            & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
            & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
            & ddt_tend_qi = ddt_tend_qi                 ,    & !< out: tendency QI
            & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
            & ddt_tend_qs = ddt_tend_qs                 ,    & !< out: tendency QS
            & idbg=msg_level/2                          ,    &
            & l_cv=.TRUE.                               ,    &
            & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) !< in: latent heat choice

        CASE(3)  ! improved ice nucleation scheme by C. Koehler based on hydci_pp

          CALL hydci_pp_ice (                                &
            & nvec   =nproma                            ,    & !> in:  actual array size
            & ke     =nlev                              ,    & !< in:  actual array size
            & ke1    =nlevp1                            ,    & !< in:  half model levels (start index is 1)
            & ivstart=i_startidx                        ,    & !< in:  start index of calculation
            & ivend  =i_endidx                          ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                  ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                     ,    & !< in:  timestep
            & qi0    =atm_phy_nwp_config(jg)%qi0        ,    & 
            & qc0    =atm_phy_nwp_config(jg)%qc0        ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)     ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)            ,    & !< in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)            ,    & !< in:  full level pres
            & w      =p_prog%w(:,:,jb)                  ,    & !< in:  vertical wind speed, half levs (m/s)
            & tke    =ptr_tke(:,:,jb)            ,    & !< in:  turbulent kinetik energy
            & rho    =p_prog%rho    (:,:,jb  )          ,    & !< in:  density
            & qv     =ptr_tracer (:,:,jb,iqv)    ,    & !< in:  spec. humidity
            & qc     =ptr_tracer (:,:,jb,iqc)    ,    & !< in:  cloud water
            & qi     =ptr_tracer (:,:,jb,iqi)    ,    & !< in:  cloud ice
            & qni    =ptr_tracer (:,:,jb,iqni)   ,    & !< in:  cloud ice number     ( 1/kg)
            & qni_nuc=ptr_tracer (:,:,jb,iqni_nuc),   & !< in:  activated ice nuclei ( 1/kg)
            & qr     =ptr_tracer (:,:,jb,iqr)    ,    & !< in:  rain water
            & qs     =ptr_tracer (:,:,jb,iqs)    ,    & !< in:  snow
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)     ,    & !< out: precipitation rate of snow
            & qrsflux= prm_diag%qrs_flux (:,:,jb)       ,    & !< out: precipitation flux
            & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
            & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
            & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
            & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
            & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
            & ddt_tend_qi = ddt_tend_qi                 ,    & !< out: tendency QI
            & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
            & ddt_tend_qs = ddt_tend_qs                 ,    & !< out: tendency QS
            & idbg=msg_level/2                          ,    &
            & l_cv=.TRUE. )

        CASE(4)  ! two-moment scheme 

          CALL two_moment_mcrph(                       &
                       isize  = nproma,                &!in: array size
                       ke     = nlev,                  &!in: end level/array size
                       is     = i_startidx,            &!in: start index
                       ie     = i_endidx,              &!in: end index
                       ks     = kstart_moist(jg),      &!in: start level
                       dt     = tcall_gscp_jg ,        &!in: time step
                       dz     = p_metrics%ddqz_z_full(:,:,jb),  &!in: vertical layer thickness
                       hhl    = p_metrics%z_ifc(:,:,jb),        &!in: height of half levels
                       rho    = p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       qv     = ptr_tracer (:,:,jb,iqv), &!inout:sp humidity
                       qc     = ptr_tracer (:,:,jb,iqc), &!inout:cloud water
                       qnc    = ptr_tracer (:,:,jb,iqnc),&!inout: cloud droplet number
                       qr     = ptr_tracer (:,:,jb,iqr), &!inout:rain
                       qnr    = ptr_tracer (:,:,jb,iqnr),&!inout:rain droplet number
                       qi     = ptr_tracer (:,:,jb,iqi), &!inout: ice
                       qni    = ptr_tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = ptr_tracer (:,:,jb,iqs), &!inout: snow
                       qns    = ptr_tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = ptr_tracer (:,:,jb,iqg), &!inout: graupel
                       qng    = ptr_tracer (:,:,jb,iqng),&!inout: graupel number
                       qh     = ptr_tracer (:,:,jb,iqh), &!inout: hail
                       qnh    = ptr_tracer (:,:,jb,iqnh),&!inout: hail number
                       ninact = ptr_tracer (:,:,jb,ininact), &!inout: IN number
                       tk     = p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),  &!inout precp rate hail
                       qrsflux= prm_diag%qrs_flux(:,:,jb),      & !inout: 3D precipitation flux for LHN
                       msg_level = msg_level,                   &
                       & l_cv=.TRUE.,                           &
                       & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) !< in: latent heat choice

        CASE(5)  ! two-moment scheme with prognostic cloud droplet number
                 ! and budget equations for CCN and IN

          CALL two_moment_mcrph(                       &
                       isize  = nproma,                &!in: array size
                       ke     = nlev,                  &!in: end level/array size
                       is     = i_startidx,            &!in: start index
                       ie     = i_endidx,              &!in: end index
                       ks     = kstart_moist(jg),      &!in: start level
                       dt     = tcall_gscp_jg ,        &!in: time step
                       dz     = p_metrics%ddqz_z_full(:,:,jb),  &!in: vertical layer thickness
                       hhl    = p_metrics%z_ifc(:,:,jb),        &!in: height of half levels
                       rho    = p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       qv     = ptr_tracer (:,:,jb,iqv), &!inout: humidity
                       qc     = ptr_tracer (:,:,jb,iqc), &!inout: cloud water
                       qnc    = ptr_tracer (:,:,jb,iqnc),&!inout: cloud droplet number
                       qr     = ptr_tracer (:,:,jb,iqr), &!inout: rain
                       qnr    = ptr_tracer (:,:,jb,iqnr),&!inout: rain drop number
                       qi     = ptr_tracer (:,:,jb,iqi), &!inout: ice
                       qni    = ptr_tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = ptr_tracer (:,:,jb,iqs), &!inout: snow
                       qns    = ptr_tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = ptr_tracer (:,:,jb,iqg), &!inout: graupel
                       qng    = ptr_tracer (:,:,jb,iqng),&!inout: graupel number
                       qh     = ptr_tracer (:,:,jb,iqh), &!inout: hail
                       qnh    = ptr_tracer (:,:,jb,iqnh),&!inout: hail number
                       nccn   = ptr_tracer (:,:,jb,inccn),&!inout: CCN number
                       ninpot = ptr_tracer (:,:,jb,ininpot), &!inout: IN number
                       ninact = ptr_tracer (:,:,jb,ininact), &!inout: IN number
                       tk     = p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
                       qrsflux= prm_diag%qrs_flux(:,:,jb),      & !inout: 3D precipitation flux for LHN
                       msg_level = msg_level                ,    &
                       & l_cv=.TRUE.                        ,    &
                       & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) !< in: latent heat choice
#ifdef __ICON_ART
        CASE(6)  ! two-moment scheme with prognostic cloud droplet number
                 ! and chemical composition taken from the ART extension

          CALL art_clouds_interface_2mom(                        &
                       isize  = nproma,                          &!in: array size
                       ke     = nlev,                            &!in: end level/array size
                       jg     = jg,                              &!in: domain index
                       jb     = jb,                              &!in: block index
                       is     = i_startidx,                      &!in: start index
                       ie     = i_endidx,                        &!in: end index
                       ks     = kstart_moist(jg),                &!in: start level
                       dt     = tcall_gscp_jg ,                  &!in: time step
                       dz     = p_metrics%ddqz_z_full(:,:,jb),   &!in: vertical layer thickness
                       rho    = p_prog%rho(:,:,jb  )       ,     &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,     &!in:  pressure
                       tke    = ptr_tke(:,:,jb),          &!in:  turbulent kinetik energy
                       p_trac = ptr_tracer (:,:,jb,:),    &!inout: all tracers
                       tk     = p_diag%temp(:,:,jb),             &!inout: temp 
                       w      = p_prog%w(:,:,jb),                &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),   &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),    &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),   &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
! not impl yet!        qrsflux= prm_diag%qrs_flux(:,:,jb),      & !inout: 3D precipitation flux for LHN
                       tkvh   = prm_diag%tkvh(:,:,jb),           &!in: turbulent diffusion coefficients for heat     (m/s2 )
                       l_cv=.TRUE.     )
#endif
        CASE(7)  ! two-moment scheme with liquid water on graupel and hail (lwf scheme)

          CALL two_moment_mcrph(                       &
                       isize  = nproma,                &!in: array size
                       ke     = nlev,                  &!in: end level/array size
                       is     = i_startidx,            &!in: start index
                       ie     = i_endidx,              &!in: end index
                       ks     = kstart_moist(jg),      &!in: start level
                       dt     = tcall_gscp_jg ,        &!in: time step
                       dz     = p_metrics%ddqz_z_full(:,:,jb),  &!in: vertical layer thickness
                       hhl    = p_metrics%z_ifc(:,:,jb),        &!in: height of half levels
                       rho    = p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       qv     = ptr_tracer (:,:,jb,iqv), &!inout:sp humidity
                       qc     = ptr_tracer (:,:,jb,iqc), &!inout:cloud water
                       qnc    = ptr_tracer (:,:,jb,iqnc),&!inout: cloud droplet number
                       qr     = ptr_tracer (:,:,jb,iqr), &!inout:rain
                       qnr    = ptr_tracer (:,:,jb,iqnr),&!inout:rain droplet number
                       qi     = ptr_tracer (:,:,jb,iqi), &!inout: ice
                       qni    = ptr_tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = ptr_tracer (:,:,jb,iqs), &!inout: snow
                       qns    = ptr_tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = ptr_tracer (:,:,jb,iqg), &!inout: graupel
                       qng    = ptr_tracer (:,:,jb,iqng),&!inout: graupel number
                       qgl    = ptr_tracer (:,:,jb,iqgl),&!inout: liquid water on graupel
                       qh     = ptr_tracer (:,:,jb,iqh), &!inout: hail
                       qnh    = ptr_tracer (:,:,jb,iqnh),&!inout: hail number
                       qhl    = ptr_tracer (:,:,jb,iqhl),&!inout: liquid water on hail
                       ninact = ptr_tracer (:,:,jb,ininact), &!inout: IN number
                       tk     = p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
                       qrsflux= prm_diag%qrs_flux  (:,:,jb)     ,    & !inout: 3D precipitation flux for LHN
                       msg_level = msg_level                ,    &
                       & l_cv=.TRUE.                        ,    &
                       & ithermo_water=atm_phy_nwp_config(jg)%ithermo_water) !< in: latent heat choice

        CASE(9)  ! Kessler scheme (warm rain scheme)

          CALL kessler (                                     &
            & nvec   =nproma                            ,    & ! in:  actual array size
            & ke     =nlev                              ,    & ! in:  actual array size
            & ivstart =i_startidx                       ,    & ! in:  start index of calculation
            & ivend   =i_endidx                         ,    & ! in:  end index of calculation
            & kstart =kstart_moist(jg)                  ,    & ! in:  vertical start index
            & zdt    =tcall_gscp_jg                     ,    & ! in:  timestep
            & qc0    = atm_phy_nwp_config(jg)%qc0       ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)     ,    & ! in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)            ,    & ! in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)            ,    & ! in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )          ,    & ! in:  density
            & qv     =ptr_tracer (:,:,jb,iqv)    ,    & ! in:  spec. humidity
            & qc     =ptr_tracer (:,:,jb,iqc)    ,    & ! in:  cloud water
            & qr     =ptr_tracer (:,:,jb,iqr)    ,    & ! in:  rain water
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & ! out: precipitation rate of rain
            & qrsflux= prm_diag%qrs_flux (:,:,jb)       ,    & !< out: precipitation flux
            & ldiag_ttend = ldiag_ttend                 ,    & !< in:  if temp. tendency shall be diagnosed
            & ldiag_qtend = ldiag_qtend                 ,    & !< in:  if moisture tendencies shall be diagnosed
            & ddt_tend_t  = ddt_tend_t                  ,    & !< out: tendency temperature
            & ddt_tend_qv = ddt_tend_qv                 ,    & !< out: tendency QV
            & ddt_tend_qc = ddt_tend_qc                 ,    & !< out: tendency QC
            & ddt_tend_qr = ddt_tend_qr                 ,    & !< out: tendency QR
            & idbg   =msg_level/2                       ,    &
            & l_cv    =.TRUE. )

          IF (ldiag_qtend) THEN
            ddt_tend_qi(:,:) = 0._wp
            ddt_tend_qs(:,:) = 0._wp
          ENDIF

        CASE DEFAULT

          CALL finish('mo_nwp_gscp_interface', 'Unknown cloud physics scheme [1-5].')

        END SELECT

        IF (timers_level > 10) CALL timer_stop(timer_phys_micro_specific) 

        IF (ldiag_ttend) THEN
          !$acc parallel default(present)
          !$acc loop gang vector collapse(2)
          DO jk = kstart_moist(jg), nlev
            DO jc = i_startidx, i_endidx
              prm_nwp_tend%ddt_temp_gscp(jc,jk,jb) = ddt_tend_t(jc,jk)   ! tendency temperature
            ENDDO
          ENDDO
          !$acc end parallel
        ENDIF
        IF (ldiag_qtend) THEN
          !$acc parallel default(present)
          !$acc loop gang vector
          DO jk = kstart_moist(jg), nlev
            DO jc = i_startidx, i_endidx
              prm_nwp_tend%ddt_tracer_gscp(jc,jk,jb,iqv) = ddt_tend_qv(jc,jk)  ! tendency QV
              prm_nwp_tend%ddt_tracer_gscp(jc,jk,jb,iqc) = ddt_tend_qc(jc,jk)  ! tendency QC
              prm_nwp_tend%ddt_tracer_gscp(jc,jk,jb,iqi) = ddt_tend_qi(jc,jk)  ! tendency QI
              prm_nwp_tend%ddt_tracer_gscp(jc,jk,jb,iqr) = ddt_tend_qr(jc,jk)  ! tendency QR
              prm_nwp_tend%ddt_tracer_gscp(jc,jk,jb,iqs) = ddt_tend_qs(jc,jk)  ! tendency QS
            ENDDO
          ENDDO
          !$acc end parallel
        ENDIF


        !-------------------------------------------------------------------------
        !>
        !! Calculate surface precipitation
        !!
        !-------------------------------------------------------------------------
      
        IF (atm_phy_nwp_config(jg)%lcalc_acc_avg) THEN
          SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
          CASE(4,5,6,7)

#ifdef _OPENACC
           CALL finish('mo_nwp_gscp_interface:','only graupel microphysics (inwp_gscp=2) is supported on GPU!')
#endif

!DIR$ IVDEP
           DO jc =  i_startidx, i_endidx
             prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)                         &
                  &                   + tcall_gscp_jg * prm_diag%rain_gsp_rate (jc,jb)
             prm_diag%ice_gsp(jc,jb)  = prm_diag%ice_gsp(jc,jb)                          & 
                  &                   + tcall_gscp_jg * prm_diag%ice_gsp_rate (jc,jb)
             prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)                         &
                  &                   + tcall_gscp_jg * prm_diag%snow_gsp_rate (jc,jb)
             prm_diag%hail_gsp(jc,jb) = prm_diag%hail_gsp(jc,jb)                         &
                  &                   + tcall_gscp_jg * prm_diag%hail_gsp_rate (jc,jb)
             prm_diag%graupel_gsp(jc,jb) = prm_diag%graupel_gsp(jc,jb)                   &
                  &                   + tcall_gscp_jg * prm_diag%graupel_gsp_rate (jc,jb)
             !
             prm_diag%prec_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)  &
               &                      + prm_diag%ice_gsp(jc,jb)   &
               &                      + prm_diag%snow_gsp(jc,jb)  &
               &                      + prm_diag%hail_gsp(jc,jb)  &
               &                      + prm_diag%graupel_gsp(jc,jb)

           ENDDO

          CASE(2)

!DIR$ IVDEP
           !$acc parallel default(present)
           !$acc loop gang vector
           DO jc =  i_startidx, i_endidx

             prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)           &
               &                      + tcall_gscp_jg                      &
               &                      * prm_diag%rain_gsp_rate (jc,jb)
             prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)           &
               &                      + tcall_gscp_jg                      &
               &                      * prm_diag%snow_gsp_rate (jc,jb)
             prm_diag%ice_gsp(jc,jb) = prm_diag%ice_gsp(jc,jb)             &
               &                      + tcall_gscp_jg                      &
               &                      * prm_diag%ice_gsp_rate (jc,jb)
             prm_diag%graupel_gsp(jc,jb) = prm_diag%graupel_gsp(jc,jb)     &
               &                      + tcall_gscp_jg                      &
               &                      * prm_diag%graupel_gsp_rate (jc,jb)

             ! note: ice is deliberately excluded here because it predominantly contains blowing snow
             prm_diag%prec_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)  &
               &                      + prm_diag%snow_gsp(jc,jb)  &
               &                      + prm_diag%graupel_gsp(jc,jb)

           ENDDO
           !$acc end parallel

          CASE DEFAULT

!DIR$ IVDEP
           !$acc parallel default(present)
           !$acc loop gang vector
           DO jc =  i_startidx, i_endidx

             prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)         & 
               &                      + tcall_gscp_jg                    &
               &                      * prm_diag%rain_gsp_rate (jc,jb)
             prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)           &
               &                      + tcall_gscp_jg                      &
               &                      * prm_diag%snow_gsp_rate (jc,jb)
             prm_diag%ice_gsp(jc,jb) = prm_diag%ice_gsp(jc,jb)             &
               &                      + tcall_gscp_jg                      &
               &                      * prm_diag%ice_gsp_rate (jc,jb)

             ! note: ice is deliberately excluded here because it predominantly contains blowing snow
             prm_diag%prec_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)  &
               &                      + prm_diag%snow_gsp(jc,jb)

           ENDDO
           !$acc end parallel

          END SELECT
        ENDIF

        ! saturation adjustment after microphysics
        ! - this is the second satad call
        ! - first satad in physics interface before microphysics

        IF (timers_level > 10) CALL timer_start(timer_phys_micro_satad) 
        IF (lsatad) THEN

          IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN   ! EDMF DUALM: no satad in PBL
 
            CALL satad_v_3d(                                 &           
               & maxiter  = 10                            ,& !> IN
               & tol      = 1.e-3_wp                      ,& !> IN
               & te       = p_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = ptr_tracer (:,:,jb,iqv),& !> INOUT
               & qce      = ptr_tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = p_prog%rho        (:,:,jb)    ,& !> IN
               & qtvar    = ptr_tracer (:,:,jb,iqtvar) ,& !> IN
               & idim     = nproma                        ,& !> IN
               & kdim     = nlev                          ,& !> IN
               & ilo      = i_startidx                    ,& !> IN
               & iup      = i_endidx                      ,& !> IN
               & klo      = kstart_moist(jg)              ,& !> IN
               & kup      = nlev                           & !> IN
               )

          ELSE

#ifdef _OPENACC
            CALL satad_v_3d_gpu(                             &
#else
            CALL satad_v_3d(                                 &
#endif
               & maxiter  = 10                            ,& !> IN
               & tol      = 1.e-3_wp                      ,& !> IN
               & te       = p_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = ptr_tracer (:,:,jb,iqv),& !> INOUT
               & qce      = ptr_tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = p_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                        ,& !> IN
               & kdim     = nlev                          ,& !> IN
               & ilo      = i_startidx                    ,& !> IN
               & iup      = i_endidx                      ,& !> IN
               & klo      = kstart_moist(jg)              ,& !> IN
               & kup      = nlev                           & !> IN
               )

          ENDIF

        ENDIF

        IF (timers_level > 10) CALL timer_stop(timer_phys_micro_satad) 

        ! Update tt_lheat to be used in LHN
        IF (lcompute_tt_lheat) THEN
            !$acc parallel default(present)
            !$acc loop gang vector collapse(2)
            DO jk=1,nlev
              DO jc=i_startidx,i_endidx
                prm_diag%tt_lheat(jc,jk,jb) = prm_diag%tt_lheat(jc,jk,jb) + p_diag%temp(jc,jk,jb)
              ENDDO
            ENDDO
            !$acc end parallel
        ENDIF

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

 
    ! Some more run time diagnostics (can also be used for other schemes)
    IF (msg_level>14 .AND. atm_phy_nwp_config(jg)%l2moment) THEN
       CALL nwp_diag_output_minmax_micro(p_patch, p_prog, p_diag, ptr_tracer)
    END IF

    !$acc end data
     
  END SUBROUTINE nwp_microphysics

END MODULE mo_nwp_gscp_interface

