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
!----------------------------

MODULE mo_nwp_gscp_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message
  USE mo_parallel_config,      ONLY: nproma

  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, iss, iorg, iso4, idu
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, iqr, iqs,       &
                                     iqni, iqni_nuc, iqg, iqh, iqnr, iqns,     &
                                     iqng, iqnh, iqnc, inccn, ininpot, ininact
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, iprog_aero
  USE gscp_kessler,            ONLY: kessler
  USE gscp_cloudice,           ONLY: cloudice
  USE gscp_graupel,            ONLY: graupel
  USE gscp_hydci_pp_ice,       ONLY: hydci_pp_ice
  USE mo_exception,            ONLY: finish
  USE mo_mcrph_sb,             ONLY: two_moment_mcrph, set_qnc, &
       &                             set_qnr,set_qni,set_qns,set_qng
  USE mo_art_clouds_interface, ONLY: art_clouds_interface_2mom
  USE mo_nwp_diagnosis,        ONLY: nwp_diag_output_minmax_micro
  USE gscp_data,               ONLY: cloud_num
  USE mo_cpl_aerosol_microphys,ONLY: specccn_segalkhain, ncn_from_tau_aerosol_speccnconst, &
                                     specccn_segalkhain_simple
  USE mo_grid_config,          ONLY: l_limited_area
  USE mo_satad,                ONLY: satad_v_3D

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
                            &   p_prog_rcf,                   & !>inout
                            &   p_diag ,                      & !>inout
                            &   prm_diag                      ) !>inout 



    TYPE(t_patch)          , INTENT(in)   :: p_patch        !!<grid/patch info.
    TYPE(t_nh_metrics)     , INTENT(in)   :: p_metrics
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog          !<the dyn prog vars
    TYPE(t_nh_prog)        , INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag)        , INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag        !<the atm phys vars

    REAL(wp)               , INTENT(in)   :: tcall_gscp_jg   !< time interval for 
                                                             !< microphysics
    LOGICAL                , INTENT(in)   :: lsatad          !< satad on/off

    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels !CK<
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_rlstart, i_rlend

    ! Local scalars:

    INTEGER :: jc,jb,jg,jk               !<block indices

    REAL(wp) :: zncn(nproma,p_patch%nlev),qnc(nproma,p_patch%nlev),qnc_s(nproma)
    LOGICAL  :: l_nest_other_micro
    LOGICAL  :: ltwomoment

    ! local variables
    !

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1 !CK<

    ! domain ID
    jg = p_patch%id

    ! logical for SB two-moment scheme
    ltwomoment = ( atm_phy_nwp_config(jg)%inwp_gscp==4 .OR. atm_phy_nwp_config(jg)%inwp_gscp==5 .OR. &
               &   atm_phy_nwp_config(jg)%inwp_gscp==6)

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

    SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
    CASE(4,5,6)

       ! Update lateral boundaries of nested domains
       IF ( (l_limited_area.AND.jg==1) .OR. l_nest_other_micro) THEN

          IF (msg_level > 10) &
               & CALL message('mo_nwp_gscp_interface: ',"lateral boundaries for number densities")

          i_rlstart  = 1
          i_rlend    = grf_bdywidth_c
          i_startblk = p_patch%cells%start_blk(i_rlstart,1)
          i_endblk   = p_patch%cells%end_blk(i_rlend,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk, i_endblk

             CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                  i_startidx, i_endidx, i_rlstart, i_rlend)
             
             DO jk = 1, nlev
                DO jc = i_startidx, i_endidx
                   p_prog_rcf%tracer(jc,jk,jb,iqnc) = set_qnc(p_prog_rcf%tracer(jc,jk,jb,iqc))
                   p_prog_rcf%tracer(jc,jk,jb,iqnr) = set_qnr(p_prog_rcf%tracer(jc,jk,jb,iqr))
                   p_prog_rcf%tracer(jc,jk,jb,iqni) = set_qni(p_prog_rcf%tracer(jc,jk,jb,iqi))
                   p_prog_rcf%tracer(jc,jk,jb,iqns) = set_qns(p_prog_rcf%tracer(jc,jk,jb,iqs))
                   p_prog_rcf%tracer(jc,jk,jb,iqng) = set_qng(p_prog_rcf%tracer(jc,jk,jb,iqg))
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
    IF (msg_level>14 .AND. ltwomoment) THEN
       CALL nwp_diag_output_minmax_micro(p_patch, p_prog, p_diag, p_prog_rcf)
    END IF
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zncn,qnc,qnc_s) ICON_OMP_GUIDED_SCHEDULE
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
            p_prog%w(:,:,jb), p_prog_rcf%tracer(:,:,jb,iqc), p_prog%rho(:,:,jb), p_metrics%z_ifc(:,:,jb), qnc)

        ELSE IF (atm_phy_nwp_config(jg)%icpl_aero_gscp == 1) THEN

          IF (iprog_aero == 0) THEN
            qnc_s(i_startidx:i_endidx) = prm_diag%cloud_num(i_startidx:i_endidx,jb)
          ELSE

            CALL ncn_from_tau_aerosol_speccnconst (nproma, nlev, i_startidx, i_endidx, nlev, nlev, &
              p_metrics%z_ifc(:,:,jb), prm_diag%aerosol(:,iss,jb), prm_diag%aerosol(:,iso4,jb),    &
              prm_diag%aerosol(:,iorg,jb), prm_diag%aerosol(:,idu,jb), zncn)

            CALL specccn_segalkhain_simple (nproma, i_startidx, i_endidx, zncn(:,nlev), prm_diag%cloud_num(:,jb))
            qnc_s(i_startidx:i_endidx) = prm_diag%cloud_num(i_startidx:i_endidx,jb)

          ENDIF

        ELSE

          qnc_s(i_startidx:i_endidx) = cloud_num

        ENDIF


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
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)   ,    & !< inout:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)   ,    & !< inout:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)   ,    & !< inout:  cloud ice
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)   ,    & !< inout:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)   ,    & !< inout:  snow
            & qnc    = qnc_s                           ,    & !< cloud number concentration
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)    ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)    ,    & !< out: precipitation rate of snow
            & idbg=msg_level/2                         ,    &
            & l_cv=.TRUE. )
          

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
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)    ,    & !< in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)    ,    & !< in:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)    ,    & !< in:  cloud ice
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)    ,    & !< in:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)    ,    & !< in:  snow
            & qg     =p_prog_rcf%tracer (:,:,jb,iqg)    ,    & !< in:  graupel
            & qnc    = qnc_s                            ,    & !< cloud number concentration
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)     ,    & !< out: precipitation rate of snow
            & prg_gsp=prm_diag%graupel_gsp_rate (:,jb)  ,    & !< out: precipitation rate of snow
            & idbg=msg_level/2                          ,    &
            & l_cv=.TRUE. )
          
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
            & tke    =p_prog_rcf%tke(:,:,jb)            ,    & !< in:  turbulent kinetik energy
            & rho    =p_prog%rho    (:,:,jb  )          ,    & !< in:  density
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)    ,    & !< in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)    ,    & !< in:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)    ,    & !< in:  cloud ice
            & qni    =p_prog_rcf%tracer (:,:,jb,iqni)   ,    & !< in:  cloud ice number     ( 1/kg)
            & qni_nuc=p_prog_rcf%tracer (:,:,jb,iqni_nuc),   & !< in:  activated ice nuclei ( 1/kg)            
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)    ,    & !< in:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)    ,    & !< in:  snow
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)     ,    & !< out: precipitation rate of snow
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
                       rho    = p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       qv     = p_prog_rcf%tracer (:,:,jb,iqv), &!inout:sp humidity
                       qc     = p_prog_rcf%tracer (:,:,jb,iqc), &!inout:cloud water
                       qnc    = p_prog_rcf%tracer (:,:,jb,iqnc),&!inout: cloud droplet number 
                       qr     = p_prog_rcf%tracer (:,:,jb,iqr), &!inout:rain
                       qnr    = p_prog_rcf%tracer (:,:,jb,iqnr),&!inout:rain droplet number 
                       qi     = p_prog_rcf%tracer (:,:,jb,iqi), &!inout: ice
                       qni    = p_prog_rcf%tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = p_prog_rcf%tracer (:,:,jb,iqs), &!inout: snow 
                       qns    = p_prog_rcf%tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = p_prog_rcf%tracer (:,:,jb,iqg), &!inout: graupel 
                       qng    = p_prog_rcf%tracer (:,:,jb,iqng),&!inout: graupel number
                       qh     = p_prog_rcf%tracer (:,:,jb,iqh), &!inout: hail 
                       qnh    = p_prog_rcf%tracer (:,:,jb,iqnh),&!inout: hail number
                       ninact = p_prog_rcf%tracer (:,:,jb,ininact), &!inout: IN number
                       tk     = p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
                       msg_level = msg_level                ,    &
                       l_cv=.TRUE.          )    

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
                       rho    = p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       qv     = p_prog_rcf%tracer (:,:,jb,iqv), &!inout: humidity
                       qc     = p_prog_rcf%tracer (:,:,jb,iqc), &!inout: cloud water
                       qnc    = p_prog_rcf%tracer (:,:,jb,iqnc),&!inout: cloud droplet number 
                       qr     = p_prog_rcf%tracer (:,:,jb,iqr), &!inout: rain
                       qnr    = p_prog_rcf%tracer (:,:,jb,iqnr),&!inout: rain drop number 
                       qi     = p_prog_rcf%tracer (:,:,jb,iqi), &!inout: ice
                       qni    = p_prog_rcf%tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = p_prog_rcf%tracer (:,:,jb,iqs), &!inout: snow 
                       qns    = p_prog_rcf%tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = p_prog_rcf%tracer (:,:,jb,iqg), &!inout: graupel 
                       qng    = p_prog_rcf%tracer (:,:,jb,iqng),&!inout: graupel number
                       qh     = p_prog_rcf%tracer (:,:,jb,iqh), &!inout: hail 
                       qnh    = p_prog_rcf%tracer (:,:,jb,iqnh),&!inout: hail number
                       nccn   = p_prog_rcf%tracer (:,:,jb,inccn),&!inout: CCN number
                       ninpot = p_prog_rcf%tracer (:,:,jb,ininpot), &!inout: IN number
                       ninact = p_prog_rcf%tracer (:,:,jb,ininact), &!inout: IN number
                       tk     = p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
                       msg_level = msg_level                ,    &
                       l_cv=.TRUE.     )
    
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
                       tke    = p_prog_rcf%tke(:,:,jb),          &!in:  turbulent kinetik energy
                       p_trac = p_prog_rcf%tracer (:,:,jb,:),    &!inout: all tracers
                       tk     = p_diag%temp(:,:,jb),             &!inout: temp 
                       w      = p_prog%w(:,:,jb),                &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),   &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),    &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),   &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
                       tkvh   = prm_diag%tkvh(:,:,jb),           &!in: turbulent diffusion coefficients for heat     (m/s2 )
                       msg_level = msg_level,                    &
                       l_cv=.TRUE.     )
    
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
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)    ,    & ! in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)    ,    & ! in:  cloud water
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)    ,    & ! in:  rain water
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & ! out: precipitation rate of rain
            & idbg   =msg_level/2                       ,    &
            & l_cv    =.TRUE. )


        CASE DEFAULT

          CALL finish('mo_nwp_gscp_interface', 'Unknown cloud physics scheme [1-5].')

        END SELECT

        !-------------------------------------------------------------------------
        !>
        !! Calculate surface precipitation
        !!
        !-------------------------------------------------------------------------
      
        IF (atm_phy_nwp_config(jg)%lcalc_acc_avg) THEN
          SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
          CASE(4,5,6)

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
              prm_diag%tot_prec(jc,jb) = prm_diag%tot_prec(jc,jb)                         &
                   &                   + tcall_gscp_jg * ( prm_diag%rain_gsp_rate (jc,jb) &
                   &                                     + prm_diag%ice_gsp_rate (jc,jb)  &
                   &                                     + prm_diag%snow_gsp_rate (jc,jb) &
                   &                                     + prm_diag%hail_gsp_rate (jc,jb) &
                   &                                     + prm_diag%graupel_gsp_rate (jc,jb) )
           ENDDO

          CASE(2)

!DIR$ IVDEP
           DO jc =  i_startidx, i_endidx

            prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)           &
   &                                 + tcall_gscp_jg                      &
   &                                 * prm_diag%rain_gsp_rate (jc,jb)
            prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)           &
   &                                 + tcall_gscp_jg                      &
   &                                 * prm_diag%snow_gsp_rate (jc,jb)
            prm_diag%graupel_gsp(jc,jb) = prm_diag%graupel_gsp(jc,jb)     &
   &                                 + tcall_gscp_jg                      &
   &                                 * prm_diag%graupel_gsp_rate (jc,jb)
            prm_diag%tot_prec(jc,jb) = prm_diag%tot_prec(jc,jb)           &
   &                                 +  tcall_gscp_jg                     &
   &                                 * ( prm_diag%rain_gsp_rate (jc,jb)   &
   &                                 +   prm_diag%snow_gsp_rate (jc,jb)   &
   &                                 +   prm_diag%graupel_gsp_rate (jc,jb) )
           ENDDO

          CASE DEFAULT

!DIR$ IVDEP
           DO jc =  i_startidx, i_endidx

            prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)         & 
   &                                 + tcall_gscp_jg                    &
   &                                 * prm_diag%rain_gsp_rate (jc,jb)
            prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)         &
   &                                 + tcall_gscp_jg                    &
   &                                 * prm_diag%snow_gsp_rate (jc,jb)
            prm_diag%tot_prec(jc,jb) = prm_diag%tot_prec(jc,jb)         &
   &                                 +  tcall_gscp_jg                   &
   &                                 * ( prm_diag%rain_gsp_rate (jc,jb) &
   &                                 +   prm_diag%snow_gsp_rate (jc,jb) )
           ENDDO

          END SELECT
        ENDIF

        ! saturation adjustment after microphysics
        ! - this is the second satad call
        ! - first satad in physics interface before microphysics

        IF (lsatad) THEN

          CALL satad_v_3d(                                 &           
               & maxiter  = 10                            ,& !> IN
               & tol      = 1.e-3_wp                      ,& !> IN
               & te       = p_diag%temp       (:,:,jb)    ,& !> INOUT
               & qve      = p_prog_rcf%tracer (:,:,jb,iqv),& !> INOUT
               & qce      = p_prog_rcf%tracer (:,:,jb,iqc),& !> INOUT
               & rhotot   = p_prog%rho        (:,:,jb)    ,& !> IN
               & idim     = nproma                        ,& !> IN
               & kdim     = nlev                          ,& !> IN
               & ilo      = i_startidx                    ,& !> IN
               & iup      = i_endidx                      ,& !> IN
               & klo      = kstart_moist(jg)              ,& !> IN
               & kup      = nlev                           & !> IN
               )

        ENDIF

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
 
    ! Some more run time diagnostics (can also be used for other schemes)
    IF (msg_level>14 .AND. ltwomoment) THEN
       CALL nwp_diag_output_minmax_micro(p_patch, p_prog, p_diag, p_prog_rcf)
    END IF

     
  END SUBROUTINE nwp_microphysics

END MODULE mo_nwp_gscp_interface

