!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! turbulence parameterisations:
!! inwp_turb == 1 == turbulence scheme by M. Raschendorfer run in COSMO
!! inwp_turb == 2 == turbulence scheme imported from the GME
!! This module handles the atmospheric turbulent diffusion, only.
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

MODULE mo_nwp_turbdiff_interface

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_model_domain,           ONLY: t_patch
  USE mo_impl_constants,         ONLY: min_rlcell_int, igme, icosmo, &
    &                                  max_ntracer
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_physical_constants,     ONLY: lh_v=>alv
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,          ONLY: phy_params 
  USE mo_nwp_lnd_types,          ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,        ONLY: nproma
  USE mo_run_config,             ONLY: msg_level, iqv, iqc, iqi, iqtke
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_nonhydrostatic_config,  ONLY: kstart_moist
  USE mo_data_turbdiff,          ONLY: get_turbdiff_param, lsflcnd
  USE src_turbdiff,              ONLY: organize_turbdiff, modvar
  USE mo_gme_turbdiff,           ONLY: partura, progimp_turb

  USE mo_art_config,             ONLY: art_config
  USE mo_advection_config,       ONLY: advection_config
  USE mo_turbdiff_config,        ONLY: turbdiff_config
  USE mo_art_turbdiff_interface, ONLY: art_turbdiff_interface

  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_turbdiff

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbdiff  ( tcall_turb_jg,                     & !>in
                          & p_patch,                          & !>in
                          & p_metrics,                        & !>in
                          & ext_data,                         & !>in
                          & p_prog,                           & !>in
                          & p_prog_now_rcf,                   & !>in
                          & p_prog_rcf,                       & !>inout
                          & p_diag ,                          & !>inout
                          & prm_diag, prm_nwp_tend,           & !>inout
                          & wtr_prog_now,                     & !>in 
                          & lnd_prog_now,                     & !>in 
                          & lnd_diag                          ) !>in


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
  TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog),      TARGET,INTENT(in)   :: p_prog_now_rcf  !<progs with red.
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
  TYPE(t_wtr_prog),            INTENT(in)   :: wtr_prog_now    !< prog vars for wtr
  TYPE(t_lnd_prog),            INTENT(in)   :: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for 
                                                               !< turbulence

  ! Local array bounds

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: i_nchdom                !< domain index

  ! Local scalars:

  INTEGER :: jc,jk,jb,jg      !loop indices

  ! local variables for turbdiff

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''

  INTEGER  :: nlev, nlevp1, nlevcm                  !< number of full, half and canopy levels

  REAL(wp) :: tke_inc_ic(nproma)                    !< TKE increment at half levels

  REAL(wp) :: l_hori(nproma)                        !< horizontal length scale

  REAL(wp) :: z_tvs(nproma,p_patch%nlevp1,1)        !< aux turbulence velocity scale [m/s]

  ! type structure to hand over additional tracers to turbdiff
  TYPE(modvar) :: ptr(max_ntracer)

  INTEGER :: ncloud_offset                          !< offset for ptr-indexing in ART 
                                                    !< interface due to additionally 
                                                    !< diffused cloud fields

!--------------------------------------------------------------


  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  nlevcm = nlevp1

  IF (msg_level >= 15) CALL message('mo_nwp_turbdiff:', 'turbulence')
    
  ! local variables related to the blocking
  
  i_nchdom  = MAX(1,p_patch%n_childdom)
  jg        = p_patch%id
  
  ! exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

  
  IF ( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) THEN
     CALL get_turbdiff_param(jg)
  ENDIF



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,ierrstat,errormsg,eroutine,tke_inc_ic,z_tvs, &
!$OMP            ncloud_offset,ptr,l_hori)  &
!$OMP ICON_OMP_GUIDED_SCHEDULE

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

   !-------------------------------------------------------------------------
   !<  turbulent diffusion
   !-------------------------------------------------------------------------

   !<  NOTE: since  turbulence is a fast process it is
   !!        allowed to do a sequential updating except for wind speed
   !!        because back-and-forth interpolation would cause too large errors
   !!  (GZ, 2011-08-29): Nevertheless, tendency fields are now passed to turbdiff
   !!        to have them available for extended diagnostic output


    IF ( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) THEN

!-------------------------------------------------------------------------
!< COSMO turbulence scheme by M. Raschendorfer  
!-------------------------------------------------------------------------


      !
      ! convert TKE to the turbulence velocity scale SQRT(2*TKE) as required by turbdiff
      ! INPUT to turbdiff is timestep now
      DO jk=1, nlevp1
        DO jc=i_startidx, i_endidx
          z_tvs(jc,jk,1) = SQRT(2._wp* (p_prog_now_rcf%tke(jc,jk,jb))) 
         ENDDO
      ENDDO



      IF (advection_config(jg)%iadv_tke > 0) THEN
        ! Interpolate advective tvs tendency from full levels to half levels
        ! Note that both the advective TKE tendency and ddt_tke actually carry time tendencies
        ! of tvs; attempts to horizontally advect TKE failed because of numerical instability
        !
        DO jk=2, nlev
          DO jc=i_startidx, i_endidx

            tke_inc_ic(jc) = p_metrics%wgtfac_c(jc,jk,jb) * p_diag%ddt_tracer_adv(jc,jk,jb,iqtke) &
              &             + (1._wp - p_metrics%wgtfac_c(jc,jk,jb)) &
              &             * p_diag%ddt_tracer_adv(jc,jk-1,jb,iqtke)

            ! add advective TKE (actually tvs) tendency to ddt_tke, which is provided to turbdiff as input
            prm_nwp_tend%ddt_tke(jc,jk,jb) = prm_nwp_tend%ddt_tke(jc,jk,jb) + tke_inc_ic(jc)

          ENDDO  ! jc
        ENDDO  ! jk
        !
        DO jc=i_startidx, i_endidx

          ! zero gradient assumption for TKE increment at model bottom (top level not needed)
          prm_nwp_tend%ddt_tke(jc,nlevp1,jb) = prm_nwp_tend%ddt_tke(jc,nlevp1,jb) + p_diag%ddt_tracer_adv(jc,nlev,jb,iqtke)

        ENDDO  ! jc
      ENDIF



      ierrstat = 0

      !KF tendencies  have to be set to zero
      !GZ: this should be replaced by an appropriate switch in turbdiff
      prm_nwp_tend%ddt_u_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_v_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_temp_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) = 0._wp


      IF (turbdiff_config(jg)%ldiff_qi) THEN
        ! offset for ptr-indexing in ART-Interface
        ncloud_offset = 1

        ! reinit
        prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi) = 0._wp

        ! register cloud ice for turbulent diffusion
        ptr(1)%av => p_prog_rcf%tracer(:,:,jb,iqi)
        ptr(1)%at => prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi)
        ptr(1)%sv => NULL()
!!$        ptr(1)%fc =  .TRUE.
      ELSE
        ! offset for ptr-indexing in ART-Interface
        ncloud_offset = 0
      ENDIF

      CALL art_turbdiff_interface( 'setup_ptr', p_patch, p_prog_rcf, prm_nwp_tend,  &
        &                          ncloud_offset=ncloud_offset,                     &
        &                          ptr=ptr(:), dt=tcall_turb_jg,                    &
        &                          p_rho=p_prog%rho(:,:,:),                         &
        &                          p_metrics=p_metrics,                             &
        &                          p_diag=p_diag, prm_diag=prm_diag,                &
        &                          jb=jb )

      !should be dependent on location in furture!
      l_hori(i_startidx:i_endidx)=phy_params(jg)%mean_charlen

      ! turbdiff
      CALL organize_turbdiff( &
        &  iini=0, lturatm=.TRUE. , ltursrf=.FALSE., lstfnct=.TRUE. , & !atmosph. turbulence and vertical diffusion
        &          lnsfdia=.FALSE., ltkeinp=.FALSE., lgz0inp=.FALSE., & !but no surface-layer turbulence (turbtran)
        &  itnd=0, lum_dif=.TRUE. , lvm_dif=.TRUE. , lscadif=.TRUE. , & !and thus (implicitly) neither surface-layer diagn.
        &          lsrflux=.FALSE., lsfluse=lsflcnd, lqvcrst=.FALSE., & !nor surface-flux calculation (both in turbtran)
!
        &  dt_var=tcall_turb_jg, dt_tke=tcall_turb_jg,                                & !in
        &  nprv=1, ntur=1, ntim=1, trop_mask=prm_diag%tropics_mask(:,jb),             & !in
        &  ie=nproma, ke=nlev, ke1=nlevp1, kcm=nlevcm,                                & !in
        &  i_st=i_startidx, i_en=i_endidx, i_stp=i_startidx, i_enp=i_endidx,          &
        &  l_hori=l_hori, hhl=p_metrics%z_ifc(:,:,jb),                                & !in
        &  dp0=p_diag%dpres_mc(:,:,jb),                                               & !in
        &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb),  & !in
        &  h_ice=wtr_prog_now%h_ice(:,jb),                                            & !in
        &  sai=ext_data%atm%sai(:,jb), d_pat = ext_data%atm%sso_stdh_raw(:,jb),       &
        &  tkred_sfc=prm_diag%tkred_sfc(:,jb),  gz0=prm_diag%gz0(:,jb),               & !inout 
        &  t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),                      & !in
        &  ps=p_diag%pres_sfc(:,jb),                                                  & !in
        &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb),                & !in
        &  t=p_diag%temp(:,:,jb), prs=p_diag%pres(:,:,jb),                            & !in
        &  rho=p_prog%rho(:,:,jb), epr=p_prog%exner(:,:,jb),                          & !in
        &  qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),        & !in
        &  ptr=ptr(:), ndtr=art_config(jg)%nturb_tracer,                              & !diffusion of additional tracer variables!
        &  tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),                            & !out
        &  tvm=prm_diag%tvm(:,jb), tvh=prm_diag%tvh(:,jb),                            & !inout
        &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb),    & !inout
        &  tke=z_tvs(:,:,:),                                                          & !inout
        &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh(:,:,jb),                    & !inout
        &  rcld=prm_diag%rcld(:,:,jb),                                                & !inout
        &  hdef2=p_diag%hdef_ic(:,:,jb), hdiv=p_diag%div_ic(:,:,jb),                  & !in
        &  dwdx=p_diag%dwdx(:,:,jb), dwdy=p_diag%dwdy(:,:,jb),                        & !in
        &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb),                                    & !inout
        &  v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb),                                    & !inout
        &  t_tens=prm_nwp_tend%ddt_temp_turb(:,:,jb),                                 & !inout
        &  qv_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),                          & !inout
        &  qc_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc),                          & !inout
        &  tketens=prm_nwp_tend%ddt_tke(:,:,jb),                                      & !inout
        &  ut_sso=REAL(prm_nwp_tend%ddt_u_sso(:,:,jb),wp),                            & !in
        &  vt_sso=REAL(prm_nwp_tend%ddt_v_sso(:,:,jb),wp),                            & !in
        &  tket_conv=prm_nwp_tend%ddt_tke_pconv(:,:,jb),                              & !in
        &  tket_hshr=prm_nwp_tend%ddt_tke_hsh(:,:,jb),                                & !out
        &  shfl_s=prm_diag%shfl_s(:,jb), qvfl_s=prm_diag%qhfl_s(:,jb),                & !in
        &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

      ! preparation for concentration boundary condition. Usually inactive for standard ICON runs.
      IF ( .NOT. lsflcnd ) THEN
        prm_diag%lhfl_s(i_startidx:i_endidx,jb) = &
          &  prm_diag%qhfl_s(i_startidx:i_endidx,jb) * lh_v
      END IF

      CALL art_turbdiff_interface( 'update_ptr', p_patch, p_prog_rcf, prm_nwp_tend,  &
        &                          ncloud_offset=ncloud_offset,                      &
        &                          ptr=ptr(:), dt=tcall_turb_jg,                     &
        &                          i_st=i_startidx, i_en=i_endidx )


      IF (ierrstat.NE.0) THEN
        CALL finish(eroutine, errormsg)
      END IF

      ! transform updated turbulent velocity scale back to TKE
      ! Note: ddt_tke is purely diagnostic and has already been added to z_tvs
      DO jk=1, nlevp1
        DO jc=i_startidx, i_endidx
          p_prog_rcf%tke(jc,jk,jb) = 0.5_wp*(z_tvs(jc,jk,1))**2
        ENDDO
      ENDDO


      ! Interpolate updated TKE (actually tvs) back to main levels
      ! Note that TKE at lowest main level is re-computed in nwp_turbtrans, after surface TKE
      ! has been updated.
      IF (advection_config(jg)%iadv_tke > 0) THEN
        DO jk=1, nlev
          DO jc=i_startidx, i_endidx
            p_prog_rcf%tracer(jc,jk,jb,iqtke) = 0.5_wp* ( z_tvs(jc,jk,1) + z_tvs(jc,jk+1,1) )
          ENDDO
        ENDDO
      ENDIF


    ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == igme ) THEN

!-------------------------------------------------------------------------
!> GME turbulence scheme 
!-------------------------------------------------------------------------

      ! turbulent diffusion coefficients in atmosphere
      CALL partura( zh=p_metrics%z_ifc(:,:,jb), zf=p_metrics%z_mc(:,:,jb),                 & !in
        &           u=p_diag%u(:,:,jb),         v=p_diag%v(:,:,jb), t=p_diag%temp(:,:,jb), & !in
        &           qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),    & !in
        &           ph=p_diag%pres_ifc(:,:,jb), pf=p_diag%pres(:,:,jb),                    & !in
        &           ie=nproma, ke=nlev, ke1=nlevp1,                                        & !in
        &           i_startidx=i_startidx, i_endidx=i_endidx,                              & !in
        &           tkvm=prm_diag%tkvm(:,2:nlev,jb), tkvh=prm_diag%tkvh(:,2:nlev,jb)       ) !inout


      ! tendencies from turbulent diffusion
      CALL progimp_turb( t=p_diag%temp(:,:,jb), qv=p_prog_rcf%tracer(:,:,jb,iqv),      & !in
        &                qc=p_prog_rcf%tracer(:,:,jb,iqc),                             & !in
        &                u=p_diag%u(:,:,jb),    v=p_diag%v(:,:,jb),                    & !in
        &                zh=p_metrics%z_ifc(:,:,jb), zf=p_metrics%z_mc(:,:,jb),        & !in
        &                rho=p_prog%rho(:,:,jb), ps=p_diag%pres_ifc(:,nlevp1,jb),      & !in
        &                tkvm=prm_diag%tkvm(:,2:nlev,jb),                              & !in
        &                tkvh=prm_diag%tkvh(:,2:nlev,jb),                              & !in
        &                t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),         & !in
        &                h_ice=wtr_prog_now%h_ice(:,jb),                               & !in
        &                tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),               & !in
        &                ie=nproma, ke=nlev, ke1=nlevp1,                               & !in
        &                i_startidx=i_startidx, i_endidx=i_endidx, dt=tcall_turb_jg,   & !in
        &                du_turb=prm_nwp_tend%ddt_u_turb(:,:,jb),                      & !out
        &                dv_turb=prm_nwp_tend%ddt_v_turb(:,:,jb),                      & !out
        &                dt_turb=prm_nwp_tend%ddt_temp_turb(:,:,jb),                   & !out
        &                dqv_turb=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),            & !out
        &                dqc_turb=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc))              !out
!       &                shfl_s=prm_diag%shfl_s_t(:,jb,1),                             & !out
!       &                lhfl_s=prm_diag%lhfl_s_t(:,jb,1),                             & !out
!       &                qhfl_s=prm_diag%qhfl_s_t(:,jb,1),                             & !out
!       &                umfl_s=prm_diag%umfl_s(:,jb), vmfl_s=prm_diag%vmfl_s(:,jb)    ) !out

!DR NOTE: computation of sensible and latent heat fluxes (over non-land points) should be 
!DR moved either to the turbtran interface, or the surface interface!!

!     DO jc = i_startidx, i_endidx
!       prm_diag%shfl_s(jc,jb) = prm_diag%shfl_s_t(jc,jb,1)
!       prm_diag%lhfl_s(jc,jb) = prm_diag%lhfl_s_t(jc,jb,1)
!       prm_diag%qhfl_s(jc,jb) = prm_diag%qhfl_s_t(jc,jb,1)
!     ENDDO

    ENDIF !inwp_turb



    ! Update wind speed, QV and temperature with turbulence tendencies
    ! Note: wind speed is updated here, in order to pass u and v at the correct time level
    ! to turbtran and the convection scheme. However, the prognostic variable vn is updated
    ! at the end of the NWP interface by first interpolating the u/v tendencies to the 
    ! velocity points (in order to minimize interpolation errors) and then adding the tendencies
    ! to vn
    DO jk = 1, nlev
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tracer(jc,jk,jb,iqv) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqv) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv))
        p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
         &  + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        p_diag%u(jc,jk,jb) = p_diag%u(jc,jk,jb) + tcall_turb_jg*prm_nwp_tend%ddt_u_turb(jc,jk,jb)
        p_diag%v(jc,jk,jb) = p_diag%v(jc,jk,jb) + tcall_turb_jg*prm_nwp_tend%ddt_v_turb(jc,jk,jb)
      ENDDO
    ENDDO
    ! QC is updated only in that part of the model domain where moisture physics is active
    DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tracer(jc,jk,jb,iqc) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
      ENDDO
    ENDDO

    IF (turbdiff_config(jg)%ldiff_qi) THEN
      ! QI is updated only in that part of the model domain where moisture physics is active
      DO jk = kstart_moist(jg), nlev
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqi) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqi) &
               &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqi))
        ENDDO
      ENDDO
    ENDIF  ! ldiff_qi

    ! VN is updated in nwp_nh_interface (for efficiency reasons)

  ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE nwp_turbdiff

END MODULE mo_nwp_turbdiff_interface
