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
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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

MODULE mo_nwp_turbdiff_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,        ONLY: phy_params 
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_data_turbdiff,        ONLY: get_turbdiff_param
  USE src_turbdiff,            ONLY: turbdiff
  USE mo_gme_turbdiff,         ONLY: partura, progimp_turb
  USE mo_sgs_turbulence,       ONLY: drive_subgrid_diffusion

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_turbdiff

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbdiff  ( tcall_turb_jg,                     & !>in
                          & p_patch,                          & !>in
                          & p_metrics,                        & !>in
                          & p_int,                            & !>in
                          & ext_data,                         & !>in
                          & p_prog,                           & !>in
                          & p_prog_now_rcf,                   & !>in
                          & p_prog_rcf,                       & !>inout
                          & p_diag ,                          & !>inout
                          & prm_diag, prm_nwp_tend,           & !>inout
                          & wtr_prog_now,                     & !>in 
                          & lnd_prog_now,                     & !>in 
                          & lnd_prog_new,                     & !>inout only for idealized LES
                          & lnd_diag                          ) !>in


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_int_state),    INTENT(in),TARGET   :: p_int          !< single interpolation state
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
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new    !< prog vars for sfc
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

  INTEGER  :: nlev, nlevp1                          !< number of full and half levels

  REAL(wp) :: z_tvs(nproma,p_patch%nlevp1,1)        !< aux turbulence velocity scale [m/s]



!--------------------------------------------------------------


  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1


  IF (msg_level >= 15) CALL message('mo_nwp_turbdiff:', 'turbulence')
    
  ! local variables related to the blocking
  
  i_nchdom  = MAX(1,p_patch%n_childdom)
  jg        = p_patch%id
  
  ! exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

  
  IF ( atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN
     CALL get_turbdiff_param(jg)
  ENDIF

 
  !For 3D turbulence the whole patch needs to be passed. Therefore, this call
  !is made outside the block loop next. However, the tendencies it calculates
  !is then used inside the block loop (see at the end) to update u,v,t,qv,qc
  IF ( atm_phy_nwp_config(jg)%inwp_turb == 5 )THEN
    CALL message('mo_nwp_turbdiff:', '3D turbulence')
    CALL drive_subgrid_diffusion(p_prog,       & !inout for w (it is updated inside)
                                 p_prog_rcf,   & !in
                                 p_diag,       & !inout
                                 p_metrics,    & !in
                                 p_patch,      & !in
                                 p_int,        & !in
                                 lnd_prog_now, & !in
                                 lnd_prog_new, & !inout only for idealized cases
                                 lnd_diag,     & !inout
                                 prm_diag,     & !inout
                                 prm_nwp_tend, & !inout
                                 tcall_turb_jg & !in
                                 )
  END IF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,ierrstat,errormsg,eroutine,z_tvs)  &
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


    IF ( atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN

!-------------------------------------------------------------------------
!< COSMO turbulence scheme by M. Raschendorfer  
!-------------------------------------------------------------------------

      ierrstat = 0

      !KF tendencies  have to be set to zero
      !GZ: this should be replaced by an appropriate switch in turbdiff
      prm_nwp_tend%ddt_u_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_v_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_temp_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) = 0._wp
      

      ! note that TKE must be converted to the turbulence velocity scale SQRT(2*TKE)
      ! for turbdiff
      ! INPUT to turbdiff is timestep now
      z_tvs(i_startidx:i_endidx,:,1) =  &
        &           SQRT(2._wp * p_prog_now_rcf%tke(i_startidx:i_endidx,:,jb))


      ! GZ, 2013-02-22: always use concentration lower boundary condition for momentum 
      ! (corresponding to imode_turb = 2) because the flux condition suffers from stability problems
      CALL turbdiff(iini=0, lstfnct=.TRUE.,                                             & !in
         &  dt_var=tcall_turb_jg, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1,         & !in
         &  ie=nproma, ke=nlev, ke1=nlevp1,  kcm=nlevp1,                                & !in
         &  istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx,   & !in
         &  l_hori=phy_params(jg)%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),            & !in
         &  dp0=p_diag%dpres_mc(:,:,jb),                                                & !in
         &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb),   & !in
         &  h_ice=wtr_prog_now%h_ice (:,jb), ps=p_diag%pres_sfc(:,jb),                  & !in
         &  t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),                       & !in
         &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb),                 & !in
         &  T=p_diag%temp(:,:,jb),                                                      & !in
         &  qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),         & !in
         &  prs=p_diag%pres(:,:,jb), rho=p_prog%rho(:,:,jb), epr=p_prog%exner(:,:,jb),  & !in
         &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),     & !inout
         &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb),     & !inout
         &  tke=z_tvs (:,:,:),                                                          & !inout
         &  tkvm=prm_diag%tkvm(:,2:nlevp1,jb), tkvh=prm_diag%tkvh(:,2:nlevp1,jb),       & !inout
         &  rcld=prm_diag%rcld(:,:,jb),                                                 & !inout
         &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb),                                     & !inout
         &  v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb),                                     & !inout
         &  t_tens=prm_nwp_tend%ddt_temp_turb(:,:,jb),                                  & !inout
         &  qv_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),                           & !inout
         &  qc_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc),                           & !inout
         &  tketens=prm_nwp_tend%ddt_tke(:,:,jb),                                       & !inout
         &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb),                                      & !in
         &  vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) ,                                     & !in
         &  shfl_s=prm_diag%shfl_s(:,jb), qhfl_s=prm_diag%qhfl_s(:,jb),                 & !in
         &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine                     ) !inout
          
      IF (ierrstat.NE.0) THEN
        CALL finish(eroutine, errormsg)
      END IF

      ! transform updated turbulent velocity scale back to TKE
      ! Note: ddt_tke is purely diagnostic and has already been added to z_tvs
      p_prog_rcf%tke(i_startidx:i_endidx,:,jb)= 0.5_wp                            &
        &                                     * (z_tvs(i_startidx:i_endidx,:,1))**2


    ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 2 ) THEN

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
    ! Note: the update of wind speed is done here in order to pass u and v at the correct time level
    ! to turbtran and the convection scheme. However, the update of the prognostic variable vn
    ! is done at the end of the NWP interface by first interpolating the u/v tendencies to the 
    ! velocity points (in order to minimize interpolation errors) and then adding the tendencies
    ! to vn
    DO jk = 1, nlev
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
      DO jc = i_startidx, i_endidx
        p_prog_rcf%tracer(jc,jk,jb,iqc) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
             &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
      ENDDO
    ENDDO
    ! VN is updated in nwp_nh_interface (for efficiency reasons)


  ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE nwp_turbdiff

END MODULE mo_nwp_turbdiff_interface
