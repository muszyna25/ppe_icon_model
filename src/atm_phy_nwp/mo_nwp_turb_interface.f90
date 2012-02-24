!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! turbulence parameterisations:
!! inwp_turb == 1 == turbulence scheme by M. Raschendorfer run in COSMO
!! inwp_turb == 2 == turbulence scheme by Brinkop and Roeckner run in ECHAM
!! inwp_turb == 3 == EDMF DUAL turbulence scheme by Koehler and Neggers from IFS
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
MODULE mo_nwp_turb_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish

  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_physical_constants,   ONLY: alv, rd_o_cpd

  USE mo_ext_data,             ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag,&
    &                                t_nh_metrics
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend, phy_params 
  USE mo_nwp_lnd_state,        ONLY: t_lnd_prog, t_lnd_diag

  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, &
    &                                iqi, iqr, iqs, nqtendphy
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_data_turbdiff,        ONLY: get_turbdiff_param
  USE src_turbdiff,            ONLY: organize_turbdiff
  USE mo_satad,                ONLY: sat_pres_water, spec_humi  
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_vdiff_config,         ONLY: vdiff_config
  USE mo_vdiff_driver,         ONLY: vdiff
  USE mo_advection_config,     ONLY: advection_config
  USE mo_vdfouter,             ONLY: vdfouter
  USE mo_run_config,           ONLY: ltestcase
  USE mo_nh_testcases,         ONLY: nh_test_name
  USE mo_nh_wk_exp,            ONLY: qv_max_wk

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_turbulence

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbulence ( tcall_turb_jg,                     & !>input
                          & p_patch,p_metrics,                 & !>input
                          & ext_data,                          & !>input
                          & p_prog,                            & !>inout
                          & p_prog_now_rcf, p_prog_rcf,        & !>in/inout
                          & p_diag ,                           & !>inout
                          & prm_diag, prm_nwp_tend,            & !>inout 
                          & lnd_prog_now, lnd_diag             ) !>inout


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
  TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog),      TARGET,INTENT(IN)   :: p_prog_now_rcf  !<progs with red.
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_nwp_phy_tend),TARGET, INTENT(inout):: prm_nwp_tend    !< atm tend vars
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for 
                                                               !< turbulence

  ! Local array bounds:

  INTEGER :: nblks_c, nblks_e        !> number of blocks for cells / edges
  INTEGER :: npromz_e, npromz_c      !> length of last block line
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: i_nchdom                !< domain index

  ! Local scalars:

  INTEGER :: jc,jk,jb,jt,jg      !block indeces

  ! local variables for turbdiff
  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''
  REAL(wp) ::                  &     !< aux turbulence velocity scale [m/s]
    &  z_tvs    (nproma,p_patch%nlevp1,p_patch%nblks_c,1) !DR, &
!DR    &  z_ddt_tke(nproma,p_patch%nlevp1,p_patch%nblks_c) ! tvs tendency [m/s2]

  ! local variables for vdiff
  INTEGER, PARAMETER :: itrac = 1
  REAL(wp) ::  &                     !< cloud water + cloud ice 
    &  z_plitot(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< fraction of land,seaice, open water in the grid box
    &  zfrc(nproma,nsfc_type,p_patch%nblks_c)         
  REAL(wp) ::  &                     !< dummy variable for input
    &  zdummy_tsfc(nproma,1:nsfc_type,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for input
    &  zdummy_qvsfc(nproma,1:nsfc_type,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    & z_dummy_shflx(nproma,1:nsfc_type,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    & z_dummy_lhflx(nproma,1:nsfc_type,p_patch%nblks_c)
  !
  REAL(wp) ::  &                     !< dummy variable for input
    &  zdummy_i(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for input
    &  zdummy_it(nproma,p_patch%nlev,itrac,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for input
    &  zdummy_ith(nproma,itrac,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o1(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o2(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o3(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o4(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o5(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o6(nproma,p_patch%nlev,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o7(nproma,p_patch%nlev,p_patch%nblks_c) 
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_o8(nproma,p_patch%nlev,p_patch%nblks_c) 
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_ot3(nproma,p_patch%nlev,itrac,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_ot2(nproma,p_patch%nlev,itrac,p_patch%nblks_c)
  REAL(wp) ::  &                     !< dummy variable for output
    &  zdummy_oh(nproma,p_patch%nblks_c) 
  INTEGER  :: idummy_oh(nproma ,p_patch%nblks_c)    !< dummy variable for output
  INTEGER  :: nlev, nlevp1           !< number of full and half levels

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1
  
  IF (msg_level >= 12) &
        & CALL message('mo_nwp_turb:', 'turbulence')
  
  IF (msg_level >= 15) THEN
    WRITE(message_text,'(a,3E15.7)') ' bottom TKE before turbulence now = ', &
         &  p_prog_now_rcf%tke(1,nlev+1,6), p_prog_now_rcf%tke(1,nlev,6),&
         &  p_prog_now_rcf%tke(1,nlev-1,6)
    CALL message('', TRIM(message_text))
    WRITE(message_text,'(a,3E15.7)') ' bottom TKE before turbulence new = ', &
         &  p_prog_rcf%tke(1,nlev+1,6), p_prog_rcf%tke(1,nlev,6),&
        &  p_prog_rcf%tke(1,nlev-1,6)
    CALL message('', TRIM(message_text))
  ENDIF
  
  ! local variables related to the blocking
  
  nblks_c   = p_patch%nblks_int_c
  npromz_c  = p_patch%npromz_int_c
  nblks_e   = p_patch%nblks_int_e
  npromz_e  = p_patch%npromz_int_e
  i_nchdom  = MAX(1,p_patch%n_childdom)
  jg        = p_patch%id
  
  !in order to account for mesh refinement
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_blk(rl_start,1)
  i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

  
  IF ( atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN
     CALL get_turbdiff_param(jg)
  ENDIF

!$OMP PARALLEL

#ifdef __xlC__
!$OMP DO PRIVATE(jb,jt,jc,jk,i_startidx,i_endidx,ierrstat,errormsg,eroutine)
#else
!$OMP DO PRIVATE(jb,jt,jc,jk,i_startidx,i_endidx,ierrstat,errormsg,eroutine), SCHEDULE(guided)
#endif
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

   !-------------------------------------------------------------------------
   !<  turbulent transfer and diffusion
   !-------------------------------------------------------------------------

   !<  NOTE: since  turbulence is a fast process it is
   !!        allowed to do a sequential updating except for wind speed
   !!        because back-and-forth interpolation would cause too large errors
   !!  (GZ, 2011-08-29): Nevertheless, tendency fields are now passed to turbdiff
   !!        to have them available for extended diagnostic output


    IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
      ! check dry case
      IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
        lnd_diag%qv_s (:,:) = 0._wp
      ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 1) THEN
        IF ( ltestcase .AND. nh_test_name == 'wk82') THEN   
         DO jc = i_startidx, i_endidx
          lnd_prog_now%t_g(jc,jb) = p_diag%temp(jc,nlev,jb)*  &
                      ((p_diag%pres_sfc(jc,jb))/p_diag%pres(jc,nlev,jb))**rd_o_cpd
          lnd_diag%qv_s (jc,jb) = &
             &         spec_humi(sat_pres_water(lnd_prog_now%t_g(jc,jb)),&
             &                                   p_diag%pres_sfc(jc,jb) )          
          lnd_diag%qv_s(jc,jb) = MIN (lnd_diag%qv_s(jc,jb) ,p_prog_rcf%tracer(jc,nlev,jb,iqv))
         END DO
        ELSE
         !
         !> adjust  humidity at water surface because of changed surface pressure
         !
         DO jc = i_startidx, i_endidx
          lnd_diag%qv_s (jc,jb) = &
             &         spec_humi(sat_pres_water(lnd_prog_now%t_g(jc,jb)),&
             &                                   p_diag%pres_sfc(jc,jb) )
         ENDDO
        END IF
      ENDIF
    ENDIF

    IF ( atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN

      ierrstat = 0

      !KF tendencies  have to be set to zero
      ! GZ: this should be replaced by an appropriate switch in turbdiff
      prm_nwp_tend%ddt_u_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_v_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_temp_turb(:,:,jb) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) = 0._wp
      

      ! note that TKE must be converted to the turbulence velocity scale SQRT(2*TKE)
      ! for turbdiff
      ! INPUT to turbdiff is timestep now
!      ! Artificial limitation of input TKE to 60 m^2/s^2
!      z_tvs(i_startidx:i_endidx,:,jb,1)=  &
!        &     SQRT(2._wp * MIN(60._wp,p_prog_now_rcf%tke(i_startidx:i_endidx,:,jb)))
      z_tvs(i_startidx:i_endidx,:,jb,1)=  &
        &           SQRT(2._wp * p_prog_now_rcf%tke(i_startidx:i_endidx,:,jb))


!-------------------------------------------------------------------------
!< COSMO version by M. Raschendorfer  
!-------------------------------------------------------------------------

      CALL organize_turbdiff(action='tran_diff', iini=0, lstfnct=.TRUE., &
!
         &  dt_var=tcall_turb_jg, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1, &
!
         &  ie=nproma, je=1, ke=nlev, ke1=nlevp1,  kcm=nlevp1,  vst=0, &
!       
         &  istart   =i_startidx, iend   =i_endidx, istartu=i_startidx, iendu=i_endidx, &
         &  istartpar=i_startidx, iendpar=i_endidx, istartv=i_startidx, iendv=i_endidx, &
!       
         &  jstart   =1,          jend   =1       , jstartu=1         , jendu=1       , &
         &  jstartpar=1         , jendpar=1       , jstartv=1         , jendv=1       , &
!       
         &  l_hori=phy_params(jg)%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),            &
         &  dp0=p_diag%dpres_mc(:,:,jb),                                                &
!
         &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
         &  sai=prm_diag%sai(:,jb), h_ice=prm_diag%h_ice (:,jb), &
!         
         &  ps=p_diag%pres_sfc(:,jb), t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb), &
!           
         &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb), T=p_diag%temp(:,:,jb), &
         &  qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc), &
!         
         &  prs=p_diag%pres(:,:,jb), rho=p_prog%rho(:,:,jb), epr=p_prog%exner(:,:,jb), &
!         
         &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
         &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
!                
         &  tke=z_tvs (:,:,jb,:) ,&!  edr =prm_diag%edr(:,:,jb),                    &
         &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh(:,:,jb), rcld=prm_diag%rcld(:,:,jb), &
!       
         &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb), &
         &  t_tens=prm_nwp_tend%ddt_temp_turb(:,:,jb), &
         &  qv_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),&
         &  qc_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc), &
         &  tketens=prm_nwp_tend%ddt_tke(:,:,jb), &
         &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb), vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) ,&
!         
         &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), td_2m=prm_diag%td_2m(:,jb), &
         &  rh_2m=prm_diag%rh_2m(:,jb), u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb), &
         &  shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb), &
!         
         &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )
          
      IF (ierrstat.NE.0) THEN
        CALL finish(eroutine, errormsg)
      END IF

      ! transform updated turbulent velocity scale back to TKE
      p_prog_rcf%tke(i_startidx:i_endidx,:,jb)= 0.5_wp                                &
        &                                     * (z_tvs(i_startidx:i_endidx,:,jb,1))**2

       ! Update QV, QC and temperature with turbulence tendencies
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqv) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqv) &
               &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv))
          p_prog_rcf%tracer(jc,jk,jb,iqc) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
               &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
          ! Artificial limitation of turbulent temperature tendency to 0.01 K/s
          p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
            &  + tcall_turb_jg*MIN(0.01_wp,MAX(-0.01_wp,prm_nwp_tend%ddt_temp_turb(jc,jk,jb)))
!          p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
!           &  + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        ENDDO
      ENDDO

    ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 2 ) THEN

!-------------------------------------------------------------------------
!> ECHAM version 
!-------------------------------------------------------------------------

      ! GZ: setting 1 instead of i_startidx in the following assignments is needed
      ! as a workaround for the missing istart-parameter in vdiff

      DO jk = 1, nlev
        DO jc =  1, i_endidx
          z_plitot (jc,jk,jb) = p_prog_rcf%tracer(jc,jk,jb,iqc) &
&                             + p_prog_rcf%tracer(jc,jk,jb,iqi) &
&                             + p_prog_rcf%tracer(jc,jk,jb,iqr) &
&                             + p_prog_rcf%tracer(jc,jk,jb,iqs) 
        ENDDO
      ENDDO

      ! initialize dummy fields with zero
      zdummy_i  (:,:,jb)   = 0.0_wp
      zdummy_it (:,:,jb,:) = 0.0_wp
      zdummy_ith(:,:,jb)   = 0.0_wp
      zdummy_ot3(:,:,jb,:) = 0.0_wp
      z_dummy_shflx(:,:,jb)   = 0.0_wp
      z_dummy_shflx(:,:,jb)   = 0.0_wp

      ! Merge three pieces of information into one array for vdiff

      ! fraction of land in the grid box. lsmask: land-sea mask, 1.= land
      IF (ilnd <= nsfc_type) &
           &  zfrc(1:i_endidx,ilnd,jb) = 0._wp !&
!           & REAL(ext_data%atm%lsm_atm_c(i_startidx:i_endidx,jb),wp)
      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box
      IF (iwtr <= nsfc_type) &
        &  zfrc(1:i_endidx,iwtr,jb) =  1._wp
!        &        (1._wp- REAL(ext_data%atm%lsm_atm_c(i_startidx:i_endidx,jb),wp))&
!        &       *(1._wp-          lnd_diag%fr_seaice(i_startidx:i_endidx,jb))
      ! fraction of sea ice in the grid box
      IF (iice <= nsfc_type) &
           &  zfrc(1:i_endidx,iice,jb)= 0._wp !&
!           &                  lnd_diag%fr_seaice(i_startidx:i_endidx,jb)
      
      !KF tendencies in vdiff are INOUT declared, so they have to be set to zero
      prm_nwp_tend%ddt_temp_turb  (1:i_endidx,:,jb)         = 0._wp
      prm_nwp_tend%ddt_tracer_turb(1:i_endidx,:,jb,iqv:iqi) = 0._wp
      prm_nwp_tend%ddt_u_turb     (1:i_endidx,:,jb)         = 0._wp
      prm_nwp_tend%ddt_v_turb     (1:i_endidx,:,jb)         = 0._wp
!      p_prog_rcf%tke              (1:i_endidx,:,jb)         = 0._wp

      ! KF as long as if nsfc_type = 1 !!!!!
      zdummy_tsfc (1:i_endidx,nsfc_type,jb) = &
           &                              lnd_prog_now%t_g (1:i_endidx,jb)
      !zdummy_qvsfc(1:i_endidx,nsfc_type,jb) = &
      !     &                                 lnd_diag%qv_s (1:i_endidx,jb)

      ! Workarounds needed because vdiff is not coded properly for use with nesting
      DO jk = 1, nlev
        p_diag%u(1:i_startidx-1,jk,jb) = 0._wp
        p_diag%v(1:i_startidx-1,jk,jb) = 0._wp
      ENDDO
         

      CALL vdiff( lsfc_mom_flux  = vdiff_config%lsfc_mom_flux                        ,&! in
            &     lsfc_heat_flux = vdiff_config%lsfc_heat_flux                       ,&! in
            & kproma = i_endidx, kbdim   = nproma                                    ,&! in
            & klev   = nlev,   klevm1    = nlev-1,  klevp1=nlevp1                    ,&! in
            & ktrac  = itrac,  ksfc_type = nsfc_type                                 ,&! in
            & idx_wtr= iwtr,   idx_ice   = iice,   idx_lnd =ilnd                     ,&! in 
            & pdtime = tcall_turb_jg,       pstep_len = tcall_turb_jg                ,&! in
            !
            & pfrc      = zfrc(:,:,jb)                                               ,&! in
            !& pqsat_sfc = zdummy_qvsfc  (:,:,jb)                                     ,&! in
            & pocu      = prm_diag%ocu    (:,jb),  pocv   = prm_diag%ocv      (:,jb) ,&! in
            & ppsfc     = p_diag%pres_sfc(:,jb),  pcoriol= p_patch%cells%f_c(:,jb)   ,&! in
            !
            & pum1      = p_diag% u  (:,:,jb),   pvm1 = p_diag% v        (:,:,jb)    ,&! in
            & ptm1      = p_diag%temp(:,:,jb),   pqm1 = p_prog_rcf%tracer(:,:,jb,iqv),&! in
            & pxlm1     = p_prog_rcf%tracer(:,:,jb,iqc)                              ,&! in
            & pxim1     = p_prog_rcf%tracer(:,:,jb,iqi)                              ,&! in
            & pxm1      = z_plitot (:,:,jb) ,    pxtm1 = zdummy_it        (:,:,:,jb) ,&! in
            !
            & paphm1 = p_diag%pres_ifc(:,:,jb),  papm1 = p_diag%pres        (:,:,jb) ,&! in
            & pdelpm1= p_diag%dpres_mc(:,:,jb),  pgeom1= p_metrics%geopot_agl(:,:,jb),&! in 
            & ptvm1  = p_diag%tempv   (:,:,jb),  paclc = prm_diag%tot_cld(:,:,jb,icc),&! in
            & ptkem1 = p_prog_now_rcf%tke (:,2:nlevp1,jb)                            ,&! in
            !
            & pxt_emis= zdummy_ith     (:,:,jb)                                      ,&! in
            & ptsfc_tile = zdummy_tsfc(:,:,jb)                                       ,&! inout
            & pxvar   = zdummy_i     (:,:,jb)                                        ,&! inout
            & pthvvar = prm_diag%thvvar(:,2:nlevp1,jb), pustar = prm_diag%ustar(:,jb),&! inout
            & pz0m_tile = prm_diag%z0m_tile(:,jb,:)                                  ,&! inout
            & pkedisp  = prm_diag%kedisp(:,jb)                                       ,&! inout
            !
            & pute    = prm_nwp_tend%ddt_u_turb     (:,:,jb)                         ,&! inout
            & pvte    = prm_nwp_tend%ddt_v_turb     (:,:,jb)                         ,&! inout
            & ptte    = prm_nwp_tend%ddt_temp_turb  (:,:,jb)                         ,&! inout
            & pqte    = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv)                     ,&! inout
            & pxlte   = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc)                     ,&! inout
            & pxite   = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi)                     ,&! inout
            & pxtte    = zdummy_ot3(:,:,:,jb)                                        ,&! inout
            !
            & pz0m    = prm_diag%z0m(:,jb)                                           ,&! out
            & pute_vdf = zdummy_o1 (:,:,jb)                                          ,&! out
            & pvte_vdf = zdummy_o2 (:,:,jb),        ptte_vdf = zdummy_o3 (:,:,jb)    ,&! out
            & pqte_vdf = zdummy_o4 (:,:,jb),        pxlte_vdf= zdummy_o5 (:,:,jb)    ,&! out
            & pxite_vdf= zdummy_o6 (:,:,jb),        pxtte_vdf= zdummy_ot2(:,:,:,jb)  ,&! out
            !
            & pqsat_tile = zdummy_qvsfc  (:,:,jb)                                    ,&! out
            & pshflx_tile = z_dummy_shflx(:,:,jb)                                    ,&! out
            & plhflx_tile = z_dummy_lhflx(:,:,jb)                                    ,&! out
            & pxvarprod= zdummy_o7(:,:,jb),         pvmixtau = zdummy_o8 (:,:,jb)    ,&! out
            & pqv_mflux_sfc=prm_diag%qhfl_s (:,jb), pthvsig  = zdummy_oh (:,jb)      ,&! out
            & ptke     = p_prog_rcf%tke (:,2:nlevp1,jb), ihpbl = idummy_oh (:,jb)    ,&! inout
            & pghpbl   = prm_diag%ghpbl (:,jb), pri =  prm_diag%ri (:,2:nlevp1,jb)   ,&! out
            & pmixlen  = prm_diag%mixlen (:,2:nlevp1,jb)                             ,&! out
            & pcfm     = prm_diag%cfm    (:,2:nlevp1,jb)                             ,&! out
            & pcfh     = prm_diag%cfh    (:,2:nlevp1,jb)                             ,&! out
            & pcfv     = prm_diag%cfv    (:,2:nlevp1,jb)                             ,&! out
            & pcfm_tile= prm_diag%cfm_tile(:,jb,:)                                   ,&! out
            & pcfh_tile= prm_diag%cfh_tile(:,jb,:)                                   ,&! out
            & pcftke   = prm_diag%cftke  (:,2:nlevp1,jb)                             ,&! out
            & pcfthv   = prm_diag%cfthv  (:,2:nlevp1,jb))                              ! out


      !-------------------------------------------------------------------------
      !> in case of ECHAM version  update temperature and moist fields
      !-------------------------------------------------------------------------

      DO jt=1,nqtendphy ! iqv,iqc,iqi
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            p_prog_rcf%tracer(jc,jk,jb,jt) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,jt)  &
                 &             + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,jt))
          ENDDO
        ENDDO
      ENDDO

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
          &  + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        ENDDO
      ENDDO
      ! In case nsfc_type >1 , the grid-box mean should be considered instead (PR)
      IF (nsfc_type == 1) THEN
       lnd_diag%qv_s(1:i_endidx,jb)=zdummy_qvsfc(1:i_endidx,nsfc_type,jb)
       !prm_diag%lhfl_s(1:i_endidx,jb)=prm_diag%qhfl_s (1:i_endidx,jb)*alv
       prm_diag%lhfl_s(1:i_endidx,jb)=z_dummy_lhflx(1:i_endidx,nsfc_type,jb)
       prm_diag%shfl_s(1:i_endidx,jb)=z_dummy_shflx(1:i_endidx,nsfc_type,jb)
      END IF

    ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 3 ) THEN

!-------------------------------------------------------------------------
!> EDMF DUALM turbulence scheme (eddy-diffusivity/mass-flux dual mass-flux)
!-------------------------------------------------------------------------

!      !-------------------------------------------------------------------------
!      !> Calculate vertical velocity in p-system
!      !-------------------------------------------------------------------------
!      !
!      DO jk = 1,nlev
!        DO jc = i_startidx,i_endidx
!          z_omega_p(jc,jk,jb)= -p_prog%w(jc,jk,jb)*p_prog%rho(jc,jk,jb)*grav
!        ENDDO
!      ENDDO
!
!      icnt = 0
!      CALL vdfouter ( &
!        & CDCONF ???                                                  ,&! (IN)    
!        & KIDIA  = i_startidx                                         ,&! (IN)   
!        & KFDIA  = i_endidx                                           ,&! (IN)   
!        & KLON   = nproma                                             ,&! (IN)   
!        & KLEV   = nlev                                               ,&! (IN)   
!        & KLEVS  = ???                                                ,&! (IN)   
!        & KSTEP  = 1 ???                                              ,&! (IN)   
!        & KTILES = ???                                                ,&! (IN)   
!        & KTRAC  = 0  ???itrac ???                                    ,&! (IN)   
!        & KLEVSN = ???                                                ,&! (IN)   
!        & KLEVI  = 0 ???                                              ,&! (IN)   
!        & KDHVTLS =                                                   ,&! (IN)   
!        & KDHFTLS =                                                   ,&! (IN)   
!        & KDHVTSS =                                                   ,&! (IN)   
!        & KDHFTSS =                                                   ,&! (IN)   
!        & KDHVTTS =                                                   ,&! (IN)   
!        & KDHFTTS =                                                   ,&! (IN)   
!        & KDHVTIS = 0                                                 ,&! (IN)   
!        & KDHFTIS = 0                                                 ,&! (IN)   
!        & PTSPHY  = tcall_turb_jg                                     ,&! (IN)   
!        & KTVL(KLON) = ??                                             ,&! (IN)   
!        & KTVH(KLON) = ??                                             ,&! (IN)   
!        & KCNT = icnt                                                 ,&! (INOUT)
!        & PCVL(KLON) =  ??                                            ,&! (IN)   
!        & PCVH(KLON) =    ??                                          ,&! (IN)   
!        & PSIGFLT(KLON) = ext_data%atm%sso_stdh(:,jb)    (needs to be passed down)   ,&! (IN)   
!        & PUM1(KLON,KLEV) = p_diag%u(:,:,jb)                          ,&! (IN)   
!        & PVM1(KLON,KLEV) = p_diag%v(:,:,jb)                          ,&! (IN)   
!        & PTM1(KLON,KLEV) = p_diag%temp(:,:,jb)                       ,&! (IN)   
!        & PQM1(KLON,KLEV) = p_prog_rcf%tracer(:,:,jb,iqv)             ,&! (IN)   
!        & PLM1(KLON,KLEV) = p_prog_rcf%tracer(:,:,jb,iqc)             ,&! (IN)   
!        & PIM1(KLON,KLEV) = p_prog_rcf%tracer(:,:,jb,iqi) pass it down (iqi)         ,&! (IN)   
!        & PAM1(KLON,KLEV) = prm_diag%tot_cld(:,:,jb,icc) pass it down (icc)          ,&! (IN)   
!        & PCM1(KLON,KLEV,KTRAC) = dummy                               ,&! (IN)   
!        & PAPHM1(KLON,0:KLEV) = p_diag%pres_ifc(:,:,jb)               ,&! (IN)   
!        & PAPM1(KLON,KLEV)    = p_diag%pres(:,:,jb)                   ,&! (IN)   
!        & PGEOM1(KLON,KLEV)   = p_metrics%geopot_agl(:,:,jb)          ,&! (IN)   
!        & PGEOH(KLON,0:KLEV)  = p_metrics%geopot_agl_ifc(:,:,jb)      ,&! (IN)   
!        & PTSKM1M(KLON)       = ???                                   ,&! (IN)   
!        & PTSAM1M(KLON,KLEVS) = ???                                   ,&! (IN) in???  
!        & PWSAM1M(KLON,KLEVS) = ???                                   ,&! (IN) in???  
!        & PSSRFL(KLON)  = prm_diag%swflxsfc (:,jb)                    ,&! (IN)   
!        & PSLRFL(KLON)  = prm_diag%lwflxsfc (:,jb)                    ,&! (IN)   
!        & PEMIS(KLON)   = ext_data%atm%emis_rad(:,jb)                 ,&! (IN)   
!        & PHRLW(KLON,KLEV) = prm_nwp_tend%ddt_temp_radlw(:,:,jb),     ,&! (IN)   
!        & PHRSW(KLON,KLEV) = prm_nwp_tend%ddt_temp_radsw(:,:,jb)      ,&! (IN)   
!        & PTSNOW(KLON) = lnd_prog_now%t_snow(:,jb,isubs)              ,&! (IN)?? 
!        & PTICE(KLON)  = ???                                          ,&! (IN)   
!        & PHLICE(KLON) = ???                                          ,&! (IN)   
!        & PTLICE(KLON) = ???                                          ,&! (IN)   
!        & PTLWML(KLON) = ???                                          ,&! (IN)   
!        & PSST(KLON)   = lnd_prog_now%t_g(:,jb)                       ,&! (IN)   
!        & KSOTY(KLON)   ???                                           ,&! (IN)   
!        & PFRTI(KLON,KTILES) = 1 :)  ???                              ,&! (IN)   
!        & PALBTI(KLON,KTILES) = ...not active...                      ,&! (IN)   
!        & PWLMX(KLON) =         ...not active...                      ,&! (IN)   
!        & PCHAR(KLON) = const = 0.018 (no wave model)                 ,&! (IN)   
!        & PUCURR(KLON) = 0                                            ,&! (IN)   
!        & PVCURR(KLON) = 0                                            ,&! (IN)   
!        & PTSKRAD(KLON) = prm_diag%tsfctrad(:,jb)  ???                ,&! (IN)   
!        & PCFLX(KLON,KTRAC)  = 0                                      ,&! (IN)   
!        & PSOTEU(KLON,KLEV)  = 0                                      ,&! (IN)   
!        & PSOTEV(KLON,KLEV)  = 0                                      ,&! (IN)   
!        & PSOBETA(KLON,KLEV) = 0                                      ,&! (IN)   
!        & PVERVEL(KLON,KLEV) = z_omega_p(:,:,jb)                      ,&! (IN)   
!        & PZ0M(KLON) = prm_diag%z0m(:,jb)         (reduced for TOFD)  ,&! (INOUT)
!        & PZ0H(KLON) = prm_diag%z0m(:,jb) * factor ????               ,&! (INOUT)
!        & PVDIS(KLON)         = ...                                   ,&! (OUT)  
!        & PVDISG(KLON)        = ...                                   ,&! (OUT)  
!        & PDISGW3D(KLON,KLEV) = ...                                   ,&! (OUT)  
!        & PAHFLEV(KLON) =  ...                                        ,&! (OUT)  
!        & PAHFLSB(KLON) =  ...                                        ,&! (OUT)  
!        & PFWSB(KLON)   = ...                                         ,&! (OUT)  
!        & PBIR(KLON)      = ...                                       ,&! (OUT)  
!        & PVAR(KLON,KLEV) = ... qt,variance > extra tracer ??? Daniel ,&! (OUT)  
!        & PU10M(KLON)   = prm_diag%u_10m(:,jb)                        ,&! (OUT)  
!        & PV10M(KLON)   = prm_diag%v_10m(:,jb)                        ,&! (OUT)  
!        & PT2M(KLON)    = prm_diag%t_2m(:,jb)                         ,&! (OUT)  
!        & PD2M(KLON)    = prm_diag%td_2m(:,jb)                        ,&! (OUT)  
!        & PQ2M(KLON)    = prm_diag%qv_2m(:,jb)                        ,&! (OUT)  
!        & PZINV(KLON)   = out                                         ,&! (OUT)  
!        & PBLH(KLON)    = out                                         ,&! (OUT)  
!        & KHPBLN(KLON)  = out                                         ,&! (OUT)  
!        & KVARTOP(KLON) =   ???                                       ,&! (OUT)  
!        & PSSRFLTI(KLON,KTILES) = out   (in????)                      ,&! (INOUT)
!        & PEVAPSNW(KLON) = out                                        ,&! (OUT)  
!        & PGUST(KLON)    = out                                        ,&! (OUT)  
!        & PWUAVG(KLON)   = out                                        ,&! (OUT)  
!        & LDNODECP(KLON) = ???                                        ,&! (OUT)  
!        & KPBLTYPE(KLON) = out                                        ,&! (OUT)  
!        & PLDIFF(KLON,KLEV)    = out                                  ,&! (OUT)  
!        & PFPLVL(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PFPLVN(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PFHPVL(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PFHPVN(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PEXTR2(KLON,KFLDX2)      =                                  ,&! (INOUT)
!        & KFLDX2                   =                                  ,&! (IN)   
!        & PEXTRA(KLON,KLEVX,KFLDX) =                                  ,&! (INOUT)
!        & KLEVX                    =                                  ,&! (IN)   
!        & KFLDX                    =                                  ,&! (IN)   
!        & LLDIAG                   =                                  ,&! (IN)   
!        & PTE(KLON,KLEV)  = prm_nwp_tend%ddt_temp_turb(:,:,jb         ,&! (INOUT)
!        & PQE(KLON,KLEV)  = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv)  ,&! (INOUT)
!        & PLE(KLON,KLEV)  = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc)  ,&! (INOUT)
!        & PIE(KLON,KLEV)  = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi)  ,&! (INOUT)
!        & PAE(KLON,KLEV)  = ???                                       ,&! (INOUT)
!        & PVOM(KLON,KLEV) =  prm_nwp_tend%ddt_v_turb(:,:,jb)          ,&! (INOUT)
!        & PVOL(KLON,KLEV) =  prm_nwp_tend%ddt_u_turb(:,:,jb)          ,&! (INOUT)
!        & PTENC(KLON,KLEV,KTRAC) = out                                ,&! (INOUT)
!        & PTSKE1(KLON)    = ???                                       ,&! (INOUT)
!        & PUSTRTI(KLON,KTILES) =                                      ,&! (INOUT)
!        & PVSTRTI(KLON,KTILES) =                                      ,&! (INOUT)
!        & PAHFSTI(KLON,KTILES) = prm_diag%shfl_s(:,jb)  (tile mean)   ,&! (INOUT)
!        & PEVAPTI(KLON,KTILES) = prm_diag%lhfl_s(:,jb)  (W/m2!!!)     ,&! (INOUT)
!        & PTSKTI(KLON,KTILES)  = lnd_prog_new%t_g(:,jb) (now or new??),&! (INOUT)
!        & PDIFTS(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PDIFTQ(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PDIFTL(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PDIFTI(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PSTRTU(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PSTRTV(KLON,0:KLEV)  = out                                  ,&! (OUT)  
!        & PTOFDU(KLON)         = out                                  ,&! (INOUT)
!        & PTOFDV(KLON)         = out                                  ,&! (INOUT)
!        & PSTRSOU(KLON,0:KLEV) = out                                  ,&! (OUT)  
!        & PSTRSOV(KLON,0:KLEV) = out                                  ,&! (OUT)  
!        & PKH(KLON,KLEV) = prm_diag%tkvh(:,:,jb)                      ,&! (OUT)  
!        & LDLAND(KLON) = convert to logical ... ext_data%atm%fr_land(:,jb)     ,&!      (IN)  
!        & PDHTLS(KLON,KTILES,KDHVTLS+KDHFTLS) = out                   ,&! (OUT)  
!        & PDHTSS(KLON,KLEVSN,KDHVTSS+KDHFTSS) = out                   ,&! (OUT)  
!        & PDHTTS(KLON,KLEVS,KDHVTTS+KDHFTTS)  = out                   ,&! (OUT)  
!        & PDHTIS(KLON,KLEVI,KDHVTIS+KDHFTIS)  = out      )              ! (OUT)  

    ENDIF !inwp_turb

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  IF (msg_level >= 15) THEN
    WRITE(message_text,'(a,3E15.7)') ' bottom TKE after turbulence now = ', &
         &  p_prog_now_rcf%tke(1,nlev+1,6), p_prog_now_rcf%tke(1,nlev,6),&
         &  p_prog_now_rcf%tke(1,nlev-1,6)
    CALL message('', TRIM(message_text))
    WRITE(message_text,'(a,3E15.7)') ' bottom TKE after turbulence new = ', &
         &  p_prog_rcf%tke(1,nlev+1,6), p_prog_rcf%tke(1,nlev,6),&
         &  p_prog_rcf%tke(1,nlev-1,6)
    CALL message('', TRIM(message_text))
  ENDIF


END SUBROUTINE nwp_turbulence


END MODULE mo_nwp_turb_interface


