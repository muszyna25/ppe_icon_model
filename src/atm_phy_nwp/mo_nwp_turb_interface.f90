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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_turb_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_physical_constants,   ONLY: alv, rd_o_cpd, grav
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend, phy_params 
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_lnd_diag
  USE mo_phyparam_soil,        ONLY: z0_lu
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, iqr, iqs, iqtvar, nqtendphy
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
  USE mo_lnd_nwp_config,       ONLY: nlev_soil, nlev_snow, nsfc_subs

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
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for 
                                                               !< turbulence

  ! Local array bounds

  INTEGER :: nblks_c, nblks_e        !> number of blocks for cells / edges
  INTEGER :: npromz_e, npromz_c      !> length of last block line
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: i_nchdom                !< domain index

  ! Local scalars:

  INTEGER :: jc,jk,jb,jt,jg      !loop indices

  ! local variables for turbdiff

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''
  REAL(wp) ::                  &     !< aux turbulence velocity scale [m/s]
    &  z_tvs    (nproma,p_patch%nlevp1,p_patch%nblks_c,1) !DR, &
!DR    &  z_ddt_tke(nproma,p_patch%nlevp1,p_patch%nblks_c) ! tvs tendency [m/s2]
  REAL(wp) :: gz0(nproma)

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
  INTEGER  :: nlev, nlevp1                          !< number of full and half levels
  INTEGER  :: lc_class                              !< land-cover class

  ! local variables for edmf

  INTEGER  :: icnt
  INTEGER, PARAMETER :: itrac_vdf = 0
  INTEGER  :: idummy_vdf_0a(nproma), idummy_vdf_0b(nproma), idummy_vdf_0c(nproma), & 
    &         idummy_vdf_0d(nproma), idummy_vdf_0e(nproma), idummy_vdf_0f(nproma)
  REAL(wp) :: zdummy_vdf_1a(nproma), zdummy_vdf_1b(nproma), zdummy_vdf_1c(nproma), &
    &         zdummy_vdf_1d(nproma), zdummy_vdf_1e(nproma), zdummy_vdf_1f(nproma), &
    &         zdummy_vdf_1g(nproma), zdummy_vdf_1h(nproma), zdummy_vdf_1i(nproma), &
    &         zdummy_vdf_1j(nproma), zdummy_vdf_1k(nproma), zdummy_vdf_1l(nproma), &
    &         zdummy_vdf_1m(nproma), zdummy_vdf_1n(nproma), zdummy_vdf_1o(nproma), &
    &         zdummy_vdf_1p(nproma), zdummy_vdf_1q(nproma), zdummy_vdf_1r(nproma), &
    &         zdummy_vdf_1s(nproma)
  REAL(wp) :: zdummy_vdf_2a(nproma,p_patch%nlev),   zdummy_vdf_2b(nproma,p_patch%nlev)
  REAL(wp) :: zdummy_vdf_3a(nproma,p_patch%nlev+1), zdummy_vdf_3b(nproma,p_patch%nlev+1), &
    &         zdummy_vdf_3c(nproma,p_patch%nlev+1), zdummy_vdf_3d(nproma,p_patch%nlev+1), &
    &         zdummy_vdf_3e(nproma,p_patch%nlev+1), zdummy_vdf_3f(nproma,p_patch%nlev+1), &
    &         zdummy_vdf_3g(nproma,p_patch%nlev+1), zdummy_vdf_3h(nproma,p_patch%nlev+1), &
    &         zdummy_vdf_3i(nproma,p_patch%nlev+1), zdummy_vdf_3j(nproma,p_patch%nlev+1), &
    &         zdummy_vdf_3k(nproma,p_patch%nlev+1), zdummy_vdf_3l(nproma,p_patch%nlev+1)
  REAL(wp) :: zdummy_vdf_4a(nproma,nlev_soil)     , zdummy_vdf_4b(nproma,nlev_soil)
  REAL(wp) :: zdummy_vdf_5a(nproma,nsfc_subs)
  REAL(wp) :: zdummy_vdf_6a(nproma,p_patch%nlev,itrac_vdf), &
    &         zdummy_vdf_6b(nproma,p_patch%nlev,itrac_vdf), &
    &         zdummy_vdf_6c(nproma,itrac_vdf)
  REAL(wp) :: zdummy_vdf_7a(nproma,0), zdummy_vdf_7b(nproma,p_patch%nlev,0)

  REAL(wp) :: z_omega_p(nproma,p_patch%nlev), zchar(nproma),                    &
    &         zucurr(nproma)                , zvcurr(nproma),                   &
    &         zsoteu(nproma,p_patch%nlev)   , zsotev(nproma,p_patch%nlev),      &
    &         zsobeta(nproma,p_patch%nlev)  , sobs_t(nproma,nsfc_subs),         &
    &         shfl_s_t(nproma,nsfc_subs)    , &
    &         evap_s_t(nproma,nsfc_subs)    , tskin_t(nproma,nsfc_subs),        &
    &         ustr_s_t(nproma,nsfc_subs)    , vstr_s_t(nproma,nsfc_subs),       &
    &         zae(nproma,p_patch%nlev)      , zvar(nproma,p_patch%nlev),        &
    &         ztice(nproma)                 , ztske1(nproma),                   &
    &         ztskm1m(nproma)               , ztskrad(nproma),                  &
    &         zsigflt(nproma)               , zfrti(nproma,nsfc_subs)
  LOGICAL  :: l_land(nproma), ldummy_vdf_a(nproma)

  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1
  
  IF (msg_level >= 15) &
        & CALL message('mo_nwp_turb:', 'turbulence')
    
  ! local variables related to the blocking
  
  nblks_c   = p_patch%nblks_int_c
  npromz_c  = p_patch%npromz_int_c
  nblks_e   = p_patch%nblks_int_e
  npromz_e  = p_patch%npromz_int_e
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

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,jk,i_startidx,i_endidx,ierrstat,errormsg, &
!$OMP eroutine, &
!$OMP icnt, lc_class, gz0, &
!$OMP idummy_vdf_0a, idummy_vdf_0b, idummy_vdf_0c, & 
!$OMP idummy_vdf_0d, idummy_vdf_0e, idummy_vdf_0f, &
!$OMP zdummy_vdf_1a, zdummy_vdf_1b, zdummy_vdf_1c, &
!$OMP zdummy_vdf_1d, zdummy_vdf_1e, zdummy_vdf_1f, &
!$OMP zdummy_vdf_1g, zdummy_vdf_1h, zdummy_vdf_1i, &
!$OMP zdummy_vdf_1j, zdummy_vdf_1k, zdummy_vdf_1l, &
!$OMP zdummy_vdf_1m, zdummy_vdf_1n, zdummy_vdf_1o, &
!$OMP zdummy_vdf_1p, zdummy_vdf_1q, zdummy_vdf_1r, &
!$OMP zdummy_vdf_1s, &
!$OMP zdummy_vdf_2a, zdummy_vdf_2b, &
!$OMP zdummy_vdf_3a, zdummy_vdf_3b, &
!$OMP zdummy_vdf_3c, zdummy_vdf_3d, &
!$OMP zdummy_vdf_3e, zdummy_vdf_3f, &
!$OMP zdummy_vdf_3g, zdummy_vdf_3h, &
!$OMP zdummy_vdf_3i, zdummy_vdf_3j, &
!$OMP zdummy_vdf_3k, zdummy_vdf_3l, &
!$OMP zdummy_vdf_4a, zdummy_vdf_4b, &
!$OMP zdummy_vdf_5a, &
!$OMP zdummy_vdf_6a, &
!$OMP zdummy_vdf_6b, &
!$OMP zdummy_vdf_6c, &
!$OMP zdummy_vdf_7a, zdummy_vdf_7b, &
!$OMP z_omega_p, zchar, &
!$OMP zucurr,    zvcurr, &
!$OMP zsoteu,    zsotev, &
!$OMP zsobeta,   sobs_t, &
!$OMP shfl_s_t, &
!$OMP evap_s_t,  tskin_t, &
!$OMP ustr_s_t,  vstr_s_t, &
!$OMP zae,       zvar, &
!$OMP ztice,     ztske1, &
!$OMP ztskm1m,   ztskrad, &
!$OMP zsigflt,   zfrti, &
!$OMP l_land,    ldummy_vdf_a &
!$OMP ) ICON_OMP_GUIDED_SCHEDULE

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
        lnd_diag%qv_s (:,jb) = 0._wp
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
    ELSE IF (atm_phy_nwp_config(jg)%itype_z0 == 2) THEN
      ! specify land-cover-related roughness length over land points
      ! note:  water points are set in turbdiff
      gz0(:) = 0._wp
      DO jt = 1, nsfc_subs
        DO jc = i_startidx, i_endidx
          IF (ext_data%atm%fr_land(jc,jb) > 0.5_wp) THEN
            lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
            gz0(jc) = gz0(jc) + ext_data%atm%lc_frac_t(jc,jb,jt) * grav * ( &
             (1._wp-lnd_diag%snowfrac_t(jc,jb,jt))*z0_lu(lc_class) +        &
              lnd_diag%snowfrac_t(jc,jb,jt)*0.5_wp*z0_lu(21) ) ! 21 = snow/ice
          ENDIF
        ENDDO
      ENDDO
      DO jc = i_startidx, i_endidx
        IF (ext_data%atm%fr_land(jc,jb) > 0.5_wp) THEN
          prm_diag%gz0(jc,jb) = gz0(jc)
        ENDIF
      ENDDO
    ENDIF

    IF ( atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN

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
!     ! Artificial limitation of input TKE to 60 m^2/s^2
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
!          ! Artificial limitation of turbulent temperature tendency to 0.025 K/s
!          p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
!            &  + tcall_turb_jg*MIN(0.025_wp,MAX(-0.025_wp,prm_nwp_tend%ddt_temp_turb(jc,jk,jb)))
          p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
           &  + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
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
!     p_prog_rcf%tke              (1:i_endidx,:,jb)         = 0._wp

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
            !& pqsat_sfc = zdummy_qvsfc  (:,:,jb)                                    ,&! in
            & pocu      = prm_diag%ocu   (:,jb),  pocv   = prm_diag%ocv     (:,jb)   ,&! in
            & ppsfc     = p_diag%pres_sfc(:,jb),  pcoriol= p_patch%cells%f_c(:,jb)   ,&! in
            !
            & pum1      = p_diag%u   (:,:,jb),   pvm1 = p_diag%v         (:,:,jb)    ,&! in
            & ptm1      = p_diag%temp(:,:,jb),   pqm1 = p_prog_rcf%tracer(:,:,jb,iqv),&! in
            & pxlm1     = p_prog_rcf%tracer(:,:,jb,iqc)                              ,&! in
            & pxim1     = p_prog_rcf%tracer(:,:,jb,iqi)                              ,&! in
            & pxm1      = z_plitot (:,:,jb) ,    pxtm1 = zdummy_it       (:,:,:,jb)  ,&! in
            !
            & paphm1 = p_diag%pres_ifc(:,:,jb),  papm1 = p_diag%pres         (:,:,jb),&! in
            & pdelpm1= p_diag%dpres_mc(:,:,jb),  pgeom1= p_metrics%geopot_agl(:,:,jb),&! in 
            & ptvm1  = p_diag%tempv   (:,:,jb),  paclc = prm_diag%tot_cld(:,:,jb,icc),&! in
            & ptkem1 = p_prog_now_rcf%tke (:,2:nlevp1,jb)                            ,&! in
            !
            & pxt_emis= zdummy_ith    (:,:,jb)                                       ,&! in
            & ptsfc_tile = zdummy_tsfc(:,:,jb)                                       ,&! inout
            & pxvar   = zdummy_i      (:,:,jb)                                       ,&! inout
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
            & pghpbl   = prm_diag%ghpbl (:,jb),     pri = prm_diag%ri (:,2:nlevp1,jb),&! out
            & pmixlen  = prm_diag%mixlen(:,2:nlevp1,jb)                              ,&! out
            & pcfm     = prm_diag%cfm   (:,2:nlevp1,jb)                              ,&! out
            & pcfh     = prm_diag%cfh   (:,2:nlevp1,jb)                              ,&! out
            & pcfv     = prm_diag%cfv   (:,2:nlevp1,jb)                              ,&! out
            & pcfm_tile= prm_diag%cfm_tile(:,jb,:)                                   ,&! out
            & pcfh_tile= prm_diag%cfh_tile(:,jb,:)                                   ,&! out
            & pcftke   = prm_diag%cftke (:,2:nlevp1,jb)                              ,&! out
            & pcfthv   = prm_diag%cfthv (:,2:nlevp1,jb))                               ! out


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

!     Calculate vertical velocity in p-system
      
      DO jk = 1,p_patch%nlev
        DO jc = i_startidx,i_endidx
          z_omega_p(jc,jk) = - p_prog%w(jc,jk,jb) * p_prog%rho(jc,jk,jb) * grav
        ENDDO
      ENDDO

      icnt = 0

      DO jc = i_startidx, i_endidx
        zchar  (jc) = 0.018_wp ! default value from IFS if no wave model
        zucurr (jc) = 0.0_wp
        zvcurr (jc) = 0.0_wp
        ztice  (jc) = 273.0_wp ! top level ice temperature ???????
        ztske1 (jc) = 0.0_wp   ! skin temperature tendency
        ztskm1m(jc) = lnd_prog_now%t_g (jc,jb) ! skin temperature
        ztskrad(jc) = lnd_prog_now%t_g (jc,jb) ! skin temperature at last radiation step ????
        zsigflt(jc) = 0.0_wp   ! just for testing (standard dev. of filtered orogrphy)
      ENDDO

      DO jk = 1,p_patch%nlev
        DO jc = i_startidx, i_endidx
          zsoteu (jc,jk) = 0.0_wp
          zsotev (jc,jk) = 0.0_wp
          zsobeta(jc,jk) = 0.0_wp
          zae    (jc,jk) = 0.0_wp   ! cloud tendency ???
          zvar   (jc,jk) = 0.0_wp   ! qt,var should be prognostic !!!
        ENDDO
      ENDDO

      DO jt = 1,nsfc_subs
        DO jc = i_startidx, i_endidx
          sobs_t  (jc,jt) = prm_diag%swflxsfc(jc,jb)   ! simple grid-mean flux (should be tile albedo specific)!!!!!
          shfl_s_t(jc,jt) = prm_diag%shfl_s  (jc,jb)   ! should be tile specific !!!
          evap_s_t(jc,jt) = prm_diag%lhfl_s  (jc,jb) / alv ! evaporation [kg/(m2 s)]  -"-
          tskin_t (jc,jt) = lnd_prog_now%t_g (jc,jb)   ! should be tile specific
          ustr_s_t(jc,jt) = 0.0_wp                     ! prognostic surface stress U  !!!
          vstr_s_t(jc,jt) = 0.0_wp                     ! prognostic surface stress V  !!!
          zfrti   (jc,jt) = 0.0_wp                     ! all zero but tile1=1.0 ... all ocean ???
          zdummy_vdf_5a(jc,jt) = 0.3_wp                ! surface albedo ????????
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        zfrti(jc,1) = 1.0_wp                           ! all zero but tile1=1.0 ... all ocean ???
      ENDDO

      DO jk = 1,nlev_soil
        DO jc = i_startidx, i_endidx
          zdummy_vdf_4a(jc,jk) = lnd_prog_now%t_so_t(jc,jk,jb,1) ! simple: take one tile #1 ???
          zdummy_vdf_4b(jc,jk) = lnd_prog_now%w_so_t(jc,jk,jb,1) ! ---
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        IF ( ext_data%atm%fr_land(jc,jb) > 0.5_wp ) THEN
          l_land(jc) = .true.
        ELSE
          l_land(jc) = .false.
        ENDIF
      ENDDO

      DO jc = i_startidx, i_endidx
        idummy_vdf_0a(jc) = 16    !KTVL  ???
        idummy_vdf_0b(jc) = 3     !KTVH  ???
        zdummy_vdf_1a(jc) = 0.5_wp!PCVL  ???
        zdummy_vdf_1b(jc) = 0.5_wp!PCVL  ???
        idummy_vdf_0c(jc) = 1     !KSOTY ??? soil type (needs to be specified)
        zdummy_vdf_1f(jc) = 1.0_wp!maximum skin reservoir capacity (~1mm = 1kg/m2) ??? needs to be done physically by TERRA
        zdummy_vdf_1c(jc) = 0.0_wp!lake ice thickness ??? (no lakes???)
        zdummy_vdf_1d(jc) = 273.0_wp !lake ice temperature ??? (no lakes???)
      ENDDO

!     Tendencies are set to include dynamics and radiation
!     ATTENTION: currently for simplicity all input tendencies = 0
!                when updated after vdfouter difference needed (see convection)

      prm_nwp_tend%ddt_u_turb     (:,:,jb)     = 0._wp
      prm_nwp_tend%ddt_v_turb     (:,:,jb)     = 0._wp
      prm_nwp_tend%ddt_temp_turb  (:,:,jb)     = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc) = 0._wp
      prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi) = 0._wp

      CALL vdfouter ( &
        & CDCONF  = ' '                                        ,&! (IN)  unused
        & KIDIA   = i_startidx                                 ,&! (IN)   
        & KFDIA   = i_endidx                                   ,&! (IN)   
        & KLON    = nproma                                     ,&! (IN)   
        & KLEV    = p_patch%nlev                               ,&! (IN)   
        & KLEVS   = nlev_soil                                  ,&! (IN)
        & KSTEP   = 0                                          ,&! (IN)  unused: current time step
        & KTILES  = nsfc_subs                                  ,&! (IN)
        & KTRAC   = itrac_vdf                                  ,&! (IN)  default 0 (itrac?)
        & KLEVSN  = nlev_snow                                  ,&! (IN)  # snow layers (1!)
        & KLEVI   = 1                                          ,&! (IN)  # sea ice layers
        & KDHVTLS = 3                                          ,&! (IN)  DDH dimensions
        & KDHFTLS = 8                                          ,&! (IN)   - " -
        & KDHVTSS = 6                                          ,&! (IN)   - " -  
        & KDHFTSS = 9                                          ,&! (IN)   - " -  
        & KDHVTTS = 4                                          ,&! (IN)   - " -  
        & KDHFTTS = 11                                         ,&! (IN)   - " -  
        & KDHVTIS = 4                                          ,&! (IN)   - " - 
        & KDHFTIS = 9                                          ,&! (IN)   - " - 
        & PTSPHY  = tcall_turb_jg                              ,&! (IN)   
        & KTVL    = idummy_vdf_0a                              ,&! (IN)  input for TESSEL   
        & KTVH    = idummy_vdf_0b                              ,&! (IN)  input for TESSEL
        & KCNT    = icnt                                       ,&! (INOUT)
        & PCVL    = zdummy_vdf_1a                              ,&! (IN)  input for TESSEL
        & PCVH    = zdummy_vdf_1b                              ,&! (IN)  input for TESSEL   
!xmk ?  & PSIGFLT = ext_data%atm%sso_stdh(:,jb)                ,&! (IN)  input for TOFD (needs to be passed down!!!)
        & PSIGFLT = zsigflt                                    ,&! (IN)  input for TOFD (needs to be passed down!!!)
        & PUM1    = p_diag%u(:,:,jb)                           ,&! (IN)   
        & PVM1    = p_diag%v(:,:,jb)                           ,&! (IN)   
        & PTM1    = p_diag%temp(:,:,jb)                        ,&! (IN)   
        & PQM1    = p_prog_rcf%tracer(:,:,jb,iqv)              ,&! (IN)   
        & PLM1    = p_prog_rcf%tracer(:,:,jb,iqc)              ,&! (IN)   
        & PIM1    = p_prog_rcf%tracer(:,:,jb,iqi)              ,&! (IN)   
        & PAM1    = prm_diag%tot_cld (:,:,jb,icc)              ,&! (IN)   
        & PCM1    = zdummy_vdf_6a                              ,&! (IN)  tracer - for VDF transport
        & PAPHM1  = p_diag%pres_ifc         (:,:,jb)           ,&! (IN)   
        & PAPM1   = p_diag%pres             (:,:,jb)           ,&! (IN)   
        & PGEOM1  = p_metrics%geopot_agl    (:,:,jb)           ,&! (IN)   
        & PGEOH   = p_metrics%geopot_agl_ifc(:,:,jb)           ,&! (IN)   
        & PTSKM1M = ztskm1m                                    ,&! (IN)  T,skin
        & PTSAM1M = zdummy_vdf_4a                              ,&! (IN)  T,soil
        & PWSAM1M = zdummy_vdf_4b                              ,&! (IN)  Q,soil
        & PSSRFL  = prm_diag%swflxsfc(:,jb)                    ,&! (IN)   
        & PSLRFL  = prm_diag%lwflxsfc(:,jb)                    ,&! (IN)   
        & PEMIS   = ext_data%atm%emis_rad(:,jb)                ,&! (IN)   
        & PHRLW   = prm_nwp_tend%ddt_temp_radlw(:,:,jb)        ,&! (IN)   
        & PHRSW   = prm_nwp_tend%ddt_temp_radsw(:,:,jb)        ,&! (IN)   
        & PTSNOW  = lnd_prog_now%t_snow_t(:,jb,1)              ,&! (IN)  T,snow - unused (attention: tile 1????)
        & PTICE   = ztice                                      ,&! (IN)  T,ice  - unused   
        & PHLICE  = zdummy_vdf_1c                              ,&! (IN)  lake ice thickness   - unused
        & PTLICE  = zdummy_vdf_1d                              ,&! (IN)  lake ice temperature - unused 
        & PTLWML  = zdummy_vdf_1e                              ,&! (IN)  lake mean water T    - unused
        & PSST    = lnd_prog_now%t_g(:,jb)                     ,&! (IN)  SST
        & KSOTY   = idummy_vdf_0c                              ,&! (IN)  soil type
!xmk ?  & PFRTI   = ext_data%atm%lc_frac_t(:,jb,:)             ,&! (IN)  tile fraction 
        & PFRTI   = zfrti                                      ,&! (IN)  tile fraction 
        & PALBTI  = zdummy_vdf_5a                              ,&! (IN)  tile albedo
        & PWLMX   = zdummy_vdf_1f                              ,&! (IN)  maximum skin reservoir capacity
        & PCHAR   = zchar                                      ,&! (IN)  Charnock parameter (for z0 over ocean)
        & PUCURR  = zucurr                                     ,&! (IN)  Ocean current x 
        & PVCURR  = zvcurr                                     ,&! (IN)  Ocean current y
        & PTSKRAD = ztskrad                                    ,&! (IN)  unused: T,skin at last radiation step 
        & PCFLX   = zdummy_vdf_6c                              ,&! (IN)  unused: surface trace flux
        & PSOTEU  = zsoteu                                     ,&! (IN)  unused: Explicit part of U-tendency from SSO
        & PSOTEV  = zsotev                                     ,&! (IN)  unused: Explicit part of V-tendency from SSO  
        & PSOBETA = zsobeta                                    ,&! (IN)  unused: Implicit part of subgrid orography
        & PVERVEL = z_omega_p                                  ,&! (IN)   
        & PZ0M    = prm_diag%z0m(:,jb)                         ,&! (INOUT) z0,m (reduced for TOFD???) 
        & PZ0H    = prm_diag%z0m(:,jb)                         ,&! (INOUT) z0,h (* factor ????)    
        & PVDIS   = zdummy_vdf_1g                              ,&! (OUT) optional out: turbulent dissipation
        & PVDISG  = zdummy_vdf_1h                              ,&! (OUT) optional out: SO dissipation
        & PDISGW3D= zdummy_vdf_2a                              ,&! (OUT) optional out: 3D stoch. phys. dissipation
        & PAHFLEV = zdummy_vdf_1i                              ,&! (OUT) optional out: latent heat flux (snow/ice free part)
        & PAHFLSB = zdummy_vdf_1j                              ,&! (OUT) optional out: latent heat flux (snow/ice covered part)
        & PFWSB   = zdummy_vdf_1k                              ,&! (OUT) optional out: evaporation of snow
        & PBIR    = zdummy_vdf_1l                              ,&! (OUT) optional out: BIR buoyancy flux integral ratio
        & PVAR    = p_prog_rcf%tracer(:,:,jb,iqtvar)           ,&! (INOUT) qt,variance - prognostic advected tracer
        & PU10M   = prm_diag%u_10m(:,jb)                       ,&! (OUT)  
        & PV10M   = prm_diag%v_10m(:,jb)                       ,&! (OUT)  
        & PT2M    = prm_diag%t_2m (:,jb)                       ,&! (OUT)  
        & PD2M    = prm_diag%td_2m(:,jb)                       ,&! (OUT)  
        & PQ2M    = prm_diag%qv_2m(:,jb)                       ,&! (OUT)  
        & PZINV   = zdummy_vdf_1m                              ,&! (OUT) optional out: PBL HEIGHT (moist parcel, not for stable PBL)
        & PBLH    = zdummy_vdf_1n                              ,&! (OUT) optional out: PBL HEIGHT (dry diagnostic based on Ri#)
        & KHPBLN  = idummy_vdf_0d                              ,&! (OUT) optional out: PBL top level 
        & KVARTOP = idummy_vdf_0e                              ,&! (OUT) optional out: top level of predictied qt,var
        & PSSRFLTI= sobs_t                                     ,&! (INOUT) net SW sfc flux for each tile (use tile ablbedo!!!)
        & PEVAPSNW= zdummy_vdf_1o                              ,&! (OUT) optional out: evaporation from snow under forest
        & PGUST   = zdummy_vdf_1p                              ,&! (OUT) optional out: 10m gust
        & PWUAVG  = zdummy_vdf_1q                              ,&! (OUT) optional out: w,up averaged
        & LDNODECP= ldummy_vdf_a                               ,&! (OUT) optional out: no decoupling allowed 
        & KPBLTYPE= idummy_vdf_0f                              ,&! (OUT) optional out: PBL type
        & PLDIFF  = zdummy_vdf_2b                              ,&! (OUT) optional out: contrib to PBL cond. by passive clouds
        & PFPLVL  = zdummy_vdf_3a                              ,&! (OUT) optional out: PBL rain flux 
        & PFPLVN  = zdummy_vdf_3b                              ,&! (OUT) optional out: PBL snow flux
        & PFHPVL  = zdummy_vdf_3c                              ,&! (OUT) optional out: PBL rain enthalpy flux
        & PFHPVN  = zdummy_vdf_3d                              ,&! (OUT) optional out: PBL snow enthalpy flux
        & PEXTR2  = zdummy_vdf_7a                              ,&! (IN)    optional out:  - " -
        & PEXTRA  = zdummy_vdf_7b                              ,&! (INOUT) optional out:  - " -
        & KLEVX   = p_patch%nlev                               ,&! (IN)    out
        & KFLDX   = 0                                          ,&! (IN)    out
        & KFLDX2  = 0                                          ,&! (IN)    out
        & LLDIAG  = .FALSE.                                    ,&! (IN)    out
        & PTE     = prm_nwp_tend%ddt_temp_turb  (:,:,jb)       ,&! (INOUT)
        & PQE     = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv)   ,&! (INOUT)
        & PLE     = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc)   ,&! (INOUT)
        & PIE     = prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqi)   ,&! (INOUT)
        & PAE     = zae                                        ,&! (INOUT)
        & PVOM    = prm_nwp_tend%ddt_u_turb(:,:,jb)            ,&! (INOUT)
        & PVOL    = prm_nwp_tend%ddt_v_turb(:,:,jb)            ,&! (INOUT)
        & PTENC   = zdummy_vdf_6b                              ,&! (INOUT) optional inout: tracer tendency
        & PTSKE1  = ztske1                                     ,&! (INOUT) unused: T,skin tendency
        & PUSTRTI = ustr_s_t                                   ,&! (INOUT) tile u stress
        & PVSTRTI = vstr_s_t                                   ,&! (INOUT) tile v stress
        & PAHFSTI = shfl_s_t                                   ,&! (INOUT) tile sensible heat flux
        & PEVAPTI = evap_s_t                                   ,&! (INOUT) tile latent heat flux
        & PTSKTI  = tskin_t                                    ,&! (INOUT) now! ???ocean???  lnd_prog_now%t_g(:,jb)
        & PDIFTS  = zdummy_vdf_3e                              ,&! (OUT)  optional out: turbulent heat flux
        & PDIFTQ  = zdummy_vdf_3f                              ,&! (OUT)  optional out: turbulent moisture flux
        & PDIFTL  = zdummy_vdf_3g                              ,&! (OUT)  optional out: turbulent liquid water flux
        & PDIFTI  = zdummy_vdf_3h                              ,&! (OUT)  optional out: turbulent ice water flux
        & PSTRTU  = zdummy_vdf_3i                              ,&! (OUT)  optional out: turbulent U flux
        & PSTRTV  = zdummy_vdf_3j                              ,&! (OUT)  optional out: turbulent V flux
        & PTOFDU  = zdummy_vdf_1r                              ,&! (OUT)  optional out: TOFD U flux
        & PTOFDV  = zdummy_vdf_1s                              ,&! (OUT)  optional out: TOFD V flux
        & PSTRSOU = zdummy_vdf_3k                              ,&! (OUT)  optional out: SSO U flux 
        & PSTRSOV = zdummy_vdf_3l                              ,&! (OUT)  optional out: SSO V flux
        & PKH     = prm_diag%tkvh(:,:,jb)                      ,&! (OUT)
        & LDLAND  = l_land                                      &! (IN)   logical for land
!       & PDHTLS  = ...                                        ,&! (OUT)  optional out: DDH
!       & PDHTSS  = ...                                        ,&! (OUT)  optional out: DDH
!       & PDHTTS  = ...                                        ,&! (OUT)  optional out: DDH
!       & PDHTIS  = ...                                         &! (OUT)  optional out: DDH
        & )

! Update QV, QC, QI and temperature with turbulence tendencies

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqv) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqv) &
               &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv))
          p_prog_rcf%tracer(jc,jk,jb,iqc) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
               &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
          p_prog_rcf%tracer(jc,jk,jb,iqi) =MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqi) &
               &           + tcall_turb_jg*prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqi))
          p_diag%temp(jc,jk,jb) =                           p_diag%temp(jc,jk,jb) &
           &               + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        ENDDO
      ENDDO

! Some diagnostic values have to be set !!!!
      DO jc = i_startidx, i_endidx
        prm_diag%rh_2m(jc,jb)  = 0.0_wp
        prm_diag%shfl_s(jc,jb) = shfl_s_t(jc,1)           ! should be tile mean !!!
        prm_diag%lhfl_s(jc,jb) = evap_s_t(jc,1) * alv     ! should be tile mean !!!
      ENDDO

    ENDIF !inwp_turb

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE nwp_turbulence


END MODULE mo_nwp_turb_interface


