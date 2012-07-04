!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! turbulence parameterisations:
!! inwp_turb == 1 == turbulence scheme by M. Raschendorfer run in COSMO
!! inwp_turb == 2 == turbulence scheme imported from the GME (to be implemented)
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
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, iqr, iqs
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_data_turbdiff,        ONLY: get_turbdiff_param
  USE src_turbdiff,            ONLY: turbtran, turbdiff
  USE mo_satad,                ONLY: sat_pres_water, spec_humi  
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_run_config,           ONLY: ltestcase
  USE mo_nh_testcases,         ONLY: nh_test_name
  USE mo_nh_wk_exp,            ONLY: qv_max_wk
  USE mo_lnd_nwp_config,       ONLY: nlev_soil, nlev_snow, nsfc_subs, lmulti_snow

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
                          & lnd_prog_now, lnd_prog_new,        & !>inout 
                          & lnd_diag                           ) !>inout


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
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new    !< prog vars for sfc
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
    &  z_tvs    (nproma,p_patch%nlevp1,p_patch%nblks_c,1)
  REAL(wp) :: gz0(nproma)

  INTEGER  :: nlev, nlevp1                          !< number of full and half levels
  INTEGER  :: lc_class, i_lc_si                     !< land-cover class
  INTEGER  :: ks, ke, ke1    ! For turbtran call

!--------------------------------------------------------------


  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  i_lc_si= ext_data%atm%i_lc_snow_ice(1)

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
!$OMP DO PRIVATE(jb,jt,jc,jk,i_startidx,i_endidx,ierrstat,errormsg,eroutine, &
!$OMP lc_class, gz0, ks, ke, ke1) ICON_OMP_GUIDED_SCHEDULE

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
             (1._wp-lnd_diag%snowfrac_t(jc,jb,jt))*ext_data%atm%z0_lcc(lc_class) +        &
              lnd_diag%snowfrac_t(jc,jb,jt)*0.5_wp*ext_data%atm%z0_lcc(i_lc_si) ) ! 21 = snow/ice
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
      z_tvs(i_startidx:i_endidx,:,jb,1)=  &
        &           SQRT(2._wp * p_prog_now_rcf%tke(i_startidx:i_endidx,:,jb))


!-------------------------------------------------------------------------
!< COSMO version by M. Raschendorfer  
!-------------------------------------------------------------------------

      ! To avoid copying in/out more field elements than needed
      ks = nlev-1
      ke = nlev
      ke1 = nlevp1

      CALL turbtran(iini=0, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1, &
!
         &  ie=nproma, ke=2, ke1=3,                                                   &
         &  istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx, &
!       
         &  l_hori=phy_params(jg)%mean_charlen, hhl=p_metrics%z_ifc(:,ks:ke1,jb),   &
!
         &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
         &  sai=prm_diag%sai(:,jb), h_ice=prm_diag%h_ice (:,jb), &
!         
         &  ps=p_diag%pres_sfc(:,jb), t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb), &
!           
         &  u=p_diag%u(:,ks:ke,jb), v=p_diag%v(:,ks:ke,jb), w=p_prog%w(:,ks:ke1,jb), &
         &  T=p_diag%temp(:,ks:ke,jb), prs=p_diag%pres(:,ks:ke,jb),  &
         &  qv=p_prog_rcf%tracer(:,ks:ke,jb,iqv), qc=p_prog_rcf%tracer(:,ks:ke,jb,iqc), &     
!         
         &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
         &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
!                
         &  tke=z_tvs (:,ks:ke1,jb,:) ,&!  edr =prm_diag%edr(:,ks:ke1,jb),                    &
         &  tkvm=prm_diag%tkvm(:,ks:ke,jb), tkvh=prm_diag%tkvh(:,ks:ke,jb), &
         &  rcld=prm_diag%rcld(:,ks:ke1,jb),                                &
!       
         &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), td_2m=prm_diag%td_2m(:,jb), &
         &  rh_2m=prm_diag%rh_2m(:,jb), u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb), &
!         
         &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )
          

      CALL turbdiff(iini=0, lstfnct=.TRUE., &
!
         &  dt_var=tcall_turb_jg, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1, &
!
         &  ie=nproma, ke=nlev, ke1=nlevp1,  kcm=nlevp1,  &
         &  istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx,  &
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
          p_diag%temp(jc,jk,jb) = p_diag%temp(jc,jk,jb)  &
           &  + tcall_turb_jg*prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        ENDDO
      ENDDO

    ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 2 ) THEN

!-------------------------------------------------------------------------
!> GME version 
!-------------------------------------------------------------------------

      CALL finish('nwp_turb_interface','GME turbulence scheme not yet implemented')

    ENDIF !inwp_turb

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE nwp_turbulence

END MODULE mo_nwp_turb_interface
