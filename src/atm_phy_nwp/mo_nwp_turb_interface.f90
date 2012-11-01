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
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_physical_constants,   ONLY: rd_o_cpd, grav
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,        ONLY: phy_params 
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_data_turbdiff,        ONLY: get_turbdiff_param
  USE src_turbdiff,            ONLY: turbtran, turbdiff
  USE mo_satad,                ONLY: sat_pres_water, spec_humi  
  USE mo_run_config,           ONLY: ltestcase
  USE mo_nh_testcases,         ONLY: nh_test_name
  USE mo_nh_wk_exp,            ONLY: qv_max_wk
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, nlists_water, ntiles_water, &
    &                                isub_water, isub_seaice

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
                          & wtr_prog_now,                      & !>in 
                          & lnd_prog_now,                      & !>inout 
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
  TYPE(t_wtr_prog),            INTENT(in)   :: wtr_prog_now    !< prog vars for wtr
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for 
                                                               !< turbulence

  ! Local array bounds

  INTEGER :: nblks_c, nblks_e        !> number of blocks for cells / edges
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: i_nchdom                !< domain index

  ! Local scalars:

  INTEGER :: jc,jk,jb,jt,jg,ic,i_count,jt1      !loop indices

  ! local variables for turbdiff

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''

  INTEGER  :: nlev, nlevp1                          !< number of full and half levels
  INTEGER  :: lc_class, i_lc_si                     !< land-cover class

  REAL(wp) :: z_tvs(nproma,p_patch%nlevp1,1)        !< aux turbulence velocity scale [m/s]

  REAL(wp) :: fr_land_t(nproma),depth_lk_t(nproma),h_ice_t(nproma),area_frac

  ! Local fields needed to reorder turbtran input/output fields for tile approach

  ! 2D fields
  REAL(wp), DIMENSION(nproma,ntiles_total+nlists_water) :: gz0t_t, tcm_t, tch_t, tfm_t, tfh_t, tfv_t, &  
   t_2m_t, qv_2m_t, td_2m_t, rh_2m_t, u_10m_t, v_10m_t, t_g_t, qv_s_t, pres_sfc_t, sai_t

  ! 3D full-level fields
  REAL(wp), DIMENSION(nproma,2,ntiles_total+nlists_water) :: u_t, v_t, temp_t, pres_t, qv_t, qc_t, tkvm_t, tkvh_t

  ! 3D half-level fields
  REAL(wp), DIMENSION(nproma,3,ntiles_total+nlists_water) :: z_ifc_t, w_t, rcld_t

  ! SQRT(2*TKE)
  REAL(wp) :: tvs_t(nproma,3,1,ntiles_total+nlists_water)

  INTEGER, POINTER :: ilist(:)  ! pointer to tile index list

!--------------------------------------------------------------


  ! number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  i_lc_si= ext_data%atm%i_lc_snow_ice

  IF (msg_level >= 15) CALL message('mo_nwp_turb:', 'turbulence')
    
  ! local variables related to the blocking
  
  nblks_c   = p_patch%nblks_int_c
  nblks_e   = p_patch%nblks_int_e
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
!$OMP DO PRIVATE(jb,jt,jc,jk,ic,jt1,ilist,i_startidx,i_endidx,i_count,ierrstat,errormsg,eroutine,&
!$OMP lc_class,z_tvs,gz0t_t,tcm_t,tch_t,tfm_t,tfh_t,tfv_t,t_g_t,qv_s_t,t_2m_t,qv_2m_t,           &  
!$OMP td_2m_t,rh_2m_t,u_10m_t,v_10m_t,tvs_t,pres_sfc_t,u_t,v_t,temp_t,pres_t,qv_t,qc_t,tkvm_t,   &
!$OMP tkvh_t,z_ifc_t,w_t,rcld_t,sai_t,fr_land_t,depth_lk_t,h_ice_t,area_frac) ICON_OMP_GUIDED_SCHEDULE

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
         !> adjust humidity at water surface because of changed surface pressure
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
      DO jt = 1, ntiles_total
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data%atm%gp_count_t(jb,jt)
          jc = ext_data%atm%idx_lst_t(ic,jb,jt)
          lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
          prm_diag%gz0_t(jc,jb,jt) = grav * (                                     &
           (1._wp-lnd_diag%snowfrac_t(jc,jb,jt))*ext_data%atm%z0_lcc(lc_class) +  &
            lnd_diag%snowfrac_t(jc,jb,jt)*0.5_wp*ext_data%atm%z0_lcc(i_lc_si) )
        ENDDO
      ENDDO
      ! water points - preliminary implementation
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, ext_data%atm%spw_count(jb)
        jc = ext_data%atm%idx_lst_spw(ic,jb)
        prm_diag%gz0_t(jc,jb,isub_water) = prm_diag%gz0(jc,jb)
      ENDDO
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, ext_data%atm%fp_count(jb)
        jc = ext_data%atm%idx_lst_fp(ic,jb)
        prm_diag%gz0_t(jc,jb,isub_water) = prm_diag%gz0(jc,jb)
      ENDDO
      ! seaice points - preliminary implementation
!CDIR NODEP,VOVERTAKE,VOB
      DO ic = 1, ext_data%atm%spi_count(jb)
        jc = ext_data%atm%idx_lst_spi(ic,jb)
        prm_diag%gz0_t(jc,jb,isub_seaice) = prm_diag%gz0(jc,jb)
      ENDDO
      IF (ntiles_total == 1) THEN
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data%atm%lp_count(jb)
          jc = ext_data%atm%idx_lst_lp(ic,jb)
          prm_diag%gz0(jc,jb) = prm_diag%gz0_t(jc,jb,1)
        ENDDO
      ENDIF
    ELSE ! uniform tile-averaged roughness length if SSO contribution is to be included
      DO jt = 1, ntiles_total + ntiles_water
        DO jc = i_startidx, i_endidx
          prm_diag%gz0_t(jc,jb,jt) = prm_diag%gz0(jc,jb)
        ENDDO
      ENDDO
    ENDIF

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


      ! First call of turbtran for all grid points (water points with > 50% water
      ! fraction and tile 1 of the land points)
      IF (ntiles_total == 1) THEN ! tile approach not used; use tile-averaged fields from extpar

        CALL turbtran(iini=0, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1,                     &
          & ie=nproma, ke=nlev, ke1=nlevp1,                                                     &
          & istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx,           &
          & l_hori=phy_params(jg)%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),                    &
          & fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb),           &
          & sai=prm_diag%sai(:,jb), h_ice=wtr_prog_now%h_ice(:,jb), ps=p_diag%pres_sfc(:,jb),   &
          & t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),                               &
          & u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb),                         &
          & T=p_diag%temp(:,:,jb), prs=p_diag%pres(:,:,jb),                                     &
          & qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),                 &
          & gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),             &
          & tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb),             &
          & tke=z_tvs(:,:,:), tkvm=prm_diag%tkvm(:,:,jb),  &!  edr =prm_diag%edr(:,:,jb),       &
          & tkvh=prm_diag%tkvh(:,:,jb), rcld=prm_diag%rcld(:,:,jb),                             &
          & t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), td_2m=prm_diag%td_2m(:,jb),   &
          & rh_2m=prm_diag%rh_2m(:,jb), u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb), &
          & ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

      ELSE ! tile approach used

        ! preset variables for land tile indices
        fr_land_t(:)  = 1._wp
        depth_lk_t(:) = 0._wp
        h_ice_t(:)    = 0._wp

        ! Loop over index lists for land tile points, sea/lake points and seaice points
        DO  jt = 1, ntiles_total + nlists_water

          IF (jt <= ntiles_total) THEN ! land tile points
            jt1 = jt
            i_count = ext_data%atm%gp_count_t(jb,jt)
            ilist => ext_data%atm%idx_lst_t(:,jb,jt)
          ELSE IF (jt == ntiles_total + 1) THEN ! sea points (seaice points excluded)
            jt1 = isub_water   !jt
            i_count = ext_data%atm%spw_count(jb)
            ilist => ext_data%atm%idx_lst_spw(:,jb)
            fr_land_t(:) = 0._wp
          ELSE IF (jt == ntiles_total + 2) THEN ! lake points
            jt1 = isub_water   !ntiles_total + 1
            i_count = ext_data%atm%fp_count(jb)
            ilist => ext_data%atm%idx_lst_fp(:,jb)
            ! depth_lk_t(:) = ...
            ! h_ice_t(:) = ...
          ELSE ! IF (jt == ntiles_total + 3) THEN ! seaice points
            jt1 = isub_seaice    !ntiles_total + 2
            i_count = ext_data%atm%spi_count(jb)
            ilist => ext_data%atm%idx_lst_spi(:,jb)
            fr_land_t (:) = 0._wp
            depth_lk_t(:) = 0._wp
            h_ice_t   (:) = 1._wp  ! prelim. implementation. Only needed for checking whether 
                                   ! ice is present or not
          ENDIF

          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

          ! Copy input fields to the local re-indexed variables
          ! It remains to be determined which of the model levels are actually needed for non-init calls
          DO ic = 1, i_count
            jc = ilist(ic)
            gz0t_t(ic,jt) = prm_diag%gz0_t(jc,jb,jt1)
            t_g_t(ic,jt)  = lnd_prog_now%t_g_t(jc,jb,jt1)
            qv_s_t(ic,jt) = lnd_diag%qv_s_t(jc,jb,jt1)
            sai_t(ic,jt)  = ext_data%atm%sai_t(jc,jb,jt1)
            z_ifc_t(ic,1:3,jt) = p_metrics%z_ifc(jc,nlev-1:nlev+1,jb)
            pres_sfc_t(ic,jt)  = p_diag%pres_sfc(jc,jb)
            u_t(ic,1:2,jt)     = p_diag%u(jc,nlev-1:nlev,jb)
            v_t(ic,1:2,jt)     = p_diag%v(jc,nlev-1:nlev,jb)
            w_t(ic,1:3,jt)     = p_prog%w(jc,nlev-1:nlev+1,jb)
            temp_t(ic,1:2,jt)  = p_diag%temp(jc,nlev-1:nlev,jb)
            pres_t(ic,1:2,jt)  = p_diag%pres(jc,nlev-1:nlev,jb)
            qv_t(ic,1:2,jt)    = p_prog_rcf%tracer(jc,nlev-1:nlev,jb,iqv)
            qc_t(ic,1:2,jt)    = p_prog_rcf%tracer(jc,nlev-1:nlev,jb,iqc)
            tvs_t(ic,1:3,1,jt) = z_tvs(jc,nlev-1:nlev+1,1)
            tkvm_t(ic,1:2,jt)  = prm_diag%tkvm(jc,nlev-1:nlev,jb)
            tkvh_t(ic,1:2,jt)  = prm_diag%tkvh(jc,nlev-1:nlev,jb)
            rcld_t(ic,1:3,jt)  = prm_diag%rcld(jc,nlev-1:nlev+1,jb)
          ENDDO

          CALL turbtran(iini=0, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1,                   &
            & ie=nproma, ke=2, ke1=3,                                                           &
            & istart=1, iend=i_count, istartpar=1, iendpar=i_count,                             &
            & l_hori=phy_params(jg)%mean_charlen, hhl=z_ifc_t(:,:,jt),                          &
            & fr_land=fr_land_t(:), depth_lk=depth_lk_t(:),                                     &
            & sai=sai_t(:,jt), h_ice=h_ice_t(:), ps=pres_sfc_t(:,jt),                           &
            & t_g=t_g_t(:,jt), qv_s=qv_s_t(:,jt), u=u_t(:,:,jt), v=v_t(:,:,jt), w=w_t(:,:,jt),  &
            & T=temp_t(:,:,jt), prs=pres_t(:,:,jt), qv=qv_t(:,:,jt), qc=qc_t(:,:,jt),           &
            & gz0=gz0t_t(:,jt), tcm=tcm_t(:,jt), tch=tch_t(:,jt),                               &
            & tfm=tfm_t(:,jt), tfh=tfh_t(:,jt), tfv=tfv_t(:,jt),                                &
            & tke=tvs_t(:,:,:,jt), tkvm=tkvm_t(:,:,jt),  &!  edr =prm_diag%edr(:,:,jb),         &
            & tkvh=tkvh_t(:,:,jt), rcld=rcld_t(:,:,jt),                                         &
            & t_2m=t_2m_t(:,jt), qv_2m=qv_2m_t(:,jt), td_2m=td_2m_t(:,jt),                      &
            & rh_2m=rh_2m_t(:,jt), u_10m=u_10m_t(:,jt), v_10m=v_10m_t(:,jt),                    &
            & ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )

        ENDDO

        ! Aggregate tile-based output fields of turbtran over tiles
        ! i) initialize fields to zero before starting the summation
        prm_diag%gz0  (:,jb) = 0._wp
        prm_diag%tcm  (:,jb) = 0._wp
        prm_diag%tch  (:,jb) = 0._wp
        prm_diag%tfm  (:,jb) = 0._wp
        prm_diag%tfh  (:,jb) = 0._wp
        prm_diag%tfv  (:,jb) = 0._wp
        prm_diag%t_2m (:,jb) = 0._wp
        prm_diag%qv_2m(:,jb) = 0._wp
        prm_diag%td_2m(:,jb) = 0._wp
        prm_diag%rh_2m(:,jb) = 0._wp
        prm_diag%u_10m(:,jb) = 0._wp
        prm_diag%v_10m(:,jb) = 0._wp

        z_tvs(:,nlevp1,1)          = 0._wp
        prm_diag%tkvm(:,nlev,jb)   = 0._wp
        prm_diag%tkvh(:,nlev,jb)   = 0._wp
        prm_diag%rcld(:,nlevp1,jb) = 0._wp

        ! ii) loop over index lists
        DO  jt = 1, ntiles_total + nlists_water

          IF (jt <= ntiles_total) THEN ! land tile points
            jt1 = jt
            i_count = ext_data%atm%gp_count_t(jb,jt)
            ilist => ext_data%atm%idx_lst_t(:,jb,jt)
          ELSE IF (jt == ntiles_total + 1) THEN ! sea points (seaice points excluded)
            jt1 = isub_water   !jt
            i_count = ext_data%atm%spw_count(jb)
            ilist => ext_data%atm%idx_lst_spw(:,jb)
          ELSE IF (jt == ntiles_total + 2) THEN ! lake points
            jt1 = isub_water   !ntiles_total + 1
            i_count = ext_data%atm%fp_count(jb)
            ilist => ext_data%atm%idx_lst_fp(:,jb)
          ELSE ! IF (jt == ntiles_total + 3) THEN ! seaice points
            jt1 = isub_seaice   !ntiles_total + 2
            i_count = ext_data%atm%spi_count(jb)
            ilist => ext_data%atm%idx_lst_spi(:,jb)
          ENDIF

          IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count
            jc = ilist(ic)
            area_frac = ext_data%atm%frac_t(jc,jb,jt1)
            prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb)+gz0t_t(ic,jt)*area_frac
            prm_diag%tcm(jc,jb) = prm_diag%tcm(jc,jb)+tcm_t(ic,jt) *area_frac
            prm_diag%tch(jc,jb) = prm_diag%tch(jc,jb)+tch_t(ic,jt) *area_frac
            prm_diag%tfm(jc,jb) = prm_diag%tfm(jc,jb)+tfm_t(ic,jt) *area_frac
            prm_diag%tfh(jc,jb) = prm_diag%tfh(jc,jb)+tfh_t(ic,jt) *area_frac
            prm_diag%tfv(jc,jb) = prm_diag%tfv(jc,jb)+tfv_t(ic,jt) *area_frac

            z_tvs(jc,nlevp1,1)          = z_tvs(jc,nlevp1,1)+tvs_t(ic,3,1,jt)        *area_frac
            prm_diag%tkvm(jc,nlev,jb)   = prm_diag%tkvm(jc,nlev,jb)+tkvm_t(ic,2,jt)  *area_frac
            prm_diag%tkvh(jc,nlev,jb)   = prm_diag%tkvh(jc,nlev,jb)+tkvh_t(ic,2,jt)  *area_frac
            prm_diag%rcld(jc,nlevp1,jb) = prm_diag%rcld(jc,nlevp1,jb)+rcld_t(ic,3,jt)*area_frac

            prm_diag%t_2m(jc,jb)  = prm_diag%t_2m(jc,jb) +t_2m_t(ic,jt) *area_frac
            prm_diag%qv_2m(jc,jb) = prm_diag%qv_2m(jc,jb)+qv_2m_t(ic,jt)*area_frac
            prm_diag%td_2m(jc,jb) = prm_diag%td_2m(jc,jb)+td_2m_t(ic,jt)*area_frac
            prm_diag%rh_2m(jc,jb) = prm_diag%rh_2m(jc,jb)+rh_2m_t(ic,jt)*area_frac
            prm_diag%u_10m(jc,jb) = prm_diag%u_10m(jc,jb)+u_10m_t(ic,jt)*area_frac
            prm_diag%v_10m(jc,jb) = prm_diag%v_10m(jc,jb)+v_10m_t(ic,jt)*area_frac

!DR Test
!            IF (area_frac==0._wp) write(0,*) "area_frac zero at jc,jb,jt",jc,jb,jt1
!            IF (prm_diag%gz0(jc,jb)==0._wp) write(0,*) "prm_diag%gz0(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%tcm(jc,jb)==0._wp) write(0,*) "prm_diag%tcm(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%tch(jc,jb)==0._wp) write(0,*) "prm_diag%tch(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%tfm(jc,jb)==0._wp) write(0,*) "prm_diag%tfm(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%tfh(jc,jb)==0._wp) write(0,*) "prm_diag%tfh(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%tfv(jc,jb)==0._wp) write(0,*) "prm_diag%tfv(jc,jb) zero at jc,jb,jt",jc,jb,jt

!            IF (prm_diag%tkvm(jc,nlev,jb)==0._wp) write(0,*) "prm_diag%tkvm(jc,nlev,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%tkvh(jc,nlev,jb)==0._wp) write(0,*) "prm_diag%tkvh(jc,nlev,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%rcld(jc,nlevp1,jb)==0._wp) write(0,*) "prm_diag%rcld(jc,nlevp1,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%t_2m(jc,jb)==0._wp) write(0,*) "prm_diag%t_2m(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%qv_2m(jc,jb)==0._wp) write(0,*) "prm_diag%qv_2m(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%td_2m(jc,jb)==0._wp) write(0,*) "prm_diag%td_2m(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%rh_2m(jc,jb)==0._wp) write(0,*) "prm_diag%rh_2m(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%u_10m(jc,jb)==0._wp) write(0,*) "prm_diag%u_10m(jc,jb) zero at jc,jb,jt",jc,jb,jt
!            IF (prm_diag%v_10m(jc,jb)==0._wp) write(0,*) "prm_diag%v_10m(jc,jb) zero at jc,jb,jt",jc,jb,jt
!DR Test
          ENDDO

        ENDDO

      ENDIF ! tiles / no tiles

      ! Aggregation of variables needed as input for turbdiff or as diagnostic output

      CALL turbdiff(iini=0, lstfnct=.TRUE.,                                                     &
         &  dt_var=tcall_turb_jg, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1,                 &
         &  ie=nproma, ke=nlev, ke1=nlevp1,  kcm=nlevp1,                                        &
         &  istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx,           &
         &  l_hori=phy_params(jg)%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),                    &
         &  dp0=p_diag%dpres_mc(:,:,jb),                                                        &
         &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb),           &
         &  sai=prm_diag%sai(:,jb), h_ice=wtr_prog_now%h_ice (:,jb),                            &
         &  ps=p_diag%pres_sfc(:,jb), t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),     &
         &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb), T=p_diag%temp(:,:,jb),  &
         &  qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),                 &
         &  prs=p_diag%pres(:,:,jb), rho=p_prog%rho(:,:,jb), epr=p_prog%exner(:,:,jb),          &
         &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),             &
         &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb),             &
         &  tke=z_tvs (:,:,:)   ,&!  edr =prm_diag%edr(:,:,jb),                                 &
         &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh(:,:,jb), rcld=prm_diag%rcld(:,:,jb), &
         &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb),     &
         &  t_tens=prm_nwp_tend%ddt_temp_turb(:,:,jb),                                          &
         &  qv_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqv),                                   &
         &  qc_tens=prm_nwp_tend%ddt_tracer_turb(:,:,jb,iqc),                                   &
         &  tketens=prm_nwp_tend%ddt_tke(:,:,jb),                                               &
         &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb), vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) ,      &
         &  shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb),                         &
         &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )
          
      IF (ierrstat.NE.0) THEN
        CALL finish(eroutine, errormsg)
      END IF

      ! transform updated turbulent velocity scale back to TKE
      p_prog_rcf%tke(i_startidx:i_endidx,:,jb)= 0.5_wp                            &
        &                                     * (z_tvs(i_startidx:i_endidx,:,1))**2

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
!> GME turbulence scheme 
!-------------------------------------------------------------------------

      CALL finish('nwp_turb_interface','GME turbulence scheme not yet implemented')

    ENDIF !inwp_turb

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE nwp_turbulence

END MODULE mo_nwp_turb_interface
