!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! turbulence parameterisations:
!! inwp_turb == 1 == turbulence scheme by M. Raschendorfer run in COSMO
!! inwp_turb == 2 == turbulence scheme imported from the GME
!! This module handles the computation of surface transfer coefficients, only.
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

MODULE mo_nwp_turbtrans_interface

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
  USE src_turbdiff,            ONLY: turbtran
  USE mo_satad,                ONLY: sat_pres_water, spec_humi  
  USE mo_gme_turbdiff,         ONLY: parturs
  USE mo_run_config,           ONLY: ltestcase
  USE mo_nh_testcases,         ONLY: nh_test_name
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, nlists_water, ntiles_water, &
    &                                isub_water, isub_seaice

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_turbtrans

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbtrans  ( tcall_turb_jg,                     & !>in
                          & p_patch,p_metrics,                 & !>in
                          & ext_data,                          & !>in
                          & p_prog,                            & !>inout
                          & p_prog_now_rcf, p_prog_rcf,        & !>in/inout
                          & p_diag ,                           & !>inout
                          & prm_diag,                          & !>inout
                          & wtr_prog_now,                      & !>in 
                          & lnd_prog_now,                      & !>inout 
                          & lnd_diag                           ) !>inout


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
  TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
  TYPE(t_nh_prog),      TARGET,INTENT(in)   :: p_prog_now_rcf  !<progs with red.
  TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
  TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
  TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
  TYPE(t_wtr_prog),            INTENT(in)   :: wtr_prog_now    !< prog vars for wtr
  TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now    !< prog vars for sfc
  TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
  REAL(wp),                    INTENT(in)   :: tcall_turb_jg   !< time interval for 
                                                               !< turbulence

  ! Local array bounds

  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices
  INTEGER :: i_nchdom                !< domain index

  ! Local scalars:

  INTEGER :: jc,jb,jt,jg,ic,i_count,jt1      !loop indices

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
  REAL(wp), DIMENSION(nproma,ntiles_total+nlists_water) :: gz0_t, tcm_t, tch_t, tfm_t, tfh_t, tfv_t, &  
   t_2m_t, qv_2m_t, td_2m_t, rh_2m_t, u_10m_t, v_10m_t, t_g_t, qv_s_t, pres_sfc_t, sai_t, shfl_s_t, lhfl_s_t

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

  IF (msg_level >= 15) CALL message('mo_nwp_turbtrans:', 'turbulence')
    
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

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,ic,jt1,ilist,i_startidx,i_endidx,i_count,ierrstat,errormsg,eroutine,&
!$OMP lc_class,z_tvs,gz0_t,tcm_t,tch_t,tfm_t,tfh_t,tfv_t,t_g_t,qv_s_t,t_2m_t,qv_2m_t,           &  
!$OMP td_2m_t,rh_2m_t,u_10m_t,v_10m_t,tvs_t,pres_sfc_t,u_t,v_t,temp_t,pres_t,qv_t,qc_t,tkvm_t,   &
!$OMP tkvh_t,z_ifc_t,w_t,rcld_t,sai_t,fr_land_t,depth_lk_t,h_ice_t,area_frac,shfl_s_t,lhfl_s_t)  &
!$OMP ICON_OMP_GUIDED_SCHEDULE

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

   !-------------------------------------------------------------------------
   !<  turbulent transfer
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
      ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 1 .OR.   &
        &       atm_phy_nwp_config(jg)%inwp_turb == 2) THEN
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
      ! NOTE:  open water, lake and sea-ice points are set in turbtran
      DO jt = 1, ntiles_total
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, ext_data%atm%gp_count_t(jb,jt)
          ! works for the following two cases
          ! 1. snow-covered and snow_free tiles treated separately
          ! 2. snow_covered and snow_free areas combined in one tile
          jc = ext_data%atm%idx_lst_t(ic,jb,jt)
          lc_class = MAX(1,ext_data%atm%lc_class_t(jc,jb,jt)) ! to avoid segfaults
          prm_diag%gz0_t(jc,jb,jt) = grav * (                                     &
           (1._wp-lnd_diag%snowfrac_t(jc,jb,jt))*ext_data%atm%z0_lcc(lc_class) +  &
            lnd_diag%snowfrac_t(jc,jb,jt)*0.5_wp*ext_data%atm%z0_lcc(i_lc_si) )
          !
          ! factor of 0.5 added by GZ, since z0_lcc(i_lc_si)=0.01 appears to be 
          ! quite large
        ENDDO
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

      ! note that TKE must be converted to the turbulence velocity scale SQRT(2*TKE)
      ! for turbdiff
      ! INPUT to turbdiff is timestep now
      z_tvs(i_startidx:i_endidx,:,1) =  &
        &           SQRT(2._wp * p_prog_now_rcf%tke(i_startidx:i_endidx,:,jb))


      ! First call of turbtran for all grid points (water points with > 50% water
      ! fraction and tile 1 of the land points)
      IF (ntiles_total == 1) THEN ! tile approach not used; use tile-averaged fields from extpar

        CALL turbtran(iini=0, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1,            & !in
          & ie=nproma, ke=nlev, ke1=nlevp1,                                            & !in
          & istart=i_startidx, iend=i_endidx, istartpar=i_startidx, iendpar=i_endidx,  & !in
          & l_hori=phy_params(jg)%mean_charlen, hhl=p_metrics%z_ifc(:,:,jb),           & !in
          & fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb),  & !in
          & sai=prm_diag%sai(:,jb), h_ice=wtr_prog_now%h_ice(:,jb),                    & !in
          & ps=p_diag%pres_sfc(:,jb), t_g=lnd_prog_now%t_g(:,jb),                      & !in
          & qv_s=lnd_diag%qv_s(:,jb),                                                  & !in
          & u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb),                & !in
          & T=p_diag%temp(:,:,jb), prs=p_diag%pres(:,:,jb),                            & !in
          & qv=p_prog_rcf%tracer(:,:,jb,iqv), qc=p_prog_rcf%tracer(:,:,jb,iqc),        & !in
          & gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm_t(:,jb,1),                        & !inout
          & tch=prm_diag%tch_t(:,jb,1), tfm=prm_diag%tfm(:,jb),                        & !inout
          & tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv_t(:,jb,1),                        & !inout
          & tke=z_tvs(:,:,:), tkvm=prm_diag%tkvm(:,2:nlevp1,jb),                       & !inout
          & tkvh=prm_diag%tkvh(:,2:nlevp1,jb), rcld=prm_diag%rcld(:,:,jb),             & !inout
          & t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb),                      & !out
          & td_2m=prm_diag%td_2m(:,jb), rh_2m=prm_diag%rh_2m(:,jb),                    & !out
          & u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m(:,jb),                    & !out
          & shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb),                & !out
          & ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine                    ) !inout


        ! copy
        prm_diag%tcm(i_startidx:i_endidx,jb) = prm_diag%tcm_t(i_startidx:i_endidx,jb,1)
        prm_diag%tch(i_startidx:i_endidx,jb) = prm_diag%tch_t(i_startidx:i_endidx,jb,1)
        prm_diag%tfv(i_startidx:i_endidx,jb) = prm_diag%tfv_t(i_startidx:i_endidx,jb,1)

      ELSE ! tile approach used

        ! preset variables for land tile indices
        fr_land_t(:)  = 1._wp
        depth_lk_t(:) = 0._wp
        h_ice_t(:)    = 0._wp

        ! Loop over index lists for land tile points, sea/lake points and seaice points
        DO  jt = 1, ntiles_total + nlists_water

          IF (jt <= ntiles_total) THEN ! land tile points
            jt1 = jt
            i_count =  ext_data%atm%gp_count_t(jb,jt)
            ilist   => ext_data%atm%idx_lst_t(:,jb,jt)
          ELSE IF (jt == ntiles_total + 1) THEN ! sea points (open water)
            jt1 = isub_water   !jt
            i_count =  ext_data%atm%spw_count(jb)
            ilist   => ext_data%atm%idx_lst_spw(:,jb)
            fr_land_t(:) = 0._wp
          ELSE IF (jt == ntiles_total + 2) THEN ! lake points
            jt1 = isub_water   !ntiles_total + 1
            i_count =  ext_data%atm%fp_count(jb)
            ilist   => ext_data%atm%idx_lst_fp(:,jb)
            fr_land_t(:)  = 0._wp
            depth_lk_t(:) = 1._wp
            ! h_ice_t(:) = ...
          ELSE ! IF (jt == ntiles_total + 3) THEN ! seaice points
            jt1 = isub_seaice    !ntiles_total + 2
            i_count =  ext_data%atm%spi_count(jb)
            ilist   => ext_data%atm%idx_lst_spi(:,jb)
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
            gz0_t(ic,jt)       = prm_diag%gz0_t(jc,jb,jt1)
            t_g_t(ic,jt)       = lnd_prog_now%t_g_t(jc,jb,jt1)
            qv_s_t(ic,jt)      = lnd_diag%qv_s_t(jc,jb,jt1)
            sai_t(ic,jt)       = ext_data%atm%sai_t(jc,jb,jt1)
            z_ifc_t(ic,1:3,jt) = p_metrics%z_ifc(jc,nlev-1:nlev+1,jb)
            pres_sfc_t(ic,jt)  = p_diag%pres_sfc(jc,jb)
!MR: Hauptflaechengroessen nur fuer level nlev
            u_t(ic,1:2,jt)     = p_diag%u(jc,nlev-1:nlev,jb)
            v_t(ic,1:2,jt)     = p_diag%v(jc,nlev-1:nlev,jb)
            w_t(ic,1:3,jt)     = p_prog%w(jc,nlev-1:nlev+1,jb)
            temp_t(ic,1:2,jt)  = p_diag%temp(jc,nlev-1:nlev,jb)
            pres_t(ic,1:2,jt)  = p_diag%pres(jc,nlev-1:nlev,jb)
            qv_t(ic,1:2,jt)    = p_prog_rcf%tracer(jc,nlev-1:nlev,jb,iqv)
            qc_t(ic,1:2,jt)    = p_prog_rcf%tracer(jc,nlev-1:nlev,jb,iqc)
!MR: tkvm, tkvh: benoetigt fuer nlev, nlev+1 als Nebenflaechengroessen (done), 
!MR: fuer nlev+1 sogar tile-spezifisch
!MR: tke: ausserhalb der Initialisierung nur fuer ke1 noetig, dort aber tilespezifisch
            tvs_t(ic,1:3,1,jt) = z_tvs(jc,nlev-1:nlev+1,1)
            tkvm_t(ic,1:2,jt)  = prm_diag%tkvm(jc,nlev:nlevp1,jb)
            tkvh_t(ic,1:2,jt)  = prm_diag%tkvh(jc,nlev:nlevp1,jb)
!MR: rcld: benoetigt nur fuer level nlev (als Nebenflaechenvariable)
            rcld_t(ic,1:3,jt)  = prm_diag%rcld(jc,nlev-1:nlev+1,jb)
          ENDDO

          CALL turbtran(iini=0, dt_tke=tcall_turb_jg, nprv=1, ntur=1, ntim=1,         & !in
            & ie=nproma, ke=2, ke1=3,                                                 & !in
            & istart=1, iend=i_count, istartpar=1, iendpar=i_count,                   & !in
            & l_hori=phy_params(jg)%mean_charlen, hhl=z_ifc_t(:,:,jt),                & !in
            & fr_land=fr_land_t(:), depth_lk=depth_lk_t(:),                           & !in
            & sai=sai_t(:,jt), h_ice=h_ice_t(:), ps=pres_sfc_t(:,jt),                 & !in
            & t_g=t_g_t(:,jt), qv_s=qv_s_t(:,jt),                                     & !in
            & u=u_t(:,:,jt), v=v_t(:,:,jt), w=w_t(:,:,jt),                            & !in
            & T=temp_t(:,:,jt), prs=pres_t(:,:,jt), qv=qv_t(:,:,jt), qc=qc_t(:,:,jt), & !in
            & gz0=gz0_t(:,jt), tcm=tcm_t(:,jt), tch=tch_t(:,jt),                      & !inout
            & tfm=tfm_t(:,jt), tfh=tfh_t(:,jt), tfv=tfv_t(:,jt),                      & !inout
            & tke=tvs_t(:,:,:,jt), tkvm=tkvm_t(:,:,jt),                               & !inout
            & tkvh=tkvh_t(:,:,jt), rcld=rcld_t(:,:,jt),                               & !inout
            & t_2m=t_2m_t(:,jt), qv_2m=qv_2m_t(:,jt), td_2m=td_2m_t(:,jt),            & !out
            & rh_2m=rh_2m_t(:,jt), u_10m=u_10m_t(:,jt), v_10m=v_10m_t(:,jt),          & !out
            & shfl_s=shfl_s_t(:,jt), lhfl_s=lhfl_s_t(:,jt),                           & !out
            & ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine                 ) !inout

        ENDDO

        ! Aggregate tile-based output fields of turbtran over tiles
        ! i) initialize fields to zero before starting the summation
        prm_diag%gz0   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tcm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tch   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfm   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfh   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%tfv   (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%t_2m  (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%qv_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%td_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%rh_2m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%u_10m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%v_10m (i_startidx:i_endidx,jb) = 0._wp
        prm_diag%shfl_s(i_startidx:i_endidx,jb) = 0._wp
        prm_diag%lhfl_s(i_startidx:i_endidx,jb) = 0._wp

        z_tvs        (i_startidx:i_endidx,nlevp1, 1) = 0._wp
        prm_diag%tkvm(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        prm_diag%tkvh(i_startidx:i_endidx,nlevp1,jb) = 0._wp
        prm_diag%rcld(i_startidx:i_endidx,nlevp1,jb) = 0._wp

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
            prm_diag%gz0(jc,jb) = prm_diag%gz0(jc,jb)+gz0_t(ic,jt) * area_frac
            prm_diag%tcm(jc,jb) = prm_diag%tcm(jc,jb)+tcm_t(ic,jt) * area_frac
            prm_diag%tch(jc,jb) = prm_diag%tch(jc,jb)+tch_t(ic,jt) * area_frac
            prm_diag%tfm(jc,jb) = prm_diag%tfm(jc,jb)+tfm_t(ic,jt) * area_frac
            prm_diag%tfh(jc,jb) = prm_diag%tfh(jc,jb)+tfh_t(ic,jt) * area_frac
            prm_diag%tfv(jc,jb) = prm_diag%tfv(jc,jb)+tfv_t(ic,jt) * area_frac

            z_tvs(jc,nlevp1,1)          = z_tvs(jc,nlevp1,1)+tvs_t(ic,3,1,jt)         * area_frac
            prm_diag%tkvm(jc,nlevp1,jb) = prm_diag%tkvm(jc,nlevp1,jb)+tkvm_t(ic,2,jt) * area_frac
            prm_diag%tkvh(jc,nlevp1,jb) = prm_diag%tkvh(jc,nlevp1,jb)+tkvh_t(ic,2,jt) * area_frac
!MR
         !  prm_diag%rcld(jc,nlevp1,jb) = prm_diag%rcld(jc,nlevp1,jb)+rcld_t(ic,3,jt)*area_frac
!MR

            prm_diag%t_2m(jc,jb)   = prm_diag%t_2m(jc,jb) +t_2m_t(ic,jt)  * area_frac
            prm_diag%qv_2m(jc,jb)  = prm_diag%qv_2m(jc,jb)+qv_2m_t(ic,jt) * area_frac
            prm_diag%td_2m(jc,jb)  = prm_diag%td_2m(jc,jb)+td_2m_t(ic,jt) * area_frac
            prm_diag%rh_2m(jc,jb)  = prm_diag%rh_2m(jc,jb)+rh_2m_t(ic,jt) * area_frac
            prm_diag%u_10m(jc,jb)  = prm_diag%u_10m(jc,jb)+u_10m_t(ic,jt) * area_frac
            prm_diag%v_10m(jc,jb)  = prm_diag%v_10m(jc,jb)+v_10m_t(ic,jt) * area_frac
            prm_diag%shfl_s(jc,jb) = prm_diag%shfl_s(jc,jb)+shfl_s_t(ic,jt) * area_frac
            prm_diag%lhfl_s(jc,jb) = prm_diag%lhfl_s(jc,jb)+lhfl_s_t(ic,jt) * area_frac

            prm_diag%shfl_s_t(jc,jb,jt1) = shfl_s_t(ic,jt)
            prm_diag%lhfl_s_t(jc,jb,jt1) = lhfl_s_t(ic,jt)
            prm_diag%tch_t(jc,jb,jt1)    = tch_t(ic,jt)    ! needed by TERRA
            prm_diag%tcm_t(jc,jb,jt1)    = tcm_t(ic,jt)    ! needed by TERRA
            prm_diag%tfv_t(jc,jb,jt1)    = tfv_t(ic,jt)    ! needed by TERRA
            prm_diag%gz0_t(jc,jb,jt1)    = gz0_t(ic,jt)    ! input for turbtran at n+1
                                                           ! over land

!MR: tke, tkvm, tkvh tilespezifisch speichern
!DR should tile specific u10m, v10m fields be stored for TERRA?

          ENDDO

        ENDDO

      ENDIF ! tiles / no tiles

      ! transform updated turbulent velocity scale back to TKE
      p_prog_now_rcf%tke(i_startidx:i_endidx,nlevp1,jb)= 0.5_wp                        &
        &                                     * (z_tvs(i_startidx:i_endidx,nlevp1,1))**2


    ELSE IF ( atm_phy_nwp_config(jg)%inwp_turb == 2 ) THEN

!-------------------------------------------------------------------------
!> GME turbulence scheme 
!-------------------------------------------------------------------------

      ! turbulent diffusion coefficients at the surface
      CALL parturs( zsurf=p_metrics%z_ifc(:,nlevp1,jb), z1=p_metrics%z_mc(:,nlev,jb),   & !in
        &           u1=p_diag%u(:,nlev,jb), v1=p_diag%v(:,nlev,jb),                     & !in
        &           t1=p_diag%temp(:,nlev,jb), qv1=p_prog_rcf%tracer(:,nlev,jb,iqv),    & !in
        &           t_g=lnd_prog_now%t_g(:,jb), qv_s=lnd_diag%qv_s(:,jb),               & !in
        &           fr_land=ext_data%atm%fr_land(:,jb), h_ice=wtr_prog_now%h_ice(:,jb), & !in
        &           ie=nproma, i_startidx=i_startidx, i_endidx=i_endidx,                & !in
        &           tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb),                     & !out
        &           gz0=prm_diag%gz0(:,jb)                                              ) !inout


      ! Copy transfer coefficients to tile-based variables, which are used in TERRA
      DO jt = 1, ntiles_total+ntiles_water
        DO jc = i_startidx, i_endidx
          prm_diag%tcm_t(jc,jb,jt) = prm_diag%tcm(jc,jb)
          prm_diag%tch_t(jc,jb,jt) = prm_diag%tch(jc,jb)
          ! the GME turbulence scheme does not have tfv. Set tfv=1
!          prm_diag%tfv_t(jc,jb,jt) = 1._wp    !   prm_diag%tfv(jc,jb)
          prm_diag%tfv_t(jc,jb,jt) = prm_diag%tfv(jc,jb)
        ENDDO
      ENDDO

    ENDIF !inwp_turb

  ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE nwp_turbtrans

END MODULE mo_nwp_turbtrans_interface
