#ifdef __xlC__
@PROCESS SPILL(4096)
#endif
!OPTION! -cont -msg o
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! turbulence parameterisations that include the call to the surface scheme
!!
!! inwp_turb == 3 == EDMF DUAL turbulence scheme by Koehler and Neggers from IFS
!! inwp_turb == 4 == turbulence scheme by Brinkop and Roeckner run in ECHAM
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

MODULE mo_nwp_turb_sfc_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_physical_constants,   ONLY: alv, rd_o_cpd, grav
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,        ONLY: phy_params 
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_lnd_diag
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, iqr, iqs, iqtvar, nqtendphy
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_data_turbdiff,        ONLY: t0_melt, zt_ice
  USE mo_satad,                ONLY: sat_pres_water, spec_humi  
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_vdiff_config,         ONLY: vdiff_config
  USE mo_vdiff_driver,         ONLY: vdiff
  USE mo_vdfouter,             ONLY: vdfouter
  USE mo_run_config,           ONLY: ltestcase
  USE mo_nh_testcases,         ONLY: nh_test_name
  USE mo_nh_wk_exp,            ONLY: qv_max_wk
  USE mo_lnd_nwp_config,       ONLY: nlev_soil, nlev_snow, ntiles_total, lmulti_snow
  USE mo_sync,                 ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, &
                                     SYNC_C, SYNC_C1, global_max, global_min, global_sum_array
  USE mo_run_config,           ONLY: ntracer

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_turbulence_sfc

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
SUBROUTINE nwp_turbulence_sfc ( tcall_turb_jg,                     & !>input
                              & p_patch,p_metrics,                 & !>input
                              & ext_data,                          & !>input
                              & p_prog,                            & !>inout
                              & p_prog_now_rcf, p_prog_rcf,        & !>in/inout
                              & p_diag ,                           & !>inout
                              & prm_diag, prm_nwp_tend,            & !>inout 
                              & lnd_prog_now, lnd_prog_new,        & !>inout 
                              & lnd_diag                           ) !>inout


  TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
  TYPE(t_external_data),       INTENT(inout):: ext_data        !< external data
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

  INTEGER :: jc,jk,jb,jt,jtile,jg    !loop indices

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

  ! local variables for edmf

  INTEGER  :: icnt
  INTEGER, PARAMETER :: itrac_vdf = 0
  INTEGER, PARAMETER :: ntiles_edmf = 8
  INTEGER  :: KTVL(nproma)         , KTVH(nproma)         , zsoty(nproma)        , & 
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
  REAL(wp) :: zalbti(nproma,ntiles_edmf)
  REAL(wp) :: zdummy_vdf_6a(nproma,p_patch%nlev,itrac_vdf), &
    &         zdummy_vdf_6b(nproma,p_patch%nlev,itrac_vdf), &
    &         zdummy_vdf_6c(nproma,itrac_vdf)
  REAL(wp) :: zdummy_vdf_7a(nproma,0), zdummy_vdf_7b(nproma,p_patch%nlev,0)

  REAL(wp) :: z_omega_p(nproma,p_patch%nlev), zchar(nproma),                    &
    &         zucurr(nproma)                , zvcurr(nproma),                   &
    &         zsoteu(nproma,p_patch%nlev)   , zsotev(nproma,p_patch%nlev),      &
    &         zsobeta(nproma,p_patch%nlev)  , zz0h(nproma),                     &
    &         shfl_s_t(nproma,ntiles_total) , evap_s_t(nproma,ntiles_total),    &
    &         tskin_t(nproma,ntiles_total)  , &
    &         ustr_s_t(nproma,ntiles_total) , vstr_s_t(nproma,ntiles_total),    &
    &         zae(nproma,p_patch%nlev)      ,                                   &
    &         ztice(nproma)                 , ztske1(nproma),                   &
    &         ztskm1m(nproma)               , ztskrad(nproma),                  &
    &         zsigflt(nproma)               , zfrti(nproma,ntiles_edmf),        &
    &         tch_ex(nproma,ntiles_total)   , tcm_ex(nproma,ntiles_total),      &
    &         tfv_ex(nproma,ntiles_total)
  LOGICAL  :: ldummy_vdf_a(nproma)

  INTEGER, SAVE :: nstep_turb = 0

! estimates of tile albedo ???
  REAL(wp) :: zalbti_est(ntiles_edmf)
  DATA        zalbti_est  / 0.06, 0.80, 0.20, 0.20, 0.70, 0.15, 0.35, 0.40 /
! conversion of soil types TERRA to TESSEL
  REAL(wp) :: soiltyp_conv(8)
  DATA        soiltyp_conv /0,0,1,2,3,4,5,6/
! conversion of vegetation types TERRA (Tab11.5) to TESSEL (cy36r1 Tab8.1) 
  REAL(wp) :: vegtyp_conv(23)
  DATA        vegtyp_conv / 6, 5,19, 3, 4,  18,20,18,18,11, &
                         & 16,17,13,11, 7,   1, 1,18, 8,15, &
                         & 12, 8,15 /

! Globcover 2009: tile index used in TESSEL (IFS) (see also mo_ext_data_state.f90)
!   Number of landcover classes provided by external parameter data
!   Needs to be changed into a variable if landcover classifications 
!   with a different number of classes become available
  INTEGER,  PARAMETER          :: num_lcc = 23
  REAL(wp), DIMENSION(num_lcc) :: jtessel_gcv2009  ! Tessel index table GlobCover2009
!                      jtessel 
  DATA jtessel_gcv2009 / 4,   & !  1 irrigated croplands                           
                     &   4,   & !  2 rainfed croplands                             
                     &   4,   & !  3 mosaic cropland (50-70%) - vegetation (20-50%)
                     &   4,   & !  4 mosaic vegetation (50-70%) - cropland (20-50%)
                     &   6,   & !  5 closed broadleaved evergreen forest           
                     &   6,   & !  6 closed broadleaved deciduous forest           
                     &   6,   & !  7 open broadleaved deciduous forest             
                     &   6,   & !  8 closed needleleaved evergreen forest          
                     &   6,   & !  9 open needleleaved deciduous forest            
                     &   6,   & ! 10 mixed broadleaved and needleleaved forest     
                     &   4,   & ! 11 mosaic shrubland (50-70%) - grassland (20-50%)
                     &   4,   & ! 12 mosaic grassland (50-70%) - shrubland (20-50%)
                     &   4,   & ! 13 closed to open shrubland                      
                     &   4,   & ! 14 closed to open herbaceous vegetation          
                     &   4,   & ! 15 sparse vegetation                             
                     &   6,   & ! 16 closed to open forest regulary flooded        
                     &   6,   & ! 17 closed forest or shrubland permanently flooded
                     &   4,   & ! 18 closed to open grassland regularly flooded    
                     &   8,   & ! 19 artificial surfaces                           
                     &   8,   & ! 20 bare areas                                    
                     &   1,   & ! 21 water bodies                                  
                     &   5,   & ! 22 permanent snow and ice                        
                     &   8    / ! 23 undefined                                  

!--------------------------------------------------------------

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


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,jk,i_startidx,i_endidx, icnt, &
!$OMP KTVL, KTVH, zsoty, & 
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
!$OMP zalbti, &
!$OMP zdummy_vdf_6a, &
!$OMP zdummy_vdf_6b, &
!$OMP zdummy_vdf_6c, &
!$OMP zdummy_vdf_7a, zdummy_vdf_7b, &
!$OMP z_omega_p, zchar, &
!$OMP zucurr,    zvcurr, &
!$OMP zsoteu,    zsotev, &
!$OMP zsobeta,   zz0h, &
!$OMP shfl_s_t, &
!$OMP evap_s_t,  tskin_t, &
!$OMP ustr_s_t,  vstr_s_t, &
!$OMP zae, &
!$OMP ztice,     ztske1, &
!$OMP ztskm1m,   ztskrad, &
!$OMP zsigflt,   zfrti, &
!$OMP tch_ex,    tcm_ex,  tfv_ex, &
!$OMP ldummy_vdf_a &
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
        ENDIF
      ENDIF
    ENDIF


    IF ( atm_phy_nwp_config(jg)%inwp_turb == 4 ) THEN

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
      ENDIF

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

!     TERRA variable initializtion

      DO jt = 1,ntiles_total
        DO jc = i_startidx, i_endidx
          lnd_prog_new%t_snow_t  (jc,jb,jt) = lnd_prog_now%t_snow_t  (jc,jb,jt) 
          lnd_prog_new%t_s_t     (jc,jb,jt) = lnd_prog_now%t_s_t     (jc,jb,jt)   
          lnd_prog_new%t_g_t     (jc,jb,jt) = lnd_prog_now%t_g_t     (jc,jb,jt)
          lnd_prog_new%w_snow_t  (jc,jb,jt) = lnd_prog_now%w_snow_t  (jc,jb,jt)    
          lnd_prog_new%rho_snow_t(jc,jb,jt) = lnd_prog_now%rho_snow_t(jc,jb,jt)  
          lnd_prog_new%w_i_t     (jc,jb,jt) = lnd_prog_now%w_i_t     (jc,jb,jt)       
          lnd_prog_new%t_so_t(jc,nlev_soil+2,jb,jt) = lnd_prog_now%t_so_t(jc,nlev_soil+2,jb,jt)
          IF(lmulti_snow) THEN
            lnd_prog_new%t_snow_mult_t(jc,nlev_snow+1,jb,jt) &
                                             = lnd_prog_now%t_snow_mult_t(jc,nlev_snow+1,jb,jt)
          ENDIF
        ENDDO
      ENDDO
      DO jt = 1,ntiles_total
        DO jc = i_startidx, i_endidx
          IF(lmulti_snow) THEN
            DO jk=1,nlev_snow
              lnd_prog_new%t_snow_mult_t  (jc,jk,jb,jt) = lnd_prog_now%t_snow_mult_t  (jc,jk,jb,jt)
              lnd_prog_new%rho_snow_mult_t(jc,jk,jb,jt) = lnd_prog_now%rho_snow_mult_t(jc,jk,jb,jt)
              lnd_prog_new%wliq_snow_t    (jc,jk,jb,jt) = lnd_prog_now%wliq_snow_t    (jc,jk,jb,jt) 
              lnd_prog_new%wtot_snow_t    (jc,jk,jb,jt) = lnd_prog_now%wtot_snow_t    (jc,jk,jb,jt)
              lnd_prog_new%dzh_snow_t     (jc,jk,jb,jt) = lnd_prog_now%dzh_snow_t     (jc,jk,jb,jt) 
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO jt = 1,ntiles_total
        DO jc = i_startidx, i_endidx
          DO jk=1,nlev_soil+1
            lnd_prog_new%t_so_t    (jc,jk,jb,jt) = lnd_prog_now%t_so_t    (jc,jk,jb,jt)         
            lnd_prog_new%w_so_t    (jc,jk,jb,jt) = lnd_prog_now%w_so_t    (jc,jk,jb,jt)    
            lnd_prog_new%w_so_ice_t(jc,jk,jb,jt) = lnd_prog_now%w_so_ice_t(jc,jk,jb,jt)
          ENDDO
        ENDDO
      ENDDO

!     Various variables for VDFOUTER
      DO jc = i_startidx, i_endidx
        zchar  (jc) = 0.018_wp ! default value from IFS if no wave model
        zucurr (jc) = 0.0_wp
        zvcurr (jc) = 0.0_wp
        ztice  (jc) = 273.0_wp ! top level ice temperature ???????
        ztske1 (jc) = 0.0_wp   ! skin temperature tendency
        ztskm1m(jc) = lnd_prog_now%t_g (jc,jb) ! skin temperature (prognostic ???)
        ztskrad(jc) = lnd_prog_now%t_g (jc,jb) ! skin temperature at last radiation step ????
        zsigflt(jc) = 0.0_wp   ! just for testing (standard dev. of filtered orogrphy)
        zz0h   (jc) = 0.0_wp   ! diagnostic z0,h - should be in diagnostic output ???
if (ztskm1m(jc) > 400.0 ) then
  write(*,*) 'turb_sfc1: ', ztskm1m(jc)
endif
      ENDDO

      DO jk = 1,p_patch%nlev
        DO jc = i_startidx, i_endidx
          zsoteu (jc,jk) = 0.0_wp
          zsotev (jc,jk) = 0.0_wp
          zsobeta(jc,jk) = 0.0_wp
          zae    (jc,jk) = 0.0_wp   ! cloud tendency ???
        ENDDO
      ENDDO

      DO jt = 1,ntiles_total
        DO jc = i_startidx, i_endidx
          shfl_s_t(jc,jt) = prm_diag%shfl_s  (jc,jb)   ! should be tile specific !!!
          evap_s_t(jc,jt) = prm_diag%lhfl_s  (jc,jb) / alv ! evaporation [kg/(m2 s)]  -"-
          tskin_t (jc,jt) = lnd_prog_now%t_g (jc,jb)   ! should be tile specific, and prognostic !!!
          ustr_s_t(jc,jt) = 0.0_wp                     ! prognostic surface stress U  !!!
          vstr_s_t(jc,jt) = 0.0_wp                     ! prognostic surface stress V  !!!
        ENDDO
      ENDDO

      DO jt = 1,ntiles_edmf
        DO jc = i_startidx, i_endidx
          zfrti   (jc,jt) = 0.0_wp                     ! all zero but tile1=1.0 ... all ocean ???
          zalbti  (jc,jt) = zalbti_est(jt)             ! surface albedo ????????
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        KTVL(jc) = 16                                  ! KTVL: dummy default ???
        KTVH(jc) = 3                                   ! KTVH: dummy default ???
      ENDDO

      DO jt = ntiles_total, 1, -1
        DO jc = i_startidx, i_endidx
          IF ( ext_data%atm%llsm_atm_c(jc,jb) ) THEN ! land
            JTILE =                 jtessel_gcv2009(ext_data%atm%lc_class_t(jc,jb,jt))
!JTILE=8 !?????????????
!???            IF (JTILE == 4)  KTVL(jc) = vegtyp_conv(ext_data%atm%lc_class_t(jc,jb,jt))
!???            IF (JTILE == 6)  KTVH(jc) = vegtyp_conv(ext_data%atm%lc_class_t(jc,jb,jt))
!write(*,*) 'turb_sfc1: ',  JTILE, KTVL(jc), KTVH(jc), ext_data%atm%lc_class_t(jc,jb,jt)
            IF ( lnd_diag%snowfrac_lc_t(jc,jb,jt) > 0.5_wp ) THEN
              SELECT CASE ( JTILE )
                CASE (4)
                  JTILE = 5   ! snow over low vegetation
                CASE (6)
                  JTILE = 7   ! snow over high vegetation
                CASE (8)
                  JTILE = 5   ! snow over bare ground
              END SELECT
            ENDIF
           !IF ( t_g(jl) > (t0_melt + zt_ice) ) THEN   ! salt water freezing temperature
           !  JTILE = 2              ! sea ice         ! this may never be active as not LAND
           !ENDIF
           !interception layer (#3) missing ???
            zfrti(jc,jtile) = zfrti(jc,jtile) + ext_data%atm%lc_frac_t(jc,jb,jt)
          ENDIF
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
        IF ( .NOT. ext_data%atm%llsm_atm_c(jc,jb) ) THEN ! ocean or sea-ice point
          zfrti(jc,:) = 0.0_wp                           ! all zero
          zfrti(jc,1) = 1.0_wp                           ! ocean one
!ELSE
!  zfrti(jc,:) = 0.0_wp
!  zfrti(jc,7) = 1.0_wp  !8,3,4,5 work, 6?
        ENDIF
       !ELSE !!! IF ( lnd_prog_now%t_g(jc,jb) > (t0_melt + zt_ice) ) THEN  ! salt water freezing temperature
       !  zfrti(jc,1) = 1.0_wp                           ! open ocean
       !!!!ELSE
       !!!!  zfrti(jc,2) = 1.0_wp                        ! sea ice (requires sea ice cover)
      ENDDO

      DO jc = i_startidx, i_endidx
        IF ( (SUM(zfrti(jc,:)) > 1.01_wp) .or. (SUM(zfrti(jc,:)) < 0.99_wp) ) THEN
          write(*,*) 'turb_sfc2: ', zfrti(jc,:)
        ENDIF
      ENDDO

      DO jk = 1,nlev_soil
        DO jc = i_startidx, i_endidx
          zdummy_vdf_4a(jc,jk) = lnd_prog_now%t_so_t(jc,jk,jb,1) ! simple: take one tile #1 ???
          zdummy_vdf_4b(jc,jk) = lnd_prog_now%w_so_t(jc,jk,jb,1) ! ---
        ENDDO
      ENDDO

      DO jc = i_startidx, i_endidx
       !KTVL(jc) = 16    !KTVL  ???
       !KTVH(jc) = 3     !KTVH  ???
        zdummy_vdf_1a(jc) = 0.5_wp!PCVL  ???
        zdummy_vdf_1b(jc) = 0.5_wp!PCVH  ???
        zsoty(jc)         = soiltyp_conv(ext_data%atm%soiltyp_t(jc,jb,1)) !KSOTY ???
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
        & CDCONF  = 'T'                                        ,&! (IN)  used in surfexcdriver
        & KIDIA   = i_startidx                                 ,&! (IN)   
        & KFDIA   = i_endidx                                   ,&! (IN)   
        & KLON    = nproma                                     ,&! (IN)   
        & KLEV    = p_patch%nlev                               ,&! (IN)   
        & KLEVS   = nlev_soil                                  ,&! (IN)
        & KSTEP   = nstep_turb                                 ,&! (IN)  ??? used in surfexdriver!!
        & KTILES  = ntiles_edmf                                ,&! (IN)
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
        & KTVL    = KTVL                                       ,&! (IN)  input for TESSEL   
        & KTVH    = KTVH                                       ,&! (IN)  input for TESSEL
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
        & PTSKM1M = lnd_prog_now%t_g(:,jb)                     ,&! (IN)  T,skin
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
        & PSST    = ztskm1m                                    ,&! (IN)  SST
        & KSOTY   = zsoty                                      ,&! (IN)  soil type
!xmk ?  & PFRTI   = ext_data%atm%lc_frac_t(:,jb,:)             ,&! (IN)  tile fraction 
        & PFRTI   = zfrti                                      ,&! (IN)  tile fraction 
        & PALBTI  = zalbti                                     ,&! (IN)  tile albedo
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
        & PZ0M    = prm_diag%z0m(:,jb)                         ,&! (INOUT) z0,m (calculated in vupdz0) 
        & PZ0H    = zz0h                                       ,&! (INOUT) z0,h (should be diagnostic output ???)    
        & PVDIS   = zdummy_vdf_1g                              ,&! (OUT) optional out: turbulent dissipation
        & PVDISG  = zdummy_vdf_1h                              ,&! (OUT) optional out: SO dissipation
        & PDISGW3D= zdummy_vdf_2a                              ,&! (OUT) optional out: 3D stoch. phys. dissipation
        & PAHFLEV = zdummy_vdf_1i                              ,&! (OUT) optional out: latent heat flux (snow/ice free part)
        & PAHFLSB = zdummy_vdf_1j                              ,&! (OUT) optional out: latent heat flux (snow/ice covered part)
        & PFWSB   = zdummy_vdf_1k                              ,&! (OUT) optional out: evaporation of snow
        & PBIR    = zdummy_vdf_1l                              ,&! (OUT) optional out: BIR buoyancy flux integral ratio
        & PVAR    = p_prog_rcf%tracer(:,:,jb,iqtvar)           ,&! (INOUT) qt,variance - prognostic advected tracer
        & PU10M   = prm_diag%u_10m(:,jb)                       ,&! (OUT) need to be initialized???
        & PV10M   = prm_diag%v_10m(:,jb)                       ,&! (OUT)  "-"
        & PT2M    = prm_diag%t_2m (:,jb)                       ,&! (OUT)  "-"
        & PD2M    = prm_diag%td_2m(:,jb)                       ,&! (OUT)  "-" 
        & PQ2M    = prm_diag%qv_2m(:,jb)                       ,&! (OUT)  "-" 
        & PZINV   = zdummy_vdf_1m                              ,&! (OUT) optional out: PBL HEIGHT (moist parcel, not for stable PBL)
        & PBLH    = zdummy_vdf_1n                              ,&! (OUT) optional out: PBL HEIGHT (dry diagnostic based on Ri#)
        & KHPBLN  = idummy_vdf_0d                              ,&! (OUT) optional out: PBL top level 
        & KVARTOP = idummy_vdf_0e                              ,&! (OUT) optional out: top level of predictied qt,var
        & PSSRFLTI= prm_diag%swflxsfc_t(:,jb,:)                ,&! (INOUT) net SW sfc flux for each tile (use tile ablbedo!!!)
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
        & LDLAND  = ext_data%atm%llsm_atm_c(:,jb)               &! (IN)   logical for land
!       & PDHTLS  = ...                                        ,&! (OUT)  optional out: DDH
!       & PDHTSS  = ...                                        ,&! (OUT)  optional out: DDH
!       & PDHTTS  = ...                                        ,&! (OUT)  optional out: DDH
!       & PDHTIS  = ...                                         &! (OUT)  optional out: DDH
! TERRA data
        & , ext_data        = ext_data                            & !in
        & , jb=jb, jg=jg                                          & ! -
        & , t_snow_ex       = lnd_prog_new%t_snow_t    (:,jb,:)   & !inout
        & , t_snow_mult_ex  = lnd_prog_new%t_snow_mult_t(:,:,jb,:)& ! -
        & , t_s_ex          = lnd_prog_new%t_s_t       (:,jb,:)   & ! -
        & , t_g_ex          = lnd_prog_new%t_g_t       (:,jb,:)   & ! - 
        & , qv_s_ex         = lnd_diag%qv_s_t          (:,jb,:)   & ! -
        & , w_snow_ex       = lnd_prog_new%w_snow_t    (:,jb,:)   & ! -
        & , w_snow_eff_ex   = lnd_diag%w_snow_eff_t    (:,jb,:)   & ! -
        & , rho_snow_ex     = lnd_prog_new%rho_snow_t  (:,jb,:)   & ! -
        & , rho_snow_mult_ex= lnd_prog_new%rho_snow_mult_t(:,:,jb,:) & ! -
        & , h_snow_ex       = lnd_diag%h_snow_t        (:,jb,:)   & ! -
        & , w_i_ex          = lnd_prog_new%w_i_t       (:,jb,:)   & ! -
        & , t_so_ex         = lnd_prog_new%t_so_t      (:,:,jb,:) & ! -
        & , w_so_ex         = lnd_prog_new%w_so_t      (:,:,jb,:) & ! -
        & , w_so_ice_ex     = lnd_prog_new%w_so_ice_t  (:,:,jb,:) & ! -
!       & , t_2m_ex         = prm_diag%t_2m            (:,jb)     & ! -
!       & , u_10m_ex        = prm_diag%u_10m           (:,jb)     & ! - 
!       & , v_10m_ex        = prm_diag%v_10m           (:,jb)     & ! -
        & , freshsnow_ex    = lnd_diag%freshsnow_t     (:,jb,:)   & ! -
        & , snowfrac_lc_ex  = lnd_diag%snowfrac_lc_t   (:,jb,:)   & ! -
        & , snowfrac_ex     = lnd_diag%snowfrac_t      (:,jb,:)   & ! -
        & , wliq_snow_ex    = lnd_prog_new%wliq_snow_t (:,:,jb,:) & ! -
        & , wtot_snow_ex    = lnd_prog_new%wtot_snow_t (:,:,jb,:) & ! -
        & , dzh_snow_ex     = lnd_prog_new%dzh_snow_t  (:,:,jb,:) & ! -
        & , prr_con_ex      = prm_diag%tracer_rate     (:,jb,3)   & !in
        & , prs_con_ex      = prm_diag%tracer_rate     (:,jb,4)   & ! -
        & , prr_gsp_ex      = prm_diag%tracer_rate     (:,jb,1)   & ! -
        & , prs_gsp_ex      = prm_diag%tracer_rate     (:,jb,2)   & ! -
        & , tch_ex          = tch_ex                   (:,:)      & !inout
        & , tcm_ex          = tcm_ex                   (:,:)      & ! -
        & , tfv_ex          = tfv_ex                   (:,:)      & ! -
        & , sobs_ex         = prm_diag%swflxsfc_t      (:,jb,:)   & !in
        & , thbs_ex         = prm_diag%lwflxsfc_t      (:,jb,:)   & ! -
        & , pabs_ex         = prm_diag%swflxsfc_t      (:,jb,:)   & ! -
        & , runoff_s_ex     = lnd_diag%runoff_s_t      (:,jb,:)   & !inout
        & , runoff_g_ex     = lnd_diag%runoff_g_t      (:,jb,:)   & ! -
        & , t_g             = lnd_prog_new%t_g         (:,jb)     & ! -
        & , qv_s            = lnd_diag%qv_s            (:,jb)     ) ! -


! Turbulence updating strategy:
! * Update T, QV, QC, QI with turbulence tendencies
! * Give U, V tendencies to dynamics

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          p_prog_rcf%tracer(jc,jk,jb,iqv) = MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqv) &
                        & + tcall_turb_jg * prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv))
          p_prog_rcf%tracer(jc,jk,jb,iqc) = MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqc) &
                        & + tcall_turb_jg * prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqc))
          p_prog_rcf%tracer(jc,jk,jb,iqi) = MAX(0._wp, p_prog_rcf%tracer(jc,jk,jb,iqi) &
                        & + tcall_turb_jg * prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqi))
          p_diag%temp(jc,jk,jb) =                          p_diag%temp(jc,jk,jb) &
                        & + tcall_turb_jg * prm_nwp_tend%ddt_temp_turb(jc,jk,jb)
        ENDDO
      ENDDO

! Diagnostic output variables:  ???? TERRA-tiles or TESSEL-tiles???

      prm_diag%rh_2m (:,:) = 0.0_wp
      prm_diag%shfl_s(:,:) = 0.0_wp
      prm_diag%lhfl_s(:,:) = 0.0_wp
      prm_diag%tch   (:,:) = 0.0_wp
      prm_diag%tcm   (:,:) = 0.0_wp
      prm_diag%tfv   (:,:) = 0.0_wp
      DO jt = 1, ntiles_total
        DO jc = i_startidx, i_endidx
          prm_diag%rh_2m (jc,jb) = prm_diag%rh_2m (jc,jb) + 0.0_wp  !???

        ! prm_diag%shfl_s(jc,jb) = prm_diag%shfl_s(jc,jb) + shfl_s_t(jc,jt)    * ext_data%atm%lc_frac_t(jc,jb,jt)
        ! prm_diag%lhfl_s(jc,jb) = prm_diag%lhfl_s(jc,jb) + evap_s_t(jc,jt)*alv* ext_data%atm%lc_frac_t(jc,jb,jt)
          prm_diag%shfl_s(jc,jb) = prm_diag%shfl_s(jc,jb) + shfl_s_t(jc,jt)    * zfrti(jc,jt)   !bad, but needed for EDMF???
          prm_diag%lhfl_s(jc,jb) = prm_diag%lhfl_s(jc,jb) + evap_s_t(jc,jt)*alv* zfrti(jc,jt)   ! ----

          prm_diag%tch   (jc,jb) = prm_diag%tch   (jc,jb) + tch_ex  (jc,jt)    * ext_data%atm%lc_frac_t(jc,jb,jt) !??land only??
          prm_diag%tcm   (jc,jb) = prm_diag%tcm   (jc,jb) + tcm_ex  (jc,jt)    * ext_data%atm%lc_frac_t(jc,jb,jt)
          prm_diag%tfv   (jc,jb) = prm_diag%tfv   (jc,jb) + tfv_ex  (jc,jt)    * ext_data%atm%lc_frac_t(jc,jb,jt)
        ENDDO           
      ENDDO

    ENDIF !inwp_turb

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

nstep_turb = nstep_turb + 1

CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=p_prog_rcf%tracer, &
                          lpart4d=.TRUE.)
CALL sync_patch_array_mult(SYNC_C, p_patch, 2, p_diag%tempv, p_diag%exner_old)
CALL sync_patch_array     (SYNC_C, p_patch, p_diag%tempv)
CALL sync_patch_array_mult(SYNC_C1,p_patch, 2, &
                          prm_nwp_tend%ddt_u_turb, prm_nwp_tend%ddt_v_turb)
CALL sync_patch_array_mult(SYNC_C1,p_patch, 2, prm_nwp_tend%ddt_u_turb, &
                          prm_nwp_tend%ddt_v_turb)

CALL sync_patch_array(SYNC_C, p_patch, prm_diag%rh_2m)
CALL sync_patch_array(SYNC_C, p_patch, prm_diag%shfl_s)
CALL sync_patch_array(SYNC_C, p_patch, prm_diag%lhfl_s)
CALL sync_patch_array(SYNC_C, p_patch, prm_diag%tch)
CALL sync_patch_array(SYNC_C, p_patch, prm_diag%tcm)
CALL sync_patch_array(SYNC_C, p_patch, prm_diag%tfv)

END SUBROUTINE nwp_turbulence_sfc

END MODULE mo_nwp_turb_sfc_interface
