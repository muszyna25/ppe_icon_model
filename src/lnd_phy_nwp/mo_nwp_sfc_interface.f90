!>
!! This module is the interface between nwp_nh_interface to the 
!! surface parameterisations:
!! inwp_sfc  == 1 == surface scheme TERRA run in COSMO
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
MODULE mo_nwp_sfc_interface

  USE mo_kind,                ONLY: wp
!  USE mo_exception,           ONLY: message, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int, zml_soil
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_ext_data,            ONLY: t_external_data !DR, nclass_lu
  USE mo_nonhydro_state,      ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_state,       ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_state,       ONLY: t_lnd_prog, t_lnd_diag !!$, t_tiles
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv !,msg_level
  USe mo_extpar_config,       ONLY: itopo
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, nztlev, nlev_snow, nsfc_subs !, &
!    &                               lseaice ,llake, lmulti_snow
  USE mo_satad,               ONLY: sat_pres_water, spec_humi  
  USE mo_soil_ml,             ONLY: terra_multlay
!  USE mo_aggregate_surface,   ONLY: subsmean,subs_disaggregate_radflux,subsmean_albedo
!  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, igbm, iwtr, iice, ilnd
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_surface, nwp_surface_init


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_surface( tcall_sfc_jg, jstep,            & !>in
                        & p_patch,                        & !>in
                        & ext_data,                       & !>in
                        & p_prog_rcf,                     & !>in/inout
                        & p_diag ,                        & !>inout
                        & prm_diag,                       & !>inout 
                        & lnd_prog_now, lnd_prog_new,     & !>inout
                        & lnd_diag                        ) !>inout

    INTEGER ,INTENT(in)          :: jstep
    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch       !<grid/patch info
    TYPE(t_external_data),       INTENT(in)   :: ext_data      !< external data
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf    !<call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag        !< diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag      !< atm phys vars
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_now  !< prog vars for sfc
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog_new  !< prog vars for sfc
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag      !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: tcall_sfc_jg  !< time interval for 
                                                               !< surface

    ! Local array bounds:
    !
    INTEGER :: nblks_c                 !> number of blocks for cells
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: nlev                    !< number of full and half levels
    INTEGER :: isubs !DR, n_lu            

    ! Local scalars:
    !
    INTEGER :: jc,jb,jg      !loop indices
 
    ! local prognostic variables
    !
    REAL(wp) :: t_snow_t       (nproma,               p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: t_snow_mult_t  (nproma, 0:nlev_snow  ,p_patch%nblks_c, nztlev, nsfc_subs) 
    REAL(wp) :: t_so_t         (nproma, 0:nlev_soil+1,p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: w_so_t         (nproma, 1:nlev_soil+1,p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: w_so_ice_t     (nproma, 1:nlev_soil+1,p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: w_i_t          (nproma,               p_patch%nblks_c, nztlev, nsfc_subs) 
    REAL(wp) :: t_s_t          (nproma,               p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: w_snow_t       (nproma,               p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: rho_snow_t     (nproma,               p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: rho_snow_mult_t(nproma, 1:nlev_snow  ,p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: wliq_snow_t    (nproma, 1:nlev_snow  ,p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: wtot_snow_t    (nproma, 1:nlev_snow  ,p_patch%nblks_c, nztlev, nsfc_subs)
    REAL(wp) :: dzh_snow_t     (nproma, 1:nlev_snow  ,p_patch%nblks_c, nztlev, nsfc_subs)


    ! local diagnostic variables
    !
    REAL(wp) :: h_snow_t(nproma, p_patch%nblks_c, nztlev, nsfc_subs) 

!DR    REAL(wp) :: pp_t(nproma, p_patch%nlev, p_patch%nblks_c, nztlev)

    REAL(wp) :: tch_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: tcm_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: tfv_t     (nproma, p_patch%nblks_c, nsfc_subs)
    INTEGER  :: soiltyp_t (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: plcov_t   (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: rootdp_t  (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: sai_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: tai_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: eai_t     (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: landmask_t(nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: t_2m_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: u_10m_t   (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: v_10m_t   (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: sobs_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: thbs_t    (nproma, p_patch%nblks_c, nsfc_subs)
    REAL(wp) :: pabs_t    (nproma, p_patch%nblks_c, nsfc_subs)

!DR    REAL(wp) :: lu_class_frac(nproma,p_patch%nblks_c,nclass_lu)

!--------------------------------------------------------------

    ! local variables related to the blocking

    nblks_c   = p_patch%nblks_int_c

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev

    ! in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    write(0,*) "SFC-DIAGNOSIS INTERFACE ",jstep



!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,isubs), SCHEDULE(guided)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      IF( atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
        ! check dry case
        IF( atm_phy_nwp_config(jg)%inwp_satad == 0) THEN
          lnd_diag%qv_s (:,:) = 0._wp
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

      IF (  atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN
          
!#ifdef __BOUNDCHECK

!Example for Debugging, Diagnosis
!!$    IF (msg_level >= 0) THEN
!!$      WRITE(message_text,'(a,3E15.7)') ' Soil Temperature = ', &
!!$           &  lnd_prog_now%t_so(1,:,:,:,12)
!!$     CALL message('', TRIM(message_text))
!!$    ENDIF

!!$    nlev_soil       = 7     ! 7 = default value for number of soil layers
!!$    nztlev          = 2     ! 2 = default value for time integration scheme
!!$    nlev_snow       = 1     ! 0 = default value for number of snow layers
!!$    nsfc_subs       = 1     ! 1 = default value for number of TILES
!!$    lseaice    = .FALSE.
!!$    llake      = .FALSE.
!!$    lmulti_snow= .FALSE.
!!$    zml_soil=(/ 0.005,0.02,0.06,0.18,0.54,1.62,4.86,14.58 /)



        !
        ! Since the data structure of TERRA differs from that in ICON, 
        ! we have to do some matching (i.e. copying)



        DO isubs = 1,nsfc_subs
          tch_t(1:i_endidx,jb,isubs)      = prm_diag%tch(1:i_endidx,jb)
          tcm_t(1:i_endidx,jb,isubs)      = prm_diag%tcm(1:i_endidx,jb)
          tfv_t(1:i_endidx,jb,isubs)      = prm_diag%tfv(1:i_endidx,jb)
          soiltyp_t(1:i_endidx,jb,isubs)  = ext_data%atm%soiltyp(1:i_endidx,jb) 
          plcov_t(1:i_endidx,jb,isubs)    = ext_data%atm%plcov_mx(1:i_endidx,jb)
          rootdp_t(1:i_endidx,jb,isubs)   = ext_data%atm%rootdp(1:i_endidx,jb) 
          sai_t(1:i_endidx,jb,isubs)      = prm_diag%sai(1:i_endidx,jb)        
          tai_t(1:i_endidx,jb,isubs)      = prm_diag%tai(1:i_endidx,jb)        
          eai_t(1:i_endidx,jb,isubs)      = prm_diag%eai(1:i_endidx,jb)        
          t_2m_t(1:i_endidx,jb,isubs)     = prm_diag%t_2m(1:i_endidx,jb) 
          u_10m_t(1:i_endidx,jb,isubs)    = prm_diag%u_10m(1:i_endidx,jb)
          v_10m_t(1:i_endidx,jb,isubs)    = prm_diag%v_10m(1:i_endidx,jb)  
          sobs_t(1:i_endidx,jb,isubs)     = prm_diag%swflxsfc(1:i_endidx,jb) 
          thbs_t(1:i_endidx,jb,isubs)     = prm_diag%lwflxsfc(1:i_endidx,jb) 
          pabs_t(1:i_endidx,jb,isubs)     = prm_diag%swflxsfc(1:i_endidx,jb) 
          landmask_t(1:i_endidx,jb,isubs) = ext_data%atm%fr_land(1:i_endidx,jb)


          ! copy diagnostic variables
          h_snow_t(1:i_endidx,jb,1,isubs) = lnd_diag%h_snow(1:i_endidx,jb,isubs)
          h_snow_t(1:i_endidx,jb,2,isubs) = lnd_diag%h_snow(1:i_endidx,jb,isubs)


          ! copy prognostic variables
          t_snow_t(1:i_endidx,jb,1,isubs) = lnd_prog_now%t_snow(1:i_endidx,jb,isubs)
          t_snow_t(1:i_endidx,jb,2,isubs) = lnd_prog_now%t_snow(1:i_endidx,jb,isubs)
          t_snow_mult_t(1:i_endidx,0:nlev_snow,jb,1,isubs) = &
            & lnd_prog_now%t_snow_mult(1:i_endidx,1:nlev_snow+1,jb,isubs)
          t_snow_mult_t(1:i_endidx,0:nlev_snow,jb,2,isubs) = &
            & lnd_prog_now%t_snow_mult(1:i_endidx,1:nlev_snow+1,jb,isubs)
          t_s_t(1:i_endidx,jb,1,isubs)    = lnd_prog_now%t_s(1:i_endidx,jb,isubs)
          t_s_t(1:i_endidx,jb,2,isubs)    = lnd_prog_now%t_s(1:i_endidx,jb,isubs)
          w_snow_t(1:i_endidx,jb,1,isubs) = lnd_prog_now%w_snow(1:i_endidx,jb,isubs)
          w_snow_t(1:i_endidx,jb,2,isubs) = lnd_prog_now%w_snow(1:i_endidx,jb,isubs)
          rho_snow_t(1:i_endidx,jb,1,isubs) = lnd_prog_now%rho_snow(1:i_endidx,jb,isubs)
          rho_snow_t(1:i_endidx,jb,2,isubs) = lnd_prog_now%rho_snow(1:i_endidx,jb,isubs)
          rho_snow_mult_t(1:i_endidx,1:nlev_snow,jb,1,isubs) = &
            & lnd_prog_now%rho_snow_mult(1:i_endidx,1:nlev_snow,jb,isubs)
          rho_snow_mult_t(1:i_endidx,1:nlev_snow,jb,2,isubs) = &
            & lnd_prog_now%rho_snow_mult(1:i_endidx,1:nlev_snow,jb,isubs)
          w_i_t(1:i_endidx,jb,1,isubs) = lnd_prog_now%w_i(1:i_endidx,jb,isubs)
          w_i_t(1:i_endidx,jb,2,isubs) = lnd_prog_now%w_i(1:i_endidx,jb,isubs)
          t_so_t(1:i_endidx,0:nlev_soil+1,jb,1,isubs) = &
            & lnd_prog_now%t_so(1:i_endidx,1:nlev_soil+2,jb,isubs)
          t_so_t(1:i_endidx,0:nlev_soil+1,jb,2,isubs) = &
            & lnd_prog_now%t_so(1:i_endidx,1:nlev_soil+2,jb,isubs)
          w_so_t(1:i_endidx,1:nlev_soil+1,jb,1,isubs) = &
            & lnd_prog_now%w_so(1:i_endidx,1:nlev_soil+1,jb,isubs)
          w_so_t(1:i_endidx,1:nlev_soil+1,jb,2,isubs) = &
            & lnd_prog_now%w_so(1:i_endidx,1:nlev_soil+1,jb,isubs)
          w_so_ice_t(1:i_endidx,1:nlev_soil+1,jb,1,isubs) = &
            & lnd_prog_now%w_so_ice(1:i_endidx,1:nlev_soil+1,jb,isubs)
          w_so_ice_t(1:i_endidx,1:nlev_soil+1,jb,2,isubs) = &
            & lnd_prog_now%w_so_ice(1:i_endidx,1:nlev_soil+1,jb,isubs)  
          wliq_snow_t(1:i_endidx,1:nlev_snow,jb,1,isubs) = &
            & lnd_prog_now%wliq_snow(1:i_endidx,1:nlev_snow,jb,isubs)
          wliq_snow_t(1:i_endidx,1:nlev_snow,jb,2,isubs) = &
            & lnd_prog_now%wliq_snow(1:i_endidx,1:nlev_snow,jb,isubs)
          wtot_snow_t(1:i_endidx,1:nlev_snow,jb,1,isubs) = &
            & lnd_prog_now%wtot_snow(1:i_endidx,1:nlev_snow,jb,isubs)
          wtot_snow_t(1:i_endidx,1:nlev_snow,jb,2,isubs) = &
            & lnd_prog_now%wtot_snow(1:i_endidx,1:nlev_snow,jb,isubs)
          dzh_snow_t(1:i_endidx,1:nlev_snow,jb,1,isubs)  = &
            & lnd_prog_now%dzh_snow(1:i_endidx,1:nlev_snow,jb,isubs)
          dzh_snow_t(1:i_endidx,1:nlev_snow,jb,2,isubs)  = &
            & lnd_prog_now%dzh_snow(1:i_endidx,1:nlev_snow,jb,isubs)
        ENDDO

!DR        DO n_lu =1, nclass_lu
!DR          lu_class_frac(1:i_endidx,jb,n_lu)= ext_data%atm%lu_class_fraction(1:i_endidx,jb,n_lu)
!DR        ENDDO


        CALL terra_multlay(                &
        &  ie=nproma, je=1                , & ! array dimensions
        &  istartpar=i_startidx, iendpar=i_endidx    , & ! optional start/end indicies
        &  jstartpar=1, jendpar=1                    , & ! optional start/end indicies
        &  ke=nlev, nsubs0=1, nsubs1=nsfc_subs       , & ! nsfc_subs
        &  ke_soil=nlev_soil, ke_snow=nlev_snow      , &
        &  czmls=zml_soil                            , & ! processing soil level structure 
        &  dt=tcall_sfc_jg                           , &
        &  soiltyp_subs  = soiltyp_t(:,jb,:), & !ext_data%atm%soiltyp(:,jb) , &! (:,1,jb) , & ! type of the soil (keys 0-9)  --
        &  plcov         = plcov_t(:,jb,:)  , & !ext_data%atm%plcov_mx(:,jb), &! ,1, fraction of plant cover                   --
        &  rootdp        = rootdp_t(:,jb,:) , & !ext_data%atm%rootdp(:,jb) , & ! ,1, depth of the roots                        ( m  )
        &  sai           = sai_t(:,jb,:)    , & !prm_diag%sai(:,jb)        , & ! ,1,surface area index                           --
        &  tai           = tai_t(:,jb,:)    , & !prm_diag%tai(:,jb)        , & ! ,1, transpiration area index                    --
        &  eai           = eai_t(:,jb,:)    , & !prm_diag%eai(:,jb)        , & ! ,1, earth area (evaporative surface area) index --
        &  landmask      = landmask_t(:,jb,:) ,& !ext_data%atm%fr_land(:,jb) , & ! ,1, landpoint mask                         --
        &  rsmin2d       = ext_data%atm%rsmin(:,jb),  & ! minimum stomata resistance    ( s/m )
           !
        &  u  =  p_diag%u(:,:,jb)             , & ! zonal wind speed                       ( m/s )
        &  v  =  p_diag%v(:,:,jb)             , & ! meridional wind speed                  ( m/s )
        &  t  =  p_diag%temp(:,:,jb)          , & ! temperature                            (  k  )
        &  qv = p_prog_rcf%tracer(:,:,jb,iqv) , & ! specific water vapor content           (kg/kg)
        &  p0 = p_diag%pres(:,:,jb)           , & !!!!JH base state pressure               ( Pa  ) 
!       &  pp = pp_t(:,1:nlev,jb,:)           , & ! deviation from the reference pressure  ( pa  )
        &  ps = p_diag%pres_sfc(:,jb)         , & ! surface pressure                       ( pa  )
           !
        &  t_snow        =  t_snow_t(:,jb,:,:)    , & ! temperature of the snow-surface              (  K  )
        &  t_snow_mult   =  t_snow_mult_t(:,0:nlev_snow,jb,:,:), & ! temperature of the snow-surface (  K  )
        &  t_s           =  t_s_t(:,jb,:,:)       , & ! temperature of the ground surface            (  K  )
        &  t_g           =  lnd_prog_now%t_gt(:,jb,:)     , & ! weighted surface temperature                 (  K  )
        &  qv_s          =  lnd_diag%qv_st(:,jb,:)  , & ! specific humidity at the surface             (kg/kg)
        &  w_snow        =  w_snow_t(:,jb,:,:)    , & ! water content of snow                        (m H2O)
        &  rho_snow      =  rho_snow_t(:,jb,:,:)  , & ! snow density                                 (kg/m**3)
        &  rho_snow_mult =  rho_snow_mult_t(:,1:nlev_snow,jb,:,:), & ! snow density                  (kg/m**3)
        &  h_snow        =  h_snow_t(:,jb,:,:)    , & ! snow height                                  (  m  ) 
        &  w_i           =  w_i_t(:,jb,:,:)       , & ! water content of interception water          (m H2O)
        &  t_so          =  t_so_t(:,0:nlev_soil+1,jb,:,:),&      ! soil temperature (main level)    (  K  )
        &  w_so          =  w_so_t(:,1:nlev_soil+1,jb,:,:)    , & ! total water conent (ice + liquid water)       (m H20)
        &  w_so_ice      =  w_so_ice_t(:,1:nlev_soil+1,jb,:,:), & ! ice content                                   (m H20)
        &  t_2m          =  t_2m_t(:,jb,:),& !prm_diag%t_2m(:,jb)      , & ! ,nsfc_subs, temperature in 2m        (  K  )
        &  u_10m         =  u_10m_t(:,jb,:),& !prm_diag%u_10m(:,jb)     , & ! ,nsfc_subs, zonal wind in 10m       ( m/s )
        &  v_10m         =  v_10m_t(:,jb,:),& !prm_diag%v_10m(:,jb)     , & ! ,nsfc_subs,  meridional wind in 10m ( m/s )
        &  freshsnow     =  lnd_diag%freshsnow(:,jb,:) , & ! indicator for age of snow in top of snow layer       (  -  )
        &  wliq_snow     =  wliq_snow_t(:,1:nlev_snow,jb,:,:) , & ! liquid water content in the snow              (m H2O)
        &  wtot_snow     =  wtot_snow_t(:,1:nlev_snow,jb,:,:) , & ! total (liquid + solid) water content of snow  (m H2O)
        &  dzh_snow      =  dzh_snow_t(:,1:nlev_snow,jb,:,:)  , & ! layer thickness between half levels in snow   (  m  )
        &  subsfrac      =  lnd_diag%subsfrac(:,jb,:)  , & ! 
           !
        &  prr_con       =  prm_diag%rain_con(:,jb)    , & ! precipitation rate of rain, convective        (kg/m2*s)
        &  prs_con       =  prm_diag%snow_con(:,jb)    , & ! precipitation rate of snow, convective        (kg/m2*s)
        &  prr_gsp       =  prm_diag%rain_gsp(:,jb)    , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
        &  prs_gsp       =  prm_diag%snow_gsp(:,jb)    , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
        &  prg_gsp       =  prm_diag%rain_gsp(:,jb)    , & !!! TEST!! JH precipitation rate of graupel, grid-scale     (kg/m2*s)
           !
        &  tch           = tch_t(:,jb,:),& !prm_diag%tch(:,jb) , & ! ,nsfc_subs, turbulent transfer coefficient for heat     ( -- )
        &  tcm           = tcm_t(:,jb,:),& !prm_diag%tcm(:,jb) , & ! ,nsfc_subs, turbulent transfer coefficient for momentum ( -- )
        &  tfv           = tfv_t(:,jb,:),& !prm_diag%tfv(:,jb) , & ! ,nsfc_subs, laminar reduction factor for evaporation    ( -- )
           ! 
        &  sobs          = sobs_t(:,jb,:) ,& !prm_diag%swflxsfc (:,jb,jb) , & ! ,nsfc_subs, solar radiation at the ground   (W/m2)
        &  thbs          = thbs_t(:,jb,:) ,& !prm_diag%lwflxsfc (:,jb,jb) , & ! ,nsfc_subs, thermal radiation at the ground (W/m2)
        &  pabs          = pabs_t(:,jb,:) ,& !prm_diag%swflxsfc (:,jb)    , & ! ,nsfc_subs, photosynthetic active radiation (W/m2)
           ! 
        &  runoff_s      =  lnd_diag%runoff_s(:,jb,:)  , & ! surface water runoff; sum over forecast      (kg/m2)
        &  runoff_g      =  lnd_diag%runoff_g(:,jb,:)  , & ! soil water runoff; sum over forecast         (kg/m2)
           !
!        &  nstart        = 0                         , & ! first time step of the forecast
        &   ntstep       = jstep-4                    , & ! actual time step
                                ! indices for permutation of three time levels
!        &   nold          = 1                         , & ! corresponds to ntstep - 1
        &   nnow         = 1                         , & ! corresponds to ntstep 
        &   nnew         = 2                           & ! corresponds to ntstep + 1
                                                       )


       ! copy updated variables back to ICON-prognostic fields
        DO isubs = 1,nsfc_subs
          lnd_prog_new%t_snow(1:i_endidx,jb,isubs) = t_snow_t(1:i_endidx,jb,2,isubs)
          lnd_prog_new%t_snow_mult(1:i_endidx,1:nlev_snow+1,jb,isubs) = &
            & t_snow_mult_t(1:i_endidx,0:nlev_snow,jb,2,isubs)
          lnd_prog_new%t_s(1:i_endidx,jb,isubs) = t_s_t(1:i_endidx,jb,2,isubs)
          lnd_prog_new%w_snow(1:i_endidx,jb,isubs) = w_snow_t(1:i_endidx,jb,2,isubs)
          lnd_prog_new%rho_snow(1:i_endidx,jb,isubs) = rho_snow_t(1:i_endidx,jb,2,isubs)
          lnd_prog_new%rho_snow_mult(1:i_endidx,1:nlev_snow,jb,isubs) = &
            & rho_snow_mult_t(1:i_endidx,1:nlev_snow,jb,2,isubs)
          lnd_prog_new%w_i(1:i_endidx,jb,isubs) = w_i_t(1:i_endidx,jb,2,isubs) 
          lnd_prog_new%t_so(1:i_endidx,1:nlev_soil+2,jb,isubs) = &
            & t_so_t(1:i_endidx,0:nlev_soil+1,jb,2,isubs)
          lnd_prog_new%w_so(1:i_endidx,1:nlev_soil+1,jb,isubs) = &
            & w_so_t(1:i_endidx,1:nlev_soil+1,jb,2,isubs)
          lnd_prog_new%w_so_ice(1:i_endidx,1:nlev_soil+1,jb,isubs)= &
            & w_so_ice_t(1:i_endidx,1:nlev_soil+1,jb,2,isubs)
          lnd_prog_new%wliq_snow(1:i_endidx,1:nlev_snow,jb,isubs) = &
            & wliq_snow_t(1:i_endidx,1:nlev_snow,jb,2,isubs)
          lnd_prog_new%wtot_snow(1:i_endidx,1:nlev_snow,jb,isubs) = &
            & wtot_snow_t(1:i_endidx,1:nlev_snow,jb,2,isubs)
          lnd_prog_new%dzh_snow(1:i_endidx,1:nlev_snow,jb,isubs)  = &
            & dzh_snow_t(1:i_endidx,1:nlev_snow,jb,2,isubs)


          lnd_diag%h_snow(1:i_endidx,jb,isubs) = h_snow_t(1:i_endidx,jb,2,isubs)

          !DR ATTENTION: only valid, if nsfc_subs=1 !!!!!
          WHERE (ext_data%atm%fr_land(i_startidx:i_endidx,jb)>0.5_wp)
            lnd_prog_new%t_gt(i_startidx:i_endidx,jb,isubs) = &
              lnd_prog_now%t_gt(i_startidx:i_endidx,jb,isubs)
            lnd_prog_new%t_g(i_startidx:i_endidx,jb) = lnd_prog_new%t_gt(i_startidx:i_endidx,jb,1)
            lnd_diag%qv_s(i_startidx:i_endidx,jb) = lnd_diag%qv_st(i_startidx:i_endidx,jb,1)
            prm_diag%tch (1:i_endidx,jb) = tch_t (1:i_endidx,jb,1)
            prm_diag%t_2m(1:i_endidx,jb) = t_2m_t(1:i_endidx,jb,1)
          END WHERE

        ENDDO


  
      ELSE IF (  atm_phy_nwp_config(jg)%inwp_surface == 2 ) THEN

          !-------------------------------------------------------------------------
          !> ECHAM version 
          !-------------------------------------------------------------------------


     
    ENDIF !inwp_sfc

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

      
  END SUBROUTINE nwp_surface

  !-------------------------------------------------------------------------
  !>
  !! Init surface model TERRA
  !!
  !! Init surface model TERRA.
  !!
  !! @par Revision History
  !! Initial revision by Ekaterina Machulskaya, DWD (2011-07-??)
  !! Modification by Daniel Reienrt, DWD (2011-07-29)
  !! - initialize climatological layer t_so(nlev_soil+2)
  !!
!!$  SUBROUTINE nwp_surface_init    ( p_patch,        & !>in
!!$                                   subsfrac,       &
!!$                                   frac_thres,     &
!!$                                   pt_tiles )
  SUBROUTINE nwp_surface_init( p_patch, ext_data, p_prog_lnd_now, &
    &                          p_prog_lnd_new )
                                   
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch       !<grid/patch info.
    TYPE(t_external_data), INTENT(IN)    :: ext_data
    TYPE(t_lnd_prog)     , INTENT(INOUT) :: p_prog_lnd_now, p_prog_lnd_new
!!$    REAL(wp)             , INTENT(IN)   :: subsfrac(nproma,1,nsfc_subs)
!!$    REAL(wp)             , INTENT(IN)   :: frac_thres     
!!$    TYPE(t_tiles), TARGET, INTENT(INOUT):: pt_tiles      !correspondence between grid & tiles
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jb,ns     
  !-------------------------------------------------------------------------

    i_nchdom  = MAX(1,p_patch%n_childdom)

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jc,ns,i_startidx,i_endidx), SCHEDULE(guided)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        IF (itopo == 1) THEN
          DO ns = 1, nsfc_subs
            DO jc = i_startidx, i_endidx

              ! initialize climatological layer (deepest layer of t_so)
              p_prog_lnd_now%t_so(jc,nlev_soil+2,jb,ns) = ext_data%atm%t_cl(jc,jb)
              p_prog_lnd_new%t_so(jc,nlev_soil+2,jb,ns) = ext_data%atm%t_cl(jc,jb)

              p_prog_lnd_now%t_gt(jc,jb,ns) = p_prog_lnd_now%t_g(jc,jb)
              p_prog_lnd_new%t_gt(jc,jb,ns) = p_prog_lnd_now%t_g(jc,jb)

            END DO
          END DO
        ENDIF


!!$        DO ns = 1, nsfc_subs
!!$
!!$          pt_tiles%length(ns,jb) = 0
!!$
!!$          DO jc = i_startidx, i_endidx
!!$
!!$            IF(subsfrac(jc,1,ns) > frac_thres) THEN
!!$              pt_tiles%length(ns,jb) = pt_tiles%length(ns,jb) + 1
!!$              pt_tiles%corrsp(pt_tiles%length(ns,jb),ns,jb) = jc
!!$            END IF
!!$          END DO
!!$
!!$        END DO

      ENDDO  ! jb loop
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_surface_init

END MODULE mo_nwp_sfc_interface

