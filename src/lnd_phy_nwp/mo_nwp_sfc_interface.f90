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

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text , finish
!  USE mo_mpi,                  ONLY: p_pe, p_nprocs
!  USE mo_parallel_nml,         ONLY: p_test_pe, p_test_run
!  USE mo_sync,                 ONLY: SYNC_C, sync_patch_array_mult

  USE mo_model_domain,         ONLY: t_patch
!  USE mo_grf_interpolation,    ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
 ! USE mo_subdivision,          ONLY: p_patch_local_parent

  USE mo_ext_data,             ONLY: t_external_data
  USE mo_nonhydro_state,       ONLY: t_nh_prog, t_nh_diag,&
   &                                 t_nh_metrics
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag,prm_diag,&
       &                             t_nwp_phy_tend
  USE mo_nwp_lnd_state,        ONLY: t_lnd_prog, t_lnd_diag, &
!<em
                                     t_tiles
!em>

  USE mo_run_nml,              ONLY: nproma, msg_level, iqv, &
    &                                 icc!, inextra_2d

  USE mo_atm_phy_nwp_nml,      ONLY: inwp_surface, inwp_satad  
  USE mo_lnd_nwp_nml,          ONLY: nlev_soil, nztlev, nlev_snow, nsfc_subs, &
    &                                lseaice, llake, lmulti_snow, zml_soil
!  USE mo_turbdiff_ras,       ONLY: organize_turbdiff
  USE mo_satad,              ONLY: sat_pres_water, spec_humi  
  USE src_turbdiff,          ONLY: organize_turbdiff
!  USE mo_icoham_sfc_indices, ONLY: nsfc_type, igbm, iwtr, iice, ilnd
  USE mo_vdiff_driver,       ONLY: vdiff
  USE mo_echam_vdiff_params, ONLY: z0m_oce
  USE mo_physical_constants, ONLY: &
     t0_melt => tmelt,& ! absolute zero for temperature
     r_v   => rv    , & ! gas constant for water vapour
     r_d   => rd    , & ! gas constant for dry air
     rvd_m_o=>vtmpc1, & ! r_v/r_d - 1
     o_m_rdv        , & ! 1 - r_d/r_v
     rdv            , & ! r_d / r_v
     lh_v  => alv   , & ! latent heat of vapourization
     lh_s  => als   , & ! latent heat of sublimation
     lh_f  => alf   , & ! latent heat of fusion
     cp_d  => cpd   , & ! specific heat of dry air at constant press
     cpdr  => rcpd  , & ! (specific heat of dry air at constant press)^-1
     g     => grav  , & ! acceleration due to gravity
     sigma => stbo  , & ! Boltzmann-constant
     rdocp => rd_o_cpd  ! r_d / cp_d
  USE mo_convect_tables,   ONLY:   &
     b1    => c1es  , & !! constants for computing the sat. vapour
     b2w   => c3les , & !! pressure over water (l) and ice (i)
     b2i   => c3ies , & !!               -- " --
     b4w   => c4les , & !!               -- " --
     b4i   => c4ies , & !!               -- " --
     b234w => c5les     !!               -- " --
  USE mo_cuparameters,       ONLY: rho_w => rhoh2o    ! density of liquid water (kg/m^3)
  USE mo_phyparam_soil
  USE mo_nwp_phy_init
  USE mo_soil_ml,            ONLY: terra_multlay
!  USE mo_aggregate_surface,  ONLY: subsmean,subs_disaggregate_radflux,subsmean_albedo

  
  IMPLICIT NONE

  PRIVATE

!<em
  PUBLIC  ::  nwp_surface, nwp_surface_init
!em>

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_surface    ( tcall_sfc_jg,                      & !>input
                            & p_patch,p_metrics,                 & !>input
                            & ext_data,                          & !>input
                            & p_prog,                            & !>inout
                            & p_prog_rcf,                        & !>in/inout
                            & p_diag ,                           & !>inout
                            & prm_diag,                          & !>inout 
                            & lnd_prog, lnd_diag              )!>inout


    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
    TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the prog vars
!    TYPE(t_nh_prog),      TARGET,INTENT(IN)   :: p_prog_now_rcf  !<progs with red.
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !< atm phys vars
!    TYPE(t_nwp_phy_tend),TARGET, INTENT(inout):: prm_nwp_tend    !< atm tend vars
    TYPE(t_lnd_prog),            INTENT(inout):: lnd_prog        !< prog vars for sfc
    TYPE(t_lnd_diag),            INTENT(inout):: lnd_diag        !< diag vars for sfc
    REAL(wp),                    INTENT(in)   :: tcall_sfc_jg    !< time interval for 
                                                                 !< surface
!    REAL(wp),                    INTENT(in)   :: mean_charlen    !< characteristic griddistance
                                                                 !< needed  by turbulence
    ! Local array bounds:

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: nlevp1           !< number of full and half levels

    ! Local scalars:

    INTEGER :: jc,jk,jb,jt,jg      !block indeces
 
    ! local variables for turbdiff
    INTEGER :: ierrstat=0
    CHARACTER (LEN=25) :: eroutine=''
    CHARACTER (LEN=80) :: errormsg=''

!<em
    REAL(wp), TARGET, ALLOCATABLE :: t_t(:,:,:,:), p_t(:,:,:,:)
    REAL(wp), TARGET, ALLOCATABLE :: u_t(:,:,:,:), v_t(:,:,:,:)
    REAL(wp), TARGET, ALLOCATABLE :: ps_t(:,:,:)
!em>

  ! local variables related to the blocking

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

  ! number of vertical levels
    nlevp1 = p_patch%nlevp1

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


    ALLOCATE(t_t(nproma, nlevp1, 2, i_startblk:i_endblk))
    ALLOCATE(p_t(nproma, nlevp1, 2, i_startblk:i_endblk))
    ALLOCATE(u_t(nproma, nlevp1, 2, i_startblk:i_endblk))
    ALLOCATE(v_t(nproma, nlevp1, 2, i_startblk:i_endblk))
    ALLOCATE(ps_t(nproma, 2, i_startblk:i_endblk))

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jt,jc,jk,i_startidx,i_endidx), SCHEDULE(guided)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        IF(inwp_surface == 0) THEN
          ! check dry case
          IF(inwp_satad == 0) THEN
            lnd_diag%qv_s (:,:) = 0._wp
          ELSE
          ! 
          !> adjust  humidity at water surface because of changed surface pressure
          !
          DO jc = i_startidx, i_endidx
            lnd_diag%qv_s (jc,jb) = &
                 &         spec_humi(sat_pres_water(lnd_prog%t_g(jc,jb)),&
                 &                                   p_diag%pres_sfc(jc,jb) )
          ENDDO
        ENDIF
        ENDIF

        IF ( inwp_surface == 1 ) THEN
          
! #ifdef __BOUNDCHECK

!Example for Debugging, Diagnosis
!!$    IF (msg_level >= 0) THEN
!!$      WRITE(message_text,'(a,3E15.7)') ' Soil Temperature = ', &
!!$           &  lnd_prog%t_so(1,:,:,:,12)
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

!<em
    t_t(:,:,1,jb) = p_diag%temp(:,:,jb)
    p_t(:,:,1,jb) = p_diag%pres(:,:,jb)
    t_t(:,:,2,jb) = p_diag%temp(:,:,jb)
    p_t(:,:,2,jb) = p_diag%pres(:,:,jb)
    u_t(:,:,1,jb) = p_diag%u(:,:,jb)
    v_t(:,:,1,jb) = p_diag%v(:,:,jb)
    u_t(:,:,2,jb) = p_diag%u(:,:,jb)
    v_t(:,:,2,jb) = p_diag%v(:,:,jb)
    ps_t(:,2,jb) = p_diag%pres_sfc(:,jb)
!em>

  write(0,*) "SFC-DIAGNOSIS ",nlevp1,tcall_sfc_jg
  write(0,*) "zml_soil: ",zml_soil
  write(0,*) "t",p_diag%temp(1,nlevp1,jb)
  write(0,*) "p",p_diag%pres(1,nlevp1,jb)
  write(0,*) "t_g",lnd_prog%t_gt(1,:,nsfc_subs,jb)
  write(0,*) "qv_s",lnd_diag%qv_st(1,:,nsfc_subs,jb)

           CALL terra_multlay(                &
                ie=nproma,je=1              , & ! array dimensions
                istartpar=i_startidx,iendpar=i_endidx         , & ! optional start/end indicies
                jstartpar=1,jendpar=1                         , & ! optional start/end indicies
                ke=nlevp1, nsubs0=1,nsubs1=nsfc_subs              , & ! nsfc_subs
                ke_soil=nlev_soil, ke_snow=nlev_snow      , &
                czmls=zml_soil                            , & ! processing soil level structure 
                dt=tcall_sfc_jg                           , &
                !
                soiltyp_subs  =  ext_data%atm%soiltyp(:,jb) , & ! (:,1,jb) , & ! type of the soil (keys 0-9)         --
                plcov         =  ext_data%atm%plcov_mx(:,jb), & ! ,1, fraction of plant cover                        --
                rootdp        =  ext_data%atm%rootdp(:,jb) , & ! ,1, depth of the roots                            ( m  )
                sai           =  prm_diag%sai(:,jb)         , & ! ,1,surface area index                              --
                tai           =  prm_diag%tai(:,jb)         , & ! ,1, transpiration area index                        --
                eai           =  prm_diag%eai(:,jb)         , & ! ,1, earth area (evaporative surface area) index     --
                rlandmask     =  ext_data%atm%fr_land(:,jb) , & ! ,1, landpoint mask                                  --
                rsmin2d       =  ext_data%atm%rsmin(:,jb),  & ! minimum stomata resistance                    ( s/m )
                !
!<em
!                u   = p_diag%u(:,:,jb)                     , & ! zonal wind speed                              ( m/s )
!                v  = p_diag%v(:,:,jb)                     , & ! meridional wind speed                         ( m/s )
!                t  = p_diag%temp(:,:,jb)                  , & ! temperature                                   (  k  )
                u  = u_t(:,:,:,jb)                     , & ! zonal wind speed                              ( m/s )
                v  = v_t(:,:,:,jb)                     , & ! meridional wind speed                         ( m/s )
                t  = t_t(:,:,:,jb)                  , & ! temperature                                   (  k  )
!em>
                qv = p_prog_rcf%tracer(:,:,jb,iqv)        , & ! specific water vapor content                  (kg/kg)
                p0 = p_diag%pres(:,:,jb)                  , & !!!!JH base state pressure                           (Pa) 
!<em
!                pp = p_diag%pres(:,:,jb)                  , & ! deviation from the reference pressure         ( pa  )
!                ps = p_diag%pres_sfc(:,jb)              , & ! surface pressure                              ( pa  )
                pp = p_t(:,:,:,jb)                  , & ! deviation from the reference pressure         ( pa  )
                ps = ps_t(:,:,jb)              , & ! surface pressure                              ( pa  )
!em>
                !
                t_snow        =  lnd_prog%t_snow(:,:,nsfc_subs,jb)    , & ! temperature of the snow-surface               (  K  )
                t_snow_mult   =  lnd_prog%t_snow_mult(:,:,:,nsfc_subs,jb), & ! temperature of the snow-surface               (  K  )
                t_s           =  lnd_prog%t_s(:,:,nsfc_subs,jb)       , & ! temperature of the ground surface             (  K  )
                t_g           =  lnd_prog%t_gt(:,:,nsfc_subs,jb)       , & ! weighted surface temperature                  (  K  )
                qv_s          =  lnd_diag%qv_st(:,:,nsfc_subs,jb)      , & ! specific humidity at the surface              (kg/kg)
                w_snow        =  lnd_prog%w_snow(:,:,nsfc_subs,jb)    , & ! water content of snow                         (m H2O)
                rho_snow      =  lnd_prog%rho_snow(:,nsfc_subs,jb)  , & ! snow density                                  (kg/m**3)
                rho_snow_mult =  lnd_prog%rho_snow_mult(:,:,:,nsfc_subs,jb), & ! snow density                                  (kg/m**3)
                h_snow        =  lnd_diag%h_snow(:,:,nsfc_subs,jb)    , & ! snow height                                   (  m  
                w_i           =  lnd_prog%w_i(:,:,nsfc_subs,jb)     , & ! water content of interception water           (m H2O)
                t_so          =  lnd_prog%t_so(:,:,:,nsfc_subs,jb)    , & ! soil temperature (main level)                 (  K  )
                w_so          =  lnd_prog%w_so(:,:,:,nsfc_subs,jb)    , & ! total water conent (ice + liquid water)       (m H20)
                w_so_ice      =  lnd_prog%w_so_ice(:,:,:,nsfc_subs,jb), & ! ice content                                   (m H20)
                t_2m          =  prm_diag%t_2m(:,jb)      , & ! ,nsfc_subs, temperature in 2m                             (  K  )
                u_10m         =  prm_diag%u_10m(:,jb)     , & ! ,nsfc_subs, zonal wind in 10m                             ( m/s )
                v_10m         =  prm_diag%v_10m(:,jb)     , & ! ,nsfc_subs,  meridional wind in 10m                        ( m/s )
                freshsnow     =  lnd_diag%freshsnow(:,nsfc_subs,jb) , & ! indicator for age of snow in top of snow layer(  -  )
                wliq_snow     =  lnd_diag%wliq_snow(:,:,nsfc_subs,jb) , & ! liquid water content in the snow              (m H2O)
                wtot_snow     =  lnd_diag%wtot_snow(:,:,nsfc_subs,jb) , & ! total (liquid + solid) water content of snow  (m H2O)
                dzh_snow      =  lnd_prog%dzh_snow(:,:,:,nsfc_subs,jb)  , & ! layer thickness between half levels in snow   (  m  )
!<em
                subsfrac      =  lnd_diag%subsfrac(:,nsfc_subs,jb)  , & ! 
!em>
                !
                prr_con       =  prm_diag%rain_con(:,jb)    , & ! precipitation rate of rain, convective        (kg/m2*s)
                prs_con       =  prm_diag%snow_con(:,jb)    , & ! precipitation rate of snow, convective        (kg/m2*s)
                prr_gsp       =  prm_diag%rain_gsp(:,jb)    , & ! precipitation rate of rain, grid-scale        (kg/m2*s)
                prs_gsp       =  prm_diag%snow_gsp(:,jb)    , & ! precipitation rate of snow, grid-scale        (kg/m2*s)
                prg_gsp       =  prm_diag%rain_gsp(:,jb)    , & !!! TEST!! JH precipitation rate of graupel, grid-scale     (kg/m2*s)
                !
                tch           =  prm_diag%tch(:,jb)       , & ! ,nsfc_subs,  turbulent transfer coefficient for heat       ( -- )
                tcm           =  prm_diag%tcm(:,jb)       , & ! ,nsfc_subs, turbulent transfer coefficient for momentum   ( -- )
                tfv           =  prm_diag%tfv(:,jb)       , & ! ,nsfc_subs, laminar reduction factor for evaporation      ( -- )
                ! 
                sobs          =  prm_diag%swflxsfc (:,jb) , & ! ,nsfc_subs, solar radiation at the ground                 ( W/m2)
                thbs          =  prm_diag%lwflxsfc (:,jb) , & ! ,nsfc_subs, thermal radiation at the ground               ( W/m2)
                pabs          =  prm_diag%swflxsfc (:,jb) , & !!!! ,nsfc_subs, photosynthetic active radiation               ( W/m2)
                ! 
                runoff_s      =  lnd_diag%runoff_s(:,nsfc_subs,jb)  , & ! surface water runoff; sum over forecast      (kg/m2)
                runoff_g      =  lnd_diag%runoff_g(:,nsfc_subs,jb)  , & ! soil water runoff; sum over forecast         (kg/m2)
                rstom         =  lnd_diag%rstom(:,nsfc_subs,jb)     , & ! stomata resistance                           ( s/m )
                lhfl_bs       =  lnd_diag%lhfl_bs(:,nsfc_subs,jb)   , & ! average latent heat flux from bare soil evap.( W/m2)
                lhfl_pl       =  lnd_diag%lhfl_pl(:,nsfc_subs,jb)   , & ! average latent heat flux from plants         ( W/m2)
                !
                nstart        = 0                         , & ! first time step of the forecast
                ntstep        = 0                         , & ! actual time step
                                ! indices for permutation of three time levels
!               nold          = 1                         , & ! corresponds to ntstep - 1
                nnow          = 1                         , & ! corresponds to ntstep 
                nnew          = 2                         , & ! corresponds to ntstep + 1
                nztlev        =2                          , & ! time step scheme 2,3
                lmelt       = .TRUE.                      , & ! soil model with melting process
                lmelt_var   = .TRUE.                      , & ! freezing temperature dependent on water content
                lmulti_snow = .FALSE.                     , & ! run the multi-layer snow model
                !
                itype_gscp=3                              , & ! type of grid-scale precipitation physics
                itype_trvg=2                              , & ! type of vegetation transpiration parameterization
                itype_evsl=2                              , & ! type of parameterization of bare soil evaporation
                itype_tran=2                              , & ! type of surface to atmospher transfer
                itype_root=1                              , & ! type of root density distribution
                itype_heatcond=1                          , & ! type of soil heat conductivity
                itype_hydbound=1                          , & ! type of hydraulic lower boundary condition
                lstomata=.true.                       , & ! map of minimum stomata resistance
                l2tls  =.true.                        , & ! forecast with 2-TL integration scheme
                lana_rho_snow=.false.                 , & ! if .TRUE., take rho_snow-values from analysis file 
                itype_subs=1                            & ! type of subscale surface treatment =1 MOSAIC, =2 TILE       
                                                      )
        IF (ierrstat.NE.0) THEN
           CALL finish(eroutine, errormsg)
        END IF

! #endif
  
      ELSE IF ( inwp_surface == 2 ) THEN

          !-------------------------------------------------------------------------
          !> ECHAM version 
          !-------------------------------------------------------------------------


     
    ENDIF !inwp_sfc

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  DEALLOCATE(t_t)
  DEALLOCATE(p_t)
  DEALLOCATE(u_t)
  DEALLOCATE(v_t)
  DEALLOCATE(ps_t)

!Example for Debugging, Diagnosis
!!$    IF (msg_level >= 15) THEN
!!$      WRITE(message_text,'(a,3E15.7)') ' bottom TKE after turbulence = ', &
!!$           &  p_prog_now_rcf%tke(30,nlev+1,12), p_prog_now_rcf%tke(30,nlev,12),&
!!$           &  p_prog_now_rcf%tke(30,nlev-11,12)
!!$     CALL message('', TRIM(message_text))
!!$      WRITE(message_text,'(a,3E15.7)') ' bottom TKE after turbulence = ', &
!!$           &  p_prog_rcf%tke(30,nlev+1,12), p_prog_rcf%tke(30,nlev,12),&
!!$          &  p_prog_rcf%tke(30,nlev-1,12)
!!$      CALL message('', TRIM(message_text))
!!$    ENDIF

      
  END SUBROUTINE nwp_surface

  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_surface_init    ( p_patch,                      & !>input
                                   subsfrac,                     &
                                   frac_thres,                   &
                                   pt_tiles )
                                   
    TYPE(t_patch),       TARGET,                        INTENT(IN)    :: p_patch         !<grid/patch info.
    REAL(KIND = ireals), DIMENSION(nproma,1,nsfc_subs), INTENT(IN)    :: subsfrac
    REAL(KIND = ireals),                                INTENT(IN)    :: frac_thres     
    TYPE(t_tiles),       TARGET,                        INTENT(INOUT) :: pt_tiles        !correspondence between grid & tiles
    
    ! Local array bounds:
    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jb,ns     

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

        DO ns = 1, nsfc_subs

          pt_tiles%length(ns,jb) = 0

          DO jc = i_startidx, i_endidx
            IF(subsfrac(jc,1,ns) > frac_thres) THEN
              pt_tiles%length(ns,jb) = pt_tiles%length(ns,jb) + 1
              pt_tiles%corrsp(pt_tiles%length(ns,jb),ns,jb) = jc
            END IF
          END DO

        END DO

      ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_surface_init

END MODULE mo_nwp_sfc_interface

