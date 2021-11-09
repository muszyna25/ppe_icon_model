!>
!! @brief Subroutine interface_cloud_two calls NWP two-moment bulk microphysics
!!
!! @author Monika Esch, MPI-M, 2020-04
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_cloud_two

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime
  USE mo_copy                ,ONLY: copy

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_cov_config    ,ONLY: echam_cov_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
       &                            t_echam_phy_tend,  prm_tend

  USE mo_cloud_two_types     ,ONLY: t_cloud_two_input, t_cloud_two_output
  USE mo_cloud_two_memory    ,ONLY:   cloud_two_input,   cloud_two_output
  USE mo_cloud_two           ,ONLY:   cloud_two

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_two

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqh,  &
       &                            iqnc, iqni, iqnr, iqns, iqng, iqnh,  &
       &                            ininact, msg_level
  ! USE mo_cloud_two_config    ,ONLY: cloud_two_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_cloud_two

CONTAINS

  SUBROUTINE interface_cloud_two(jg,jb,jcs,jce        ,&
       &                         nproma,nlev,ntracer  ,& 
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         pdtime               )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev,ntracer
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Pointers
    !
    ! to echam_phy_memory
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend
    !
    ! to cloud_two_memory
    TYPE(t_cloud_two_input ),POINTER    :: input
    TYPE(t_cloud_two_output),POINTER    :: output

    ! Local variables
    !
    LOGICAL  :: lparamcpl
    INTEGER  :: fc_two
    !
    REAL(wp) :: tend_ta_two(nproma,nlev)
    REAL(wp) :: tend_qv_two(nproma,nlev)
    REAL(wp) :: tend_qc_two(nproma,nlev)
    REAL(wp) :: tend_qnc_two(nproma,nlev)
    REAL(wp) :: tend_qi_two(nproma,nlev)
    REAL(wp) :: tend_qni_two(nproma,nlev)
    REAL(wp) :: tend_qr_two(nproma,nlev)
    REAL(wp) :: tend_qnr_two(nproma,nlev)
    REAL(wp) :: tend_qs_two(nproma,nlev)
    REAL(wp) :: tend_qns_two(nproma,nlev)
    REAL(wp) :: tend_qg_two(nproma,nlev)
    REAL(wp) :: tend_qng_two(nproma,nlev)
    REAL(wp) :: tend_qh_two(nproma,nlev)
    REAL(wp) :: tend_qnh_two(nproma,nlev)
    REAL(wp) :: tend_ninact_two(nproma,nlev)
    !
    ! Dummy variable for 2mom scheme. Used only when
    !  assimilation of radar data using latent heat nudging
    !  is switched on (ldass_lhn=true). 
    REAL(wp) :: zqrsflux(nproma,nlev)
    !
    INTEGER  :: jk, jks, jke, jl
    INTEGER  :: jc

    ! associate pointers
    !
    ! to echam_phy_memory
    field  => prm_field(jg)
    tend   => prm_tend (jg)
    !
    ! to cloud_two memory
    input  => cloud_two_input (jg)
    output => cloud_two_output(jg)

    ! associate pointers
    lparamcpl =  echam_phy_config(jg)%lparamcpl
    fc_two    =  echam_phy_config(jg)%fc_two

    jks       =  echam_cov_config(jg)% jkscov
    jke       =  nlev

    ! security
    tend_ta_two(:,:)   = 0.0_wp
    tend_qv_two(:,:)   = 0.0_wp
    tend_qc_two(:,:)   = 0.0_wp
    tend_qnc_two(:,:)   = 0.0_wp
    tend_qi_two(:,:)   = 0.0_wp
    tend_qni_two(:,:)   = 0.0_wp
    tend_qr_two(:,:)   = 0.0_wp
    tend_qnr_two(:,:)   = 0.0_wp
    tend_qs_two(:,:)   = 0.0_wp
    tend_qns_two(:,:)   = 0.0_wp
    tend_qg_two(:,:)   = 0.0_wp
    tend_qng_two(:,:)   = 0.0_wp
    tend_qh_two(:,:)   = 0.0_wp
    tend_qnh_two(:,:)   = 0.0_wp
    tend_ninact_two(:,:)   = 0.0_wp

    zqrsflux(:,:)  = 0.0_wp

    ! store input in memory for output
    !
    ! input parameters
    !
    IF (ASSOCIATED(input% jg        )) CALL copy(jcs,jce, jg       , input% jg (:,jb))
    IF (ASSOCIATED(input% jcs       )) CALL copy(jcs,jce, jcs      , input% jcs (:,jb))
    IF (ASSOCIATED(input% jce       )) CALL copy(jcs,jce, jce      , input% jce (:,jb))
    IF (ASSOCIATED(input% msg_level )) CALL copy(jcs,jce, msg_level, input% msg_level (:,jb))
    IF (ASSOCIATED(input% pdtime    )) CALL copy(jcs,jce, pdtime   , input% pdtime    (:,jb))
    !
    ! input fields, in
    !
    IF (ASSOCIATED(input% dz    )) CALL copy(jcs,jce, jks,jke, field% dz        (:,:,jb)    , input% dz    (:,:,jb))
    IF (ASSOCIATED(input% zh    )) CALL copy(jcs,jce, jks,jke, field% zh        (:,:,jb)    , input% zh    (:,:,jb))
    IF (ASSOCIATED(input% rho   )) CALL copy(jcs,jce, jks,jke, field% rho       (:,:,jb)    , input% rho   (:,:,jb))
    IF (ASSOCIATED(input% pf    )) CALL copy(jcs,jce, jks,jke, field% pfull     (:,:,jb)    , input% pf    (:,:,jb))
    IF (ASSOCIATED(input% cpair )) CALL copy(jcs,jce, jks,jke, field% cpair     (:,:,jb)    , input% cpair (:,:,jb))
    !
    ! input fields, inout
    !
    IF (ASSOCIATED(input% ta    )) CALL copy(jcs,jce, jks,jke, field% ta   (:,:,jb)        , input% ta    (:,:,jb))
    IF (ASSOCIATED(input% qv    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqv)    , input% qv    (:,:,jb))
    IF (ASSOCIATED(input% qc    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqc)    , input% qc    (:,:,jb))
    IF (ASSOCIATED(input% qnc   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqnc)   , input% qnc   (:,:,jb))
    IF (ASSOCIATED(input% qi    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqi)    , input% qi    (:,:,jb))
    IF (ASSOCIATED(input% qni   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqni)   , input% qni   (:,:,jb))
    IF (ASSOCIATED(input% qr    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqr)    , input% qr    (:,:,jb))
    IF (ASSOCIATED(input% qnr   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqnr)   , input% qnr   (:,:,jb))
    IF (ASSOCIATED(input% qs    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqs)    , input% qs    (:,:,jb))
    IF (ASSOCIATED(input% qns   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqns)   , input% qns   (:,:,jb))
    IF (ASSOCIATED(input% qg    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqg)    , input% qg    (:,:,jb))
    IF (ASSOCIATED(input% qng   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqng)   , input% qng   (:,:,jb))
    IF (ASSOCIATED(input% qh    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqh)    , input% qh    (:,:,jb))
    IF (ASSOCIATED(input% qnh   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqnh)   , input% qnh   (:,:,jb))
    IF (ASSOCIATED(input% ninact)) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,ininact), input% ninact(:,:,jb))
    IF (ASSOCIATED(input% w     )) CALL copy(jcs,jce, jks,jke, field% wa(:,:,jb)           , input% w     (:,:,jb))
    !
    IF (ASSOCIATED(input% pr_rain)) CALL copy(jcs,jce, field% rain_gsp_rate (:,jb), input% pr_rain   (:,jb))
    IF (ASSOCIATED(input% pr_ice )) CALL copy(jcs,jce, field%  ice_gsp_rate (:,jb), input% pr_ice    (:,jb))
    IF (ASSOCIATED(input% pr_snow)) CALL copy(jcs,jce, field% snow_gsp_rate (:,jb), input% pr_snow   (:,jb))
    IF (ASSOCIATED(input% pr_grpl)) CALL copy(jcs,jce, field% graupel_gsp_rate (:,jb), input% pr_grpl   (:,jb))
    IF (ASSOCIATED(input% pr_hail)) CALL copy(jcs,jce, field% hail_gsp_rate (:,jb), input% pr_hail   (:,jb))

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN

          IF (ltimer) CALL timer_start(timer_two)
          !

          CALL cloud_two( jg,  nproma, nlev                     ,& !< in : grid index
               &          jcs, jce                              ,& !< in : column index range
               &          jks, jke                              ,& !< in : column index range
               &          msg_level                             ,& !< in : message level 
               &          pdtime                                ,& !< in : timestep
               &          field% dz        (:,:,jb)       ,& !< in : vertical layer thickness
               &          field% zh        (:,:,jb)       ,& !< in : height of half levels
               &          field% rho       (:,:,jb)       ,& !< in : density
               &          field% pfull     (:,:,jb)       ,& !< in : pressure
               &          field% cpair     (:,:,jb)       ,& !< in : specific heat of air
               &          field% ta        (:,:,jb)       ,& !< inout : temperature
               &          field% qtrc      (:,:,jb,iqv)   ,& !< inout : sp humidity
               &          field% qtrc      (:,:,jb,iqc)   ,& !< inout : cloud water
               &          field% qtrc      (:,:,jb,iqnc)  ,& !< inout : cloud water number
               &          field% qtrc      (:,:,jb,iqi)   ,& !< inout : ice
               &          field% qtrc      (:,:,jb,iqni)  ,& !< inout : ice number
               &          field% qtrc      (:,:,jb,iqr)   ,& !< inout : rain
               &          field% qtrc      (:,:,jb,iqnr)  ,& !< inout : rain number
               &          field% qtrc      (:,:,jb,iqs)   ,& !< inout : snow
               &          field% qtrc      (:,:,jb,iqns)  ,& !< inout : snow number
               &          field% qtrc      (:,:,jb,iqg)   ,& !< inout : graupel
               &          field% qtrc      (:,:,jb,iqng)  ,& !< inout : graupel number
               &          field% qtrc      (:,:,jb,iqh)   ,& !< inout : hail
               &          field% qtrc      (:,:,jb,iqnh)  ,& !< inout : hail number
               &          field% qtrc      (:,:,jb,ininact) ,& !< inout : activated ice nuclei
               &          field% wa        (:,:,jb)     ,& !< in : vertical velocity
               &          tend_ta_two      (:,:)          ,& !< out: tendency of temperature
               &          tend_qv_two      (:,:)          ,& !< out: tendency of water vapor
               &          tend_qc_two      (:,:)          ,& !< out: tendency of cloud water
               &          tend_qnc_two     (:,:)          ,& !< out: tendency of cloud water
               &          tend_qi_two      (:,:)          ,& !< out: tendency of cloud ice
               &          tend_qni_two     (:,:)          ,& !< out: tendency of cloud ice
               &          tend_qr_two      (:,:)          ,& !< out: tendency of rain
               &          tend_qnr_two     (:,:)          ,& !< out: tendency of rain droplets
               &          tend_qs_two      (:,:)          ,& !< out: tendency of snow
               &          tend_qns_two     (:,:)          ,& !< out: tendency of snow droplets
               &          tend_qg_two      (:,:)          ,& !< out: tendency of graupel 
               &          tend_qng_two     (:,:)          ,& !< out: tendency of graupel droplets
               &          tend_qh_two      (:,:)          ,& !< out: tendency of hail
               &          tend_qnh_two     (:,:)          ,& !< out: tendency of hail droplets
               &          tend_ninact_two  (:,:)          ,& !< out: tendency of activated ice
               &          field% rain_gsp_rate (:,jb)           ,& !& inout: precip rate rain
               &          field%  ice_gsp_rate (:,jb)           ,& !& inout: precip rate ice
               &          field% snow_gsp_rate (:,jb)           ,& !& inout: precip rate snow
               &          field% graupel_gsp_rate (:,jb)        ,& !& inout: precip rate graupel
               &          field% hail_gsp_rate (:,jb)           ,& !& inout: precip rate hail
               &          zqrsflux         (:,:)          )  !< inout: unused 3d total prec. rate
          !
          IF (ltimer) CALL timer_stop(timer_two)

          !
          ! Calculate rain and snow
          !
          DO jc = jcs, jce
            field% rsfl(jc,jb) = field% rain_gsp_rate (jc,jb)
            field% ssfl(jc,jb) = field% snow_gsp_rate (jc,jb) + field% graupel_gsp_rate (jc,jb) &
            &                  + field%  ice_gsp_rate (jc,jb) + field% hail_gsp_rate (jc,jb)
          END DO
          !
          !
          ! store output in memory for output or recycling
          !
          IF (ASSOCIATED(output% tend_ta_two   )) CALL copy(jcs,jce, jks,jke, tend_ta_two   (:,:), output% tend_ta_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qv_two   )) CALL copy(jcs,jce, jks,jke, tend_qv_two   (:,:), output% tend_qv_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qc_two   )) CALL copy(jcs,jce, jks,jke, tend_qc_two   (:,:), output% tend_qc_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qnc_two  )) CALL copy(jcs,jce, jks,jke, tend_qnc_two  (:,:), output% tend_qnc_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qi_two   )) CALL copy(jcs,jce, jks,jke, tend_qi_two   (:,:), output% tend_qi_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qni_two  )) CALL copy(jcs,jce, jks,jke, tend_qni_two  (:,:), output% tend_qni_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qr_two   )) CALL copy(jcs,jce, jks,jke, tend_qr_two   (:,:), output% tend_qr_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qnr_two  )) CALL copy(jcs,jce, jks,jke, tend_qnr_two  (:,:), output% tend_qnr_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qs_two   )) CALL copy(jcs,jce, jks,jke, tend_qs_two   (:,:), output% tend_qs_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qns_two  )) CALL copy(jcs,jce, jks,jke, tend_qns_two  (:,:), output% tend_qns_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qg_two   )) CALL copy(jcs,jce, jks,jke, tend_qg_two   (:,:), output% tend_qg_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qng_two  )) CALL copy(jcs,jce, jks,jke, tend_qng_two  (:,:), output% tend_qng_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_qh_two   )) CALL copy(jcs,jce, jks,jke, tend_qh_two   (:,:), output% tend_qh_two   (:,:,jb))
          IF (ASSOCIATED(output% tend_qnh_two  )) CALL copy(jcs,jce, jks,jke, tend_qnh_two  (:,:), output% tend_qnh_two  (:,:,jb))
          IF (ASSOCIATED(output% tend_ninact_two)) CALL copy(jcs,jce, jks,jke, tend_ninact_two(:,:), output% tend_ninact_two(:,:,jb))
          !
          IF (ASSOCIATED(output% pr_rain     )) CALL copy(jcs,jce, field% rain_gsp_rate (:,jb), output% pr_rain (:,jb))
          IF (ASSOCIATED(output% pr_ice      )) CALL copy(jcs,jce, field%  ice_gsp_rate (:,jb), output% pr_ice  (:,jb))
          IF (ASSOCIATED(output% pr_snow     )) CALL copy(jcs,jce, field% snow_gsp_rate (:,jb), output% pr_snow (:,jb))
          IF (ASSOCIATED(output% pr_grpl     )) CALL copy(jcs,jce, field% graupel_gsp_rate (:,jb), output% pr_grpl (:,jb))
          IF (ASSOCIATED(output% pr_hail     )) CALL copy(jcs,jce, field% hail_gsp_rate (:,jb), output% pr_hail (:,jb))
          !
       ELSE    ! is_active
          !
          ! retrieve output from memory for recycling
          !
          IF (ASSOCIATED(output% tend_ta_two   )) CALL copy(jcs,jce, jks,jke, output% tend_ta_two   (:,:,jb), tend_ta_two   (:,:))
          IF (ASSOCIATED(output% tend_qv_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qv_two   (:,:,jb), tend_qv_two   (:,:))
          IF (ASSOCIATED(output% tend_qc_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qc_two   (:,:,jb), tend_qc_two   (:,:))
          IF (ASSOCIATED(output% tend_qnc_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qnc_two  (:,:,jb), tend_qnc_two  (:,:))
          IF (ASSOCIATED(output% tend_qi_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qi_two   (:,:,jb), tend_qi_two   (:,:))
          IF (ASSOCIATED(output% tend_qni_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qni_two  (:,:,jb), tend_qni_two  (:,:))
          IF (ASSOCIATED(output% tend_qr_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qr_two   (:,:,jb), tend_qr_two   (:,:))
          IF (ASSOCIATED(output% tend_qnr_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qnr_two  (:,:,jb), tend_qnr_two  (:,:))
          IF (ASSOCIATED(output% tend_qs_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qs_two   (:,:,jb), tend_qs_two   (:,:))
          IF (ASSOCIATED(output% tend_qns_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qns_two  (:,:,jb), tend_qns_two  (:,:))
          IF (ASSOCIATED(output% tend_qg_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qg_two   (:,:,jb), tend_qg_two   (:,:))
          IF (ASSOCIATED(output% tend_qng_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qng_two  (:,:,jb), tend_qng_two  (:,:))
          IF (ASSOCIATED(output% tend_qh_two   )) CALL copy(jcs,jce, jks,jke, output% tend_qh_two   (:,:,jb), tend_qh_two   (:,:))
          IF (ASSOCIATED(output% tend_qnh_two  )) CALL copy(jcs,jce, jks,jke, output% tend_qnh_two  (:,:,jb), tend_qnh_two  (:,:))
          IF (ASSOCIATED(output% tend_ninact_two)) CALL copy(jcs,jce, jks,jke, output% tend_ninact_two(:,:,jb), tend_ninact_two(:,:))
          !
       END IF  ! is_active

       !

       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_two)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          DO jk = jks, jke
            DO jl = jcs, jce
              ! use tendency to update the model state
              tend%   ta_phy(jl,jk,jb)        = tend%   ta_phy(jl,jk,jb)         + tend_ta_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqv)    = tend% qtrc_phy(jl,jk,jb,iqv)     + tend_qv_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqc)    = tend% qtrc_phy(jl,jk,jb,iqc)     + tend_qc_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqnc)   = tend% qtrc_phy(jl,jk,jb,iqnc)    + tend_qnc_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqi)    = tend% qtrc_phy(jl,jk,jb,iqi)     + tend_qi_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqni)   = tend% qtrc_phy(jl,jk,jb,iqni)    + tend_qni_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqr)    = tend% qtrc_phy(jl,jk,jb,iqr)     + tend_qr_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqnr)   = tend% qtrc_phy(jl,jk,jb,iqnr)    + tend_qnr_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqs)    = tend% qtrc_phy(jl,jk,jb,iqs)     + tend_qs_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqns)   = tend% qtrc_phy(jl,jk,jb,iqns)    + tend_qns_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqg)    = tend% qtrc_phy(jl,jk,jb,iqg)     + tend_qg_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqng)   = tend% qtrc_phy(jl,jk,jb,iqng)    + tend_qng_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqh)    = tend% qtrc_phy(jl,jk,jb,iqh)     + tend_qh_two    (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqnh)   = tend% qtrc_phy(jl,jk,jb,iqnh)    + tend_qnh_two   (jl,jk)
              tend% qtrc_phy(jl,jk,jb,ininact)= tend% qtrc_phy(jl,jk,jb,ininact) + tend_ninact_two(jl,jk)
            END DO
          END DO
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
         SELECT CASE(fc_two)
         CASE(0)
            ! diagnostic, do not use tendency
         CASE(1)
             DO jk = jks, jke
               DO jl = jcs, jce
                 field%   ta(jl,jk,jb)        = field%   ta(jl,jk,jb)         + tend_ta_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqv)    = field% qtrc(jl,jk,jb,iqv)     + tend_qv_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqc)    = field% qtrc(jl,jk,jb,iqc)     + tend_qc_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqnc)   = field% qtrc(jl,jk,jb,iqnc)    + tend_qnc_two(jl,jk)   *pdtime
                 field% qtrc(jl,jk,jb,iqi)    = field% qtrc(jl,jk,jb,iqi)     + tend_qi_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqni)   = field% qtrc(jl,jk,jb,iqni)    + tend_qni_two(jl,jk)   *pdtime
                 field% qtrc(jl,jk,jb,iqr)    = field% qtrc(jl,jk,jb,iqr)     + tend_qr_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqnr)   = field% qtrc(jl,jk,jb,iqnr)    + tend_qnr_two(jl,jk)   *pdtime
                 field% qtrc(jl,jk,jb,iqs)    = field% qtrc(jl,jk,jb,iqs)     + tend_qs_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqns)   = field% qtrc(jl,jk,jb,iqns)    + tend_qns_two(jl,jk)   *pdtime
                 field% qtrc(jl,jk,jb,iqg)    = field% qtrc(jl,jk,jb,iqg)     + tend_qg_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqng)   = field% qtrc(jl,jk,jb,iqng)    + tend_qng_two(jl,jk)   *pdtime
                 field% qtrc(jl,jk,jb,iqh)    = field% qtrc(jl,jk,jb,iqh)     + tend_qh_two(jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqnh)   = field% qtrc(jl,jk,jb,iqnh)    + tend_qnh_two(jl,jk)   *pdtime
                 field% qtrc(jl,jk,jb,ininact)= field% qtrc(jl,jk,jb,ininact) + tend_ninact_two(jl,jk)*pdtime
               END DO
             END DO
             !
    ! output fields, inout
    !
    IF (ASSOCIATED(output% ta    )) CALL copy(jcs,jce, jks,jke, field% ta   (:,:,jb)        , output% ta    (:,:,jb))
    IF (ASSOCIATED(output% qv    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqv)    , output% qv    (:,:,jb))
    IF (ASSOCIATED(output% qc    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqc)    , output% qc    (:,:,jb))
    IF (ASSOCIATED(output% qnc   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqnc)   , output% qnc   (:,:,jb))
    IF (ASSOCIATED(output% qi    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqi)    , output% qi    (:,:,jb))
    IF (ASSOCIATED(output% qni   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqni)   , output% qni   (:,:,jb))
    IF (ASSOCIATED(output% qr    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqr)    , output% qr    (:,:,jb))
    IF (ASSOCIATED(output% qnr   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqnr)   , output% qnr   (:,:,jb))
    IF (ASSOCIATED(output% qs    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqs)    , output% qs    (:,:,jb))
    IF (ASSOCIATED(output% qns   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqns)   , output% qns   (:,:,jb))
    IF (ASSOCIATED(output% qg    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqg)    , output% qg    (:,:,jb))
    IF (ASSOCIATED(output% qng   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqng)   , output% qng   (:,:,jb))
    IF (ASSOCIATED(output% qh    )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqh)    , output% qh    (:,:,jb))
    IF (ASSOCIATED(output% qnh   )) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,iqnh)   , output% qnh   (:,:,jb))
    IF (ASSOCIATED(output% ninact)) CALL copy(jcs,jce, jks,jke, field% qtrc (:,:,jb,ininact), output% ninact(:,:,jb))
    IF (ASSOCIATED(output% w     )) CALL copy(jcs,jce, jks,jke, field% wa(:,:,jb)           , output% w     (:,:,jb))
    !
!!$         CASE(2)
!!$            ! use tendency as forcing in the dynamics
!!$            ...
         END SELECT
       END IF
       !
    ELSE       ! is_in_sd_ed_interval : bypass area
       !
       IF (ASSOCIATED(output% tend_ta_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_ta_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qv_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qv_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qc_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qc_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qnc_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qnc_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qi_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qi_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qni_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qni_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qr_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qr_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qnr_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qnr_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qs_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qs_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qns_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qns_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qg_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qg_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qng_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qng_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_qh_two     )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qh_two     (:,:,jb))
       IF (ASSOCIATED(output% tend_qnh_two    )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_qnh_two    (:,:,jb))
       IF (ASSOCIATED(output% tend_ninact_two )) CALL copy(jcs,jce, jks,jke, 0._wp, output% tend_ninact_two (:,:,jb))
       !
       IF (ASSOCIATED(output% pr_rain )) CALL copy(jcs,jce, 0._wp, output% pr_rain (:,jb))
       IF (ASSOCIATED(output% pr_ice  )) CALL copy(jcs,jce, 0._wp, output% pr_ice  (:,jb))
       IF (ASSOCIATED(output% pr_snow )) CALL copy(jcs,jce, 0._wp, output% pr_snow (:,jb))
       IF (ASSOCIATED(output% pr_grpl )) CALL copy(jcs,jce, 0._wp, output% pr_grpl (:,jb))
       IF (ASSOCIATED(output% pr_hail )) CALL copy(jcs,jce, 0._wp, output% pr_hail (:,jb))
       !
    END IF     ! is_in_sd_ed_interval

    ! disassociate pointers

    NULLIFY(field)
    NULLIFY(tend)

  END SUBROUTINE interface_cloud_two

END MODULE mo_interface_cloud_two
