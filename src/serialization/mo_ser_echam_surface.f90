!--------------------------------------------------------------------
!
! Serialization routine for ECHAM surface
!
!--------------------------------------------------------------------

MODULE mo_ser_echam_surface

  USE mo_kind,        ONLY: vp, wp
  USE mo_ser_common,  ONLY: init
  IMPLICIT NONE

  LOGICAL :: writeIn = .TRUE.
  LOGICAL :: writeOut = .TRUE.

  PUBLIC :: serialize_input
  PUBLIC :: serialize_output

  CONTAINS

  SUBROUTINE serialize_input(                                   &
                           & pfrc,                              &
                           & pcfh_tile, pcfm_tile,              &
                           & pfac_sfc, pocu, pocv,              &
                           & aa, aa_btm, bb, bb_btm,            &
                           & pcpt_tile, pqsat_tile,             &
                           & ptsfc_tile,                        &
                           & plhflx_tile, pshflx_tile,          &
                           !! optional
                           & lsm,                               &
                           & alake,                             &
                           & pu,                                &
                           & pv,                                &
                           & ptemp,                             &
                           & pq,                                &
                           & prsfl,                             &
                           & prsfc,                             &
                           & pssfl,                             &
                           & pssfc,                             &
                           & rlds,                              &
                           & rlus,                              &
                           & rsds,                              &
                           & rsus,                              &
                           !
                           & rvds_dir,                          &
                           & rpds_dir,                          &
                           & rnds_dir,                          &
                           & rvds_dif,                          &
                           & rpds_dif,                          &
                           & rnds_dif,                          &
                           !
                           & ps,                                &
                           & pcosmu0,                           &
                           & pch_tile,                          &
                           !! for JSBACH
                           & pcsat,                             &
                           & pcair,                             &
                           !
                           & z0m_tile, z0h_lnd,                 &
                           & albvisdir, albnirdir, albvisdif, albnirdif, &
                           & albvisdir_tile,                    &
                           & albnirdir_tile,                    &
                           & albvisdif_tile,                    &
                           & albnirdif_tile,                    &
                           & albedo, albedo_tile,               &
                           & rsns_tile, rlns_tile,              &
                           !! Sea ice
                           & Tsurf,                             &
                           & T1,                                &
                           & T2,                                &
                           & hi,                                &
                           & hs,                                &
                           & Qtop,                              &
                           & Qbot,                              &
                           & conc,                              &
                           & albvisdir_ice, albvisdif_ice,      &
                           & albnirdir_ice, albnirdif_ice)       
    REAL(wp),INTENT(IN) :: pfrc      (:,:)
    REAL(wp),INTENT(IN) :: pcfh_tile (:,:)
    REAL(wp),INTENT(IN) :: pcfm_tile (:,:)
    REAL(wp),INTENT(IN) :: pfac_sfc  (:)
    REAL(wp),INTENT(IN) :: pocu      (:)
    REAL(wp),INTENT(IN) :: pocv      (:)
    REAL(wp),INTENT(INOUT) :: aa     (:,:,:,:)
    REAL(wp),INTENT(INOUT) :: aa_btm (:,:,:,:)
    REAL(wp),INTENT(INOUT) :: bb     (:,:,:)
    REAL(wp),INTENT(INOUT) :: bb_btm (:,:,:)
    REAL(wp),INTENT(INOUT) :: pcpt_tile (:,:)
    REAL(wp),INTENT(INOUT) :: pqsat_tile(:,:)
    REAL(wp),INTENT(INOUT) :: ptsfc_tile (:,:)

    REAL(wp),INTENT(INOUT) :: plhflx_tile (:,:)
    REAL(wp),INTENT(INOUT) :: pshflx_tile (:,:)
    !! JSBACH input
    REAL(wp),OPTIONAL,INTENT(IN) :: lsm(:)
    REAL(wp),OPTIONAL,INTENT(IN) :: alake(:)
    REAL(wp),OPTIONAL,INTENT(IN) :: pu        (:)              ! zonal wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: pv        (:)              ! meridional wind lowest level
    REAL(wp),OPTIONAL,INTENT(IN) :: ptemp     (:)              ! temperature of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: pq        (:)              ! humidity of lowest atmospheric level
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfl     (:)              ! rain large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: prsfc     (:)              ! rain convective
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfl     (:)              ! snow large scale
    REAL(wp),OPTIONAL,INTENT(IN) :: pssfc     (:)              ! snow convective
    REAL(wp),OPTIONAL,INTENT(IN) :: rlds      (:)              ! downward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN) :: rsds      (:)              ! downward surface shortwave flux [W/m2]
    
    REAL(wp),INTENT(IN) :: rvds_dir(:)        ! all-sky   vis. dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rpds_dir(:)        ! all-sky   par  dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rnds_dir(:)        ! all-sky   nir  dir. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rvds_dif(:)        ! all-sky   vis. dif. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rpds_dif(:)        ! all-sky   par  dif. downward flux at current   time [W/m2]
    REAL(wp),INTENT(IN) :: rnds_dif(:)        ! all-sky   nir  dif. downward flux at current   time [W/m2]

    REAL(wp),OPTIONAL,INTENT(IN) :: ps        (:)              ! surface pressure
    REAL(wp),OPTIONAL,INTENT(IN) :: pcosmu0   (:)              ! cos of zenith angle
    REAL(wp),OPTIONAL,INTENT(IN) :: pch_tile  (:,:)
    !! JSBACH output
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: z0h_lnd(:), z0m_tile(:,:)  ! OUT
    !
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albedo(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir(:), albvisdif(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir(:), albnirdif(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albedo_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlus     (:)           ! INOUT upward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(IN)    :: rsus     (:)           ! IN upward surface shortwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rsns_tile(:,:) ! shortwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlns_tile(:,:) ! longwave net flux at surface on tiles
    !! Sea ice
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Tsurf(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T1   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T2   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hi   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hs   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qtop (:,:) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qbot (:,:) ! OUT
    REAL(wp),OPTIONAL,INTENT(IN)    :: conc (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_ice(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_ice(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_ice(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_ice(:,:)

    !$ser verbatim IF (writeIn) THEN
    !$ser verbatim call init()
    !$ser savepoint echam_surface-input
#if defined SERIALIZE_CREATE_REFERENCE 
    !$ser mode write
#elif defined SERIALIZE_PERTURB_REFERENCE
    !$ser mode read-perturb
#elif defined SERIALIZE_READ_REFERENCE
    !$ser mode read
#else
#error SERIALIZATION MODE IS NOT SET
#endif 
    !$ser data pfrc=pfrc pcfh_tile=pcfh_tile pcfm_tile=pcfm_tile               &
    !$ser&     pfac_sfc=pfac_sfc pocu=pocu pocv=pocv aa=aa aa_btm=aa_btm bb_bb &
    !$ser&     bb_btm=bb_btm pcpt_tile=pcpt_tile pqsat_tile=pqsat_tile         &
    !$ser&     ptsfc_tile=ptsfc_tile plhflx_tile=plhflx_tile                   &
    !$ser&     pshflx_tile=pshflx_tile lsm=lsm alake=alake pu=pu pv=pv         &
    !$ser&     ptemp=ptemp pq=pq prsfl=prsfl prsfc=prsfc pssfl=pssfl           &
    !$ser&     pssfc=pssfc rlds=rlds rlus=rlus rsds=rsds rsus=rsus             &
    !$ser&     rvds_dir=rvds_dir rpds_dir=rpds_dir rnds_dir=rnds_dir           &
    !$ser&     rvds_dif=rvds_dif rpds_dif=rpds_dif rnds_dif=rnds_dif ps=ps     &
    !$ser&     pcosmu0=pcosmu0 pch_tile=pch_tile pcsat=pcsat pcair=pcair       &
    !$ser&     z0m_tile=z0m_tile z0h_lnd=z0h_lnd                               &
    !$ser&     albvisdir=albvisdir albnirdir=albnirdir albvisdif=albvisdif     &
    !$ser&     albnirdif=albnirdif albvisdir_tile=albvisdir_tile               &
    !$ser&     albnirdir_tile=albnirdir_tile albvisdif_tile=albvisdif_tile     &
    !$ser&     albnirdif_tile=albnirdif_tile albedo=albedo                     &
    !$ser&     albedo_tile=albedo_tile Tsurf=Tsurf T1=T1 T2=T2 hi=hi hs=hs     &
    !$ser&     Qtop=Qtop Qbot=Qbot conc=conc albvisdir_ice=albvisdir_ice       &
    !$ser&     albvisdif_ice=albvisdif_ice albnirdir_ice=albnirdir_ice         &
    !$ser&     albnirdif_ice=albnirdif_ice
    !$ser verbatim writeIn = .FALSE.
    !$ser verbatim ENDIF

  END SUBROUTINE serialize_input

  SUBROUTINE serialize_output(                                  &
                           & aa, aa_btm, bb, bb_btm,            &
                           & pcpt_tile, pqsat_tile,             &
                           & ptsfc_tile,                        &
                           & pu_stress_gbm, pv_stress_gbm,      &
                           & plhflx_gbm, pshflx_gbm,            &
                           & pevap_gbm,                         &
                           & pu_stress_tile,   pv_stress_tile,  &
                           & plhflx_tile, pshflx_tile,          &
                           & pevap_tile,                        &
                           !! optional
                           & rlus,                              &
                           !! for JSBACH
                           & pcsat,                             &
                           & pcair,                             &
                           & q_snocpymlt,                       &
                           !
                           & z0m_tile, z0h_lnd,                 &
                           & albvisdir, albnirdir, albvisdif, albnirdif, &
                           & albvisdir_tile,                    &
                           & albnirdir_tile,                    &
                           & albvisdif_tile,                    &
                           & albnirdif_tile,                    &
                           & albedo, albedo_tile,               &
                           & ptsfc,                             &
                           & ptsfc_rad,                         &
                           & rsns_tile, rlns_tile,              &
                           & lake_ice_frc,                      &
                           !! Sea ice
                           & Tsurf,                             &
                           & T1,                                &
                           & T2,                                &
                           & hi,                                &
                           & hs,                                &
                           & Qtop,                              &
                           & Qbot,                              &
                           & albvisdir_ice, albvisdif_ice,      &
                           & albnirdir_ice, albnirdif_ice)       
    REAL(wp),INTENT(INOUT) :: aa     (:,:,:,:)
    REAL(wp),INTENT(INOUT) :: aa_btm (:,:,:,:)
    REAL(wp),INTENT(INOUT) :: bb     (:,:,:)
    REAL(wp),INTENT(INOUT) :: bb_btm (:,:,:)
    REAL(wp),INTENT(INOUT) :: pcpt_tile (:,:)
    REAL(wp),INTENT(INOUT) :: pqsat_tile(:,:)
    REAL(wp),INTENT(INOUT) :: ptsfc_tile (:,:)

    REAL(wp),INTENT(OUT)   :: pu_stress_gbm (:)
    REAL(wp),INTENT(OUT)   :: pv_stress_gbm (:)
    REAL(wp),INTENT(OUT)   ::    plhflx_gbm (:)
    REAL(wp),INTENT(OUT)   ::    pshflx_gbm (:)
    REAL(wp),INTENT(OUT)   ::     pevap_gbm (:)

    REAL(wp),INTENT(OUT)   :: pu_stress_tile (:,:)
    REAL(wp),INTENT(OUT)   :: pv_stress_tile (:,:)
    REAL(wp),INTENT(INOUT) :: plhflx_tile (:,:)   ! OUT
    REAL(wp),INTENT(INOUT) :: pshflx_tile (:,:)   ! OUT
    REAL(wp),INTENT(OUT)   :: pevap_tile (:,:)

    !! JSBACH output
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcsat(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: pcair(:)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: q_snocpymlt(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: z0h_lnd(:), z0m_tile(:,:)  ! OUT
    !
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albedo(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir(:), albvisdif(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir(:), albnirdif(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albedo_tile(:,:)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc    (:)
    REAL(wp),OPTIONAL,INTENT(OUT)   :: ptsfc_rad(:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlus     (:)           ! INOUT upward surface  longwave flux [W/m2]
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rsns_tile(:,:) ! shortwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(INOUT) :: rlns_tile(:,:) ! longwave net flux at surface on tiles
    REAL(wp),OPTIONAL,INTENT(OUT)   :: lake_ice_frc(:)        ! fraction of ice on lakes
    !! Sea ice
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Tsurf(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T1   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: T2   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hi   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: hs   (:,:) ! for coupled ocean only
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qtop (:,:) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: Qbot (:,:) ! OUT
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdir_ice(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albvisdif_ice(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdir_ice(:,:)
    REAL(wp),OPTIONAL,INTENT(INOUT) :: albnirdif_ice(:,:)

    !$ser verbatim if (writeOut) then
    !$ser verbatim call init()
    !$ser savepoint echam_surface-output
    !$ser mode write
    !$ser data aa=aa aa_btm=aa_btm bb=bb bb_btm=bb_btm pcpt_tile=pcpt_tile     &
    !$ser&     pqsat_tile=pqsat_tile ptsfc_tile=ptsfc_tile                     &
    !$ser&     pu_stress_gbm=pu_stress_gbm pv_stress_gbm=pv_stress_gbm         &
    !$ser&     plhflx_tile=plhflx_tile pshflx_tile=pshflx_tile                 &
    !$ser&     pevap_tile=pevap_tile rlus=rlus pcsat=pcsat pcair=pcair         &
    !$ser&     q_snocpymlt=q_snocpymlt z0m_tile=z0m_tile z0h_lnd=z0h_lnd       &
    !$ser&     albvisdir=albvisdir albnirdir=albnirdir albvisdif=albvisdif     &
    !$ser&     albnirdif=albnirdif albvisdir_tile=albvisdir_tile               &
    !$ser&     albnirdir_tile=albnirdir_tile albvisdif_tile=albvisdif_tile     &
    !$ser&     albnirdif_tile=albnirdif_tile albedo=albedo                     &
    !$ser&     albedo_tile=albedo_tile ptsfc=ptsfc ptsfc_rad=ptsfc_rad         &
    !$ser&     rsns_tile=rsns_tile rlns_tile=rlns_tile                         &
    !$ser&     lake_ice_frc=lake_ice_frc Tsurf=Tsurf T1=T1 T2=T2 hi=hi hs=hs   &
    !$ser&     Qtop=Qtop Qbot=Qbot albvisdir_ice=albvisdir_ice                 &
    !$ser&     albvisdif_ice=albvisdif_ice albnirdir_ice=albnirdir_ice         &
    !$ser&     albnirdif_ice=albnirdif_ice
    !$ser verbatim writeOut = .FALSE.
    !$ser verbatim endif

  END SUBROUTINE serialize_output

END MODULE mo_ser_echam_surface
