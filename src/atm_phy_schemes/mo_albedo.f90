!>
!! Calculation of surface albedo
!!
!! Calculation of surface albedo taking soil type, vegetation 
!! and snow/ice conditions into account.
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial Revision by Daniel Reinert, DWD (2012-03-19)
!! Moved to a central place from mo_nwp_rad_interface and 
!! mo_nwp_rrtm_interface.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_albedo

  USE mo_kind,                 ONLY: wp
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_radiation_config,     ONLY: rad_csalbw
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, ntiles_water, ntiles_lnd,  &
    &                                lseaice, llake, isub_water, isub_lake,   &
    &                                isub_seaice
  USE mo_phyparam_soil,        ONLY: csalb, csalb_snow_fe, csalb_snow_fd,     &
    &                                csalb_snow_min, csalb_snow_max, csalb_p, csalb_snow
  USE mo_physical_constants,   ONLY: tmelt, tf_salt
  USE mo_data_flake,           ONLY: albedo_whiteice_ref, albedo_blueice_ref, &
    &                                c_albice_MR, tpl_T_f, h_Ice_min_flk
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_impl_constants,       ONLY: min_rlcell_int

  IMPLICIT NONE

  PRIVATE



  PUBLIC  :: sfc_albedo
  PUBLIC  :: sfc_albedo_modis


CONTAINS


  !>
  !! Calculation of surface albedo
  !!
  !! Calculation of surface albedo based on tabulated shortwave bare soil 
  !! albedo data. In addition, soil type, vegetation and snow/ice conditions 
  !! are taken into account
  !!
  !! @par Revision History
  !! Initial Revision by Thorsten Reinhardt, AGeoBw, Offenbach
  !! Modification by Daniel Reinert, DWD (2012-03-19)
  !! - Moved here from mo_nwp_rad_interface and mo_nwp_rrtm_interface.
  !!   Adaption to TERRA-tile approach.
  !! - Modification by Daniel Reinert, DWD (2013-07-03)
  !!   Albedo for lake-ice points based on an empirical formula proposed by 
  !!   proposed by Mironov and Ritter (2004)
  !! - Modification by Daniel Reinert, DWD (2013-08-07)
  !!   Added albedo for direct radiation (VIS and NIR spectral bands)
  !!
  SUBROUTINE sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)

    TYPE(t_patch),          INTENT(   in):: pt_patch  !< grid/patch info.

    TYPE(t_external_data),  INTENT(   in):: ext_data  !< external data

    TYPE(t_lnd_prog),       INTENT(   in):: lnd_prog  !< land prognostic state (new)

    TYPE(t_wtr_prog),       INTENT(   in):: wtr_prog  !< water prognostic state (new)

    TYPE(t_lnd_diag),       INTENT(   in):: lnd_diag  !< land diagnostic state

    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag

    ! Local scalars:
    REAL(wp):: zvege                   !< plant cover fraction
    REAL(wp):: zsnow                   !< snow cover fraction
    REAL(wp):: zminsnow_alb            !< temperature-dependent minimum snow albedo
    REAL(wp):: t_fac                   !< factor for temperature dependency of zminsnow_alb over glaciers
    REAL(wp):: zsalb_snow              !< snow albedo (predictor)
    REAL(wp):: zsnow_alb               !< snow albedo (corrector)

    INTEGER :: jg                      !< patch ID
    INTEGER :: jb, jc, ic, jt          !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: ist
    INTEGER :: ilu                     !< land cover class
    INTEGER :: i_count_lnd             !< number of land points
    INTEGER :: i_count_sea             !< number of sea points
    INTEGER :: i_count_flk             !< number of lake points
    INTEGER :: i_count_seaice          !< number of seaice points

    !-----------------------------------------------------------------------

    jg = pt_patch%id
    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,ic,i_startidx,i_endidx,ist,zvege,zsnow,  &
!$OMP            zsalb_snow,zsnow_alb,ilu,i_count_lnd,i_count_sea, &
!$OMP            i_count_flk,i_count_seaice,zminsnow_alb,t_fac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk


      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                 i_startidx, i_endidx, rl_start, rl_end)


      !------------------------------------------------------------------------------
      ! Calculation of surface albedo taking soil type,              
      ! vegetation and snow/ice conditions into account
      !------------------------------------------------------------------------------
      
      IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN


        !
        ! 1. Consider land points (may have tiles)
        !
        ! - loop over surface tiles
        ! - note that different grid points may have different numbers 
        !   of active tiles (1<=ntiles<=ntiles_total). Therefore each tile has a 
        !   separate index list.
        ! 
        DO jt = 1, ntiles_total

          i_count_lnd = ext_data%atm%gp_count_t(jb,jt)

          IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty

!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count_lnd

            jc = ext_data%atm%idx_lst_t(ic,jb,jt)

            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included

            ! surface albedo including moisture correction
            prm_diag%albdif_t(jc,jb,jt) = csalb(ist)&
              &                         - rad_csalbw(ist)*lnd_prog%w_so_t(jc,1,jb,jt)

          ENDDO  ! ic


          ! Account for Snow cover and vegetation
          !---------------------------------------
          IF (ntiles_lnd > 1) THEN     ! tile approach


!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count_lnd

              jc = ext_data%atm%idx_lst_t(ic,jb,jt)

              zsnow= 0.0_wp

              ! temperature-dependent minimum snow albedo over glaciers
              ! this is needed to prevent unrealistically low albedos over Antarctica
              IF (ext_data%atm%lc_class_t(jc,jb,jt) == ext_data%atm%i_lc_snow_ice) THEN
                t_fac = MIN(1._wp,0.1_wp*(tmelt-lnd_prog%t_snow_t(jc,jb,jt)))
                zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*csalb_snow
              ELSE
                zminsnow_alb = csalb_snow_min
              ENDIF
              !
              ! 1. Consider effects of aging on solar snow albedo
              !
              zsalb_snow = zminsnow_alb + lnd_diag%freshsnow_t(jc,jb,jt)*(csalb_snow_max-zminsnow_alb)

              ! special treatment for forests and artificial surfaces
              ! - aging of snow not considered
              ! - instead, snow albedo is limited to land-class specific value

              ! get land cover class
              ilu = ext_data%atm%lc_class_t(jc,jb,jt) 
              zsnow_alb = MIN(zsalb_snow, ABS(ext_data%atm%snowalb_lcc(ilu)))


              ! 2. account for plant cover and snow cover
              !
              ! plant cover
              zvege = ext_data%atm%plcov_t(jc,jb,jt)
              ! snow cover fraction
              zsnow = lnd_diag%snowfrac_t(jc,jb,jt)

              ! 3. compute final solar snow albedo
              !
              prm_diag%albdif_t(jc,jb,jt) = zsnow * zsnow_alb           &  ! snow covered
                &  + (1.0_wp - zsnow) * (zvege * csalb_p                &  ! snow-free with vege
                &  + (1.0_wp - zvege) * prm_diag%albdif_t(jc,jb,jt))       ! snow-free bare

            ENDDO  ! ic


          ELSE  ! without tiles

!CDIR NODEP,VOVERTAKE,VOB
            DO ic = 1, i_count_lnd

              jc = ext_data%atm%idx_lst_t(ic,jb,jt)

              zsnow= 0.0_wp

              ! temperature-dependent minimum snow albedo over glaciers
              ! this is needed to prevent unrealistically low albedos over Antarctica
              IF (ext_data%atm%lc_class_t(jc,jb,jt) == ext_data%atm%i_lc_snow_ice) THEN
                t_fac = MIN(1._wp,0.1_wp*(tmelt-lnd_prog%t_snow_t(jc,jb,jt)))
                zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*csalb_snow
              ELSE
                zminsnow_alb = csalb_snow_min
              ENDIF
              !
              ! 1. Consider effects of aging on solar snow albedo
              !
              zsalb_snow = zminsnow_alb + lnd_diag%freshsnow_t(jc,jb,jt)*(csalb_snow_max-zminsnow_alb)

              ! special treatment for forests
              zsnow_alb = zsalb_snow*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb)) &
                + csalb_snow_fe * ext_data%atm%for_e(jc,jb)                       &
                + csalb_snow_fd * ext_data%atm%for_d(jc,jb)

              ! 2. account for plant cover and snow cover
              !
              ! plant cover
              zvege = ext_data%atm%plcov_t(jc,jb,jt)
              ! snow cover fraction
              zsnow = lnd_diag%snowfrac_t(jc,jb,jt)

              ! 3. compute final solar snow albedo
              !
              prm_diag%albdif_t(jc,jb,jt) = zsnow * zsnow_alb           &  ! snow covered
                &  + (1.0_wp - zsnow) * (zvege * csalb_p                &  ! snow-free with vege
                &  + (1.0_wp - zvege) * prm_diag%albdif_t(jc,jb,jt))       ! snow-free bare

            ENDDO  ! ic

          ENDIF ! ntiles_lnd

        ENDDO  !ntiles





        ! 2. Consider water points with/without seaice model
        !
        IF ( (lseaice) ) THEN  ! seaice model switched on

          !
          ! Open sea points
          !
          ! - loop over sea points (points which are at least partly ice-free) 
          !
          i_count_sea = ext_data%atm%spw_count(jb)

          DO ic = 1, i_count_sea
            jc = ext_data%atm%idx_lst_spw(ic,jb)

            ist = 9  ! sea water

            prm_diag%albdif_t(jc,jb,isub_water) = csalb(ist)
          ENDDO


          !
          ! Sea-ice points
          !
          ! - loop over sea-ice points (points which are at least partly ice-covered)
          !
          i_count_seaice = ext_data%atm%spi_count(jb)

          DO ic = 1, i_count_seaice
            jc = ext_data%atm%idx_lst_spi(ic,jb)

            ist = 10   ! seaice

            ! In case the sea ice model is used, compute ice albedo for seaice 
            ! points with an empirical formula taken from GME.
            ! The ice albedo is the lower the warmer, and therefore wetter 
            ! the ice is. Use ice temperature at time level nnew 
            ! (2-time level scheme in sea ice model).
            prm_diag%albdif_t(jc,jb,isub_seaice) = csalb(ist) * ( 1.0_wp - 0.3846_wp    &
              &                         * EXP(-0.35_wp*(tmelt-wtr_prog%t_ice(jc,jb))))
            ! gives alb_max = 0.70
            !       alb_min = 0.43
            ! compare with Mironov et. al (2012), Tellus
            !       alb_max = 0.65
            !       alb_min = 0.40
          ENDDO

        ELSE   ! no seaice model

          !
          ! Sea points
          !
          ! - loop over sea points (including sea-ice points)
          !
          i_count_sea = ext_data%atm%spw_count(jb)

          DO ic = 1, i_count_sea
            jc = ext_data%atm%idx_lst_spw(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_water) < tf_salt) THEN 
              ist = 10  ! sea ice
            ELSE
              ist = 9   ! sea water 
            ENDIF

            prm_diag%albdif_t(jc,jb,isub_water) = csalb(ist)
          ENDDO

        ENDIF  ! seaice



        ! 3. Consider lake points with/without lake model
        !
        IF (llake) THEN  ! Lake model switched on

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%fp_count(jb)

!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count_flk
            jc = ext_data%atm%idx_lst_fp(ic,jb)


            !  In case ice is present at lake points, compute ice albedo 
            !  for lake points with an empirical formulation 
            !  proposed by Mironov and Ritter (2004) for use in GME 
            !  [ice_albedo=function(ice_surface_temperature)].
            !  Use surface temperature at time level "nnew".

            ! special handling for ice-covered lake points
            IF (wtr_prog%h_ice(jc,jb) > h_Ice_min_flk) THEN 

              prm_diag%albdif_t(jc,jb,isub_lake) = albedo_whiteice_ref                      &
                &              - (albedo_whiteice_ref - albedo_blueice_ref)                 &
                &              * EXP(-c_albice_MR*(tpl_T_f-lnd_prog%t_g_t(jc,jb,isub_lake)) &
                &              /tpl_T_f)
            ELSE
              ist = 9   ! water
              prm_diag%albdif_t(jc,jb,isub_lake) = csalb(ist)
            ENDIF
          ENDDO

        ELSE  ! Lake model switched off

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%fp_count(jb)

          DO ic = 1, i_count_flk
            jc = ext_data%atm%idx_lst_fp(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_lake) < tpl_T_f) THEN 
              ist = 10  ! sea ice
            ELSE
              ist = 9   ! water
            ENDIF

            prm_diag%albdif_t(jc,jb,isub_lake) = csalb(ist)
          ENDDO

        ENDIF  ! llake





        !*****************************!
        !                             !
        !  Aggregate surface albedo   !
        !                             !
        !*****************************!

        ! Loop over ALL grid points


        !
        ! Aggregate surface albedo on all points
        !
        IF (ntiles_total == 1) THEN
 
          DO jc = i_startidx, i_endidx
            prm_diag%albdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
            ! For RRTM copy albdif to albnirdif and albvisdif 
            ! i.e. no distiction is made between shortwave, vis and nir albedo
            prm_diag%albvisdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
            prm_diag%albnirdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
          ENDDO

        ELSE ! aggregate fields over tiles

          prm_diag%albdif(i_startidx:i_endidx,jb)    = 0._wp

          DO jt = 1, ntiles_total+ntiles_water
            DO jc = i_startidx, i_endidx
              prm_diag%albdif(jc,jb) = prm_diag%albdif(jc,jb)      &
                &                    + prm_diag%albdif_t(jc,jb,jt) &
                &                    * ext_data%atm%frac_t(jc,jb,jt)
            ENDDO
          ENDDO

          DO jc = i_startidx, i_endidx
            ! For RRTM copy albdif to albnirdif and albvisdif 
            ! i.e. no distiction is made between shortwave, vis and nir albedo
            prm_diag%albvisdif(jc,jb) = prm_diag%albdif(jc,jb)
            prm_diag%albnirdif(jc,jb) = prm_diag%albdif(jc,jb)
          ENDDO
        ENDIF  ! ntiles_total = 1


      ELSE  ! surface model switched OFF


        DO jc = i_startidx, i_endidx
          
          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog%t_g(jc,jb) >= tf_salt ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF

          prm_diag%albdif(jc,jb) = csalb(ist)
          ! For RRTM copy albdif to albnirdif and albvisdif
          ! i.e. no distiction is made between shortwave, vis and nir albedo
          prm_diag%albvisdif(jc,jb) = csalb(ist)
          prm_diag%albnirdif(jc,jb) = csalb(ist)
          
        ENDDO

      ENDIF ! inwp_surface=1


      ! albvisdir, albnirdir only needed for RRTM 
      !
      ! compute black sky albedo from white sky albedo and solar zenith angle formula 
      ! as in Ritter-Geleyn's fesft. So far we do not distinguish between  
      ! visible and NIR spectral bands.
      DO jc = i_startidx, i_endidx

        prm_diag%albvisdir(jc,jb) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),prm_diag%albvisdif(jc,jb))

        ! no need to do the computation twice, since albvisdif=albnirdif=albdif
        ! Thus: just copy
        prm_diag%albnirdir(jc,jb) = prm_diag%albvisdir(jc,jb)
      ENDDO

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE sfc_albedo



  !>
  !! Calculation of surface albedo based on MODIS data
  !!
  !! Calculation of surface albedo based on snow-free MODIS data. In addition, 
  !! snow/ice conditions are taken into account. The snow-free MODIS albedo is updated 
  !! on a daily basis.
  !! We distinguish between
  !! - shortwave broadband albedo  (diffuse, 0.3-5.0�m): albdif
  !! - UV visible broadband albedo (diffuse, 0.3-0.7�m): albvisdif
  !! - near IR broadband albedo    (diffuse, 0.7-5.0�m): albnirdif 
  !!
  !! albvisdif and albnirdif are exclusively used by the RRTM scheme
  !!
  !! land points (snow-free): MODIS albedo is used with separate values for visible and 
  !!                          nera-IR spectral bands. No distinction for tiles
  !! land points (snow covered): based on land-class specific tabulated values with
  !!                             consideration of aging snow. No distinction between 
  !!                             visible and near-IR spectral bands, yet
  !! sea points: tabulated value. No distinction yet between visible and near-IR spectral 
  !!             bands.
  !! lake points: tabulated value. No distinction yet between visible and near-IR spectral 
  !!              bands.
  !! sea-ice points: based on an empirical formula taken from GME. No distinction yet 
  !!                 between visible and near-IR spectral bands.
  !! lake-ice points: based on an empirical formula taken from GME. No distinction yet 
  !!                 between visible and near-IR spectral bands.
  !!
  !! Possible improvements:
  !!=======================
  !! snow-covered land points: separate values for visible and near-IR spectral bands.
  !! sea-ice points: separate values for visible and near-IR spectral bands.
  !! lake-ice points: separate values for visible and near-IR spectral bands.
  !! sea/lake points: separate values for direct and diffuse radiation (see IFS)
  !!
  !! Note that when using separate values for visible and near-IR spectral bands, the 
  !! shortwave albedo albdif_t must be derived by spectral integration of the visible 
  !! and near-IR albedo.                         
  !! 
  !! @par Revision History
  !! Initial Revision by Daniel Reinert, DWD (2013-05-15)
  !! - Modification by Daniel Reinert, DWD (2013-07-03)
  !!   Albedo for lake-ice points based on an empirical formula proposed by 
  !!   proposed by Mironov and Ritter (2004)
  !! - Modification by Daniel Reinert, DWD (2013-08-07)
  !!   Added albedo for direct radiation (VIS and NIR spectral bands)
  !! - Modification by Daniel Reinert, DWD (2013-08-08)
  !!   Added albedo for direct radiation (now separate computation for VIS and NIR spectral bands)
  !!
  SUBROUTINE sfc_albedo_modis(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)

    TYPE(t_patch),          INTENT(   in):: pt_patch  !< grid/patch info.

    TYPE(t_external_data),  INTENT(   in):: ext_data  !< external data

    TYPE(t_lnd_prog),       INTENT(   in):: lnd_prog  !< land prognostic state (new)

    TYPE(t_wtr_prog),       INTENT(   in):: wtr_prog  !< water prognostic state (new)

    TYPE(t_lnd_diag),       INTENT(   in):: lnd_diag  !< land diagnostic state

    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag

    ! Local scalars:
    REAL(wp):: snow_frac               !< snow cover fraction
    REAL(wp):: zminsnow_alb            !< temperature-dependent minimum snow albedo
    REAL(wp):: zmaxsnow_alb            !< maximum snow albedo depending on landuse class
    REAL(wp):: zlimsnow_alb            !< upper limit snow albedo depending on snow depth and roughness length
    REAL(wp):: zsnowalb_lu             !< maximum snow albedo specified in landuse table
    REAL(wp):: t_fac                   !< factor for temperature dependency of zminsnow_alb over glaciers
    REAL(wp):: zsnow_alb               !< snow albedo

    ! Auxiliaries for tile-specific calculation of direct beam albedo
    REAL(wp):: zalbvisdir_t, zalbnirdir_t, zalbdirfac, zalbdirlim 

    INTEGER :: jg                      !< patch ID
    INTEGER :: jb, jc, ic, jt          !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: ist
    INTEGER :: ilu                     !< land cover class
    INTEGER :: i_count_lnd             !< number of land points
    INTEGER :: i_count_sea             !< number of sea points
    INTEGER :: i_count_flk             !< number of lake points
    INTEGER :: i_count_seaice          !< number of seaice points

    !-----------------------------------------------------------------------

    jg = pt_patch%id
    i_nchdom  = MAX(1,pt_patch%n_childdom)

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,jc,ic,i_startidx,i_endidx,ist,snow_frac,t_fac,               &
!$OMP            zsnow_alb,ilu,i_count_lnd,i_count_sea,i_count_flk,                 &
!$OMP            i_count_seaice,zminsnow_alb,zmaxsnow_alb,zlimsnow_alb,zsnowalb_lu, &
!$OMP            zalbvisdir_t,zalbnirdir_t,zalbdirfac,zalbdirlim) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk


      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                 i_startidx, i_endidx, rl_start, rl_end)


      !------------------------------------------------------------------------------
      ! Calculation of land surface albedo based on MODIS input data
      !------------------------------------------------------------------------------
      
      IF ( atm_phy_nwp_config(jg)%inwp_surface == 1 ) THEN


        !
        ! 1. Consider land points (may have tiles)
        !
        ! - loop over surface tiles
        ! - note that different grid points may have different numbers 
        !   of active tiles (1<=ntiles<=ntiles_total). Therefore each tile has a 
        !   separate index list.
        ! 
        DO jt = 1, ntiles_total

          i_count_lnd = ext_data%atm%gp_count_t(jb,jt)

          IF (i_count_lnd == 0) CYCLE ! skip loop if the index list for the given tile is empty

!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count_lnd

            jc = ext_data%atm%idx_lst_t(ic,jb,jt)

            ! current land-cover class
            ilu = ext_data%atm%lc_class_t(jc,jb,jt)

            ! maximum snow albedo specified in landuse table
            zsnowalb_lu = ABS(ext_data%atm%snowalb_lcc(ilu))

            IF (ext_data%atm%alb_dif(jc,jb) > csalb_snow_min .AND. ilu == ext_data%atm%i_lc_snow_ice) THEN
              ! temperature-dependent minimum snow albedo over glaciers (i.e. alb_dif = 0.7)
              ! this is needed to prevent unrealistically low albedos over Antarctica
              t_fac = MIN(1._wp,0.1_wp*(tmelt-lnd_prog%t_snow_t(jc,jb,jt)))
              zminsnow_alb = (1._wp-t_fac)*csalb_snow_min + t_fac*ext_data%atm%alb_dif(jc,jb)
            ELSE
              ! otherwise take minimum snow albedo as 60% of the landuse-class specific maximum snow albedo
              zminsnow_alb = MAX(0.4_wp*csalb_snow_min,MIN(csalb_snow_min,0.6_wp*zsnowalb_lu))
            ENDIF

            ! maximum snow albedo depending on landuse class if tile approach is used
            ! (without tiles, this would conflict with the snow albedo limit depending on the forest fraction)
            IF (ntiles_lnd > 1) THEN
              zmaxsnow_alb = MIN(csalb_snow_max,zsnowalb_lu)
              zlimsnow_alb = MIN(csalb_snow_max,SQRT(zsnowalb_lu))
            ELSE
              zmaxsnow_alb = csalb_snow_max
              zlimsnow_alb = csalb_snow_max
            ENDIF

            ! Further limitation of maximum snow albedo in case of thin snow cover depending on landuse 
            ! roughness length and SSO standard deviation; without tiles, it is assumed that non-forest
            ! vegetation is fully covered by snow if the snow depth exceeds three times the roughness length.
            ! The snow depth is taken to be at least 5 cm here because the effect of partial snow
            ! cover is considered below.
            !
            zmaxsnow_alb = MIN(zmaxsnow_alb, zlimsnow_alb *                                         &
              MIN(1._wp,SQRT(0.25_wp+0.25_wp*MAX(0.05_wp,lnd_diag%h_snow_t(jc,jb,jt))/              &
              MAX(MIN(0.5_wp,ext_data%atm%z0_lcc(ilu)),1.e-3_wp*ext_data%atm%sso_stdh_raw(jc,jb)) )))

            ! Consider effects of aging on solar snow albedo
            !
            zsnow_alb = zminsnow_alb + lnd_diag%freshsnow_t(jc,jb,jt)*(zmaxsnow_alb-zminsnow_alb)

            IF (ntiles_lnd == 1) THEN
              ! special treatment for forests
              ! - no landuse-class specific limitation of snow albedo
              ! - instead, snow albedo is limited as a function of the forest fraction
              zsnow_alb = zsnow_alb*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb))       &
                + csalb_snow_fe * ext_data%atm%for_e(jc,jb) + csalb_snow_fd * ext_data%atm%for_d(jc,jb)
            ENDIF

            ! snow cover fraction
            snow_frac = lnd_diag%snowfrac_t(jc,jb,jt)

            ! shortwave broadband surface albedo (white sky)
            prm_diag%albdif_t(jc,jb,jt) = snow_frac * zsnow_alb  &
              &                  + (1._wp - snow_frac)* ext_data%atm%alb_dif(jc,jb)

            ! UV visible broadband surface albedo (white sky)
            prm_diag%albvisdif_t(jc,jb,jt) = snow_frac * zsnow_alb  &
              &                  + (1._wp - snow_frac)* ext_data%atm%albuv_dif(jc,jb)

            ! near IR broadband surface albedo (white sky)
            prm_diag%albnirdif_t(jc,jb,jt) = snow_frac * zsnow_alb  &
              &                  + (1._wp - snow_frac)* ext_data%atm%albni_dif(jc,jb)


          ENDDO  ! ic

        ENDDO  !ntiles





        ! 2. Consider water points with/without seaice model
        !
        IF ( (lseaice) ) THEN  ! seaice model switched on

          !
          ! Open sea points
          !
          ! - loop over sea points  (points which are at least partly ice-free) 
          !
          i_count_sea = ext_data%atm%spw_count(jb)

          DO ic = 1, i_count_sea
            jc = ext_data%atm%idx_lst_spw(ic,jb)

            ist = 9  ! sea water

            prm_diag%albdif_t   (jc,jb,isub_water) = csalb(ist)
            prm_diag%albvisdif_t(jc,jb,isub_water) = csalb(ist)
            prm_diag%albnirdif_t(jc,jb,isub_water) = csalb(ist)
          ENDDO



          !
          ! Sea-ice points
          !
          ! - loop over sea-ice points (points which are at least partly ice-covered)
          !
          i_count_seaice = ext_data%atm%spi_count(jb)

          DO ic = 1, i_count_seaice
            jc = ext_data%atm%idx_lst_spi(ic,jb)

            ist = 10   ! sea-ice

            ! In case the sea ice model is used, compute ice albedo for seaice 
            ! points with an empirical formula taken from GME.
            ! The ice albedo is the lower the warmer, and therefore wetter 
            ! the ice is. Use ice temperature at time level nnew 
            ! (2-time level scheme in sea ice model).
            prm_diag%albdif_t(jc,jb,isub_seaice) = csalb(ist) * ( 1.0_wp - 0.3846_wp  &
              &                         * EXP(-0.35_wp*(tmelt-wtr_prog%t_ice(jc,jb))))
            ! gives alb_max = 0.70
            !       alb_min = 0.43
            ! compare with Mironov et. al (2012), Tellus
            !       alb_max = 0.65
            !       alb_min = 0.40

            prm_diag%albvisdif_t(jc,jb,isub_seaice) = prm_diag%albdif_t(jc,jb,isub_seaice)
            prm_diag%albnirdif_t(jc,jb,isub_seaice) = prm_diag%albdif_t(jc,jb,isub_seaice)
          ENDDO


        ELSE

          !
          ! Consider sea points (including sea-ice points)
          !
          ! - loop over sea points
          !
          i_count_sea = ext_data%atm%spw_count(jb)

          DO ic = 1, i_count_sea
            jc = ext_data%atm%idx_lst_spw(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_water) < tf_salt) THEN 
              ist = 10  ! sea ice
            ELSE
              ist = 9   ! sea water
            ENDIF

            prm_diag%albdif_t   (jc,jb,isub_water) = csalb(ist)
            prm_diag%albvisdif_t(jc,jb,isub_water) = csalb(ist)
            prm_diag%albnirdif_t(jc,jb,isub_water) = csalb(ist)
          ENDDO

        ENDIF



        ! 3. Consider lake points with/without lake model
        !
        IF (llake) THEN

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%fp_count(jb)

!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, i_count_flk
            jc = ext_data%atm%idx_lst_fp(ic,jb)


            !  In case ice is present at lake points, compute ice albedo 
            !  for lake points with an empirical formulation 
            !  proposed by Mironov and Ritter (2004) for use in GME 
            !  [ice_albedo=function(ice_surface_temperature)].
            !  Use surface temperature at time level "nnew".

            ! special handling for ice-covered lake points
            IF (wtr_prog%h_ice(jc,jb) > h_Ice_min_flk) THEN 

              prm_diag%albdif_t(jc,jb,isub_lake) = albedo_whiteice_ref                      &
                &              - (albedo_whiteice_ref - albedo_blueice_ref)                 &
                &              * EXP(-c_albice_MR*(tpl_T_f-lnd_prog%t_g_t(jc,jb,isub_lake)) &
                &              /tpl_T_f)
            ELSE
              ist = 9   ! water
              prm_diag%albdif_t(jc,jb,isub_lake) = csalb(ist)
            ENDIF

            prm_diag%albvisdif_t(jc,jb,isub_lake) = prm_diag%albdif_t(jc,jb,isub_lake)
            prm_diag%albnirdif_t(jc,jb,isub_lake) = prm_diag%albdif_t(jc,jb,isub_lake) 
          ENDDO


        ELSE

          !
          ! Lake points
          !
          ! - loop over lake points
          !
          i_count_flk = ext_data%atm%fp_count(jb)

          DO ic = 1, i_count_flk
            jc = ext_data%atm%idx_lst_fp(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_lake) < tpl_T_f) THEN 
              ist = 10 ! sea ice
            ELSE
              ist = 9  ! water
            ENDIF

            prm_diag%albdif_t   (jc,jb,isub_lake) = csalb(ist)
            prm_diag%albvisdif_t(jc,jb,isub_lake) = csalb(ist)
            prm_diag%albnirdif_t(jc,jb,isub_lake) = csalb(ist)
          ENDDO


        ENDIF





        !*****************************!
        !                             !
        !  Aggregate surface albedo   !
        !                             !
        !*****************************!

        ! Loop over ALL grid points


        !
        ! Aggregate surface albedo on all points
        !
        IF (ntiles_total == 1) THEN
 
          DO jc = i_startidx, i_endidx
            prm_diag%albdif(jc,jb) = prm_diag%albdif_t(jc,jb,1)
            ! albvisdif, albnirdif only needed for RRTM 
            prm_diag%albvisdif(jc,jb) = prm_diag%albvisdif_t(jc,jb,1)
            prm_diag%albnirdif(jc,jb) = prm_diag%albnirdif_t(jc,jb,1)
          ENDDO

        ELSE ! aggregate fields over tiles

          prm_diag%albdif   (i_startidx:i_endidx,jb) = 0._wp
          prm_diag%albvisdif(i_startidx:i_endidx,jb) = 0._wp
          prm_diag%albnirdif(i_startidx:i_endidx,jb) = 0._wp
          prm_diag%albvisdir(i_startidx:i_endidx,jb) = 0._wp
          prm_diag%albnirdir(i_startidx:i_endidx,jb) = 0._wp

          DO jt = 1, ntiles_total+ntiles_water
            DO jc = i_startidx, i_endidx

              ! Tile-specific direct beam albedo
              IF (ext_data%atm%frac_t(jc,jb,jt) > 0._wp) THEN

                zalbvisdir_t = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),prm_diag%albvisdif_t(jc,jb,jt))

                zalbnirdir_t = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),prm_diag%albnirdif_t(jc,jb,jt))

                ! Limit direct beam albedo to diffuse albedo for landuse classes with "rough" vegetation
                ! and in mountainous regions
                ilu = ext_data%atm%lc_class_t(jc,jb,jt)
                IF (jt <= ntiles_total) THEN

                  ! Full limitation is applied for a roughness length >= 15 cm or an SSO standard deviation >= 150 m
                  zalbdirfac = MAX(0.01_wp*(ext_data%atm%sso_stdh_raw(jc,jb) - 50._wp), &
                                   10._wp *(ext_data%atm%z0_lcc(ilu) - 0.05_wp)         )
                  zalbdirfac = MIN(1._wp, MAX(0._wp, zalbdirfac))

                  zalbdirlim = zalbdirfac*prm_diag%albvisdif_t(jc,jb,jt) + (1._wp-zalbdirfac)*zalbvisdir_t
                  zalbvisdir_t = MIN(zalbvisdir_t,zalbdirlim)
                  zalbdirlim = zalbdirfac*prm_diag%albnirdif_t(jc,jb,jt) + (1._wp-zalbdirfac)*zalbnirdir_t
                  zalbnirdir_t = MIN(zalbnirdir_t,zalbdirlim)

                ENDIF

              ELSE
                zalbvisdir_t = 0._wp
                zalbnirdir_t = 0._wp
              ENDIF


              prm_diag%albdif(jc,jb) = prm_diag%albdif(jc,jb)         &
                &                    + prm_diag%albdif_t(jc,jb,jt)    &
                &                    * ext_data%atm%frac_t(jc,jb,jt)

              ! albvisdif, albnirdif only needed for RRTM 
              prm_diag%albvisdif(jc,jb) = prm_diag%albvisdif(jc,jb)   &
                &                    + prm_diag%albvisdif_t(jc,jb,jt) &
                &                    * ext_data%atm%frac_t(jc,jb,jt)

              prm_diag%albnirdif(jc,jb) = prm_diag%albnirdif(jc,jb)   &
                &                    + prm_diag%albnirdif_t(jc,jb,jt) &
                &                    * ext_data%atm%frac_t(jc,jb,jt)

              prm_diag%albvisdir(jc,jb) = prm_diag%albvisdir(jc,jb)                  &
                &                       + zalbvisdir_t * ext_data%atm%frac_t(jc,jb,jt)

              prm_diag%albnirdir(jc,jb) = prm_diag%albnirdir(jc,jb)                  &
                &                       + zalbnirdir_t * ext_data%atm%frac_t(jc,jb,jt)

            ENDDO
          ENDDO

        ENDIF  ! ntiles_total = 1


      ELSE  ! surface model switched OFF


        DO jc = i_startidx, i_endidx
          
          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog%t_g(jc,jb) >= tf_salt ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF

          prm_diag%albdif(jc,jb) = csalb(ist)
          ! albvisdif, albnirdif only needed for RRTM 
          prm_diag%albvisdif(jc,jb) = csalb(ist)
          prm_diag%albnirdif(jc,jb) = csalb(ist)
          
        ENDDO

      ENDIF ! inwp_surface=1


      ! albvisdir, albnirdir only needed for RRTM 
      !
      ! compute black sky albedo from white sky albedo and solar zenith angle formula 
      ! as in Ritter-Geleyn's fesft.
      IF (ntiles_total == 1 .OR. atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
        DO jc = i_startidx, i_endidx

          prm_diag%albvisdir(jc,jb) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),prm_diag%albvisdif(jc,jb))

          prm_diag%albnirdir(jc,jb) = sfc_albedo_dir_rg(prm_diag%cosmu0(jc,jb),prm_diag%albnirdif(jc,jb))

        ENDDO
      ENDIF


    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE sfc_albedo_modis




  !>
  !! Surface albedo (SA) for direct beam
  !!
  !! Surface albedo for direct beam, according to Ritter-Geleyn 
  !! radiation scheme.
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2015-06-22)
  !!
  FUNCTION sfc_albedo_dir_rg (cosmu0, alb_dif)  RESULT (alb_dir)
    !
    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)

    REAL(wp), INTENT(IN) :: alb_dif          !< diffuse albedo (NIR or VIS or broadband)

    REAL(wp) :: alb_dir
  !------------------------------------------------------------------------------

     alb_dir = MIN(0.999_wp,( 1.0_wp                             &
       &     +  0.5_wp * (cosmu0 * (1.0_wp/alb_dif - 1.0_wp)))   &
       &     / (1.0_wp + (cosmu0 * (1.0_wp/alb_dif - 1.0_wp)))**2)

  END FUNCTION sfc_albedo_dir_rg


!!$  !>
!!$  !! Surface albedo (SA) for direct beam over land after Zaengl
!!$  !!
!!$  !! Surface albedo (SA) for direct beam over land after Zaengl.
!!$  !! It builds upon the parameterization used in the 
!!$  !! Ritter-Geleyn radiation scheme. However, the derived direct 
!!$  !! albedo is not allowed to exceed the corresponding diffuse 
!!$  !! albedo over 'rough' surfaces. Surfaces are regarded as 'rough', 
!!$  !! if z0>= 15 cm or SSO standard deviation >= 150 m. 
!!$  !!
!!$  !!
!!$  !! @par Revision History
!!$  !! Initial revision by Daniel Reinert, DWD (2015-06-22)
!!$  !!
!!$  FUNCTION sfc_albedo_dir_zaengl (cosmu0, alb_dif, z0, sso_stdh)  RESULT (alb_dir)
!!$    !
!!$    REAL(wp), INTENT(IN) :: cosmu0           !< cosine of solar zenith angle (SZA)
!!$
!!$    REAL(wp), INTENT(IN) :: alb_dif          !< diffuse albedo (NIR or VIS or broadband)
!!$
!!$    REAL(wp), INTENT(IN) :: z0               !< surface roughness length
!!$
!!$    REAL(wp), INTENT(IN) :: sso_stdh         !< Standard deviation of sub-grid scale orography
!!$
!!$    REAL(wp) :: alb_dir
!!$
!!$    REAL(wp) :: zalb_dir                     !< unlimited direct albedo
!!$    REAL(wp) :: zalbdirlim                   !< limit, which the computed albedo must not exceed
!!$    REAL(wp) :: zalbdirfac                   !< factor for limit computation
!!$
!!$  !------------------------------------------------------------------------------
!!$
!!$    ! unlimited direct albedo according to Ritter-Geleyn's formulation
!!$    zalb_dir = sfc_albedo_dir_rg(cosmu0, alb_dif)
!!$
!!$
!!$    ! Full limitation is applied for a roughness length >= 15 cm or an SSO standard deviation >= 150 m
!!$    zalbdirfac = MAX(0.01_wp*(sso_stdh - 50._wp), 10._wp *(z0 - 0.05_wp) )
!!$    zalbdirfac = MIN(1._wp, MAX(0._wp, zalbdirfac))
!!$
!!$    zalbdirlim = zalbdirfac*alb_dif + (1._wp-zalbdirfac)*zalb_dir
!!$    alb_dir     = MIN(zalb_dir, zalbdirlim)
!!$
!!$  END FUNCTION sfc_albedo_dir_zaengl

END MODULE mo_albedo

