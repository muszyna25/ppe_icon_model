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
!! Initial Revision by daniel Reinert, DWD (2012-03-19)
!! Moved to a central place from mo_nwp_rad_interface and 
!! mo_nwp_rrtm_interface.
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

MODULE mo_albedo

  USE mo_kind,                 ONLY: wp
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_radiation_config,     ONLY: rad_csalbw
  USE mo_lnd_nwp_config,       ONLY: ntiles_total, ntiles_water, lseaice,     &
    &                                isub_water, isub_seaice
  USE mo_phyparam_soil,        ONLY: csalb, csalb_snow_fe, csalb_snow_fd,     &
    &                                csalb_snow_min, csalb_snow_max, cf_snow, &
    &                                csalb_p
  USE mo_physical_constants,   ONLY: tmelt, tf_salt
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_impl_constants,       ONLY: min_rlcell_int

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: sfc_albedo

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS


  !>
  !! Calculation of surface albedo
  !!
  !! Calculation of surface albedo taking soil type, vegetation 
  !! and snow/ice conditions into account
  !!
  !! @par Revision History
  !! Initial Revision by Thorsten Reinhardt, AGeoBw, Offenbach
  !! Modification by Daniel Reinert, DWD (2012-03-19)
  !! - Moved here from mo_nwp_rad_interface and mo_nwp_rrtm_interface.
  !!   Adaption to TERRA-tile approach.
  !!
  SUBROUTINE sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)

    TYPE(t_patch),          INTENT(   in):: pt_patch  !< grid/patch info.

    TYPE(t_external_data),  INTENT(   in):: ext_data  !< external data

    TYPE(t_lnd_prog),       INTENT(   in):: lnd_prog  !< land prognostic state (new)

    TYPE(t_wtr_prog),       INTENT(   in):: wtr_prog  !< water prognostic state (new)

    TYPE(t_lnd_diag),       INTENT(   in):: lnd_diag  !< land diagnostic state

    TYPE(t_nwp_phy_diag),   INTENT(inout):: prm_diag

    ! Local scalars:
    REAL(wp):: zvege, zsnow, zsalb_snow, zsnow_alb

    INTEGER :: jg                      !< patch ID
    INTEGER :: jb, jc, ic, jt          !< loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    INTEGER :: ist
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
!$OMP DO PRIVATE(jb,jt,jc,ic,i_startidx,i_endidx,ist,zvege,zsnow, &
!$OMP            zsalb_snow,zsnow_alb,i_count_lnd,i_count_sea,    &
!$OMP            i_count_flk,i_count_seaice) ICON_OMP_DEFAULT_SCHEDULE
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
            prm_diag%albvisdif_t(jc,jb,jt) = csalb(ist)&
              &                         - rad_csalbw(ist)*lnd_prog%w_so_t(jc,1,jb,jt)


            ! Account for Snow cover and vegetation
            !---------------------------------------
            zvege= 0.0_wp
            zsnow= 0.0_wp

            ! consider effects of aging on solar snow albedo
            !
            zsalb_snow = csalb_snow_min + &
              & lnd_diag%freshsnow_t(jc,jb,jt)*(csalb_snow_max-csalb_snow_min)
            zsnow_alb = zsalb_snow*(1._wp-ext_data%atm%for_e(jc,jb)-ext_data%atm%for_d(jc,jb)) &
              + csalb_snow_fe * ext_data%atm%for_e(jc,jb)                       &
              + csalb_snow_fd * ext_data%atm%for_d(jc,jb)

            ! account for snow cover and plant cover and compute final solar
            ! snow albedo
            !
            zvege = ext_data%atm%plcov_t(jc,jb,jt)
            IF (lnd_prog%w_snow_t(jc,jb,jt) > 0.0_wp) THEN
              zsnow = MIN(1.0_wp, lnd_prog%w_snow_t(jc,jb,jt) / cf_snow)
            ENDIF

            prm_diag%albvisdif_t(jc,jb,jt) = zsnow * zsnow_alb            &
              &  + (1.0_wp - zsnow) * (zvege * csalb_p + (1.0_wp - zvege) &
              &  * prm_diag%albvisdif_t(jc,jb,jt))

          ENDDO

        ENDDO  !ntiles



        ! Different treatment of water points with/without seaice model
        !
        IF ( (lseaice) ) THEN  ! seaice model switched on

          !
          ! 2a. Consider (open) sea points
          !
          ! - loop over sea points
          !
          i_count_sea = ext_data%atm%spw_count(jb)
          jt = isub_water

          DO ic = 1, i_count_sea
            jc = ext_data%atm%idx_lst_spw(ic,jb)

            ist = 9  ! sea water

            prm_diag%albvisdif_t(jc,jb,jt) = csalb(ist)
          ENDDO



          !
          ! 3a. Consider lake points (no tiles)
          !
          ! - loop over lake points (same jt as water points)
          !
          i_count_flk = ext_data%atm%fp_count(jb)

          DO ic = 1, i_count_flk
            jc = ext_data%atm%idx_lst_fp(ic,jb)

            ist = 9  ! sea water

            prm_diag%albvisdif_t(jc,jb,jt) = csalb(ist)
          ENDDO



          !
          ! 4. Consider sea-ice points
          !
          ! - loop over sea-ice points
          !
          i_count_seaice = ext_data%atm%spi_count(jb)
          jt = isub_seaice

          DO ic = 1, i_count_seaice
            jc = ext_data%atm%idx_lst_spi(ic,jb)

            ist = 10   ! seaice

            ! In case the sea ice model is used, compute ice albedo for seaice 
            ! points with an empirical formula taken from GME.
            ! The ice albedo is the lower the warmer, and therefore wetter 
            ! the ice is. Use ice temperature at time level nnew 
            ! (2-time level scheme in sea ice model).
            prm_diag%albvisdif_t(jc,jb,jt) = csalb(ist) * ( 1.0_wp - 0.3846_wp    &
              &                            * EXP(-0.35_wp*(tmelt-wtr_prog%t_ice(jc,jb))))
            ! gives alb_max = 0.70
            !       alb_min = 0.43
            ! compare with Mironov et. al (2012), Tellus
            !       alb_max = 0.65
            !       alb_min = 0.40
          ENDDO


        ELSE   ! no seaice model

          !
          ! 2b. Consider sea points
          !
          ! - loop over sea points
          !
          i_count_sea = ext_data%atm%spw_count(jb)
          jt = isub_water

          DO ic = 1, i_count_sea
            jc = ext_data%atm%idx_lst_spw(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_water) < tf_salt) THEN 
              ist = 10  ! sea ice
            ELSE
              ist = ext_data%atm%soiltyp(jc,jb)
            ENDIF

            prm_diag%albvisdif_t(jc,jb,jt) = csalb(ist)
          ENDDO



          !
          ! 3b. Consider lake points (no tiles)
          !
          ! - loop over lake points (same jt as water points)
          !
          i_count_flk = ext_data%atm%fp_count(jb)

          DO ic = 1, i_count_flk
            jc = ext_data%atm%idx_lst_fp(ic,jb)

            ! special handling of sea ice points
            IF (lnd_prog%t_g_t(jc,jb,isub_water) < tf_salt) THEN ! sea ice
              ist = 10
            ELSE
              ist = 9 ! water
            ENDIF

            prm_diag%albvisdif_t(jc,jb,jt) = csalb(ist)
          ENDDO

        ENDIF  ! lseaice





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
            prm_diag%albvisdif(jc,jb) = prm_diag%albvisdif_t(jc,jb,1)
          ENDDO

        ELSE ! aggregate fields over tiles

          prm_diag%albvisdif(i_startidx:i_endidx,jb) = 0._wp

          DO jt = 1, ntiles_total+ntiles_water
            DO jc = i_startidx, i_endidx
              prm_diag%albvisdif(jc,jb) = prm_diag%albvisdif(jc,jb)      &
                &                       + prm_diag%albvisdif_t(jc,jb,jt) &
                &                       * ext_data%atm%frac_t(jc,jb,jt)
            ENDDO
          ENDDO

        ENDIF  ! ntiles_total = 1


      ELSE  ! surface model switched OFF


        DO jc = i_startidx, i_endidx
          
          ist = 10

          IF (ext_data%atm%llsm_atm_c(jc,jb) .OR. lnd_prog%t_g(jc,jb) >= tf_salt ) THEN
            ist = ext_data%atm%soiltyp(jc,jb) ! water (ist=9) and sea ice (ist=10) included
          ENDIF

          prm_diag%albvisdif(jc,jb) = csalb(ist)
          
        ENDDO

      ENDIF ! inwp_surface=1




!      IF (atm_phy_nwp_config(jg)%lseaice) THEN
!        DO jc = i_startidx,i_endidx
!          ! In case the sea ice model is used AND water point AND ice is present,
!          ! compute ice albedo for water points with an empirical formula taken from GME.
!          ! The ice albedo is the lower the warmer, and therefore wetter, the ice is.
!          ! Use ice temperature at time level nnow (2-time level scheme in sea ice model).
!!          IF ((.NOT. llandmask(i,j)) .AND. (h_ice(i,j,nnow) > 0.0_ireals))               &
!!            zalso(i,j) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(tmelt-t_ice(i,j,nnow)))) &
!!            * csalb(10)
!          IF (( .NOT. ext_data%atm%llsm_atm_c(jc,jb)) .AND. &
!            & (prm_diag%h_ice(i,j,nnow) > 0.0_wp))          &
!            prm_diag%albvisdif(jc,jb) = (1.0_wp-0.3846_wp*EXP(-0.35_wp*(tmelt-t_ice(i,j,nnow)))) &
!            * csalb(10)
!        ENDDO
!      ENDIF

      !lake model not yet implemented
!      IF (atm_phy_nwp_config(jg)%llake) THEN
!        DO jc = i_startidx,i_endidx
!          IF((ext_data%atm%depth_lk(jc,jb)      >  0.0_wp) .AND.    &
!            (prm_diag%h_ice(jc,jb) >= h_Ice_min_flk) ) THEN
!            !  In case the lake model FLake is used AND lake point AND ice is present,
!            !  compute ice albedo for lake points with an empirical formulation 
!            !  proposed by Mironov and Ritter (2004) for use in GME 
!            !  [ice_albedo=function(ice_surface_temperature)].
!            !  Use surface temperature at time level "nnow".
!
!            prm_diag%albvisdif(jc,jb) = EXP(-c_albice_MR*(tpl_T_f-t_s(i,j,nnow))/tpl_T_f)
!            prm_diag%albvisdif(jc,jb) = albedo_whiteice_ref * (1._ireals-zalso(i,j)) +      &
!              albedo_blueice_ref  * prm_diag%albvisdif(jc,jb)
!          ENDIF
!        ENDDO
!      ENDIF


    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  END SUBROUTINE sfc_albedo

END MODULE mo_albedo

