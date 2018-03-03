!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
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
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_echam_phy_main

  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message
  USE mo_physical_constants  ,ONLY: cpd, cpv, cvd, cvv, Tf, tmelt
  USE mo_run_config          ,ONLY: nlev, nlevp1,    &
    &                               iqv, iqc, iqi
  USE mo_echam_phy_config    ,ONLY: echam_phy_tc, dt_zero
  USE mo_echam_cld_config    ,ONLY: echam_cld_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend,      &
    &                               cdimissval
  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop,  &
    &                               timer_cover
  USE mtime                  ,ONLY: datetime, isCurrentEventActive, &
    &                               OPERATOR(==), OPERATOR(<=), OPERATOR(>)
  USE mo_echam_sfc_indices   ,ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_cover               ,ONLY: cover

  USE mo_interface_echam_rad ,ONLY: interface_echam_rad
  USE mo_interface_echam_par ,ONLY: interface_echam_par
  USE mo_interface_echam_rht ,ONLY: echam_rht
  USE mo_interface_echam_vdf ,ONLY: echam_vdf
  USE mo_interface_echam_car ,ONLY: interface_echam_car
  USE mo_interface_echam_cnv ,ONLY: echam_cnv
  USE mo_interface_echam_gwd ,ONLY: echam_gwd
  USE mo_interface_echam_sso ,ONLY: echam_sso
  USE mo_interface_echam_cld ,ONLY: echam_cld
  USE mo_interface_echam_mox ,ONLY: echam_mox

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch

#ifdef __ICON_ART
  USE mo_ext_data_state         ,ONLY: ext_data
  USE mo_linked_list            ,ONLY: t_var_list
  USE mo_nonhydro_types         ,ONLY: t_nh_diag
  USE mo_nonhydro_types         ,ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_art_reaction_interface ,ONLY: art_reaction_interface
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_main

CONTAINS

  !>
  !!
  SUBROUTINE echam_phy_main( patch             &
    &                       ,rl_start, rl_end  &
    &                       ,datetime_old      &
    &                       ,pdtime            &
#ifdef __ICON_ART
    &                       ,p_prog_list       &
    &                       ,pt_prog_new       &
    &                       ,p_metrics         &
    &                       ,pt_diag           &
#endif
    &                                          )


    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    INTEGER                ,INTENT(in) :: rl_start, rl_end
    TYPE(datetime)         ,POINTER    :: datetime_old     !< old date and time (at start of this step)
    REAL(wp)               ,INTENT(in) :: pdtime           !< time step

#ifdef __ICON_ART
    TYPE(t_var_list)       ,INTENT(in) :: p_prog_list      !< current prognostic state list
    TYPE(t_nh_prog)        ,INTENT(in) :: pt_prog_new
    TYPE(t_nh_metrics)     ,INTENT(in) :: p_metrics        !< NH metrics state
    TYPE(t_nh_diag)        ,INTENT(in) :: pt_diag          !< list of diagnostic fields
#endif
    
    ! Local variables
    TYPE(t_echam_phy_field),POINTER    :: field
    TYPE(t_echam_phy_tend) ,POINTER    :: tend

    REAL(wp) :: zq_rlw_impl (nproma)                       !< additional heating by LW rad. due to implicit coupling
                                                           !  in surface energy balance [W/m2]

    INTEGER  :: ictop (nproma,patch%nblks_c)               !< from massflux

    INTEGER  :: jg                                         !< grid level/domain index
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb                                         !< block index
    INTEGER  :: jcs, jce                                   !< start/end column index within this block

    LOGICAL  :: is_in_sd_ed_interval                       !< time is in process interval [sd,ed[
    LOGICAL  :: is_active                                  !< process is active

    jg = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    ! 1. Associate pointers
    field => prm_field(jg)
    tend  => prm_tend (jg)

    field%ictop(:,:)   = nlev-1

    ! initialize physics accumulated heating
    field% q_phy(:,:,:) = 0._wp

    ! initialize output fields of unused parameterizations,
    ! which may contain values from the restart file
    ! (this block needs to be rearrranged later)
    !
    IF ( echam_phy_tc(jg)%dt_rad == dt_zero ) THEN
       field% rsd_rt      (:,:,:) = 0.0_wp
       field% rsu_rt      (:,:,:) = 0.0_wp
       field% rsdcs_rt    (:,:,:) = 0.0_wp
       field% rsucs_rt    (:,:,:) = 0.0_wp
       field% rvds_dir_rt (:,  :) = 0.0_wp
       field% rpds_dir_rt (:,  :) = 0.0_wp
       field% rnds_dir_rt (:,  :) = 0.0_wp
       field% rvds_dif_rt (:,  :) = 0.0_wp
       field% rpds_dif_rt (:,  :) = 0.0_wp
       field% rnds_dif_rt (:,  :) = 0.0_wp
       field% rvus_rt     (:,  :) = 0.0_wp
       field% rpus_rt     (:,  :) = 0.0_wp
       field% rnus_rt     (:,  :) = 0.0_wp
       field% rld_rt      (:,:,:) = 0.0_wp
       field% rlu_rt      (:,:,:) = 0.0_wp
       field% rldcs_rt    (:,:,:) = 0.0_wp
       field% rlucs_rt    (:,:,:) = 0.0_wp
    END IF
    !
    IF ( echam_phy_tc(jg)%dt_cnv == dt_zero ) THEN
       field% rsfc  (:,:) = 0.0_wp
       field% ssfc  (:,:) = 0.0_wp
       field% rtype (:,:) = 0.0_wp
    END IF
    !
    IF ( echam_phy_tc(jg)%dt_cld == dt_zero ) THEN
       field% rsfl (:,:) = 0.0_wp
       field% ssfl (:,:) = 0.0_wp
    END IF

    !------------------------------------------------------------
    ! 3. COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
    !------------------------------------------------------------

    ! 3.2b Specific heat of moist air
    !
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_cpair_cvair_qconv(jb,jcs,jce, field)
    END DO
!$OMP END PARALLEL DO 

    !-------------------------------------------------------------------
    ! Calculate surface fraction 
    !------------------------------------------------------------------- 
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL surface_fractions( jb,jcs,jce, nproma, field)
    END DO
!$OMP END PARALLEL DO 
 
    !-------------------------------------------------------------------
    ! Cloud cover (diagnostic)
    !------------------------------------------------------------------- 
    IF ( echam_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       IF   ( (echam_phy_tc(jg)%sd_cld <= datetime_old) .AND. &
            & (echam_phy_tc(jg)%ed_cld >  datetime_old) .AND. &
            & (echam_phy_tc(jg)%dt_cld >  dt_zero     )       ) THEN
          !
!$OMP PARALLEL DO PRIVATE(jcs,jce)
          DO jb = i_startblk,i_endblk
             CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
             CALL cloud_cover(jg,jb,jcs,jce, nproma, field)
          END DO
!$OMP END PARALLEL DO 
          !
       END IF
       !
    END IF

    !---------------------------------------------------------------------
    ! 3.9 DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
    !     (Needed only for sub-models. Note: sequence of arguments
    !      different from the original ECHAM6)
    !---------------------------------------------------------------------
    !
    !CALL WMO_tropopause( jce, nproma, nlev,         &! in
    !                   & ncctop, nccbot, lresum,   &! in
    !                   & field% ta(:,:,jb),        &! in
    !                   & field% presm_old(:,:,jb), &! in
    !                   & field% tropo(:,jb),       &! out for diagnostics
    !                   & itrpwmo, itrpwmop1        )! out for submodel
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! 3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
    !      (1/M**3) USED IN RADLSW AND CLOUD
    !---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_droplet_number(jg, jb, jcs,jce, nproma, field)
    END DO
!$OMP END PARALLEL DO 


    !-------------------------------------------------------------------
    ! Radiation (one interface for LW+SW)
    !-------------------------------------------------------------------
    !
    IF ( (echam_phy_tc(jg)%dt_rad > dt_zero) ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_rad <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_rad >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_rad,   datetime_old)
       !
       CALL message_forcing_action('LW and SW radiation (rad)',   &
            &                      is_in_sd_ed_interval, is_active)
       !
       ! radiative fluxes
       CALL interface_echam_rad(is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old                     )
       !
       ! radiative heating
       !
       is_active = .TRUE. !< always adjust rad heating
       !
       CALL interface_echam_par(echam_rht,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
       !
    END IF


    !-------------------------------------------------------------------
    ! Vertical diffusion, boundary layer and surface
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_vdf >  dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_vdf <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_vdf >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_vdf,   datetime_old)
       !
       CALL message_forcing_action('vertical diffusion (vdf)',    &
            &                      is_in_sd_ed_interval, is_active)
       !
       CALL interface_echam_par(echam_vdf,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
       !
    END IF
    !
    !-----------------------

    !-----------------------
    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
!$OMP PARALLEL DO PRIVATE(jcs,jce, zq_rlw_impl)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        ! Heating due to the fact that surface model only used part of longwave radiation to compute new surface temperature
        zq_rlw_impl(jcs:jce) =                                                  &
          &   ( (field%rld_rt(jcs:jce,nlev,jb)-field%rlu_rt(jcs:jce,nlev,jb))   & ! (( rlns from "radiation"
          &    -(field%rlds  (jcs:jce,jb)     -field%rlus  (jcs:jce,jb)     ))  & !   -rlns from "radheating" and "update_surface")
          &  -field%q_rlw(jcs:jce,nlev,jb)                                                ! old heating from radheat

        ! Heating accumulated
        field% q_phy(jcs:jce,nlev,jb) = field% q_phy(jcs:jce,nlev,jb) + zq_rlw_impl(jcs:jce)

        ! Tendency
        tend%ta_rlw_impl(jcs:jce,jb) = zq_rlw_impl(jcs:jce) * field% qconv(jcs:jce,nlev,jb)

        ! Tendencies accumulated
        tend%ta_phy(jcs:jce,nlev,jb) = tend%ta_phy(jcs:jce,nlev,jb) + tend%ta_rlw_impl(jcs:jce,jb)
      END DO
!$OMP END PARALLEL DO 

    END IF
    !---------------------


    !-------------------------------------------------------------------
    ! Linearized ozone chemistry of Cariolle
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_car > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_car <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_car >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_car,   datetime_old)
       !
       CALL message_forcing_action('lin. Cariolle ozone chem. (car)', &
            &                      is_in_sd_ed_interval, is_active    )
       !
       CALL interface_echam_car(is_in_sd_ed_interval, is_active, &
            &                   patch, rl_start, rl_end,         &
            &                   datetime_old, pdtime             )
       !
    END IF

#ifdef __ICON_ART
    IF (echam_phy_tc(jg)%dt_art > dt_zero) THEN

      is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_art <= datetime_old) .AND. &
           &                          (echam_phy_tc(jg)%ed_art >  datetime_old)
      is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_art,   datetime_old)

      CALL message_forcing_action('ART (rad)',                   &
            &                     is_in_sd_ed_interval, is_active)

       CALL art_reaction_interface(ext_data(jg),           & !> in
            &                      patch,                  & !> in
            &                      datetime_old,           & !> in
            &                      pdtime,                 & !> in
            &                      p_prog_list,            & !> in
            &                      pt_prog_new,            &
            &                      p_metrics,              & !> in
            &                      pt_diag,                & !> inout
            &                      field%qtrc,             &
            &                      tend = tend)
    ENDIF
#endif
    
    !-------------------------------------------------------------------
    ! Atmospheric gravity wave drag
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_gwd > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_gwd <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_gwd >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_gwd,   datetime_old)
       !
       CALL message_forcing_action('Hines gravity wave drag (gwd)', &
            &                      is_in_sd_ed_interval, is_active  )
       !
       CALL interface_echam_par(echam_gwd,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
       !
    END IF

    !-------------------------------------------------------------------
    ! Sub grid scale orographic effects: blocking and orog. gravit waves
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_sso <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_sso >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_sso,   datetime_old)
       !
       CALL message_forcing_action('sub grid scale orographic effects (sso)', &
            &                      is_in_sd_ed_interval, is_active            )
       !
       CALL interface_echam_par(echam_sso,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
       !
    END IF

    !-------------------------------------------------------------------
    ! Cumulus convection
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_cnv <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_cnv >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_cnv,   datetime_old)
       !
       CALL message_forcing_action('cumulus convection (cnv)',    &
            &                      is_in_sd_ed_interval, is_active)
       !
       CALL interface_echam_par(echam_cnv,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
       !
    END IF

    !-------------------------------------------------------------
    ! Update provisional physics state
    field% qtrc(:,:,:,iqc) = field% qtrc(:,:,:,iqc) + tend% qtrc_phy(:,:,:,iqc)*pdtime
    field% qtrc(:,:,:,iqi) = field% qtrc(:,:,:,iqi) + tend% qtrc_phy(:,:,:,iqi)*pdtime

    !-------------------------------------------------------------------
    ! Cloud processes
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_cld <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_cld >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_cld,   datetime_old)
       !
       CALL message_forcing_action('cloud microphysics (cld)',    &
            &                      is_in_sd_ed_interval, is_active)
       !
       CALL interface_echam_par(echam_cld,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
       !
    END IF

    !-------------------------------------------------------------------
    IF ( (echam_phy_tc(jg)%dt_mox > dt_zero) ) THEN
      !
      is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_mox <= datetime_old) .AND. &
           &                          (echam_phy_tc(jg)%ed_mox >  datetime_old)
      is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_mox,   datetime_old)
       !
       CALL message_forcing_action('parameterized methane oxidation (mox)', &
            &                      is_in_sd_ed_interval, is_active          )
       !
       CALL interface_echam_par(echam_mox,                       &
            &                   is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old, pdtime             )
    END IF


!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      !  total precipitation flux
      field% totprec (jcs:jce,jb)     =  field% rsfl (jcs:jce,jb) & ! rain large scale
            &                           +field% ssfl (jcs:jce,jb) & ! snow large scale
            &                           +field% rsfc (jcs:jce,jb) & ! rain convection
            &                           +field% ssfc (jcs:jce,jb)   ! snow convection

      ! vertical integral
      field% q_phy_vi(jcs:jce,jb) = SUM(field% q_phy(jcs:jce,:,jb),DIM=2)

      ! now convert the temperature tendency from physics, as computed for constant pressure conditions,
      ! to constant volume conditions, as needed for the coupling to the dynamics
      tend% ta_phy (jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb)*field%cpair(jcs:jce,:,jb)/field%cvair(jcs:jce,:,jb)
    END DO
!$OMP END PARALLEL DO 

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE echam_phy_main
  !----------------------------------------------------------------
  
  !----------------------------------------------------------------
  SUBROUTINE surface_fractions( jb,jcs,jce, nbdim, field)

    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box

    INTEGER  :: jc
 
 
    ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    DO jc=jcs,jce

      ! fraction of land in the grid box.
      ! lsmask: land-sea mask, depends on input data, either:
      ! fractional, including or excluding lakes in the land part or
      ! non-fractional, each grid cell is either land, sea, or sea-ice.
      ! See mo_echam_phy_init or input data set for details.

      zfrl(jc) = MAX(field% lsmask(jc,jb),0._wp)

      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box

      zfrw(jc) = MAX(1._wp-zfrl(jc),0._wp)*MAX(1._wp-(field%seaice(jc,jb)+field%lake_ice_frc(jc,jb)),0._wp)

      ! fraction of sea ice in the grid box
      zfri(jc) = MAX(1._wp-zfrl(jc)-zfrw(jc),0._wp)
      ! security for ice temperature with changing ice mask
      !
      IF(zfri(jc) > 0._wp .AND. field%ts_tile(jc,jb,iice) == cdimissval ) THEN
         field% ts_tile(jc,jb,iice)  = tmelt + Tf    ! = 271.35 K
      END IF
    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) field%frac_tile(jcs:jce,jb,ilnd) = zfrl(jcs:jce)
    IF (iwtr.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iwtr) = zfrw(jcs:jce)
    IF (iice.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iice) = zfri(jcs:jce)

  END SUBROUTINE surface_fractions
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  SUBROUTINE cloud_cover(jg,jb,jcs,jce, nbdim, field)

    INTEGER         ,INTENT(IN) :: jg             !< grid  index
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    INTEGER  :: itype(nbdim)              !< type of convection

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !-------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_cover)

    itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))
      
    CALL cover(    jg,                        &! in
         &         jce, nbdim,                &! in
         &         nlev, nlevp1,              &! in
         &         itype,                     &! zfrw, zfri,       &! in
         &         field%frac_tile(:,jb,iwtr),&
         &         field%frac_tile(:,jb,iice),&
         &         field% zf(:,:,jb),         &! in
         &         field% presi_old(:,:,jb),  &! in
         &         field% presm_old(:,:,jb),  &! in
         &         field%  ta(:,:,jb),        &! in    tm1
         &         field%  qtrc(:,:,jb,iqv),  &! in    qm1
         &         field%  qtrc(:,:,jb,iqi),  &! in    xim1
         &         field%  aclc(:,:,jb),      &! out   (for "radiation" and "vdiff_down")
         &         field% rintop(:,  jb)     ) ! out   (for output)

    IF (ltimer) CALL timer_stop(timer_cover)

  END SUBROUTINE cloud_cover
  !----------------------------------------------------------------

  !---------------------------------------------------------------------
  ! 3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
  !      (1/M**3) USED IN RADLSW AND CLOUD
  !---------------------------------------------------------------------
  SUBROUTINE calculate_droplet_number(jg, jb, jcs, jce, nbdim, field)
    INTEGER         ,INTENT(IN) :: jg             !< grid  index
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    INTEGER  :: jc, jk
    LOGICAL  :: lland(nbdim), lglac(nbdim)
    REAL(wp) :: zprat, zn1, zn2, zcdnc

    ! Shortcuts to components of echam_cld_config
    !
    REAL(wp), POINTER :: cn1lnd, cn2lnd, cn1sea, cn2sea
    !
    cn1lnd => echam_cld_config(jg)% cn1lnd
    cn2lnd => echam_cld_config(jg)% cn2lnd
    cn1sea => echam_cld_config(jg)% cn1sea
    cn2sea => echam_cld_config(jg)% cn2sea

    DO jc=jcs,jce
      lland(jc) = field%lfland(jc,jb)
      lglac(jc) = lland(jc).AND.field%glac(jc,jb).GT.0._wp
    END DO

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        zprat=(MIN(8._wp,80000._wp/field%presm_old(jc,jk,jb)))**2

        IF (lland(jc).AND.(.NOT.lglac(jc))) THEN
          zn1= cn1lnd
          zn2= cn2lnd
        ELSE
          zn1= cn1sea
          zn2= cn2sea
        END IF
        IF (field%presm_old(jc,jk,jb).LT.80000._wp) THEN
          zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_wp
        END IF
        field% acdnc(jc,jk,jb) = zcdnc
        !
      END DO
    END DO

  END SUBROUTINE calculate_droplet_number
  !----------------------------------------------------------------

 
  !---------------------------------------------------------------------
  SUBROUTINE calculate_cpair_cvair_qconv( jb,jcs,jce, field)
    INTEGER                 ,INTENT(in) :: jb
    INTEGER                 ,INTENT(in) :: jcs, jce
    TYPE(t_echam_phy_field) ,POINTER    :: field

    INTEGER  :: jc, jk

    DO jk = 1,nlev
       DO jc = jcs,jce
          !
          field%cpair(jc,jk,jb) = cpd+(cpv-cpd)*field%qtrc(jc,jk,jb,iqv)
          field%cvair(jc,jk,jb) = cvd+(cvv-cvd)*field%qtrc(jc,jk,jb,iqv)
          !
          field%qconv(jc,jk,jb) = 1._wp/(field%mair(jc,jk,jb)*field%cpair(jc,jk,jb))
          !
       END DO
    END DO
    
  END SUBROUTINE calculate_cpair_cvair_qconv
  !----------------------------------------------------------------

 
  !---------------------------------------------------------------------
  SUBROUTINE message_forcing_action(process, is_in_sd_ed_interval, is_active)
    CHARACTER(LEN=*) ,INTENT(in) :: process
    LOGICAL          ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL          ,INTENT(in) :: is_active

    IF (is_in_sd_ed_interval) THEN
       IF (is_active) THEN
          CALL message('echam_phy_main','compute forcing by '//process)
       ELSE
          CALL message('echam_phy_main','recycle forcing by '//process)
       END IF
    ELSE
       CALL    message('echam_phy_main','no      forcing by '//process)
    END IF

  END SUBROUTINE message_forcing_action
  !---------------------------------------------------------------------



END MODULE mo_echam_phy_main
