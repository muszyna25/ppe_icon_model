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
  USE mtime                  ,ONLY: datetime, isCurrentEventActive, &
       &                            OPERATOR(<=), OPERATOR(>)

  USE mo_model_domain        ,ONLY: t_patch

  USE mo_omp_loop            ,ONLY: omp_loop_cell, &
    &                               omp_loop_cell_tc

  USE mo_echam_phy_config    ,ONLY: echam_phy_tc, dt_zero, echam_phy_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend
  USE mo_echam_phy_diag      ,ONLY: surface_fractions, &
    &                               droplet_number,    &
    &                               cpair_cvair_qconv, &
    &                               initialize,        &
    &                               finalize

  USE mo_run_config          ,ONLY: iqc, iqi

  USE mo_interface_echam_cov ,ONLY: interface_echam_cov
  USE mo_interface_echam_rad ,ONLY: interface_echam_rad
  USE mo_interface_echam_rht ,ONLY: interface_echam_rht
  USE mo_interface_echam_vdf ,ONLY: interface_echam_vdf
  USE mo_interface_echam_car ,ONLY: interface_echam_car
  USE mo_interface_echam_art ,ONLY: interface_echam_art
  USE mo_interface_echam_cnv ,ONLY: interface_echam_cnv
  USE mo_interface_echam_gwd ,ONLY: interface_echam_gwd
  USE mo_interface_echam_sso ,ONLY: interface_echam_sso
  USE mo_interface_echam_cld ,ONLY: interface_echam_cld
  USE mo_interface_echam_mox ,ONLY: interface_echam_mox

  USE mo_run_config          ,ONLY: nlev
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_impl_constants      ,ONLY: min_rlcell_int, grf_bdywidth_c

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: echam_phy_main

CONTAINS

  !>
  !!
  SUBROUTINE echam_phy_main(patch         ,&
    &                       datetime_old  ,&
    &                       pdtime        )


    ! Arguments
    !
    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    TYPE(datetime)         ,POINTER    :: datetime_old
    REAL(wp)               ,INTENT(in) :: pdtime

    ! Pointers
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    INTEGER  :: jg                                         !< grid level/domain index

    LOGICAL  :: is_in_sd_ed_interval                       !< time is in process interval [sd,ed[
    LOGICAL  :: is_active                                  !< process is active

    REAL(wp) :: zq_rlw_impl (nproma)                       !< additional heating by LW rad. due to implicit coupling
                                                           !  in surface energy balance [W/m2]

    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb                                         !< block index
    INTEGER  :: jcs, jce                                   !< start/end column index within this block

    jg = patch%id

    ! associate pointers
    field     => prm_field(jg)
    tend      => prm_tend (jg)
    
    !-------------------------------------------------------------------
    ! Initialize (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,initialize)

    !-------------------------------------------------------------------
    ! Specific heat of moist air (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,cpair_cvair_qconv)

    !-------------------------------------------------------------------
    ! Calculate surface fraction (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,surface_fractions)
 
    !-------------------------------------------------------------------
    ! Cloud cover (diagnostic)
    !-------------------------------------------------------------------
    !
    IF (.NOT.echam_phy_config(jg)%lnew) THEN
    !
    IF ( echam_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       IF   ( (echam_phy_tc(jg)%sd_cld <= datetime_old) .AND. &
            & (echam_phy_tc(jg)%ed_cld >  datetime_old) .AND. &
            & (echam_phy_tc(jg)%dt_cld >  dt_zero     )       ) THEN
          !
          CALL omp_loop_cell(patch,interface_echam_cov)
          !
       END IF
       !
    END IF
    !
    ELSE
    !
    CALL omp_loop_cell(patch,interface_echam_cov)
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
    ! Cloud droplet number concentration (diagnostic)
    ! used in radiation and cloud
    !---------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,droplet_number)


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
       CALL message_forcing_action('LW and SW radiation (rad)'     ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       ! radiative fluxes
       CALL omp_loop_cell_tc(patch, interface_echam_rad      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
       !
       ! always compute radiative heating
       is_active = .TRUE.
       !
       ! radiative heating
       CALL omp_loop_cell_tc(patch, interface_echam_rht      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
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
       CALL message_forcing_action('vertical diffusion (vdf)'      ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       CALL omp_loop_cell_tc(patch, interface_echam_vdf      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
       !
    END IF

    IF (.NOT.echam_phy_config(jg)%lnew) THEN
    !
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    !
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)
    !
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
    !
    END IF

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
       CALL message_forcing_action('lin. Cariolle ozone chem. (car)' ,&
            &                      is_in_sd_ed_interval, is_active   )
       !
       CALL omp_loop_cell_tc(patch, interface_echam_car      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
       !
    END IF


    !-------------------------------------------------------------------
    ! Atmospheric chemistry of ART
    !-------------------------------------------------------------------
    !
    IF (echam_phy_tc(jg)%dt_art > dt_zero) THEN
      !
      is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_art <= datetime_old) .AND. &
           &                          (echam_phy_tc(jg)%ed_art >  datetime_old)
      is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_art,   datetime_old)

      CALL message_forcing_action('ART (rad)'                     ,&
            &                     is_in_sd_ed_interval, is_active )
      !
      ! OMP loops are hidden inside the ART routines. Hence the full patch needs
      ! to be passed to the ART routines and is it not possible to call the
      ! ART reaction interface inside the standard omp block loop.
      ! This should be reprogrammed.
      !
      CALL interface_echam_art(patch                           ,&
           &                   is_in_sd_ed_interval, is_active ,&
           &                   datetime_old, pdtime            )
      !
    END IF
    !


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
       CALL message_forcing_action('Hines gravity wave drag (gwd)' ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       CALL omp_loop_cell_tc(patch, interface_echam_gwd      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
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
       CALL message_forcing_action('sub grid scale orographic effects (sso)' ,&
            &                      is_in_sd_ed_interval, is_active           )
       !
       CALL omp_loop_cell_tc(patch, interface_echam_sso      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
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
       CALL message_forcing_action('cumulus convection (cnv)'      ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       CALL omp_loop_cell_tc(patch, interface_echam_cnv      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
       !
    END IF

    IF (.NOT.echam_phy_config(jg)%lnew) THEN
    ! Update provisional physics state
    field% qtrc(:,:,:,iqc) = field% qtrc(:,:,:,iqc) + tend% qtrc_phy(:,:,:,iqc)*pdtime
    field% qtrc(:,:,:,iqi) = field% qtrc(:,:,:,iqi) + tend% qtrc_phy(:,:,:,iqi)*pdtime
    END IF

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
       CALL omp_loop_cell_tc(patch, interface_echam_cld      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
       !
    END IF


    !-------------------------------------------------------------------
    ! Methane oxidation + H2O photolysis
    !-------------------------------------------------------------------
    !
    IF ( (echam_phy_tc(jg)%dt_mox > dt_zero) ) THEN
      !
      is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_mox <= datetime_old) .AND. &
           &                          (echam_phy_tc(jg)%ed_mox >  datetime_old)
      is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_mox,   datetime_old)
       !
       CALL message_forcing_action('ch4 oxidation & h2o photolysis (mox)' ,&
            &                      is_in_sd_ed_interval, is_active        )
       !
       CALL omp_loop_cell_tc(patch, interface_echam_mox      ,&
            &                is_in_sd_ed_interval, is_active ,&
            &                datetime_old, pdtime            )
    END IF


    !-------------------------------------------------------------------
    ! Finalize (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell(patch,finalize)

    ! disassociate pointers
    NULLIFY(field)
    NULLIFY(tend)

  END SUBROUTINE echam_phy_main
  !---------------------------------------------------------------------

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
