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

#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_echam_phy_main

  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message, warning
  USE mtime                  ,ONLY: datetime, isCurrentEventActive, &
       &                            OPERATOR(<=), OPERATOR(>)

  USE mo_model_domain        ,ONLY: t_patch

  USE mo_omp_loop            ,ONLY: omp_loop_cell_diag, &
    &                               omp_loop_cell_prog

  USE mo_echam_phy_config    ,ONLY: echam_phy_tc, dt_zero
  USE mo_echam_phy_diag      ,ONLY: surface_fractions, &
    &                               droplet_number,    &
    &                               cpair_cvair_qconv, &
    &                               initialize,        &
    &                               finalize

  USE mo_echam_diagnostics   ,ONLY: echam_global_diagnostics
#if defined( _OPENACC )
  USE mo_var_list_gpu        ,ONLY: gpu_h2d_var_list, gpu_d2h_var_list
#endif

  USE mo_interface_echam_cov ,ONLY: interface_echam_cov
  USE mo_interface_echam_wmo ,ONLY: interface_echam_wmo
  USE mo_interface_echam_rad ,ONLY: interface_echam_rad
  USE mo_interface_echam_rht ,ONLY: interface_echam_rht
  USE mo_interface_echam_vdf ,ONLY: interface_echam_vdf
  USE mo_interface_echam_car ,ONLY: interface_echam_car
  USE mo_interface_echam_art ,ONLY: interface_echam_art
  USE mo_interface_echam_cnv ,ONLY: interface_echam_cnv
  USE mo_interface_echam_gwd ,ONLY: interface_echam_gwd
  USE mo_interface_echam_sso ,ONLY: interface_echam_sso
  USE mo_interface_echam_cld ,ONLY: interface_echam_cld
  USE mo_interface_echam_mig ,ONLY: interface_echam_mig
  USE mo_interface_echam_mox ,ONLY: interface_echam_mox

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

    ! Local variables
    !
    INTEGER  :: jg                                         !< grid level/domain index

    LOGICAL  :: is_in_sd_ed_interval                       !< time is in process interval [sd,ed[
    LOGICAL  :: is_active                                  !< process is active

    jg = patch%id

    !-------------------------------------------------------------------
    ! Initialize (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,initialize)

    !-------------------------------------------------------------------
    ! Specific heat of moist air (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,cpair_cvair_qconv)

    !-------------------------------------------------------------------
    ! Calculate surface fraction (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,surface_fractions)
 
    !-------------------------------------------------------------------
    ! Cloud cover (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,interface_echam_cov)

    !---------------------------------------------------------------------
    ! 3.9 Determine tropopause height (diagnostic)
    !---------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,interface_echam_wmo)

    !---------------------------------------------------------------------
    ! Cloud droplet number concentration (diagnostic)
    ! used in radiation and cloud
    !---------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,droplet_number)

    !-------------------------------------------------------------------
    ! Radiation (one interface for LW+SW)
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_rad <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_rad >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_rad,   datetime_old)
       !
#if defined( _OPENACC )
       IF ( is_active ) THEN
          CALL warning('GPU:echam_rad_main','GPU host synchronization should be removed when port is done!')
          CALL gpu_d2h_var_list('prm_field_D', jg)
          CALL gpu_d2h_var_list('prm_tend_D', jg)
       END IF
#endif
       !
       CALL message_forcing_action('LW and SW radiation (rad:fluxes )' ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       ! radiative fluxes
       CALL interface_echam_rad(is_in_sd_ed_interval, is_active, &
            &                   patch,                           &
            &                   datetime_old                     )
#if defined( _OPENACC )
       IF ( is_active ) THEN
          CALL warning('GPU:echam_rad_main','GPU device synchronization should be removed when port is done!')
          CALL gpu_h2d_var_list('prm_field_D', jg)
          CALL gpu_h2d_var_list('prm_tend_D', jg)
       END IF
#endif
       !
       ! always compute radiative heating
       is_active = .TRUE.
       !
       CALL message_forcing_action('LW and SW radiation (rht:heating)' ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       ! radiative heating
       CALL omp_loop_cell_prog(patch, interface_echam_rht      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
    END IF


    !-------------------------------------------------------------------
    ! Vertical diffusion, boundary layer and surface
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_vdf <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_vdf >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_vdf,   datetime_old)
       !
       CALL message_forcing_action('vertical diffusion (vdf)'      ,&
            &                      is_in_sd_ed_interval, is_active )
       !
       CALL omp_loop_cell_prog(patch, interface_echam_vdf      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
    END IF

    !-------------------------------------------------------------------
    ! Linearized ozone chemistry of Cariolle
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_car > dt_zero ) THEN
#if defined( _OPENACC )
       CALL warning('GPU:echam_car_main','GPU host synchronization should be removed when port is done!')
       CALL gpu_d2h_var_list('prm_field_D', jg)
       CALL gpu_d2h_var_list('prm_tend_D', jg)
#endif
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_car <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_car >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_car,   datetime_old)
       !
       CALL message_forcing_action('lin. Cariolle ozone chem. (car)' ,&
            &                      is_in_sd_ed_interval, is_active   )
       !
       CALL omp_loop_cell_prog(patch, interface_echam_car      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
#if defined( _OPENACC )
       CALL warning('GPU:echam_car_main','GPU device synchronization should be removed when port is done!')
       CALL gpu_h2d_var_list('prm_field_D', jg)
       CALL gpu_h2d_var_list('prm_tend_D', jg)
#endif
    END IF


    !-------------------------------------------------------------------
    ! Atmospheric chemistry of ART
    !-------------------------------------------------------------------
    !
    IF (echam_phy_tc(jg)%dt_art > dt_zero) THEN
#if defined( _OPENACC )
       CALL warning('GPU:echam_art_main','GPU host synchronization should be removed when port is done!')
       CALL gpu_d2h_var_list('prm_field_D', jg)
       CALL gpu_d2h_var_list('prm_tend_D', jg)
#endif
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
#if defined( _OPENACC )
       CALL warning('GPU:echam_art_main','GPU device synchronization should be removed when port is done!')
       CALL gpu_h2d_var_list('prm_field_D', jg)
       CALL gpu_h2d_var_list('prm_tend_D', jg)
#endif
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
       CALL omp_loop_cell_prog(patch, interface_echam_gwd      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
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
       CALL omp_loop_cell_prog(patch, interface_echam_sso      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
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
       CALL omp_loop_cell_prog(patch, interface_echam_cnv      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
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
       CALL omp_loop_cell_prog(patch, interface_echam_cld      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
    END IF

    !-------------------------------------------------------------------
    ! Graupel (microphysics) processes
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_mig > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_mig <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_mig >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_mig,   datetime_old)
       !
       CALL message_forcing_action('graupel microphysics (mig)',    &
            &                      is_in_sd_ed_interval, is_active)
       !
       CALL omp_loop_cell_prog(patch, interface_echam_mig      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
    END IF

    !-------------------------------------------------------------------
    ! Methane oxidation + H2O photolysis
    !-------------------------------------------------------------------
    !
    IF ( echam_phy_tc(jg)%dt_mox > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (echam_phy_tc(jg)%sd_mox <= datetime_old) .AND. &
            &                          (echam_phy_tc(jg)%ed_mox >  datetime_old)
       is_active = isCurrentEventActive(echam_phy_tc(jg)%ev_mox,   datetime_old)
       !
       CALL message_forcing_action('ch4 oxidation & h2o photolysis (mox)' ,&
            &                      is_in_sd_ed_interval, is_active        )
       !
       CALL omp_loop_cell_prog(patch, interface_echam_mox      ,&
            &                  is_in_sd_ed_interval, is_active ,&
            &                  datetime_old, pdtime            )
       !
    END IF


    !-------------------------------------------------------------------
    ! Global output diagnostics
    !-------------------------------------------------------------------
    !
    CALL echam_global_diagnostics(patch)


    !-------------------------------------------------------------------
    ! Finalize (diagnostic)
    !-------------------------------------------------------------------
    !
    CALL omp_loop_cell_diag(patch,finalize)

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
