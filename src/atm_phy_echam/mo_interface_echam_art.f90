!>
!! @brief Subroutine interface_echam_art calls the art reaction interface.
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

MODULE mo_interface_echam_art

  USE mo_kind                   ,ONLY: wp
  USE mtime                     ,ONLY: datetime

  USE mo_dynamics_config        ,ONLY: nnew_rcf
  USE mo_model_domain           ,ONLY: t_patch
  USE mo_nonhydro_state         ,ONLY: p_nh_state_lists
  
!!$  USE mo_echam_phy_config       ,ONLY: echam_phy_config
  USE mo_echam_phy_memory       ,ONLY: t_echam_phy_field, prm_field, &
    &                                  t_echam_phy_tend,  prm_tend
#ifdef __ICON_ART
  USE mo_art_reaction_interface ,ONLY: art_reaction_interface
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_art

CONTAINS

  SUBROUTINE interface_echam_art(patch                ,&
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         pdtime               )

    ! Arguments
    !
    TYPE(t_patch)   ,TARGET ,INTENT(in) :: patch
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Pointers
    !
!!$    LOGICAL                 ,POINTER    :: lparamcpl
!!$    INTEGER                 ,POINTER    :: fc_art
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    INTEGER  :: jg

    jg  = patch%id

    ! associate pointers
!!$    lparamcpl => echam_phy_config(jg)%lparamcpl
!!$    fc_art    => echam_phy_config(jg)%fc_art
    field     => prm_field(jg)
    tend      => prm_tend (jg)  


    IF ( is_in_sd_ed_interval ) THEN
       !
#ifdef __ICON_ART
       IF ( is_active ) THEN
          !
          CALL art_reaction_interface(jg,                                           & !> in
               &                      datetime_old,                                 & !> in
               &                      pdtime,                                       & !> in
               &                      p_nh_state_lists(jg)%prog_list(nnew_rcf(jg)), & !> in
               &                      field%qtrc)
          !
       END IF
#endif
       !
!!$       ! accumulate tendencies for later updating the model state
!!$       SELECT CASE(fc_art)
!!$       CASE(0)
!!$          ! diagnostic, do not use tendency
!!$       CASE(1)
!!$          ! use tendencies to update the model state
!!$          ...
!!$       END SELECT
!!$       !
!!$       ! update physics state for input to the next physics process
!!$       SELECT CASE(fc_art)
!!$       CASE(0)
!!$          ! diagnostic, do not use tendency
!!$       CASE(1)
!!$          ! use tendency to update the physics state
!!$          IF (lparamcpl) THEN
!!$             ...
!!$          END IF
!!$       END SELECT
!!$       !
!!$    ELSE
!!$       !
!!$       ! reset output array of ART to default values valid for times
!!$       ! before and after the interval [start date,end date[.
       !
    END IF

    ! disassociate pointers
!!$    NULLIFY(lparamcpl)
!!$    NULLIFY(fc_art)
    NULLIFY(field)
    NULLIFY(tend)

  END SUBROUTINE interface_echam_art

END MODULE mo_interface_echam_art
