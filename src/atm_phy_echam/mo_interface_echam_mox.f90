!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author S. Rast, MPI-M
!!
!! @par Revision History
!!   original source S. Rast, MPI-M (03-24-2017)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_mox
  
  USE mo_kind                ,ONLY: wp

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_run_config          ,ONLY: nlev, iqv

  USE mtime                  ,ONLY: datetime
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
    &                               t_echam_phy_tend,  prm_tend
  
  USE mo_methox              ,ONLY: methox

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_mox

CONTAINS
  
  SUBROUTINE echam_mox(is_in_sd_ed_interval, &  
       &               is_active,            &
       &               jg, jb,jcs,jce,       &
       &               datetime_old,         &
       &               pdtime                )

    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    INTEGER                 ,INTENT(in) :: jg                  !< grid  index
    INTEGER                 ,INTENT(in) :: jb                  !< block index
    INTEGER                 ,INTENT(in) :: jcs, jce            !< start/end column index within this block
    TYPE(datetime)          ,POINTER    :: datetime_old        !< generic input, not used in echam_mox
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Local variables
    !
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! associate pointers
    field => prm_field(jg)
    tend  => prm_tend (jg)

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN
          !
          CALL methox( jcs, jce,                  &
               &       nproma,                    &
               &       nlev,                      &
               &       field%presm_old(:,:,jb),   &
               &       field%qtrc(:,:,jb,iqv),    &
               &       tend%qtrc_mox(:,:,jb,iqv)  )
          !
       END IF
       !
       tend% qtrc_phy(jcs:jce,:,jb,iqv) = tend% qtrc_phy(jcs:jce,:,jb,iqv) + tend%qtrc_mox(jcs:jce,:,jb,iqv)
       !
    ELSE
       !
       tend% qtrc_mox(jcs:jce,:,jb,iqv) = 0.0_wp
       !
    END IF
    
  END SUBROUTINE echam_mox

END MODULE mo_interface_echam_mox
