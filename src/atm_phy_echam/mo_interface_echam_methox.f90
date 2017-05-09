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

!----------------------------
MODULE mo_interface_echam_methox
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: nlev, iqv
  USE mo_methox,              ONLY: methox
CONTAINS
  SUBROUTINE interface_echam_methox(patch,rl_start,rl_end,field,tend)
    TYPE(t_patch)   ,INTENT(IN), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)         :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER  :: field    
    TYPE(t_echam_phy_tend) ,   POINTER  :: tend
    REAL(wp)                            :: dqdt(nproma,nlev)
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(patch, jb,   i_startblk,   i_endblk,    &
            &             jcs,   jce,  rl_start,     rl_end       )
       dqdt(:,:)=0._wp
       CALL methox(jcs, jce, nproma, nlev, field%qtrc(:,:,jb,iqv), dqdt, field%presm_old(:,:,jb))
      tend% qtrc(jcs:jce,:,jb,iqv) = tend% qtrc(jcs:jce,:,jb,iqv) + dqdt(jcs:jce,:)
    END DO
  END SUBROUTINE INTERFACE_ECHAM_METHOX
  END MODULE mo_interface_echam_methox
