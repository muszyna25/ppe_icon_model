!>
!! Updates dynamical fields with slow physics tendencies in the special case
!! that dynamics are switched off (ldynamics=F).
!!
!! Anurag Dipankar, MPIM (24 March 2014): Fixed some bugs and removed unnecessary 
!! sync for tracers as tracers are synced within the nwp_interface. 
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2013-11-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.


MODULE mo_update_dyn

  USE mo_kind,               ONLY: wp
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_state
  USE mo_physical_constants, ONLY: p0ref, rd, cvd_o_rd
  USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int
  USE mo_sync,               ONLY: SYNC_E, SYNC_C, sync_patch_array
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
#ifdef _OPENACC
  USE mo_mpi,                ONLY: i_am_accel_node
#endif

  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: add_slowphys

#if defined( _OPENACC )
#define ACC_DEBUG NOACC
#if defined(__NH_SUPERVISE_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.     !  THIS SHOULD BE .FALSE. AFTER VALIDATION PHASE!
#endif

CONTAINS

!-------------------------------------------------------------------------
!
!
!>
!! Updates dynamical fields with slow physics tendencies
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2013-11-28)
!! 
!!
  SUBROUTINE add_slowphys(p_nh, p_patch, nnow, nnew, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_patch),     TARGET, INTENT(INOUT) :: p_patch

    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    INTEGER :: nlev                  ! number of vertical (full) levels

    INTEGER :: jc, je, jk, jb        ! loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk  ! start and end block
    INTEGER :: i_startidx, i_endidx  ! start and end indices
    INTEGER :: i_nchdom

!-----------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev

    i_nchdom = MAX(1,p_patch%n_childdom)

!$ACC DATA PRESENT( p_nh )  IF ( i_am_accel_node .AND. acc_on )

!$ACC UPDATE DEVICE ( p_nh%prog(nnow)%exner, p_nh%prog(nnow)%vn,                           &
!$ACC                 p_nh%diag%ddt_exner_phy, p_nh%diag%ddt_vn_phy, p_nh%prog(nnew)%rho ) &
!$ACC        IF ( i_am_accel_node .AND. acc_on .AND. acc_validate )

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    rl_start = grf_bdywidth_c+1 
    rl_end   = min_rlcell_int 

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! exner (update)
          p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)  &
             &                            + dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)

          ! diagnose theta_v from updated exner
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = (p0ref/(rd*p_nh%prog(nnew)%rho(jc,jk,jb)))  &
             &                              * p_nh%prog(nnew)%exner(jc,jk,jb)**cvd_o_rd

        ENDDO
      ENDDO
!$ACC END PARALLEL

    ENDDO
!$OMP ENDDO NOWAIT


    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk  (rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
!$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          ! vn (update)
          p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)  &
            &                + dtime * p_nh%diag%ddt_vn_phy(je,jk,jb)

        ENDDO  ! je
      ENDDO  ! jk
!$ACC END PARALLEL

    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL


    ! Synchronize updated prognostic variables
    CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%exner,opt_varname="prog(nnew)%exner")
    CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn,opt_varname="prog(nnew)%vn")

!$ACC UPDATE HOST( p_nh%prog(nnew)%exner, p_nh%prog(nnew)%vn, p_nh%prog(nnew)%theta_v ) &
!$ACC        IF ( i_am_accel_node .AND. acc_on .AND. acc_validate )

!$ACC END DATA

  END SUBROUTINE add_slowphys 


END MODULE mo_update_dyn



