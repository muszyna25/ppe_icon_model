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
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_physical_constants, ONLY: p0ref, rd, cvd_o_rd
  USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int
  USE mo_sync,               ONLY: SYNC_E, SYNC_C, sync_patch_array
  USE mo_run_config,         ONLY: ntracer
  USE mo_dynamics_config,    ONLY: lcoriolis
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e

  IMPLICIT NONE

  PRIVATE


  PUBLIC  :: add_slowphys


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
  SUBROUTINE add_slowphys(p_nh, p_patch, p_int, nnow, nnew, dtime, &
                          lstep_adv, nnow_rcf, nnew_rcf)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_int_state), TARGET, INTENT(IN)    :: p_int
    TYPE(t_patch),     TARGET, INTENT(IN)    :: p_patch

    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew, nnow_rcf, nnew_rcf
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime
    LOGICAL,                   INTENT(IN)    :: lstep_adv

    INTEGER :: nlev                  ! number of vertical (full) levels

    INTEGER :: jc, je, jk, jb, jt    ! loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk  ! start and end block
    INTEGER :: i_startidx, i_endidx  ! start and end indices
    INTEGER :: i_nchdom
    REAL(wp) :: p_diag_vt, ddt_corio_vn
    ! Pointers
    INTEGER, POINTER :: iqidx(:,:,:), iqblk(:,:,:) ! to quad edge indices


!-----------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

    i_nchdom = MAX(1,p_patch%n_childdom)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    rl_start = grf_bdywidth_c+1 
    rl_end   = min_rlcell_int 

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jt,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          ! rho (simply copy)
          p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)

          ! w (simply copy)
          p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb)

          ! exner (update)
          p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)  &
             &                            + dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)

          ! diagnose theta_v from updated exner
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = (p0ref/(rd*p_nh%prog(nnew)%rho(jc,jk,jb)))  &
             &                              * p_nh%prog(nnew)%exner(jc,jk,jb)**cvd_o_rd
        ENDDO
      ENDDO

      IF ( lstep_adv ) THEN
        DO jt = 1, ntracer
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              ! tracer (simply copy)
              p_nh%prog(nnew_rcf)%tracer(jc,jk,jb,jt) = p_nh%prog(nnow_rcf)%tracer(jc,jk,jb,jt)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
!$OMP ENDDO NOWAIT


    rl_start = grf_bdywidth_e+1
    rl_end   = min_rledge_int

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk  (rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx, p_diag_vt,ddt_corio_vn)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx


          ! add Coriolis force
          IF ( lcoriolis ) THEN
            ! RBF reconstruction of tangential wind component
            p_diag_vt = &
              p_int%rbf_vec_coeff_e(1,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,1),jk,iqblk(je,jb,1)) + &
              p_int%rbf_vec_coeff_e(2,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,2),jk,iqblk(je,jb,2)) + &
              p_int%rbf_vec_coeff_e(3,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,3),jk,iqblk(je,jb,3)) + &
              p_int%rbf_vec_coeff_e(4,je,jb) * p_nh%prog(nnow)%vn(iqidx(je,jb,4),jk,iqblk(je,jb,4))

            !Coriolis tendency
            ddt_corio_vn  = - p_diag_vt * p_patch%edges%f_e(je,jb)
          ELSE
            ddt_corio_vn = 0._wp
          END IF

          ! vn (update)
          p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)  &
            &                + dtime * (ddt_corio_vn+p_nh%diag%ddt_vn_phy(je,jk,jb))

        ENDDO  ! je
      ENDDO  ! jk
    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL


    ! Synchronize updated prognostic variables
    CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%exner)
    CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn)

  END SUBROUTINE add_slowphys 


END MODULE mo_update_dyn



