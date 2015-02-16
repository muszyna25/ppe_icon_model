!! mo_les_utilities
!!
!! This module has usefull routines for LES runs
!!
!! @author Anurag Dipankar, MPI-M
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-June-01)
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

MODULE mo_les_utilities

  USE mo_kind,                ONLY: wp
  USE mo_nonhydro_types,      ONLY: t_nh_metrics
  USE mo_nonhydro_state,      ONLY: p_nh_state
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp,                ONLY: cells2edges_scalar, cells2verts_scalar, &
                                    edges2cells_scalar, edges2verts_scalar
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_sync,                ONLY: global_sum_array, SYNC_C, SYNC_E, SYNC_V, &
                                    sync_patch_array, sync_patch_array_mult
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell_int
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, rcpd, alv, alvdcp, &
                                    rd_o_cpd, rd
  USE mo_satad,               ONLY: qsat_rho
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_grid_config,         ONLY: n_dom
  USE mo_les_config,          ONLY: les_config
  IMPLICIT NONE

  PRIVATE

  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
    ddxt_z_full, &
    ddxt_z_full_c, &
    ddxt_z_full_v,&
    ddxt_z_half_e,&
    ddxt_z_half_v,&
    ddxt_z_half_c,&
    ddxn_z_full,&
    ddxn_z_full_c,&
    ddxn_z_full_v,&
    ddxn_z_half_e,&
    ddxn_z_half_c

  PUBLIC :: vert_intp_full2half_cell_3d, vert_intp_linear_1d, global_hor_mean
  PUBLIC :: vertical_derivative, brunt_vaisala_freq
  PUBLIC :: &
    ddxt_z_full, &
    ddxt_z_full_c, &
    ddxt_z_full_v,&
    ddxt_z_half_e,&
    ddxt_z_half_v,&
    ddxt_z_half_c,&
    ddxn_z_full,&
    ddxn_z_full_c,&
    ddxn_z_full_v,&
    ddxn_z_half_e,&
    ddxn_z_half_c

  PUBLIC :: init_vertical_grid_for_les

  CONTAINS


  !>
  !! init_vertical_grid_for_les
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Slavko Brdar, DWD (2014-08-29)
  SUBROUTINE init_vertical_grid_for_les(jg, p_patch, p_int, p_metrics)
    INTEGER,                   INTENT(in)     :: jg
    TYPE(t_patch),     TARGET, INTENT(in)     :: p_patch
    TYPE(t_int_state),         INTENT(in)     :: p_int
    TYPE(t_nh_metrics),        INTENT(inout)  :: p_metrics

    ! local variables
    INTEGER :: i_startblk, i_endblk
    INTEGER :: nlevp1, i_nchdom, nlev, nlen
    INTEGER :: jk, je, jc, jb, jkm1

    nlev        = p_patch%nlev
    nlevp1      = p_patch%nlevp1

    IF(.NOT.les_config(jg)%les_metric) &
      RETURN

    ALLOCATE(&
    ddxt_z_full(nproma,p_patch%nlev,p_patch%nblks_e), &
    ddxt_z_full_c(nproma,p_patch%nlev,p_patch%nblks_c), &
    ddxt_z_full_v(nproma,p_patch%nlev,p_patch%nblks_v),&
    ddxt_z_half_e(nproma,p_patch%nlevp1,p_patch%nblks_e),&
    ddxt_z_half_v(nproma,p_patch%nlevp1,p_patch%nblks_v),&
    ddxt_z_half_c(nproma,p_patch%nlevp1,p_patch%nblks_c),&
    ddxn_z_full(nproma,p_patch%nlev,p_patch%nblks_e),&
    ddxn_z_full_c(nproma,p_patch%nlev,p_patch%nblks_c),&
    ddxn_z_full_v(nproma,p_patch%nlev,p_patch%nblks_v),&
    ddxn_z_half_e(nproma,p_patch%nlevp1,p_patch%nblks_e),&
    ddxn_z_half_c(nproma,p_patch%nlevp1,p_patch%nblks_c))

    ddxt_z_full(:,:,:)   = p_metrics%ddxt_z_full(:,:,:)
    ddxt_z_half_e(:,:,:) = p_metrics%ddxt_z_half_e(:,:,:)
    ddxn_z_full(:,:,:)   = p_metrics%ddxn_z_full(:,:,:)
    ddxn_z_half_e(:,:,:) = p_metrics%ddxn_z_half_e(:,:,:)

    CALL sync_patch_array(SYNC_E, p_patch, ddxn_z_full)
    CALL sync_patch_array(SYNC_E, p_patch, ddxn_z_half_e)
    CALL sync_patch_array(SYNC_E, p_patch, ddxt_z_full)
    CALL sync_patch_array(SYNC_E, p_patch, ddxt_z_half_e)

    IF (p_test_run) THEN
      ddxt_z_half_v(:,:,:) = 0._wp
      ddxn_z_half_c(:,:,:) = 0._wp
      ddxn_z_full_c(:,:,:) = 0._wp
      ddxn_z_full_v(:,:,:) = 0._wp
      ddxt_z_half_c(:,:,:) = 0._wp
      ddxt_z_full_c(:,:,:) = 0._wp
      ddxt_z_full_v(:,:,:) = 0._wp
      p_metrics%inv_ddqz_z_full_v(:,:,:) = 0._wp
!      p_metrics%inv_ddqz_z_full_e(:,:,:) = 0._wp
    END IF

    !
    ! normal metric derivative
    !
    CALL edges2cells_scalar(ddxn_z_half_e, p_patch, p_int%e_bln_c_s, &
      ddxn_z_half_c)
    CALL sync_patch_array(SYNC_C, p_patch, ddxn_z_half_c)

    CALL edges2cells_scalar(ddxn_z_full, p_patch, p_int%e_bln_c_s, &
      ddxn_z_full_c)
    CALL sync_patch_array(SYNC_C, p_patch, ddxn_z_full_c)

    CALL cells2verts_scalar(ddxn_z_full_c, p_patch, p_int%cells_aw_verts, &
      ddxn_z_full_v)
    CALL sync_patch_array(SYNC_V, p_patch, ddxn_z_full_v)

    !
    ! tangential metric derivative
    !
    CALL edges2cells_scalar(ddxt_z_half_e, p_patch, p_int%e_bln_c_s, &
      ddxt_z_half_c)
    CALL sync_patch_array(SYNC_C, p_patch, ddxt_z_half_c)

    CALL cells2verts_scalar(ddxt_z_half_c, p_patch, p_int%cells_aw_verts, &
      ddxt_z_half_v)
    CALL sync_patch_array(SYNC_V, p_patch, ddxt_z_half_v)

    CALL edges2cells_scalar(ddxt_z_full, p_patch, p_int%e_bln_c_s, &
      ddxt_z_full_c)
    CALL sync_patch_array(SYNC_C, p_patch, ddxt_z_full_c)

    CALL cells2verts_scalar(ddxt_z_full_c, p_patch, p_int%cells_aw_verts, &
      ddxt_z_full_v)
    CALL sync_patch_array(SYNC_V, p_patch, ddxt_z_full_v)

    CALL cells2verts_scalar(p_metrics%inv_ddqz_z_full, p_patch, p_int%cells_aw_verts, &
      p_metrics%inv_ddqz_z_full_v)
    CALL sync_patch_array(SYNC_V, p_patch, p_metrics%inv_ddqz_z_full_v)
    CALL sync_patch_array(SYNC_E, p_patch, p_metrics%inv_ddqz_z_full_e)

    ! SB for testing
    IF (.false.) THEN
      ddxt_z_full(:,:,:)   = 0._wp
      ddxt_z_full_c(:,:,:) = 0._wp
      ddxt_z_full_v(:,:,:) = 0._wp
      ddxt_z_half_e(:,:,:) = 0._wp
      ddxt_z_half_v(:,:,:) = 0._wp
      ddxt_z_half_c(:,:,:) = 0._wp

      ddxn_z_full(:,:,:)   = 0._wp
      ddxn_z_full_c(:,:,:) = 0._wp
      ddxn_z_full_v(:,:,:) = 0._wp
      ddxn_z_half_e(:,:,:) = 0._wp
      ddxn_z_half_c(:,:,:) = 0._wp
    ELSE IF (.false.) THEN
      ddxt_z_full   = p_metrics%ddxt_z_full
      ddxt_z_full_c = p_metrics%ddxt_z_full_c
      ddxt_z_full_v = p_metrics%ddxt_z_full_v
      ddxt_z_half_e = p_metrics%ddxt_z_half_e
      ddxt_z_half_v = p_metrics%ddxt_z_half_v
      ddxt_z_half_c = p_metrics%ddxt_z_half_c

      ddxn_z_full   = p_metrics%ddxn_z_full
      ddxn_z_full_c = p_metrics%ddxn_z_full_c
      ddxn_z_full_v = p_metrics%ddxn_z_full_v
      ddxn_z_half_e = p_metrics%ddxn_z_half_e
      ddxn_z_half_c = p_metrics%ddxn_z_half_c
    END IF

  END SUBROUTINE init_vertical_grid_for_les


  !>
  !! vert_intp_full2half_3d
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  SUBROUTINE vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end)

    TYPE(t_nh_metrics),INTENT(in), TARGET :: p_metrics 
    TYPE(t_patch),     INTENT(in), TARGET :: p_patch
    REAL(wp), INTENT(in)                  :: varin(:,:,:)
    INTEGER,  INTENT(in)                  :: rl_start, rl_end 
    REAL(wp), INTENT(out)                 :: varout(:,:,:)                     

    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_endidx, i_startidx, nlevp1, i_nchdom, nlev
    INTEGER :: jk, jc, jb

    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlev+1
    i_nchdom  = MAX(1,p_patch%n_childdom)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx , i_endidx
         DO jk = 2 , nlev  
#else
        DO jk = 2 , nlev  
         DO jc = i_startidx , i_endidx
#endif
          varout(jc,jk,jb) = p_metrics%wgtfac_c(jc,jk,jb)*varin(jc,jk,jb) + &
                        (1._wp-p_metrics%wgtfac_c(jc,jk,jb))*varin(jc,jk-1,jb)
         END DO
        END DO
        DO jc = i_startidx, i_endidx
           varout(jc,1,jb) =                                &
             p_metrics%wgtfacq1_c(jc,1,jb)*varin(jc,1,jb) + &
             p_metrics%wgtfacq1_c(jc,2,jb)*varin(jc,2,jb) + &
             p_metrics%wgtfacq1_c(jc,3,jb)*varin(jc,3,jb)

           varout(jc,nlevp1,jb) =                               &
             p_metrics%wgtfacq_c(jc,1,jb)*varin(jc,nlev,jb)   + &
             p_metrics%wgtfacq_c(jc,2,jb)*varin(jc,nlev-1,jb) + &
             p_metrics%wgtfacq_c(jc,3,jb)*varin(jc,nlev-2,jb)
        END DO     
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL 


  END SUBROUTINE vert_intp_full2half_cell_3d


  !> vert_intp_linear_1d
  !!
  !! linear vertical interpolation: grid za to zb
  !! - It extrapolates if no data given for zb > za
  !! Taken from UCLA-LES
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  SUBROUTINE vert_intp_linear_1d(za, xa, zb, xb) 
     REAL(wp), INTENT(IN)  :: za(:), zb(:), xa(:)
     REAL(wp), INTENT(OUT) :: xb(:)
  
     REAL(wp) :: wt
     INTEGER  :: l, k, na, nb

     na = SIZE(za)
     nb = SIZE(zb)

     l = na
     DO k = nb, 1, -1
       IF (zb(k) <= za(1)) THEN
          DO WHILE ( zb(k) > za(l-1) .AND. l > 1)
             l = l-1
          END DO
          wt=(zb(k)-za(l))/(za(l-1)-za(l))
          xb(k)=xa(l)+(xa(l-1)-xa(l))*wt    
       ELSE
          wt=(zb(k)-za(1))/(za(2)-za(1))
          xb(k)=xa(1)+(xa(2)-xa(1))*wt
       END IF
    END DO

  END SUBROUTINE vert_intp_linear_1d

  !>
  !! global_hor_mean: only called for interior points
  !! Calculates horizontally averaged vertically varying quantaties 
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  SUBROUTINE global_hor_mean(p_patch, var, varout, inv_no_cells, nchdom)

    TYPE(t_patch),     INTENT(in), TARGET :: p_patch
    REAL(wp), INTENT(in)                  :: var(:,:,:), inv_no_cells
    INTEGER,  INTENT(in)                  :: nchdom
    REAL(wp), INTENT(out)                 :: varout(:)                     

    REAL(wp) :: var_aux(SIZE(var,1),SIZE(var,2),SIZE(var,3))
    INTEGER  :: i_startblk, i_endblk, rl_start
    INTEGER  :: i_endidx, i_startidx
    INTEGER  :: jk, jc, jb, nz

    !Put all fields to 0
    var_aux(:,:,:) = 0._wp

    rl_start   = grf_bdywidth_c+1
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell_int,nchdom)
    nz         = SIZE(var,2)

   !Now put values in interior nodes
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, min_rlcell_int)
#ifdef __LOOP_EXCHANGE
       DO jc = i_startidx , i_endidx
         DO jk = 1 , nz
#else
       DO jk = 1 , nz
         DO jc = i_startidx , i_endidx
#endif
             var_aux(jc,jk,jb) = var(jc,jk,jb)
         END DO
       END DO
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   DO jk = 1 , nz
    varout(jk) = global_sum_array(var_aux(:,jk,:)) * inv_no_cells
   END DO

  END SUBROUTINE global_hor_mean

  !>
  !! vertical_derivative
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  FUNCTION vertical_derivative (var, inv_dz) RESULT(dvardz)

    REAL(wp), INTENT(in) :: var(:), inv_dz(:)
                     
    REAL(wp) :: dvardz(SIZE(inv_dz))                     
    INTEGER  :: jk

!$OMP PARALLEL
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jk = 1 , SIZE(inv_dz)
      dvardz(jk) = ( var(jk) - var(jk+1) ) * inv_dz(jk)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END FUNCTION vertical_derivative

  !>
  !! Brunt Vaisala Frequency: 
  !! Calculates BVF for unsaturated and saturated case based on Durran & Klemp 1982
  !! Eq. 4. and using moist lapse rate expression from Marshall and Plumb
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2014-July-07)
  SUBROUTINE brunt_vaisala_freq(p_patch, p_metrics, thetav, bru_vais)

    TYPE(t_patch),     INTENT(in), TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in), TARGET :: p_metrics 
    REAL(wp), DIMENSION(:,:,:), INTENT(in):: thetav
    REAL(wp), INTENT(OUT)                 :: bru_vais(:,:,:)

    REAL(wp) :: thetav_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: term1, qs, temp_ic
    INTEGER  :: i_startblk, i_endblk, rl_start, rl_end
    INTEGER  :: i_endidx, i_startidx, nlev, nlevp1, i_nchdom
    INTEGER  :: jk, jc, jb

    !To be calculated at all cells at interface levels, except top/bottom 
    !boundaries

    nlev      = p_patch%nlev
    nlevp1    = nlev+1
    i_nchdom  = MAX(1,p_patch%n_childdom)

    rl_start   = 2
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, thetav, thetav_ic, rl_start, rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx , i_endidx
        DO jk = 2 , nlev
#else
      DO jk = 2 , nlev
        DO jc = i_startidx , i_endidx
#endif
          bru_vais(jc,jk,jb) = grav * ( thetav(jc,jk-1,jb) - thetav(jc,jk,jb) ) * &
                               p_metrics%inv_ddqz_z_half(jc,jk,jb)/thetav_ic(jc,jk,jb)    
        END DO
      END DO     
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL     
   
  END SUBROUTINE brunt_vaisala_freq


   
!-------------------------------------------------------------------------------

END MODULE mo_les_utilities


