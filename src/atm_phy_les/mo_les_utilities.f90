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
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_sync,                ONLY: global_sum_array
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_parallel_config,     ONLY: nproma
  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell_int
  USE mo_physical_constants,  ONLY: cpd, rcvd, p0ref, grav, rcpd, alv, alvdcp, &
                                    rd_o_cpd, rd
  USE mo_satad,               ONLY: qsat_rho
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vert_intp_full2half_cell_3d, vert_intp_linear_1d, global_hor_mean
  PUBLIC :: vertical_derivative, brunt_vaisala_freq

  CONTAINS


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

        DO jk = 2 , nlev  
         DO jc = i_startidx , i_endidx
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
       DO jk = 1 , nz
         DO jc = i_startidx , i_endidx
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
  SUBROUTINE brunt_vaisala_freq(p_patch, p_metrics, t_g, pres_sfc, theta, qc, rho, &
                                temp, bru_vais)

    TYPE(t_patch),     INTENT(in), TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in), TARGET :: p_metrics 
    REAL(wp), DIMENSION(:,:,:), INTENT(in):: theta, qc, rho, temp
    REAL(wp), DIMENSION(:,:),   INTENT(in):: t_g, pres_sfc
    REAL(wp), INTENT(OUT)                 :: bru_vais(:,:,:)

    REAL(wp) :: theta_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: term1, qs
    INTEGER  :: i_startblk, i_endblk, rl_start, rl_end
    INTEGER  :: i_endidx, i_startidx, nlev, nlevp1, i_nchdom
    INTEGER  :: jk, jc, jb

    !To be calculated at all cells at full levels

    nlev      = p_patch%nlev
    nlevp1    = nlev+1
    i_nchdom  = MAX(1,p_patch%n_childdom)

    rl_start   = 1
    rl_end     = min_rlcell_int-2
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, theta, theta_ic, rl_start, rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, qs, term1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx , i_endidx
        theta_ic(jc,nlevp1,jb) = t_g(jc,jb) * EXP(rd_o_cpd*LOG(pres_sfc(jc,jb)/p0ref))
      
        !bru_vais(jc,nlev,jb) = grav * ( 0.375_wp*theta(jc,nlev,jb) + 0.75_wp*theta(jc,nlev-1,jb)- &
        !                                0.125_wp*theta(jc,nlev-2,jb) - theta_ic(jc,nlevp1,jb) ) * &
        !                                p_metrics%inv_ddqz_z_full(jc,jk,jb)/theta(jc,nlev,jb)    
      END DO

      DO jk = 1 , nlev
        DO jc = i_startidx , i_endidx
 
          bru_vais(jc,jk,jb) = grav * ( theta_ic(jc,jk,jb) - theta_ic(jc,jk+1,jb) ) * &
                               p_metrics%inv_ddqz_z_full(jc,jk,jb)/theta(jc,jk,jb)    

          !Now for saturated case
          IF(qc(jc,jk,jb) > 0._wp)THEN 

            qs    = qsat_rho(temp(jc,jk,jb),rho(jc,jk,jb))  
            term1 = 1._wp + alv*qs/(rd*temp(jc,jk,jb))

            bru_vais(jc,jk,jb) = ( bru_vais(jc,jk,jb) + grav**2 * rcpd * ( term1 / &
                    (1._wp + 0.067_wp*alvdcp*qs) - 1._wp )/temp(jc,jk,jb) ) * term1 

          END IF

        END DO
      END DO     
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL     
   
  END SUBROUTINE brunt_vaisala_freq


   
!-------------------------------------------------------------------------------

END MODULE mo_les_utilities


