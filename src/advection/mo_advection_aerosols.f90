!>
!! Routines for idealized transport of 2D aerosol fields
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_advection_aerosols

  USE mo_kind,                ONLY: wp, vp
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_extra1
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ntracer, ltimer
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, min_rlcell, nclass_aero
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_advection_config,    ONLY: advection_config
  ! USE mo_advection_limiter,   ONLY: hflx_limiter_sm
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_intp_rbf,            ONLY: rbf_vec_interpol_edge
  USE mo_math_gradients,      ONLY: grad_green_gauss_cell
  USE mo_math_divrot,         ONLY: recon_lsq_cell_l_svd
  USE mo_advection_traj,      ONLY: btraj
  USE mo_exception,           ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: aerosol_2D_advection, setup_aerosol_advection


CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Computes indices for vertically averaged fluxes used for transport of 2D aerosol fields
  !!
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2015-11-05)
  !!
  !!
  SUBROUTINE setup_aerosol_advection(p_patch)

    TYPE(t_patch),     TARGET, INTENT(IN) :: p_patch  !< patch of current domain

    INTEGER :: nlev, jk, jk1, jg, i

    REAL(wp) :: zf, zstart(2), zend(2)

    nlev = p_patch%nlev
    jg   = p_patch%id

    ! Start and end heights for vertical averaging of mass fluxes
    ! 1. all aerosol classes except dust
    zstart(1) =  400._wp
    zend(1)   = 1200._wp
    ! 2. dust (has a larger decay scale than the others)
    zstart(2) = 1000._wp
    zend(2)   = 3000._wp

    ! Compute corresponding level indices
    DO jk = nlev, 1, -1
      jk1 = jk + p_patch%nshift_total
      zf = 0.5_wp*(vct_a(jk1)+vct_a(jk1+1))
      DO i = 1, 2
        IF (zf >= zstart(i) .AND. advection_config(jg)%kend_aero(i) == 0) THEN
          advection_config(jg)%kend_aero(i) = jk
        ENDIF
        IF (zf >= zend(i) .AND. advection_config(jg)%kstart_aero(i) == 0) THEN
          advection_config(jg)%kstart_aero(i) = jk
        ENDIF
      ENDDO
    ENDDO

    WRITE(message_text,'(a,i2)') '2D aerosol advection: start and end levels for flux averaging, domain ',jg
    CALL message('',message_text)
    WRITE(message_text,'(2(a,2i5))') 'kstart: ',advection_config(jg)%kstart_aero(1:2), &
                                     ', kend: ',advection_config(jg)%kend_aero(1:2)
    CALL message('',message_text)

  END SUBROUTINE setup_aerosol_advection

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Driver routine for idealized transport of 2D aerosol fields
  !!
  !!
  !! @par Revision History
  !! Developed by Guenther Zaengl, DWD (2015-11-05)
  !!
  !!
  SUBROUTINE aerosol_2D_advection(p_patch, p_int, dtime, aerosol, vn_traj, mflx_h, mflx_v, &
                                  deltaz_e, rhodz_now, rhodz_new)


    TYPE(t_patch),     TARGET, INTENT(IN) :: p_patch  ! patch of current domain

    TYPE(t_int_state), TARGET, INTENT(IN) :: p_int    ! interpolation state

    REAL(wp), INTENT(IN) :: dtime           ! advection time step

    REAL(wp), INTENT(INOUT) :: aerosol(:,:,:)   ! 2D aerosol optical depth fields (middle index = aerosol class)

    REAL(wp), INTENT(IN)    :: vn_traj(:,:,:)   ! edge-normal velocity for back-trajectory calculation
    REAL(wp), INTENT(IN)    :: mflx_h(:,:,:)    ! horizontal mass flux at edges (vn*rho*deltaz)
    REAL(wp), INTENT(IN)    :: mflx_v(:,:,:)    ! vertical mass flux at cell interfaces (w*rho)
    REAL(vp), INTENT(IN)    :: deltaz_e(:,:,:)  ! layer thickness at edges
    REAL(wp), INTENT(IN)    :: rhodz_now(:,:,:) ! air mass per unit area in cell at time level now
    REAL(wp), INTENT(IN)    :: rhodz_new(:,:,:) ! air mass per unit area in cell at time level new


    ! Vertically averaged edge-normal and tangential velocities
    REAL(wp) :: vn_traj_avg(nproma,2,p_patch%nblks_e), vt_traj_avg(nproma,2,p_patch%nblks_e)

    REAL(wp) :: mflx_h_int(nproma,2,p_patch%nblks_e)  !< vertically integrated horizontal mass flux

    ! Vertically integrated air masses per unit area at both time levels
    REAL(wp) :: rhodz_now_int(nproma,nclass_aero,p_patch%nblks_c), rhodz_new_int(nproma,nclass_aero,p_patch%nblks_c)

    ! Distance vector from upwind cell barycenter to end of back trajectory
    REAL(vp) :: distv_bary(nproma,2,p_patch%nblks_e,2)

    ! Horizontal gradient field of aerosols
    REAL(vp) :: grad_aero(2,nproma,nclass_aero,p_patch%nblks_c)
    !
    REAL(wp) :: flx_aero(nproma,nclass_aero,p_patch%nblks_e)

    ! Upwind cell line and block indices
    INTEGER :: cell_idx(nproma,2,p_patch%nblks_e), cell_blk(nproma,2,p_patch%nblks_e)

    REAL(vp) ::  fluxdiv_c(nproma,nclass_aero), dz(nproma)
    REAL(wp) ::  dthalf

    INTEGER  :: nlev, nlevp1
    INTEGER  :: jb, jk, jt, jc, je, jg, ilc, ibc, kst, kend
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend

    ! Mapping array between 5 aerosol types and 2 layer thickness classes
    INTEGER, PARAMETER :: ji(nclass_aero) = (/1,1,1,1,2/)

    ! Pointer to index fields
    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

   !-----------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    jg  = p_patch%id

    dthalf = 0.5_wp*dtime

    ! line and block indices of edges as seen from cells
    iidx => p_patch%cells%edge_idx
    iblk => p_patch%cells%edge_blk

    ! Compute vertically averaged mass fluxes, back-trajectory velocities, and air masses

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart  = 5
    i_rlend    = min_rledge_int-2
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jk,jt,je,i_startidx,i_endidx,dz)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, i_rlstart, i_rlend )

      DO jt = 1, 2
        vn_traj_avg(:,jt,jb) = 0._wp
        mflx_h_int(:,jt,jb)  = 0._wp
        dz(:) = 0._vp

        DO jk = advection_config(jg)%kstart_aero(jt), advection_config(jg)%kend_aero(jt)
          DO je = i_startidx, i_endidx

            vn_traj_avg(je,jt,jb) = vn_traj_avg(je,jt,jb) + deltaz_e(je,jk,jb)*vn_traj(je,jk,jb)
            mflx_h_int(je,jt,jb)  = mflx_h_int(je,jt,jb)  + mflx_h(je,jk,jb)
            dz(je) = dz(je) + deltaz_e(je,jk,jb)

          ENDDO
        ENDDO
        DO je = i_startidx, i_endidx
          vn_traj_avg(je,jt,jb) = vn_traj_avg(je,jt,jb) / dz(je)
          ! ensure that averaged velocity and averaged mass flux have the same sign
          IF (vn_traj_avg(je,jt,jb)*mflx_h_int(je,jt,jb) < 0._wp) THEN
            vn_traj_avg(je,jt,jb) = SIGN(0.01_wp,mflx_h_int(je,jt,jb))
          ENDIF
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO

    i_rlstart  = grf_bdywidth_c-1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jk,jt,jc,i_startidx,i_endidx,kst,kend)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, i_rlstart, i_rlend )

      ! Decay scales for all aerosol types but dust are the same. Therefore, the computation is done
      ! only for types 1 and 5, and the remaining ones are filled afterwards
      DO jt = 1, nclass_aero, 4
        rhodz_now_int(:,jt,jb) = 0._wp
        rhodz_new_int(:,jt,jb) = 0._wp

        kst  = advection_config(jg)%kstart_aero(ji(jt))
        kend = advection_config(jg)%kend_aero(ji(jt))

        DO jk = kst, kend
          DO jc = i_startidx, i_endidx

            rhodz_now_int(jc,jt,jb) = rhodz_now_int(jc,jt,jb) + rhodz_now(jc,jk,jb)
            rhodz_new_int(jc,jt,jb) = rhodz_new_int(jc,jt,jb) + rhodz_new(jc,jk,jb)

          ENDDO
        ENDDO

        ! Account for vertical fluxes across top and bottom levels
        DO jc = i_startidx, i_endidx
          rhodz_new_int(jc,jt,jb) = rhodz_new_int(jc,jt,jb) - dtime * &
            (mflx_v(jc,kend+1,jb) - mflx_v(jc,kst,jb))
        ENDDO

      ENDDO

      ! Copy integrated values from type 1 to types 2-4
      DO jc = i_startidx, i_endidx
        rhodz_now_int(jc,2:4,jb) = rhodz_now_int(jc,1,jb)
        rhodz_new_int(jc,2:4,jb) = rhodz_new_int(jc,1,jb)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Compute tangential component of vertically averaged back-trajectory winds
    CALL rbf_vec_interpol_edge(vn_traj_avg, p_patch, p_int, vt_traj_avg,               &
                               opt_rlstart=grf_bdywidth_e-2, opt_rlend=min_rledge_int-1)

    ! Compute back trajectories
    CALL btraj(p_patch, p_int, vn_traj_avg, vt_traj_avg, dthalf, cell_idx, cell_blk, &
               distv_bary, opt_rlstart=grf_bdywidth_e-2, opt_rlend=min_rledge_int-1, opt_elev=2 )

    ! Reconstruct 2D gradient fields of aerosol
    IF (advection_config(jg)%igrad_c_miura == 1 .AND. advection_config(jg)%llsq_svd) THEN
      CALL recon_lsq_cell_l_svd(aerosol, p_patch, p_int%lsq_lin, grad_aero, opt_rlend=min_rlcell_int-1)
    ELSE
      CALL grad_green_gauss_cell(aerosol, p_patch, p_int, grad_aero, opt_rlend=min_rlcell_int-1)
    ENDIF


!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart  = grf_bdywidth_e-2
    i_rlend    = min_rledge_int-1
    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jt,je,i_startidx,i_endidx,ilc,ibc)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


      DO jt = 1, nclass_aero
        DO je = i_startidx, i_endidx

          ilc = cell_idx(je,ji(jt),jb)
          ibc = cell_blk(je,ji(jt),jb)
          flx_aero(je,jt,jb) = ( aerosol(ilc,jt,ibc) + distv_bary(je,ji(jt),jb,1)*grad_aero(1,ilc,jt,ibc) + &
            distv_bary(je,ji(jt),jb,2)*grad_aero(2,ilc,jt,ibc) ) * mflx_h_int(je,ji(jt),jb)

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ! Does not seem to be needed because undershoots happen practically never
  !  CALL hflx_limiter_sm(p_patch, p_int, dtime, aerosol, flx_aero,      &
  !                       1, 5, rhodz_now_int, opt_rlend=min_rlcell_int-1)

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP DO PRIVATE(jb,jt,jc,i_startidx,i_endidx,fluxdiv_c)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, i_rlstart, i_rlend )


      DO jc = i_startidx, i_endidx
        DO jt = 1, nclass_aero

          fluxdiv_c(jc,jt) =                                                     &
            flx_aero(iidx(jc,jb,1),jt,iblk(jc,jb,1))*p_int%geofac_div(jc,1,jb) + &
            flx_aero(iidx(jc,jb,2),jt,iblk(jc,jb,2))*p_int%geofac_div(jc,2,jb) + &
            flx_aero(iidx(jc,jb,3),jt,iblk(jc,jb,3))*p_int%geofac_div(jc,3,jb)

        ENDDO
      ENDDO

      DO jt = 1, nclass_aero
        DO jc = i_startidx, i_endidx
          aerosol(jc,jt,jb) = MAX(0._wp, ( aerosol(jc,jt,jb)*rhodz_now_int(jc,jt,jb) - &
            dtime*fluxdiv_c(jc,jt) ) / rhodz_new_int(jc,jt,jb))
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE aerosol_2D_advection


END MODULE mo_advection_aerosols

