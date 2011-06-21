!>
!! mo_vector_operations
!!
!! Basket for subroutines that perform transformation between orthogonal,
!! covariant and contravariant vectors. Also vorticities are computed here.
!!
!! @author Almut Gassmann, MPI-M
!!
!! @par Revision History
!! Initial release by Almut Gassmann (2009-04-17)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_vector_operations

  USE mo_kind,             ONLY: wp
  USE mo_io_nml,           ONLY: l_outputtime
  USE mo_run_nml,          ONLY: nproma, ltimer, dtime
  USE mo_vertical_grid,    ONLY: nflat
  USE mo_model_domain,     ONLY: t_patch
  USE mo_interpolation,    ONLY: t_int_state,          edges2edges_scalar, &
                                 edges2cells_scalar, cells2edges_scalar, &
                                 edges2verts_scalar, verts2edges_scalar, &
                                 verts2cells_scalar, cells2verts_scalar, &
                                 i_cori_method, sick_a, sick_o, l_corner_vort
  USE mo_nonhydro_state,   ONLY: t_nh_diag, t_nh_metrics
  USE mo_math_operators,   ONLY: rot_vertex, grad_fd_norm
  USE mo_sync,             ONLY: SYNC_E, SYNC_V, SYNC_C, &
                                 sync_patch_array, sync_patch_array_mult
  USE mo_exception,        ONLY: finish
  USE mo_timer,            ONLY: timer_start, timer_stop, timer_corio
  USE mo_parallel_nml,     ONLY: p_test_run


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: covariant_velocities, orthogonal_vorticities, &
            contravariant_vorticities, vorticity_tendencies, &
            kinetic_energy, contravariant_vert_mass_flux, &
            impl_vert_adv_vn

  CONTAINS

  !----------------------------------------------------------------------------
  !>
  !! impl_vert_adv_vn 
  !!
  !! Implicit vertical advection for horizontal velocity
  !! 3 diagonal matrix is solved. (Durran,1999:pp.440-441) 
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2011-01-11)
  !!
  SUBROUTINE impl_vert_adv_vn (p_vn,p_w_e,p_patch,p_metrics,p_vi)

    ! Passed variables
    REAL(wp), INTENT(in)              :: p_vn(:,:,:)   !< horizontal velocity
    REAL(wp), INTENT(in)              :: p_w_e(:,:,:)  !< vertical contravariant veloc
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch       !< patch
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics     !< metrics
    REAL(wp), INTENT(inout)           :: p_vi(:,:,:)  &
    &  ;  !< implicit horizontal velocity: (vn(nnow)+vn(nnew))/2

    ! Local variables
    REAL(wp):: z_a(nproma,p_patch%nlev), &
      &        z_b(nproma,p_patch%nlev), &
      &        z_c(nproma,p_patch%nlev)  !< matrix coefficients
    REAL(wp):: z_fac(nproma), z_q(nproma,p_patch%nlev)      !< auxiliary vectors
    INTEGER :: nblks_e, npromz_e, jb, jk, nlen
    INTEGER :: nlev, nlevp1              !< number of full and half levels
    !--------------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, z_a, z_b, z_c, z_q, z_fac)
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1, nlev
        p_vi(1:nlen,jk,jb) = p_vn(1:nlen,jk,jb)
        z_fac(1:nlen)  = 0.25_wp*dtime/p_metrics%ddqz_z_full_e(1:nlen,jk,jb)
        z_b(1:nlen,jk) = 1.0_wp-z_fac(1:nlen)*(p_w_e(1:nlen,jk,jb)-p_w_e(1:nlen,jk+1,jb))
        z_a(1:nlen,jk) =  z_fac(1:nlen)*p_w_e(1:nlen,jk  ,jb)
        z_c(1:nlen,jk) = -z_fac(1:nlen)*p_w_e(1:nlen,jk+1,jb)
      ENDDO
      z_q (1:nlen,1   ) = -z_c(1:nlen,1)/z_b(1:nlen,1)
      p_vi(1:nlen,1,jb) = p_vi(1:nlen,1,jb)/z_b(1:nlen,1)
      DO jk = 2, nlev
        z_fac(1:nlen)      = 1.0_wp/(z_b(1:nlen,jk) + z_a(1:nlen,jk)*z_q(1:nlen,jk-1))
        z_q(1:nlen,jk)     = -z_c(1:nlen,jk)*z_fac(1:nlen)
        p_vi(1:nlen,jk,jb) = (p_vi(1:nlen,jk,jb)  &
        &                    -z_a(1:nlen,jk)*p_vi(1:nlen,jk-1,jb))*z_fac(1:nlen)
      ENDDO
      DO jk=nlev-1,1,-1
        p_vi(1:nlen,jk,jb) = p_vi(1:nlen,jk,jb)+z_q(1:nlen,jk)*p_vi(1:nlen,jk+1,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE impl_vert_adv_vn

  !----------------------------------------------------------------------------
  !>
  !! covariant_velocities
  !!
  !! The model prognoses local orthogonal velocity components.
  !! This subroutine derives the covariant velocity components needed as
  !! entries for the vorticity computations.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-17)
  !!
  SUBROUTINE covariant_velocities (p_vn,p_w,p_patch,p_int,p_metrics,p_diag)

    ! Passed variables
    REAL(wp), INTENT(in)              :: p_vn(:,:,:)
    REAL(wp), INTENT(in)              :: p_w(:,:,:)
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN)     :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag), INTENT(inout)    :: p_diag

    ! Local variables
    INTEGER :: jb, jk, nlen, npromz_e, nblks_e, npromz_c, nblks_c
    REAL(wp):: z_w_k_e(nproma,p_patch%nlev,p_patch%nblks_e)  ! w interpolated to edges
    REAL(wp):: z_w_k(nproma,p_patch%nlev,p_patch%nblks_c)
    INTEGER :: nlev, nlevp1              !< number of full and half levels
    !--------------------------------------------------------------------------

    !> @par Covariant vertical component
    !! @f$ \dot{q}_z = w \sqrt{g} @f$
    !! This includes rigid upper lid if w at the top is given with zero. Note
    !! the different ddqz_z_half at the lower boundary, which should be
    !! consistent with further use. w is supposed to be already given at the
    !! lower boundary, which is ok, as the covariant components are needed in
    !! the slow (rotational) part of the model.

    !> @par Covariant horizontal component
    !! @f$ \dot{q}_n = \dot{x}_n + w J_n @f$
    !! For flat levels covariant and orthogonal component are the same.
    !! For non flat levels, w has to be averaged first vertically,
    !! afterwards it is averaged horizontally. Last, it is multiplyied with
    !! the slope

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk)
    DO jb = 1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1, nlevp1
        p_diag%w_cov(1:nlen,jk,jb)=p_w(1:nlen,jk,jb)*p_metrics%ddqz_z_half(1:nlen,jk,jb)
      ENDDO
      DO jk = nflat+1, nlev
        z_w_k(1:nlen,jk,jb)= 0.5_wp/p_metrics%ddqz_z_full(1:nlen,jk,jb)   &
        &*(p_diag%w_cov(1:nlen,jk,jb)+p_diag%w_cov(1:nlen,jk+1,jb))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    CALL cells2edges_scalar(z_w_k, p_patch, p_int%c_lin_e, &
         &                  z_w_k_e, nflat+1, nlev)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk)
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1, nflat
        p_diag%vn_cov(1:nlen,jk,jb) = p_vn(1:nlen,jk,jb)
      ENDDO
      DO jk = nflat+1,nlev
        p_diag%vn_cov(1:nlen,jk,jb) = p_vn(1:nlen,jk,jb) +  &
        & p_metrics%ddxn_z_full(1:nlen,jk,jb)*z_w_k_e(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE covariant_velocities
  !----------------------------------------------------------------------------
  !>
  !! contravariant_vert_mass_flux
  !!
  !! The model prognosed orthogonal velocity components.
  !! This subroutine derives the contravariant vertical mass flux components,
  !! needed as entry for the implicit vertical advection
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-17)
  !! Modification by Almut Gassmann (2011-04-11)
  !! - rather compute vertical mass fluxes
  !!
  SUBROUTINE contravariant_vert_mass_flux (p_vn,p_rho_e,p_wmflx,p_patch,p_int, &
  &                                        p_metrics,p_mflx_w_con,l_mps)

    ! Passed variables
    REAL(wp), INTENT(in)               :: p_vn(:,:,:)    !< horizontal velocity
    REAL(wp), INTENT(in)               :: p_rho_e(:,:,:) !< density at edges
    REAL(wp), INTENT(in)               :: p_wmflx(:,:,:) !< orth. vert. mass flux
    TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
    TYPE(t_int_state), INTENT(IN)      :: p_int
    TYPE(t_nh_metrics), INTENT(IN)     :: p_metrics
    REAL(wp), INTENT(inout)            :: p_mflx_w_con(:,:,:) ! contravariant vertical mass flux
    LOGICAL, INTENT(in)                :: l_mps !< if true, compute in meters / seconds

    ! Local variables
    INTEGER  :: jk, jb, nblks_c, npromz_c, nblks_e, npromz_e, nlen, nstart
    REAL(wp) :: z_jn_vn_e    (nproma,p_patch%nlev,p_patch%nblks_e), &
                z_metric_corr(nproma,p_patch%nlev,p_patch%nblks_c)
    INTEGER  :: nlev, nlevp1           !< number of full levels
    !--------------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    nblks_c  = p_patch%nblks_int_c
    npromz_c = p_patch%npromz_int_c
    nblks_e  = p_patch%nblks_int_e
    npromz_e = p_patch%npromz_int_e

    !> @par Contravariant vertical component
    !! @f$ \dot{q}^z = (w-\dot{x}_nJ_n)/\sqrt(\gamma) @f$
    !! For the metric correction term the slopes at full levels are first
    !! multiplied with the horizontal velocity by applying the inner product
    !! at the primal grid. The vertical weighting of the metric correction is
    !! done by direct distance weighting. Upper and lower boundary conditions
    !! are zero (rigid upper lid, bottom free slip).
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk)
    DO jb=1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = nflat+1,nlev
        z_jn_vn_e(1:nlen,jk,jb)= p_vn(1:nlen,jk,jb) &
        &                    *p_rho_e(1:nlen,jk,jb) & 
        &      *p_metrics%ddxn_z_full(1:nlen,jk,jb) &
        &    *p_metrics%ddqz_z_full_e(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL edges2cells_scalar(z_jn_vn_e,p_patch,p_int%e_inn_c,z_metric_corr,nflat+1,nlev)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, nstart)
    DO jb = 1,nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = nflat+1,nlev
        z_metric_corr(1:nlen,jk,jb) = z_metric_corr(1:nlen,jk,jb) &
        &                   / p_metrics%ddqz_z_full(1:nlen,jk,jb)
      ENDDO
      p_mflx_w_con(1:nlen,1,jb)    = 0.0_wp ! rigid upper lid
      DO jk = 2, nflat
        p_mflx_w_con(1:nlen,jk,jb) = p_wmflx(1:nlen,jk,jb)
      ENDDO
      nstart = MAX(2,nflat+1)
      z_metric_corr(1:nlen,nstart-1,jb) = 0.0_wp
      DO jk = nstart, nlev
        p_mflx_w_con(1:nlen,jk,jb) = p_wmflx(1:nlen,jk,jb)  &
        & -0.5_wp*(z_metric_corr(1:nlen,jk-1,jb)+z_metric_corr(1:nlen,jk,jb))
      ENDDO
      p_mflx_w_con(1:nlen,nlevp1,jb) = 0.0_wp ! lower boundary cond. (free slip)
      IF(.NOT.l_mps) THEN
        DO jk = 2,nlev
          p_mflx_w_con(1:nlen,jk,jb)= p_mflx_w_con(1:nlen,jk,jb)&
          &                 /p_metrics%ddqz_z_half(1:nlen,jk,jb)
        ENDDO
      ENDIF
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE contravariant_vert_mass_flux
  !----------------------------------------------------------------------------
  !>
  !! orthogonal_vorticities
  !!
  !! The helicity bracket works with local orthogonal components.
  !! This subroutine derives the orthogonal vorticity components from the given
  !! contravariant components.
  !! The general formula for obtaining the vertical orthogonal vorticity is
  !! for a quadrilateral grid:
  !! @f$\xi_z=(\omega^z\sqrt(\gamma)+\omega^xJ_x+\omega^yJ_y) @f$
  !! There and in any triangular or hexagonal grid, the slopes of the terrain
  !! @f$J_n @f$ are only given in normal direction on the edges, but the
  !! corresponding @f$\omega^t @f$ are tangential at this position. Thus,
  !! we must first interpolate both to the same central point and take the
  !! inner product afterwards.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-21)
  !!
  SUBROUTINE orthogonal_vorticities (p_patch,p_int,p_metrics,p_diag)

    ! Passed variables
    TYPE(t_patch), TARGET, INTENT(IN)     :: p_patch
    TYPE(t_int_state), INTENT(IN)         :: p_int
    TYPE(t_nh_metrics), INTENT(IN)        :: p_metrics
    TYPE(t_nh_diag), TARGET, INTENT(inout):: p_diag

    ! Local variables
    REAL(wp) :: z_tmp_e (nproma,p_patch%nlev,p_patch%nblks_e)
    INTEGER  :: jk, jb, jv, nlen, nblks_e, npromz_e, nblks_v, npromz_v
    INTEGER  :: nlev              !< number of full levels
    INTEGER,  DIMENSION(:,:,:), POINTER :: iidx, iblk
    !----------------------------

    iidx    => p_patch%verts%edge_idx
    iblk    => p_patch%verts%edge_blk
    nblks_e =  p_patch%nblks_int_e
    npromz_e=  p_patch%npromz_int_e
    nblks_v =  p_patch%nblks_int_v
    npromz_v=  p_patch%npromz_int_v

    ! number of vertical levels
    nlev = p_patch%nlev

    ! no change in horizontal components
    p_diag%omega_t => p_diag%omega_t_con

!$OMP PARALLEL
    ! vertical component
    ! because I omit dividing and multiplying by sqrt(gamma) consecutively
    ! (see subroutine contravariant_vorticities), we have
    ! ... where the levels are flat
!$OMP WORKSHARE
    p_diag%omega_z(:,1:nflat,:) = p_diag%omega_z_con(:,1:nflat,:)
!$OMP END WORKSHARE

!$OMP DO PRIVATE(jb, nlen, jk)
    ! Metric correction term
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = nflat+1, nlev
        ! First, average the contravariant horizontal components to the main level
        z_tmp_e(1:nlen,jk,jb) = 0.5_wp*&
        & (p_metrics%ddqz_z_half_e(1:nlen,jk  ,jb)*p_diag%omega_t_con(1:nlen,jk  ,jb) &
        & +p_metrics%ddqz_z_half_e(1:nlen,jk+1,jb)*p_diag%omega_t_con(1:nlen,jk+1,jb))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!start a new parallel section here, because the next loop wants the z_tmp_e as neighbors
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk, jv)
    DO jb = 1,nblks_v
      IF (jb /= nblks_v) THEN
        nlen = nproma
      ELSE
        nlen = npromz_v
      ENDIF
      ! Second, determine North and East direction
      ! tria_north gives minus east, tria_east gives north for horizontal
      ! vorticities
      DO jv = 1,nlen
        DO jk = nflat+1, nlev
          p_diag%omega_z(jv,jk,jb) = p_diag%omega_z_con(jv,jk,jb) &
          & + p_metrics%ddnorth_z(jv,jk,jb)/p_metrics%ddqz_z_full_v(jv,jk,jb) &
          & * (z_tmp_e(iidx(jv,jb,1),jk,iblk(jv,jb,1))*p_int%tria_east(1,jv,jb) &
          &   +z_tmp_e(iidx(jv,jb,2),jk,iblk(jv,jb,2))*p_int%tria_east(2,jv,jb) &
          &   +z_tmp_e(iidx(jv,jb,3),jk,iblk(jv,jb,3))*p_int%tria_east(3,jv,jb))&
          & + p_metrics%ddeast_z (jv,jk,jb)/p_metrics%ddqz_z_full_v(jv,jk,jb) &
          & *(-z_tmp_e(iidx(jv,jb,1),jk,iblk(jv,jb,1))*p_int%tria_north(1,jv,jb) &
          &   -z_tmp_e(iidx(jv,jb,2),jk,iblk(jv,jb,2))*p_int%tria_north(2,jv,jb) &
          &   -z_tmp_e(iidx(jv,jb,3),jk,iblk(jv,jb,3))*p_int%tria_north(3,jv,jb))
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE orthogonal_vorticities
  !----------------------------------------------------------------------------
  !>
  !! contravariant_vorticities
  !!
  !! Naturally the vorticity components occur as contravariant components.
  !! They are computed here out of the given covariant velocity components.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-21)
  !!
  SUBROUTINE contravariant_vorticities (p_vn,p_vi,p_patch, p_int, p_metrics, p_diag)

    ! Passed variables
    REAL(wp), INTENT(IN)              :: p_vn(:,:,:) !< horiz. vel. (Runge Kutta) 
    REAL(wp), INTENT(INOUT)           :: p_vi(:,:,:) !< horiz. vel. (vert. impl.)
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN)     :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT)    :: p_diag

    ! Local variables
    INTEGER  :: jb, jk, nlen, nblks_e, npromz_e
    REAL(wp) :: z_wcov_grad(nproma,p_patch%nlevp1,p_patch%nblks_e)
    INTEGER  :: nlev, nlevp1              !< number of full and half levels

    !--------------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !> @par Contravariant vertical vorticity
    !! @remark
    !! Correctly we have to divide by sqrt(gamma) here, but because this
    !! contravariant components are afterwards converted into orthogonal
    !! components, which would again multiply with sqrt(gamma), we omit this in
    !! both cases for efficiency.
    CALL rot_vertex(p_diag%vn_cov,p_patch,p_int,p_diag%omega_z_con)

    ! Contravariant horizontal vorticity.
    !------------------------------------

    CALL grad_fd_norm (p_diag%w_cov,p_patch,z_wcov_grad,1,nlevp1)

    CALL sync_patch_array(SYNC_E,p_patch,p_vi)

    nblks_e  = p_patch%nblks_int_e
    npromz_e = p_patch%npromz_int_e
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk)
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      !vn_cov is no longer needed for something else
      p_diag%vn_cov(1:nlen,:,jb) = p_diag%vn_cov(1:nlen,:,jb) &
      &                 + p_vi(1:nlen,:,jb)-p_vn(1:nlen,:,jb)
      ! upper boundary: w is zero and vn has no vertical shear, level is flat
      p_diag%omega_t_con(1:nlen,1,jb) = 0.0_wp
      ! all inner levels
      DO jk = 2, nlev
        p_diag%omega_t_con(1:nlen,jk,jb) =       ( -z_wcov_grad(1:nlen,jk,jb) &
              +  p_diag%vn_cov(1:nlen,jk-1,jb) - p_diag%vn_cov(1:nlen,jk,jb))&
              /p_metrics%ddqz_z_half_e(1:nlen,jk,jb)
      ENDDO
      ! lower boundary condition: no vertical shear should be seen here,
      p_diag%omega_t_con(1:nlen,nlevp1,jb)= -z_wcov_grad(1:nlen,nlevp1,jb) &
      &                          /p_metrics%ddqz_z_half_e(1:nlen,nlevp1,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE contravariant_vorticities
  !----------------------------------------------------------------------------
  ! >
  !! kinetic_energy
  !!
  !! computes kinetic energy
  SUBROUTINE kinetic_energy (p_vn,p_vi,p_w, p_patch, p_int, p_metrics, p_diag, &
    &                        nloops,l_predictor)

    ! Passed variables
    REAL(wp), INTENT(IN),TARGET       :: p_vn(:,:,:) !< horiz. vel. (Runge Kutta) 
    REAL(wp), INTENT(IN),TARGET       :: p_vi(:,:,:) !< horiz. vel. (vert. impl.)
    REAL(wp), INTENT(IN)              :: p_w(:,:,:)  !< vert. vel.
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN)     :: p_int
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag), INTENT(inout)    :: p_diag  !< single nh diagnostic state
    INTEGER, INTENT(IN)               :: nloops !< during the model run nloops=2
                                                !< but for the total integrals nloops=1
    LOGICAL, INTENT(IN)               :: l_predictor

    ! Local variables
    REAL(wp):: z_kin_hor_c(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp):: z_kin_hor_a(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp):: z_kin_hor_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp):: z_kin_hor_v(nproma,p_patch%nlev,p_patch%nblks_v)
    REAL(wp), POINTER :: z_vn(:,:,:)
    REAL(wp):: z_weight

    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v, &
               nlen, jb, jk, jm
    INTEGER :: nlev              !< number of full levels
    !----------------------------------------------------------------
    ! number of vertical levels
    nlev = p_patch%nlev

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
    nblks_v   = p_patch%nblks_int_v
    npromz_v  = p_patch%npromz_int_v

    IF (l_predictor) THEN
!$OMP PARALLEL
!$OMP WORKSHARE
      p_diag%e_kin(:,:,:) = 0.0_wp
      p_diag%e_kinh(:,:,:) = 0.0_wp
!$OMP END WORKSHARE
!$OMP END PARALLEL
      z_weight = 1.0_wp
    ELSEIF(.NOT.l_predictor) THEN
      z_weight = 0.5_wp
    ENDIF

    DO jm = 1, nloops

      IF (jm==1) THEN
        z_vn => p_vn
      ELSE ! jm==2
        z_vn => p_vi
      ENDIF

      ! Horizontal part
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk)
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO jk = 1, nlev
          z_kin_hor_e(1:nlen,jk,jb) = 0.5_wp*z_vn(1:nlen,jk,jb)**2 &
          &              *p_metrics%ddqz_z_full_e(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      CALL edges2cells_scalar(z_kin_hor_e,p_patch,p_int%e_inn_c,z_kin_hor_c)

      CALL edges2verts_scalar(z_kin_hor_e,p_patch,p_int%e_inn_v,z_kin_hor_v)
      CALL sync_patch_array(SYNC_V,p_patch,z_kin_hor_v)
      CALL verts2cells_scalar(z_kin_hor_v,p_patch,p_int%verts_aw_cells,z_kin_hor_a)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jk = 1,nlev
          z_kin_hor_c(1:nlen,jk,jb) = &
          &  sick_o*z_kin_hor_c(1:nlen,jk,jb)+sick_a*z_kin_hor_a(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO

      ! Merge horizontal and vertical part

!$OMP DO PRIVATE(jb, nlen, jk)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        IF (jm==1) THEN
          DO jk = 1, nlev
            ! 3-dimensional kinetic energy, to be used in the vn-equation...
            ! ... + diagonostic
            p_diag%e_kin(1:nlen,jk,jb) = z_weight*(p_diag%e_kin(1:nlen,jk,jb) + &
            & (z_kin_hor_c(1:nlen,jk,jb)&
            & +0.25_wp*(p_w(1:nlen,jk  ,jb)**2*p_metrics%ddqz_z_half(1:nlen,jk  ,jb) &
            &          +p_w(1:nlen,jk+1,jb)**2*p_metrics%ddqz_z_half(1:nlen,jk+1,jb)))&
            &          /p_metrics%ddqz_z_full(1:nlen,jk,jb))
          ENDDO
        ELSE ! jm==2, 2-dimensional kinetic energy, to be used in the w-equation
          DO jk = 1, nlev
            p_diag%e_kinh(1:nlen,jk,jb) = z_weight*(p_diag%e_kinh(1:nlen,jk,jb) + &
            & z_kin_hor_c(1:nlen,jk,jb)/p_metrics%ddqz_z_full(1:nlen,jk,jb))
          ENDDO
        ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ENDDO

  END SUBROUTINE kinetic_energy
  !----------------------------------------------------------------------------
  ! >
  !! vorticity_tendencies
  !!
  !! Computes the tendencies of the vorticity term in all horizontal and
  !! vertical directions. It needs the orthogonal ingredients of the
  !! velocity and the vorticity vector components. The heli_coeffs computed in
  !! mo_interpolation may be reused here, but additionally the vertical
  !! grid variation has to be taken into account, which may be achieved
  !! by premultiplying some variables with sqrt(g) and afterwards dividing
  !! by sqrt(g).
  !!
  !! @par Revision History
  !! Initial release after a previous odd version by Almut Gassmann (2009-11-08)
  !! Complete rewriting by Almut Gassmann (2011-04-14)
  !!
  SUBROUTINE vorticity_tendencies (p_vn,p_vi,p_wmflx_ort_e,p_rho,p_rho_e,p_rho_e_l,&
  &                                p_patch,p_int,p_metrics,p_diag,l_predictor)

    ! Passed variables
    REAL(wp), INTENT(IN)                :: p_vn(:,:,:) ! horizontal velocity
    REAL(wp), INTENT(IN)                :: p_vi(:,:,:) ! horizontal velocity (implicit)
    REAL(wp), INTENT(IN)                :: p_wmflx_ort_e(:,:,:)  ! vertical orth mass flux
    REAL(wp), INTENT(IN)                :: p_rho(:,:,:)
    REAL(wp), INTENT(IN)                :: p_rho_e(:,:,:)
    REAL(wp), INTENT(IN)                :: p_rho_e_l(:,:,:)
    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state), TARGET,INTENT(IN):: p_int
    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics
    TYPE(t_nh_diag), INTENT(INOUT)      :: p_diag  !< single nh diagnostic state
    LOGICAL, INTENT(IN)                 :: l_predictor

    ! Local variables
    REAL(wp):: z_rho_v            (nproma,p_patch%nlev  ,p_patch%nblks_v)
    REAL(wp):: z_potvort_e_l      (nproma,p_patch%nlevp1,p_patch%nblks_e)
    REAL(wp):: z_hmfl_e           (nproma,p_patch%nlev  ,p_patch%nblks_e)
    REAL(wp):: z_hmfl_ew          (nproma,p_patch%nlev  ,p_patch%nblks_e)
    REAL(wp):: z_ddt_w_e_l        (nproma,p_patch%nlevp1,p_patch%nblks_e)
    REAL(wp):: z_ddt_w_vort       (nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp):: z_potvort_e        (nproma,p_patch%nlev  ,p_patch%nblks_e)
    REAL(wp):: z_potvort_c        (nproma,p_patch%nlev  ,p_patch%nblks_c)
    REAL(wp):: z_potvort_v        (nproma,p_patch%nlev  ,p_patch%nblks_v)
    REAL(wp):: z_nor_o_c          (nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp):: z_tan_o_c          (nproma,p_patch%nlevp1,p_patch%nblks_c)
    REAL(wp):: z_u_c              (nproma,p_patch%nlev  ,p_patch%nblks_c)
    REAL(wp):: z_v_c              (nproma,p_patch%nlev  ,p_patch%nblks_c)
    REAL(wp):: z_v_e              (nproma,p_patch%nlev  ,p_patch%nblks_e)
    REAL(wp):: z_u_e              (nproma,p_patch%nlev  ,p_patch%nblks_e)
    REAL(wp):: z_tmp_v            (nproma,p_patch%nlev  ,p_patch%nblks_v)
    REAL(wp):: z_tmp_v_l          (nproma,p_patch%nlevp1,p_patch%nblks_v)
    REAL(wp):: z_tmp_c_l          (nproma,p_patch%nlevp1,p_patch%nblks_c)

    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v, &
               nlen, nincr, nlev, nlevp1
    INTEGER, DIMENSION(:,:,:), POINTER :: ieidx, ieblk, icidx, icblk
    INTEGER :: jb,je,jk,ji
    !----------------------------------------------------------------

    IF(p_test_run) z_ddt_w_e_l(:,:,:)=0.0_wp

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
    nblks_v   = p_patch%nblks_int_v
    npromz_v  = p_patch%npromz_int_v

    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlevp1

    ! rho at vertex points
    CALL cells2verts_scalar(p_rho,p_patch,p_int%cells_aw_verts,z_rho_v)


    IF (l_predictor) THEN
!$OMP PARALLEL
!$OMP WORKSHARE
     p_diag%ddt_vn_vort(:,:,:)=0.0_wp
     p_diag%ddt_w_vort(:,:,:)=0.0_wp
!$OMP END WORKSHARE
!$OMP END PARALLEL
    ENDIF

 
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1, nlevp1
        ! the horizontal potential vorticity uses still the tilde averaged rho_e
        z_potvort_e_l(1:nlen,jk,jb) = p_diag%omega_t(1:nlen,jk,jb) &
        &                                 /p_rho_e_l(1:nlen,jk,jb) &
        &                   *p_metrics%ddqz_z_half_e(1:nlen,jk,jb)
      ENDDO
      DO jk = 1, nlev
        ! tendency for horizontal equation from horizontal vorticity
        p_diag%ddt_vn_vort(1:nlen,jk,jb) = p_diag%ddt_vn_vort(1:nlen,jk,jb) &
        &- 0.5_wp                    &
        &*(p_wmflx_ort_e(1:nlen,jk  ,jb)*z_potvort_e_l(1:nlen,jk  ,jb) &
        &+ p_wmflx_ort_e(1:nlen,jk+1,jb)*z_potvort_e_l(1:nlen,jk+1,jb)) &
        & /p_metrics%ddqz_z_full_e(1:nlen,jk,jb)
        ! implicit vn to be used together with the horizontal potential vorticity 
        z_hmfl_ew(1:nlen,jk,jb) = p_rho_e(1:nlen,jk,jb)*p_vi(1:nlen,jk,jb)
      ENDDO
      DO jk = 2, nlev
        ! tendency for vertical equation (*sqrtg which goes into e_inn_v/c)
        z_ddt_w_e_l(1:nlen,jk,jb) =      z_potvort_e_l(1:nlen,jk,jb)  &
        & *0.5_wp*(z_hmfl_ew(1:nlen,jk-1,jb)+z_hmfl_ew(1:nlen,jk,jb))
      ENDDO
      IF (i_cori_method <= 2) THEN
        DO jk = 1, nlev
          ! explicit vn to be used in the horizontal equation (for TRISK)
          z_hmfl_e (1:nlen,jk,jb) = p_rho_e(1:nlen,jk,jb)*p_vn(1:nlen,jk,jb)
        ENDDO
      ELSE
        DO jk = 1, nlev
          ! explicit vn to be used in the horizontal equation 
          !(*sqrtg which goes into vector reconstruction)
          z_hmfl_e (1:nlen,jk,jb) = p_rho_e(1:nlen,jk,jb)*p_vn(1:nlen,jk,jb) &
          &                          * p_metrics%ddqz_z_full_e(1:nlen,jk,jb)
        ENDDO
      ENDIF
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_E,p_patch,z_ddt_w_e_l)

    ! Hollingsworth correction for vertical velocity equation
    CALL edges2cells_scalar(z_ddt_w_e_l,p_patch,p_int%e_inn_c,z_ddt_w_vort,2,nlev)
    CALL edges2verts_scalar(z_ddt_w_e_l,p_patch,p_int%e_inn_v,z_tmp_v_l,2,nlev)
    CALL verts2cells_scalar(z_tmp_v_l,p_patch,p_int%verts_aw_cells,z_tmp_c_l,2,nlev)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 2,nlev
        p_diag%ddt_w_vort(1:nlen,jk,jb) = (sick_a*z_tmp_c_l(1:nlen,jk,jb) &
        &                             + sick_o*z_ddt_w_vort(1:nlen,jk,jb))&
        &                           / p_metrics%ddqz_z_half(1:nlen,jk,jb) &
        &                               + p_diag%ddt_w_vort(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,nlen,jk)
    DO jb = 1, nblks_v
      IF (jb /= nblks_v) THEN
        nlen = nproma
      ELSE
        nlen = npromz_v
      ENDIF
      DO jk = 1, nlev
        ! vertical potential vorticity
        z_potvort_v(1:nlen,jk,jb) = &
        &  (p_diag%omega_z(1:nlen,jk,jb)+p_patch%verts%f_v(1:nlen,jb))&
        &  /z_rho_v(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CALL verts2edges_scalar(z_potvort_v,p_patch,p_int%tria_aw_rhom,z_potvort_e)
    CALL sync_patch_array_mult(SYNC_E,p_patch,2,z_potvort_e,z_hmfl_e)

    IF (ltimer) CALL timer_start(timer_corio)

    SELECT CASE (i_cori_method)

    CASE(1,2) ! TRSK (2) or mixTRSK(1)

      ieidx => p_int%heli_vn_idx
      ieblk => p_int%heli_vn_blk

      SELECT CASE (i_cori_method)
      CASE(1)
        nincr = 14
      CASE(2)
        nincr = 10
      END SELECT
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,ji,je)
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
#ifdef __LOOP_EXCHANGE
        DO je = 1, nlen
          DO ji = 1, nincr
            DO jk = 1, nlev
#else
        DO jk = 1, nlev
!CDIR NOLOOPCHG
          DO ji = 1, nincr
            DO je = 1, nlen
#endif
              p_diag%ddt_vn_vort(je,jk,jb) = p_diag%ddt_vn_vort(je,jk,jb)&
              & +p_int%heli_coeff    (ji,je,jb) &
              & *z_hmfl_e    (ieidx(ji,je,jb),jk,ieblk(ji,je,jb))  &
              & *(z_potvort_e(ieidx(ji,je,jb),jk,ieblk(ji,je,jb))  &
              &  +z_potvort_e(je,jk,jb))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    CASE (3) ! my method for Coriolis term

      ieidx => p_patch%edges%quad_idx
      ieblk => p_patch%edges%quad_blk
      icidx => p_patch%edges%cell_idx
      icblk => p_patch%edges%cell_blk

      ! first, reconstruct mass flux vectors at centers of rhombi and hexagons
      CALL edges2cells_scalar(z_hmfl_e,p_patch,p_int%hex_east  ,z_u_c)
      CALL edges2cells_scalar(z_hmfl_e,p_patch,p_int%hex_north ,z_v_c)
      CALL edges2edges_scalar(z_hmfl_e,p_patch,p_int%quad_east ,z_u_e)
      CALL edges2edges_scalar(z_hmfl_e,p_patch,p_int%quad_north,z_v_e)

      ! second, average absolute potential vorticity from rhombi to centers
      IF (l_corner_vort) THEN
        CALL edges2verts_scalar(z_potvort_e,p_patch,p_int%e_1o3_v,z_tmp_v)
        CALL sync_patch_array(SYNC_V,p_patch,z_tmp_v)
        CALL verts2cells_scalar(z_tmp_v,p_patch,p_int%verts_aw_cells,z_potvort_c)
      ELSE
        CALL edges2cells_scalar(z_potvort_e,p_patch,p_int%e_aw_c,z_potvort_c)
      ENDIF

      ! third, multiply the absolute vorticities with the velocities,
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO jk = 1,nlev
          z_u_e(1:nlen,jk,jb) = z_u_e(1:nlen,jk,jb)*z_potvort_e(1:nlen,jk,jb) &
          &                           / p_metrics%ddqz_z_full_r(1:nlen,jk,jb)
          z_v_e(1:nlen,jk,jb) = z_v_e(1:nlen,jk,jb)*z_potvort_e(1:nlen,jk,jb) &
          &                           / p_metrics%ddqz_z_full_r(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO PRIVATE(jb,nlen,jk)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jk = 1,nlev
          z_u_c(1:nlen,jk,jb) = z_u_c(1:nlen,jk,jb)*z_potvort_c(1:nlen,jk,jb) &
          &                             / p_metrics%ddqz_z_full(1:nlen,jk,jb)
          z_v_c(1:nlen,jk,jb) = z_v_c(1:nlen,jk,jb)*z_potvort_c(1:nlen,jk,jb) &
          &                             / p_metrics%ddqz_z_full(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
      ! fourth, compute vorticity flux term
!$OMP DO PRIVATE(jb,nlen,jk,je)
      DO jb = 1, nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
#ifdef __LOOP_EXCHANGE
        DO je = 1, nlen
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = 1, nlen
#endif
            p_diag%ddt_vn_vort(je,jk,jb) = p_diag%ddt_vn_vort(je,jk,jb)&
            & + p_int%heli_coeff( 1,je,jb)*z_v_c(icidx(je,jb,1),jk,icblk(je,jb,1))&
            & + p_int%heli_coeff( 2,je,jb)*z_u_c(icidx(je,jb,1),jk,icblk(je,jb,1))&
            & + p_int%heli_coeff( 3,je,jb)*z_v_c(icidx(je,jb,2),jk,icblk(je,jb,2))&
            & + p_int%heli_coeff( 4,je,jb)*z_u_c(icidx(je,jb,2),jk,icblk(je,jb,2))&
            & + p_int%heli_coeff( 5,je,jb)*z_v_e(je            ,jk,jb            )&
            & + p_int%heli_coeff( 6,je,jb)*z_u_e(je            ,jk,jb            )&
            & + p_int%heli_coeff( 7,je,jb)*z_v_e(ieidx(je,jb,1),jk,ieblk(je,jb,1))&
            & + p_int%heli_coeff( 8,je,jb)*z_u_e(ieidx(je,jb,1),jk,ieblk(je,jb,1))&
            & + p_int%heli_coeff( 9,je,jb)*z_v_e(ieidx(je,jb,2),jk,ieblk(je,jb,2))&
            & + p_int%heli_coeff(10,je,jb)*z_u_e(ieidx(je,jb,2),jk,ieblk(je,jb,2))&
            & + p_int%heli_coeff(11,je,jb)*z_v_e(ieidx(je,jb,3),jk,ieblk(je,jb,3))&
            & + p_int%heli_coeff(12,je,jb)*z_u_e(ieidx(je,jb,3),jk,ieblk(je,jb,3))&
            & + p_int%heli_coeff(13,je,jb)*z_v_e(ieidx(je,jb,4),jk,ieblk(je,jb,4))&
            & + p_int%heli_coeff(14,je,jb)*z_u_e(ieidx(je,jb,4),jk,ieblk(je,jb,4))
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END SELECT ! i_cori_method

    IF (.NOT.l_predictor) THEN
!$OMP PARALLEL
!$OMP WORKSHARE
     p_diag%ddt_vn_vort(:,:,:)=0.5_wp*p_diag%ddt_vn_vort(:,:,:)
     p_diag%ddt_w_vort(:,:,:)=0.5_wp*p_diag%ddt_w_vort(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
    ENDIF

    IF (ltimer) CALL timer_stop(timer_corio)

    ! This part should later be moved somewhere else at the end of the time step
    IF (l_outputtime ) THEN
      CALL edges2cells_scalar(p_diag%omega_t,p_patch,p_int%hex_north,&
                              z_nor_o_c,1,nlevp1)
      CALL edges2cells_scalar(p_diag%omega_t,p_patch,p_int%hex_east,&
                              z_tan_o_c,1,nlevp1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk)
      DO jb = 1, nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        DO jk = 1, nlev
          p_diag%omega_x(1:nlen,jk,jb)=-0.5_wp*(z_nor_o_c(1:nlen,jk  ,jb)+&
          &                                     z_nor_o_c(1:nlen,jk+1,jb))
          p_diag%omega_y(1:nlen,jk,jb)= 0.5_wp*(z_tan_o_c(1:nlen,jk  ,jb)+&
          &                                     z_tan_o_c(1:nlen,jk+1,jb))
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF

  END SUBROUTINE vorticity_tendencies

END MODULE mo_vector_operations

