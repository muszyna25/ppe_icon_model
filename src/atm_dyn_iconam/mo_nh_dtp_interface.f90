!>
!! This module contains subroutines that are used at the
!! dynamics-tracer-physics interface in the non-hydrostatic dynamical core.
!!
!! @author Hui Wan, MPI-M
!! @author Daniel Reinert, DWD
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-02-02)
!! First version for non-hydrostatic core by Daniel Reinert, DWD (2010-04-14)
!! Modification by Daniel Reinert, DWD (2011-02-14)
!! - included computation of nonzero fluxes at the upper boundary 
!!   (necessery if vertical nesting is swithced on)
!! 
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_nh_dtp_interface

  USE mo_kind,               ONLY: wp
  USE mo_dynamics_config,    ONLY: idiv_method
  USE mo_parallel_config,  ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: lvert_nest, ntracer
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_state,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_interpolation,      ONLY: t_int_state
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,     ONLY: min_rledge_int, min_rlcell_int, min_rlcell
  USE mo_sync,               ONLY: SYNC_C, sync_patch_array
  USE mo_advection_config,   ONLY: advection_config

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: prepare_tracer

CONTAINS

  !--------------------------------------------------------------------------
  !>
  !! Diagnose some pressure- and velocity-related quantities
  !! for the tracer transport scheme, under the assumption that
  !! a two time level time stepping scheme is used by the
  !! dynamical core.
  !! Note that boundary fluxes (when using horizontal/vertical nesting) 
  !! are set in mo_nh_nest_utilities/boundary_interpolation.
  !!
  !! @par Revision History
  !! First version by Daniel Reinert, DWD (2010-04-14)
  !! Modification by Daniel Reinert, DWD (2010-07-23)
  !! - adaption to reduced calling frequency
  !!
  SUBROUTINE prepare_tracer( p_patch, p_now, p_new, p_metrics, p_int,  &!in
    &                        iadv_rcf, lstep_advphy, lclean_mflx,      &!in
    &                        p_nh_diag,                                &!inout
    &                        p_vn_traj, p_mass_flx_me,                 &!inout
    &                        p_w_traj, p_mass_flx_ic,                  &!inout
    &                        p_rhodz_mc_now, p_rhodz_mc_new, p_rho_ic, &!inout
    &                        p_topflx_tra                              )!out

    TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch

    TYPE(t_nh_prog),INTENT(IN)    :: p_now, p_new
    TYPE(t_nh_metrics),INTENT(IN) :: p_metrics
    TYPE(t_int_state), INTENT(IN) :: p_int
    TYPE(t_nh_diag),INTENT(INOUT) :: p_nh_diag


    INTEGER :: iadv_rcf     !< used here to switch off time-averaging
    LOGICAL :: lstep_advphy !< used here to switch on time-averaging
    LOGICAL :: lclean_mflx  !< switch for re-initializing time integrated
                            !< mass fluxes and trajectory-velocities

    REAL(wp),INTENT(INOUT)         :: p_vn_traj(:,:,:)      ! (nproma,  nlev,p_patch%nblks_e)
    REAL(wp),INTENT(INOUT)         :: p_w_traj(:,:,:)       ! (nproma,nlevp1,p_patch%nblks_c)
    REAL(wp),INTENT(INOUT), TARGET :: p_mass_flx_me(:,:,:)  ! (nproma,  nlev,p_patch%nblks_e)
    REAL(wp),INTENT(INOUT)         :: p_mass_flx_ic(:,:,:)  ! (nproma,nlevp1,p_patch%nblks_c)
    REAL(wp),INTENT(INOUT)         :: p_rhodz_mc_now(:,:,:) ! (nproma,  nlev,p_patch%nblks_c)
    REAL(wp),INTENT(INOUT)         :: p_rhodz_mc_new(:,:,:) ! (nproma,  nlev,p_patch%nblks_c)
    REAL(wp),INTENT(INOUT)         :: p_rho_ic(:,:,:)       ! (nproma,nlevp1,p_patch%nblks_c)
    REAL(wp),INTENT(OUT)           :: p_topflx_tra(:,:,:)   ! (nproma,p_patch%nblks_c,ntracer)

    ! local variables
    REAL(wp) :: r_iadv_rcf           !< reciprocal of iadv_rcf
    REAL(wp) :: z_mass_flx_me(nproma,p_patch%nlev, p_patch%nblks_e)
    REAL(wp) :: w_tavg               !< contravariant vertical velocity at n+\alpha 

    ! Pointers to quad edge indices
    INTEGER,  POINTER :: iqidx(:,:,:), iqblk(:,:,:)

    INTEGER  :: je, jc, jk, jb, jg, jt    !< loop indices and domain ID
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c, i_nchdom
    INTEGER  :: nlev, nlevp1       !< number of full and half levels
  !--------------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! refinement control start/end level for cells
    i_rlstart_c = 1
    i_rlend_c   = min_rlcell_int

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! domain ID
    jg = p_patch%id

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

!$OMP PARALLEL PRIVATE(i_rlstart_e,i_rlend_e,i_startblk,i_endblk)

    i_rlstart_e = 1
!DR Note that this is not correct: itype_hlimit has the dimension 
!DR MAX_TRACER and not MAX_DOM !!
    IF ( advection_config(jg)%itype_hlimit(jg) == 1 .OR. &
      &  advection_config(jg)%itype_hlimit(jg) == 2 .OR. &
      &  advection_config(jg)%iord_backtraj == 2 ) THEN
      i_rlend_e   = min_rledge_int - 3
    ELSE
      i_rlend_e   = min_rledge_int - 2
    ENDIF
    i_startblk  = p_patch%edges%start_blk(i_rlstart_e,1)
    i_endblk    = p_patch%edges%end_blk(i_rlend_e,i_nchdom)

    !
    ! contravariant normal velocites at n+1/2
    ! necessary for computation of backward trajectories.
    !
    ! horizontal (contravariant) mass flux at full level edges
    ! Taken from dynamical core (corrector-step)
    !
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart_e, i_rlend_e )

      ! reset mass fluxes and trajectory-velocities to start new integration sweep
      IF (lclean_mflx) THEN
        p_vn_traj    (:,:,jb) = 0._wp
        p_mass_flx_me(:,:,jb) = 0._wp
      ENDIF

      DO jk = 1,nlev
        DO je = i_startidx, i_endidx

          ! trajectory-velocity
          p_vn_traj(je,jk,jb) = p_vn_traj(je,jk,jb)                             &
            &             + 0.5_wp * ( p_now%vn(je,jk,jb) + p_new%vn(je,jk,jb) )

          ! mass flux
          p_mass_flx_me(je,jk,jb) = p_mass_flx_me(je,jk,jb)       &
            &                       + p_nh_diag%mass_fl_e(je,jk,jb)

          ! Note that p_mass_flx_me must be used for summation instead of
          ! z_mass_flx_me. The latter looses its information after leaving this
          ! routine. Therefore the following copy-command is necessary. A SAVE-attribute
          ! would nicely do the job, but is not allowed for this kind of dynamic array.
          z_mass_flx_me(je,jk,jb) = p_mass_flx_me(je,jk,jb)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

    IF (p_test_run .AND. lclean_mflx) THEN ! Reset also halo points to zero
!$OMP WORKSHARE
        p_vn_traj    (:,:,i_endblk+1:p_patch%nblks_e) = 0._wp
        p_mass_flx_me(:,:,i_endblk+1:p_patch%nblks_e) = 0._wp
!$OMP END WORKSHARE
    ENDIF


    i_startblk   = p_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

    !
    ! Time averaged contravariant vertical velocity
    ! and vertical mass flux at half level centers
    !
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,w_tavg)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

      ! reset mass fluxes and trajectory-velocities to start new integration sweep
      IF (lclean_mflx) THEN
        p_mass_flx_ic(:,:,jb) = 0._wp
        p_w_traj     (:,:,jb) = 0._wp
        p_rho_ic     (:,:,jb) = 0._wp
      ENDIF

      DO jk = 1, nlevp1
        DO jc = i_startidx, i_endidx

          p_rho_ic(jc,jk,jb) = p_rho_ic(jc,jk,jb) + p_nh_diag%rho_ic(jc,jk,jb)

! Note(DR): This is somewhat inconsistent since for horizontal trajectories
! v_n at n+1/2 is used.
          w_tavg = p_metrics%vwind_expl_wgt(jc,jb)*p_now%w(jc,jk,jb)  &
            &    + p_metrics%vwind_impl_wgt(jc,jb)*p_new%w(jc,jk,jb)  &
            &    - p_nh_diag%w_concorr_c(jc,jk,jb)

          p_w_traj(jc,jk,jb) = p_w_traj(jc,jk,jb) + w_tavg

          p_mass_flx_ic(jc,jk,jb) = p_mass_flx_ic(jc,jk,jb)              &
            &                     + p_nh_diag%rho_ic(jc,jk,jb) * w_tavg

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

    IF (p_test_run .AND. lclean_mflx) THEN ! Reset also halo points to zero
!$OMP WORKSHARE
        p_mass_flx_ic(:,:,i_endblk+1:p_patch%nblks_c) = 0._wp
        p_w_traj     (:,:,i_endblk+1:p_patch%nblks_c) = 0._wp
        p_rho_ic     (:,:,i_endblk+1:p_patch%nblks_c) = 0._wp
!$OMP END WORKSHARE
    ENDIF

    !
    ! density multiplied by layer thickness at cell center for
    ! timesteps n and n+1 (need to include the halo points!)
    !
    IF (lclean_mflx) THEN   ! first time step after call of transport

      i_endblk = p_patch%cells%end_blk(min_rlcell,i_nchdom)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, min_rlcell )

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            p_rhodz_mc_now(jc,jk,jb) =                                          &
              &           p_now%rho(jc,jk,jb) * p_metrics%ddqz_z_full(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF

    IF (lstep_advphy) THEN

      i_endblk = p_patch%cells%end_blk(min_rlcell,i_nchdom)
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, min_rlcell )

        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx

            p_rhodz_mc_new(jc,jk,jb) =                                          &
              &           p_new%rho(jc,jk,jb) * p_metrics%ddqz_z_full(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

    ENDIF



    !
    ! compute time averaged mass fluxes and trajectory-velocities
    !
    IF ( lstep_advphy .AND. iadv_rcf > 1 ) THEN


      r_iadv_rcf = 1._wp/REAL(iadv_rcf,wp)

      i_rlstart_e  = 2
!DR Note that this is not correct: itype_hlimit has the dimension 
!DR MAX_TRACER and not MAX_DOM !!
      IF ( advection_config(jg)%itype_hlimit(jg) == 1 .OR. &
        &  advection_config(jg)%itype_hlimit(jg) == 2 .OR. &
        &  advection_config(jg)%iord_backtraj == 2) THEN
        i_rlend_e   = min_rledge_int - 3
      ELSE
        i_rlend_e   = min_rledge_int - 2
      ENDIF
      i_startblk   = p_patch%edges%start_blk(i_rlstart_e,1)
      i_endblk     = p_patch%edges%end_blk(i_rlend_e,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_e, i_rlend_e )

        IF (idiv_method == 1) THEN
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx

              p_mass_flx_me(je,jk,jb) = r_iadv_rcf * p_mass_flx_me(je,jk,jb)
              p_vn_traj(je,jk,jb)     = r_iadv_rcf * p_vn_traj(je,jk,jb)

            ENDDO
          ENDDO
        ELSE ! use averaged mass fluxes for approximate consistency with averaged divergence
#ifdef __LOOP_EXCHANGE
          DO je = i_startidx, i_endidx
            DO jk = 1, nlev
#else
!CDIR UNROLL=3
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx
#endif

              p_mass_flx_me(je,jk,jb) = r_iadv_rcf * (                                       &
                &   z_mass_flx_me(je,jk,jb)*p_int%e_flx_avg(je,1,jb)                         &
                & + z_mass_flx_me(iqidx(je,jb,1),jk,iqblk(je,jb,1))*p_int%e_flx_avg(je,2,jb) &
                & + z_mass_flx_me(iqidx(je,jb,2),jk,iqblk(je,jb,2))*p_int%e_flx_avg(je,3,jb) &
                & + z_mass_flx_me(iqidx(je,jb,3),jk,iqblk(je,jb,3))*p_int%e_flx_avg(je,4,jb) &
                & + z_mass_flx_me(iqidx(je,jb,4),jk,iqblk(je,jb,4))*p_int%e_flx_avg(je,5,jb) )

              p_vn_traj(je,jk,jb)     = r_iadv_rcf * p_vn_traj(je,jk,jb)

            ENDDO
          ENDDO
        ENDIF
      ENDDO
!$OMP END DO NOWAIT


      i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

        DO jk = 1, nlevp1
          DO jc = i_startidx, i_endidx

            p_mass_flx_ic(jc,jk,jb) = r_iadv_rcf * p_mass_flx_ic(jc,jk,jb)

            p_w_traj(jc,jk,jb)      = r_iadv_rcf * p_w_traj(jc,jk,jb)

            p_rho_ic(jc,jk,jb)      = r_iadv_rcf * p_rho_ic(jc,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

    ELSE IF ( lstep_advphy .AND. iadv_rcf == 1 .AND. idiv_method == 2 ) THEN
      ! Compute only averaged mass flux


      i_rlstart_e = 2
      i_rlend_e   = min_rledge_int - 2
      i_startblk  = p_patch%edges%start_blk(i_rlstart_e,1)
      i_endblk    = p_patch%edges%end_blk(i_rlend_e,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_e, i_rlend_e )

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
!CDIR UNROLL=3
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            p_mass_flx_me(je,jk,jb) =     z_mass_flx_me(je,jk,jb)*p_int%e_flx_avg(je,1,jb) &
              & + z_mass_flx_me(iqidx(je,jb,1),jk,iqblk(je,jb,1))*p_int%e_flx_avg(je,2,jb) &
              & + z_mass_flx_me(iqidx(je,jb,2),jk,iqblk(je,jb,2))*p_int%e_flx_avg(je,3,jb) &
              & + z_mass_flx_me(iqidx(je,jb,3),jk,iqblk(je,jb,3))*p_int%e_flx_avg(je,4,jb) &
              & + z_mass_flx_me(iqidx(je,jb,4),jk,iqblk(je,jb,4))*p_int%e_flx_avg(je,5,jb)

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

    ENDIF

!$OMP END PARALLEL

   IF (lstep_advphy) CALL sync_patch_array(SYNC_C,p_patch,p_mass_flx_ic)


      !  
      ! diagnose vertical tracer fluxes at top margin, i.e. multiply horizontally 
      ! interpolated face value q_ubc with time averaged mass flux at nested 
      ! domain top. Since we make direct use of the dycore mass flux, this procedure 
      ! ensures tracer and air mass consistency.
      !
   IF (lvert_nest .AND. (p_patch%nshift > 0)) THEN ! vertical nesting

     i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
     i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

     DO jb = i_startblk, i_endblk
       CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
         &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

       DO jt = 1, ntracer
         DO jc = i_startidx, i_endidx
           p_topflx_tra(jc,jb,jt) = p_nh_diag%q_ubc(jc,jb,jt)           &
             &                    * p_mass_flx_ic(jc,1,jb)
         ENDDO
       ENDDO
     ENDDO

   ELSE                 ! no vertical nesting
     p_topflx_tra(:,:,:) = 0._wp
   ENDIF

  END SUBROUTINE prepare_tracer

END MODULE mo_nh_dtp_interface
