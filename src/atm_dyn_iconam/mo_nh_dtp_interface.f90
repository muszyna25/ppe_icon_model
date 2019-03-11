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
!!   (necessary if vertical nesting is switched on)
!! Modification by William Sawyer, CSCS (2015-02-06)
!! - OpenACC implementation
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_dtp_interface

  USE mo_kind,               ONLY: wp
  USE mo_dynamics_config,    ONLY: idiv_method
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: lvert_nest, ntracer
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants,     ONLY: min_rledge_int, min_rlcell_int, min_rlcell
  USE mo_sync,               ONLY: SYNC_C, sync_patch_array, sync_patch_array_mult
  USE mo_advection_config,   ONLY: advection_config, t_trList
  USE mo_initicon_config,    ONLY: is_iau_active, iau_wgt_adv
  USE mo_timer,              ONLY: timers_level, timer_start, timer_stop, timer_prep_tracer
  USE mo_fortran_tools,      ONLY: init
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: prepare_tracer
  PUBLIC :: compute_airmass

#if defined( _OPENACC )
#if defined(__NH_DTP_DIFFUSSION_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.     !  THIS SHOULD BE .FALSE. AFTER VALIDATION PHASE!
#endif

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
  !! Modification by Daniel Reinert, DWD (2013-05-06)
  !! - removed rho_ic which became obsolete after removing the second order
  !!   MUSCL scheme for vertical transport
  !!
  SUBROUTINE prepare_tracer( p_patch, p_now, p_new, p_metrics, p_int,         &!in
    &                        ndyn_substeps, lstep_advphy, lclean_mflx,        &!in
    &                        lfull_comp,                                      &!in
    &                        p_nh_diag,                                       &!inout
    &                        p_vn_traj, p_mass_flx_me,                        &!inout
    &                        p_mass_flx_ic,                                   &!inout
    &                        p_topflx_tra                                     )!out

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch

    TYPE(t_nh_prog),INTENT(IN)    :: p_now, p_new
    TYPE(t_nh_metrics),INTENT(IN) :: p_metrics
    TYPE(t_int_state), INTENT(IN) :: p_int
    TYPE(t_nh_diag),INTENT(INOUT) :: p_nh_diag


    INTEGER, INTENT(IN) :: ndyn_substeps !< used here to switch off time-averaging
    LOGICAL, INTENT(IN) :: lstep_advphy  !< used here to switch on time-averaging
    LOGICAL, INTENT(IN) :: lclean_mflx   !< switch for re-initializing time integrated
                                         !< mass fluxes and trajectory-velocities
    LOGICAL, INTENT(IN) :: lfull_comp    !< perform full amount of computations (comes in as .FALSE. if
                                         !< part of the precomputations has already been done in
                                         !< solve_nh and only standard settings are used)

    REAL(wp),INTENT(INOUT) :: p_vn_traj(:,:,:)      ! (nproma,  nlev,p_patch%nblks_e)
    REAL(wp),INTENT(INOUT) :: p_mass_flx_me(:,:,:)  ! (nproma,  nlev,p_patch%nblks_e)
    REAL(wp),INTENT(INOUT) :: p_mass_flx_ic(:,:,:)  ! (nproma,nlevp1,p_patch%nblks_c)
    REAL(wp),INTENT(OUT)   :: p_topflx_tra(:,:,:)   ! (nproma,p_patch%nblks_c,ntracer)

    ! local variables
    REAL(wp) :: r_ndyn_substeps           !< reciprocal of ndyn_substeps
    REAL(wp) :: z_mass_flx_me(nproma,p_patch%nlev, p_patch%nblks_e)
    REAL(wp) :: z_topflx_tra (nproma,ntracer, p_patch%nblks_c)
    REAL(wp) :: w_tavg               !< contravariant vertical velocity at n+\alpha

    ! Pointers to quad edge indices
    INTEGER,  POINTER :: iqidx(:,:,:), iqblk(:,:,:)

    INTEGER  :: je, jc, jk, jb, jg, jt, nt    !< loop indices and domain ID
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart_e, i_rlend_e, i_rlstart_c, i_rlend_c, i_nchdom
    INTEGER  :: nlev, nlevp1       !< number of full and half levels

    TYPE(t_trList), POINTER :: trAdvect      !< Pointer to tracer sublist
  !--------------------------------------------------------------------------
    IF (timers_level > 2) CALL timer_start(timer_prep_tracer)

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

    ! tracer fields which are advected
    trAdvect => advection_config(jg)%trAdvect

    ! Set pointers to quad edges
    iqidx => p_patch%edges%quad_idx
    iqblk => p_patch%edges%quad_blk

!!$    ! The full set of setup computations is NOT executed in prepare_tracer 
!!$    ! when the tracer advection is running together with the dynmical core 
!!$    ! (solve_nh) and only standard namelist settings are chosen (i.e. flux limiter,
!!$    ! first-order backward trajectory computation, CFL-safe vertical advection, idiv_method = 1)
!!$    !
!!$    ! lfull_comp is only used by the nonhydrostatic core.
!!$    lfull_computations = lfull_comp
!!$    IF ( ANY( advection_config(jg)%itype_hlimit(1:ntracer) == 1 )     .OR. &
!!$      &  ANY( advection_config(jg)%itype_hlimit(1:ntracer) == 2 )     .OR. &
!!$      &  advection_config(jg)%iord_backtraj == 2                      .OR. &
!!$      &  idiv_method  == 2                                            .OR. &
!!$      &  itime_scheme == TRACER_ONLY                                       ) THEN
!!$      lfull_computations = .TRUE.
!!$    ENDIF


!$OMP PARALLEL PRIVATE(i_rlstart_e,i_rlend_e,i_startblk,i_endblk)

    i_rlstart_e = 1
    IF ( ANY( advection_config(jg)%itype_hlimit(1:ntracer) == 1 ) .OR. &
      &  ANY( advection_config(jg)%itype_hlimit(1:ntracer) == 2 ) .OR. &
      &  advection_config(jg)%iord_backtraj == 2 .OR. idiv_method == 2) THEN
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
    IF (lfull_comp) THEN
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_e, i_rlend_e )

        ! reset mass fluxes and trajectory-velocities to start new integration sweep
        IF (lclean_mflx) THEN
          p_vn_traj    (:,:,jb) = 0._wp
          p_mass_flx_me(:,:,jb) = 0._wp
        ENDIF

        DO jk = 1,nlev
!DIR$ IVDEP
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
    ENDIF

    IF (lfull_comp .AND. p_test_run .AND. lclean_mflx) THEN ! Reset also halo points to zero
      CALL init(p_vn_traj    (:,:,i_endblk+1:p_patch%nblks_e))
      CALL init(p_mass_flx_me(:,:,i_endblk+1:p_patch%nblks_e))
!$OMP BARRIER
    ENDIF


    i_startblk   = p_patch%cells%start_blk(i_rlstart_c,1)
    i_endblk     = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

    !
    ! Time averaged contravariant vertical velocity
    ! and vertical mass flux at half level centers
    !
    IF (lfull_comp) THEN
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,w_tavg) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

        ! reset mass fluxes and trajectory-velocities to start new integration sweep
        IF (lclean_mflx) THEN
          p_mass_flx_ic(:,:,jb) = 0._wp
        ENDIF

        DO jk = 1, nlevp1
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

! Note(DR): This is somewhat inconsistent since for horizontal trajectories
! v_n at n+1/2 is used.
! w_concorr_c at TL n+1
            w_tavg = p_metrics%vwind_expl_wgt(jc,jb)*p_now%w(jc,jk,jb)  &
              &    + p_metrics%vwind_impl_wgt(jc,jb)*p_new%w(jc,jk,jb)  &
              &    - p_nh_diag%w_concorr_c(jc,jk,jb)

            p_mass_flx_ic(jc,jk,jb) = p_mass_flx_ic(jc,jk,jb)              &
              &                     + p_nh_diag%rho_ic(jc,jk,jb) * w_tavg

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

    IF (lfull_comp .AND. p_test_run .AND. lclean_mflx) THEN ! Reset also halo points to zero
      CALL init(p_mass_flx_ic(:,:,i_endblk+1:p_patch%nblks_c))
!$OMP BARRIER
    ENDIF



    !
    ! compute time averaged mass fluxes and trajectory-velocities
    !
    IF (lfull_comp .AND. lstep_advphy .AND. ndyn_substeps > 1 ) THEN


      r_ndyn_substeps = 1._wp/REAL(ndyn_substeps,wp)

      i_rlstart_e  = 2
      IF ( ANY( advection_config(jg)%itype_hlimit(1:ntracer) == 1) .OR. &
        &  ANY( advection_config(jg)%itype_hlimit(1:ntracer) == 2) .OR. &
        &  advection_config(jg)%iord_backtraj == 2) THEN
        i_rlend_e   = min_rledge_int - 3
      ELSE
        i_rlend_e   = min_rledge_int - 2
      ENDIF
      i_startblk   = p_patch%edges%start_blk(i_rlstart_e,1)
      i_endblk     = p_patch%edges%end_blk(i_rlend_e,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_e, i_rlend_e )

        IF (idiv_method == 1) THEN
          DO jk = 1, nlev
            DO je = i_startidx, i_endidx

              p_mass_flx_me(je,jk,jb) = r_ndyn_substeps * p_mass_flx_me(je,jk,jb)
              p_vn_traj(je,jk,jb)     = r_ndyn_substeps * p_vn_traj(je,jk,jb)

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

              p_mass_flx_me(je,jk,jb) = r_ndyn_substeps * (                                  &
                &   z_mass_flx_me(je,jk,jb)*p_int%e_flx_avg(je,1,jb)                         &
                & + z_mass_flx_me(iqidx(je,jb,1),jk,iqblk(je,jb,1))*p_int%e_flx_avg(je,2,jb) &
                & + z_mass_flx_me(iqidx(je,jb,2),jk,iqblk(je,jb,2))*p_int%e_flx_avg(je,3,jb) &
                & + z_mass_flx_me(iqidx(je,jb,3),jk,iqblk(je,jb,3))*p_int%e_flx_avg(je,4,jb) &
                & + z_mass_flx_me(iqidx(je,jb,4),jk,iqblk(je,jb,4))*p_int%e_flx_avg(je,5,jb) )

              p_vn_traj(je,jk,jb)     = r_ndyn_substeps * p_vn_traj(je,jk,jb)

            ENDDO
          ENDDO
        ENDIF
      ENDDO
!$OMP END DO NOWAIT


      i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

        DO jk = 1, nlevp1
          DO jc = i_startidx, i_endidx

            p_mass_flx_ic(jc,jk,jb) = r_ndyn_substeps * p_mass_flx_ic(jc,jk,jb)

          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

    ELSE IF ( lstep_advphy .AND. ndyn_substeps == 1 .AND. idiv_method == 2 ) THEN
      ! Compute only averaged mass flux


      i_rlstart_e = 2
      i_rlend_e   = min_rledge_int - 2
      i_startblk  = p_patch%edges%start_blk(i_rlstart_e,1)
      i_endblk    = p_patch%edges%end_blk(i_rlend_e,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
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
!$OMP END DO NOWAIT

    ENDIF

!$OMP END PARALLEL

    !
    ! diagnose vertical tracer fluxes at top margin, i.e. multiply horizontally
    ! interpolated face value q_ubc with time averaged mass flux at nested
    ! domain top. Since we make direct use of the dycore mass flux, this procedure
    ! ensures tracer and air mass consistency.
    !
    IF (lstep_advphy .AND. lvert_nest .AND. p_patch%nshift > 0) THEN ! vertical nesting

      IF (p_test_run) z_topflx_tra(:,:,:) = 0._wp

      i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend_c,i_nchdom)

!$OMP PARALLEL DO PRIVATE(jb,jt,jc,nt,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

        DO nt = 1, trAdvect%len ! Tracer loop

          jt = trAdvect%list(nt)
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            z_topflx_tra(jc,jt,jb) = p_nh_diag%q_ubc(jc,jb,jt)           &
              &                    * p_mass_flx_ic(jc,1,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL sync_patch_array_mult(SYNC_C,p_patch,2,p_mass_flx_ic,z_topflx_tra)

      i_startblk = p_patch%cells%start_blk(i_rlstart_c,1)
      i_endblk   = p_patch%cells%end_blk(min_rlcell,i_nchdom)

!$OMP PARALLEL DO PRIVATE(jb,jt,jc,nt,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, min_rlcell)

        DO nt = 1, trAdvect%len ! Tracer loop

          jt = trAdvect%list(nt)

          DO jc = i_startidx, i_endidx
            p_topflx_tra(jc,jb,jt) = z_topflx_tra(jc,jt,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

    ELSE IF (lstep_advphy) THEN               ! no vertical nesting

      p_topflx_tra(:,:,:) = 0._wp
      CALL sync_patch_array(SYNC_C,p_patch,p_mass_flx_ic)

    ENDIF


    IF (timers_level > 2) CALL timer_stop(timer_prep_tracer)

  END SUBROUTINE prepare_tracer



  !>
  !! Compute updated air mass within grid cell
  !!
  !! Compute updated air mass within grid cell. Note that here, the air mass is defined
  !! as \rho*\Delta z [kg m-2]. Thus, in order to get the true grid cell air mass one has to
  !! multiply with the grid cell area.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2014-03-15)
  !!
  SUBROUTINE compute_airmass (p_patch, p_metrics, p_prog, p_nh_diag, itlev)

    TYPE(t_patch)   ,  INTENT(IN)     :: p_patch
    TYPE(t_nh_metrics),INTENT(IN)     :: p_metrics
    TYPE(t_nh_prog) ,  INTENT(IN)     :: p_prog
    TYPE(t_nh_diag),   INTENT(INOUT)  :: p_nh_diag
    INTEGER         ,  INTENT(IN)     :: itlev   ! 1: computation of time level now
                                                 ! 2: computation of time level new

    INTEGER :: nlev                  ! number of vertical levels
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: i_nchdom
    INTEGER :: jc,jk,jb
  !---------------------------------------------------------!

    ! number of vertical levels
    nlev = p_patch%nlev

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    ! halo points must be included !
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$ACC DATA PRESENT( p_patch, p_nh_diag, p_prog, p_metrics ), IF( i_am_accel_node .AND. acc_on )

!TODO: add debugging mode

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

       CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
         &                 i_startidx, i_endidx, i_rlstart, i_rlend)

       IF (itlev == 1) THEN
!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
         !$ACC LOOP GANG
         DO jk = 1, nlev
!DIR$ IVDEP
           !$ACC LOOP VECTOR
           DO jc= i_startidx, i_endidx
             p_nh_diag%airmass_now(jc,jk,jb) = p_prog%rho(jc,jk,jb)*p_metrics%ddqz_z_full(jc,jk,jb)
           ENDDO  ! jc
         ENDDO  ! jk
!$ACC END PARALLEL
       ELSE !  itlev = 2
!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
         !$ACC LOOP GANG
         DO jk = 1, nlev
!DIR$ IVDEP
           !$ACC LOOP VECTOR
           DO jc= i_startidx, i_endidx
             p_nh_diag%airmass_new(jc,jk,jb) = p_prog%rho(jc,jk,jb)*p_metrics%ddqz_z_full(jc,jk,jb)
           ENDDO  ! jc
         ENDDO  ! jk
!$ACC END PARALLEL
         IF (is_iau_active) THEN ! Correct 'old' air mass for IAU density increments
!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
           !$ACC LOOP GANG
           DO jk = 1, nlev
!DIR$ IVDEP
             !$ACC LOOP VECTOR
             DO jc= i_startidx, i_endidx
               p_nh_diag%airmass_now(jc,jk,jb) = p_nh_diag%airmass_now(jc,jk,jb) + &
                 iau_wgt_adv*p_metrics%ddqz_z_full(jc,jk,jb)*p_nh_diag%rho_incr(jc,jk,jb)
             ENDDO  ! jc
           ENDDO  ! jk
!$ACC END PARALLEL
         ENDIF
       ENDIF

    ENDDO ! jb
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

!TODO: add debugging mode

!$ACC END DATA

  END SUBROUTINE compute_airmass


END MODULE mo_nh_dtp_interface

