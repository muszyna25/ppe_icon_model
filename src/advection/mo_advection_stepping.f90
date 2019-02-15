!>
!! Contains implementation of the advection time stepping for the hydrostatic
!! and non-hydrostatic dynamical core.
!!
!! Performs time integration of tracer continuity equations in flux form
!! using either Marchuk or Strang splitting between the horizontal and
!! vertical direction.
!!
!! @author Jochen Foerstner, DWD
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Jochen Foerstner, DWD (2008-05-15)
!! Modification by Jochen Foerstner, DWD (2008-07-23)
!! - Implementation of a first version of the MPDATA scheme.
!! Modified by Marco Giorgetta, MPI-M (2009-02-26)
!! - renamed tracer_ctl to transport_ctl
!! - renamed setup_tracer to setup_transport
!! Modification by Daniel Reinert, DWD (2009-06-29)
!! - deleted MPDATA related obsolete options
!! - OpenMP-parallelization of MPDATA
!! Modification by A. Gassmann, MPI-M (2009-10)
!! - remove luse_rbf_vec_int switches
!! Modification by Daniel Reinert (2009-12-16)
!! - removed subroutine init\_rk\_tracer and option for Runge-Kutta type
!!   time integration. Only explicit Euler scheme is retained.
!! Modification by Daniel Reinert (2010-02-09)
!! - restructuring of code. Subroutines for the calculation of horizontal
!!   and vertical fluxes have been transfered to separate modules.
!! Modification by Daniel Reinert (2010-03-05)
!! - moved NAMELIST and 'setup_transport' to 'mo_advection_utils'
!! Modification by Daniel Reinert (2010-04-22)
!! - added field p_w_contra_traj which is necessary, when running the
!!   NH-core.
!! Modification by Will Sawyer, CSCS (2016-07-15)
!! - added OpenACC support
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
MODULE mo_advection_stepping

  USE mo_kind,                ONLY: wp, vp
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_transport
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ntracer, ltimer, iforcing, iqv, iqtke
  USE mo_advection_hflux,     ONLY: hor_upwind_flux
  USE mo_advection_vflux,     ONLY: vert_upwind_flux
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, min_rlcell, &
    &                               inwp
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array_mult
  USE mo_advection_config,    ONLY: advection_config, t_trList
  USE mo_advection_utils,     ONLY: ptr_delp_mc_now, ptr_delp_mc_new
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_fortran_tools,       ONLY: negative2zero
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node, my_process_is_work
  USE mo_sync,                ONLY: SYNC_E, SYNC_C, check_patch_array
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: step_advection

#if defined( _OPENACC )
#if defined(__ADVECTION_STEPPING_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
  LOGICAL, PARAMETER ::  acc_validate = .FALSE.   ! ONLY SET TO .TRUE. FOR VALIDATION PHASE
#endif

CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Time stepping for tracer transport
  !!
  !!
  !! @par Revision History
  !! Initial revision by Jochen Foerstner, DWD (2008-05-15)
  !! Modification by Jochen Foerstner, DWD (2008-07-09)
  !! - included Runge-Kutta type time integration for tracer transport.
  !! Modification by Daniel Reinert, DWD (2009-06-29)
  !! - deleted MPDATA related obsolete options
  !! - OpenMP-parallelization
  !! Modification by Daniel Reinert, DWD (2009-08-18)
  !! - MPDATA and Godunov-Type methods are now treated separately. This was
  !!   mainly done in order to improve readability.
  !! Modification by Daniel Reinert, DWD (2009-12-16)
  !! - removed Runge-Kutta type time integration for tracer transport.
  !! Modification by Daniel Reinert, DWD (2010-01-29)
  !! - removed MPDATA
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - removed non-conservative option 'nonCnsv'
  !! Modification by Daniel Reinert, DWD (2010-02-22)
  !! - generalization of interface: avoid passing types related to the
  !!   hydrostatic or nonhydrostatic model
  !! Modification by Daniel Reinert, DWD (2010-02-25)
  !! - re-integration of mass continuity equation in order to achieve
  !!   tracer-mass consistency
  !! Modification by Daniel Reinert, DWD (2010-10-25)
  !! - removed boundary flux correction
  !! Modification by Daniel Reinert, DWD (2010-11-02)
  !! - splitting of tracer loop. -> reduced OpenMP/MPI overhead
  !! Modification by Daniel Reinert, DWD (2011-02-15)
  !! - new field providing the upper margin tracer flux (optional)
  !! Modification by Daniel Reinert, DWD (2013-05-06)
  !! - simplification of step_advection and vert_upwind_flux interfaces 
  !!   after removal of the MUSCL vertical advection scheme
  !! 
  !!
  SUBROUTINE step_advection( p_patch, p_int_state, p_dtime, k_step, p_tracer_now, &
    &                        p_mflx_contra_h, p_vn_contra_traj, p_mflx_contra_v,  &
    &                        p_w_contra_traj, p_cellhgt_mc_now, p_delp_mc_new,    &
    &                        p_delp_mc_now, p_grf_tend_tracer, p_tracer_new,      &
    &                        p_mflx_tracer_h, p_mflx_tracer_v, opt_topflx_tra,    &
    &                        opt_q_int, opt_ddt_tracer_adv )
  !
    TYPE(t_patch), TARGET, INTENT(INOUT) ::  &  !< patch on which computation
      &  p_patch                             !< is performed
                                             

    !> interpolation state
    TYPE(t_int_state), INTENT(IN) :: p_int_state

    REAL(wp), TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     , INTENT(IN) ::  &  !< tracer mixing ratios (specific concentrations)
      &  p_tracer_now(:,:,:,:)          !< at current time level n (before transport)
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(IN)  ::         &  !< horizontal mass flux (contravariant)
      &  p_mflx_contra_h(:,:,:)         !< NH: v_n*delta_z*\rho  [kg/m/s]
                                        !< HA: v_n*delta_p       [kg/s**3]
                                        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN)  ::         &  !< horizontal velocity component at n+1/2
      &  p_vn_contra_traj(:,:,:)        !< for calculation of backward trajectories
                                        !< [m/s] 
                                        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(INOUT)  ::      &  !< vertical mass flux (contravariant)
      &  p_mflx_contra_v(:,:,:)         !< NH: \rho*w     [kg/m**2/s]
                                        !< HA: eta_dot \partial p/\partial eta [Pa/s]
                                        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN)  ::         &  !< vertical velocity (contravariant) at n+1/2
      &  p_w_contra_traj(:,:,:)         !< for calculation of backward trajectories
                                        !< NH: w [m/s]
                                        !< HA: eta_dot \partial p/\partial eta [Pa/s]
                                        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::          &  !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)        !< time step n 
                                        !< NH: \Delta z       [m]
                                        !< HA: \Delta p       [Pa]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), TARGET, INTENT(IN) ::  &  !< NH: density weighted cell height at full levels 
      &  p_delp_mc_new(:,:,:)           !< at n+1 [kg/m**2]
                                        !< HA: pressure thickness for full levels 
                                        !< at n+1 [Pa]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), TARGET, INTENT(IN) ::  &  !< NH: density weighted cell height at full levels
      &  p_delp_mc_now(:,:,:)           !< at time step n [kg/m**2]
                                        !< HA: pressure thickness for full levels 
                                        !< at time step n [Pa]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::       &  !< interpolated tracer time tendencies for
      &  p_grf_tend_tracer(:,:,:,:)     !< updating the lateral boundaries of nested domains
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), TARGET &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     , INTENT(INOUT) ::  &  !< tracer mixing ratios (specific concentrations)
      &  p_tracer_new(:,:,:,:)             !< at time level n+1 (after transport)
                                           !< [kg/kg]  
                                           !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT)  ::  &   !< horizontal tracer mass flux at full level edges
      &  p_mflx_tracer_h(:,:,:,:)    !< NH: [kg/m/s]
                                     !< HA: [kg/s**3]
                                     !< dim: (nproma,nlev,nblks_e,ntracer)

    REAL(wp), INTENT(INOUT)  ::  &   !< vertical tracer mass flux at half level centers
      &  p_mflx_tracer_v(:,:,:,:)    !< NH: [kg/m**2/s]
                                     !< HA: [Pa/s]
                                     !< dim: (nproma,nlevp1,nblks_c,ntracer)

    REAL(wp), INTENT(IN), OPTIONAL:: &  !< vertical tracer flux at upper boundary 
      &  opt_topflx_tra(:,:,:)          !< NH: [kg/m**2/s]
                                        !< HA: [Pa/s]
                                        !< dim: (nproma,nblks_c,ntracer)

    REAL(wp), INTENT(OUT), OPTIONAL :: & !< tracer value at upper boundary of child nest 
      &  opt_q_int(:,:,:)               !< NH: [kg/kg]
                                        !< HA: [kg/kg]
                                        !< dim: (nproma,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT), OPTIONAL :: & !< advective tendency    [kg/kg/s]
      &  opt_ddt_tracer_adv(:,:,:,:)     !< dim: (nproma,nlev,nblks_c,ntracer)


    REAL(wp), INTENT(IN) :: &           !< advective time step [s]
      &  p_dtime  

    INTEGER,  INTENT(INOUT) :: &        !< time step counter [1]
      &  k_step                         !< necessary for Marchuk Splitting

    REAL(wp), TARGET ::  &              !< NH: density weighted cell height at full levels
      & z_delp_mc1(nproma,p_patch%nlev,p_patch%nblks_c) 
                                        !< at first intermediate time step [kg/m**2]
                                        !< HA: pressure thickness for full levels 
                                        !< at first intermediate time step [Pa]
                                        !< dim: (nproma,nlev,nblks_c) 

    REAL(wp), TARGET ::  &              !< NH: density weighted cell height at full levels
      & z_delp_mc2(nproma,p_patch%nlev,p_patch%nblks_c) 
                                        !< at second intermediate time step [kg/m**2]
                                        !< HA: pressure thickness for full levels 
                                        !< at first intermediate time step [Pa]
                                        !< dim: (nproma,nlev,nblks_c) 

    REAL(vp) ::  &                      !< flux divergence at cell center
      &  z_fluxdiv_c(nproma,p_patch%nlev) 

    REAL(wp), POINTER   &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS        &
#endif
      & :: ptr_current_tracer(:,:,:,:)  !< pointer to tracer field

    REAL(wp) :: pdtime_mod        !< modified time step
                                  !< (multiplied by cSTR * coeff_grid)

    INTEGER  :: nlev              !< number of full and half levels
    INTEGER  :: jb, jk, jt, jc, jg, nt            !< loop indices
    INTEGER  :: ikp1                              !< vertical level + 1
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

    INTEGER, DIMENSION(:,:,:), POINTER :: &  !< Pointer to line and block indices (array)
      &  iidx, iblk                          !< of edges

    TYPE(t_trList), POINTER :: trAdvect      !< Pointer to tracer sublist
    TYPE(t_trList), POINTER :: trNotAdvect   !< Pointer to tracer sublist

    LOGICAL  :: is_present_opt_ddt_tracer_adv


   !-----------------------------------------------------------------------

    IF(ltimer) CALL timer_start(timer_transport)

    ! number of vertical levels
    nlev   = p_patch%nlev

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg  = p_patch%id

    ! precompute modified timestep (timestep multiplied by Strang-splitting
    ! coefficient and a second coefficient to account for either a
    ! vertical z- or p-system)
    pdtime_mod = advection_config(jg)%cSTR * advection_config(jg)%coeff_grid &
      &        * p_dtime

    is_present_opt_ddt_tracer_adv = PRESENT( opt_ddt_tracer_adv )

    ! line and block indices of edges as seen from cells
    iidx => p_patch%cells%edge_idx
    iblk => p_patch%cells%edge_blk

    ! tracer fields which are advected
    trAdvect => advection_config(jg)%trAdvect       ! 2018-06-05: cray bug do not add to PRESENT list
    !
    ! tracer fields which are not advected
    trNotAdvect => advection_config(jg)%trNotAdvect ! 2018-06-05: cray bug do not add to PRESENT list

    !---------------------------------------------------!
    !                                                   !
    !  time integration of tracer continuity-equation   !
    !                                                   !
    !---------------------------------------------------!

    !
    ! Godunov-type advection with two-time level scheme
    !

    ptr_current_tracer => p_tracer_now
    ptr_delp_mc_now    => p_delp_mc_now

!$ACC DATA  PCOPYIN( p_tracer_now, p_mflx_contra_h, p_mflx_contra_v,    &
!$ACC                p_vn_contra_traj, p_w_contra_traj,                 &
!$ACC                p_cellhgt_mc_now, p_delp_mc_now, p_delp_mc_new),   &
!$ACC       PCOPYOUT( p_tracer_new, p_mflx_tracer_h, p_mflx_tracer_v ), &
!$ACC       CREATE( z_delp_mc1, z_delp_mc2 ),              &
!$ACC       PRESENT( p_int_state, advection_config, iidx, iblk ),       &
!$ACC       IF( i_am_accel_node .AND. acc_on )
!$ACC UPDATE DEVICE( p_tracer_now, p_mflx_contra_h, p_mflx_contra_v,     & 
!$ACC                p_vn_contra_traj, p_w_contra_traj,                  &
!$ACC                p_cellhgt_mc_now, p_delp_mc_now, p_delp_mc_new),    &
!$ACC        IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    !*********************************!
    ! Vertical advection              !
    !*********************************!
    ! calculate vertical tracer flux
    ! subroutine vertical_upwind_flux will be called if either the
    ! current time step is even or complete Strang-splitting is chosen.
    !
    IF ( advection_config(jg)%lvadv_tracer ) THEN

      IF ( MOD( k_step, 2 ) == 0 .OR. advection_config(jg)%lstrang ) THEN
 
        ! vertical transport and subsequent calculations include all halo points in order
        ! to avoid an additional synchronization step
        ! nest boundary points are needed as well because the subsequent horizontal transport
        ! accesses part of them
        !
        i_rlstart  = 2
        i_rlend    = min_rlcell
        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


        ! calculation of intermediate layer thickness (density)
        ! necessary for tracer-mass consistency.
        ptr_delp_mc_new  => z_delp_mc1

        ! integration of tracer continuity equation in vertical
        ! direction. Must be computed prior to the vertical tracer flux,
        ! since it is required for FCT (to be implemented).
        ! This is independent of the tracer and thus must be
        ! computed only once.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = i_startidx, i_endidx
          ! integration of mass continuity equation
              ptr_delp_mc_new(jc,jk,jb) =                      &
            &              ptr_delp_mc_now(jc,jk,jb)       &
            &              - pdtime_mod                                         &
            &              * ( p_mflx_contra_v(jc,jk+1,jb) &
            &              - p_mflx_contra_v(jc,jk,jb) )
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( ptr_delp_mc_new ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

        ! compute vertical tracer flux
        CALL vert_upwind_flux( p_patch, ptr_current_tracer,          &! in
          &              p_mflx_contra_v,                            &! inout
          &              p_w_contra_traj,                            &! in
          &              advection_config(jg)%cSTR*p_dtime,          &! in
          &              p_cellhgt_mc_now,                           &! in
          &              ptr_delp_mc_now,                            &! in
          &              advection_config(jg)%ivadv_tracer,          &! in
          &              advection_config(jg)%itype_vlimit,          &! in
          &              advection_config(jg)%iubc_adv,              &! in
          &              advection_config(jg)%iadv_slev(:),          &! in
          &              .TRUE.,                                     &! print CFL number
          &              p_mflx_tracer_v,                            &! out
          &              opt_topflx_tra=opt_topflx_tra,              &! in
          &              opt_q_int=opt_q_int,                        &! out
          &              opt_rlstart=i_rlstart,                      &! in
          &              opt_rlend=i_rlend                           )! in


        ! compute vertical flux divergence for each tracer
        !

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jt,jc,nt,ikp1,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, i_rlstart, i_rlend)


!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG PRIVATE(jt)

          ! computation of vertical flux divergences
          DO nt = 1, trAdvect%len ! Tracer loop

            jt = trAdvect%list(nt)

            !$ACC LOOP WORKER PRIVATE(ikp1)
            DO jk = advection_config(jg)%iadv_slev(jt), nlev

              ! index of top half level
              ikp1 = jk + 1

!DIR$ IVDEP
              !$ACC LOOP VECTOR
              DO jc = i_startidx, i_endidx

                p_tracer_new(jc,jk,jb,jt) =                                         &
                  &  ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb)  &
                  &  - pdtime_mod * ( p_mflx_tracer_v(jc,ikp1,jb,jt)                &
                  &               -   p_mflx_tracer_v(jc,jk  ,jb,jt) ) )            &
                  &  / ptr_delp_mc_new(jc,jk,jb)

              END DO
            END DO
          END DO  ! Tracer loop
!$ACC END PARALLEL
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_tracer_new ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

        ptr_current_tracer => p_tracer_new
        ptr_delp_mc_now    => z_delp_mc1

        IF ( advection_config(jg)%lstrang ) THEN

          ! prepare new intermediate densities for horizontal advection,
          ! by integrating the continuity equation in horizontal direction.

          ! calculation of delp_mc2 (intermediate layer thickness, necessary for
          ! tracer-mass consistency)
          ptr_delp_mc_new   => z_delp_mc2

          ! integration of tracer continuity equation in vertical
          ! direction. This is independent of the tracer and thus must be
          ! computed only once.

          ! Note that intermediate densities are needed by the horizontal 
          ! limiter as well. That is, why we have to extend the following 
          ! computation to halo points.
          i_rlstart  = 1
          i_rlend    = min_rlcell
          i_startblk = p_patch%cells%start_blk(i_rlstart,1)
          i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk, i_endblk

            CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
              &                 i_startidx, i_endidx, i_rlstart, i_rlend )

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
            !$ACC LOOP GANG
            DO jk = 1, nlev
              !$ACC LOOP VECTOR
              DO jc = i_startidx, i_endidx
                ptr_delp_mc_new(jc,jk,jb) =                   &
              &          p_delp_mc_new(jc,jk,jb)          &
              &          + pdtime_mod                                          &
              &          * ( p_mflx_contra_v(jc,jk+1,jb)  &
              &          - p_mflx_contra_v(jc,jk,jb) )
              ENDDO
            ENDDO
!$ACC END PARALLEL
          ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( ptr_delp_mc_new ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

        ELSE !(no lstrang, but (lvadv_tracer .AND. ( MOD( k_step, 2 ) == 0))

          ! the following integration of the continuity equation in the
          ! horizontal direction should result in $\delp^{n+1}$
          ptr_delp_mc_new  => p_delp_mc_new
        ENDIF


      ELSE
        ! If vertical advection does not proceed horizontal advection and
        ! no complete strang splitting is used, then prepare new intermediate
        ! densities for horizontal advection, by integrating the continuity
        ! equation in horizontal direction.

        ! calculation of delp_mc2 (intermediate layer thickness, necessary for
        ! tracer-mass consistency)
        ptr_delp_mc_new   => z_delp_mc2

        ! integration of tracer continuity equation in vertical
        ! direction. This is independent of the tracer and thus must be
        ! computed only once.

        i_rlstart  = grf_bdywidth_c-1
        ! Note that delp is needed for the horizontal limiter
        i_rlend    = min_rlcell
        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )
 
          ! The intermediate cell mass is limited to 10% of the original value in order
          ! to prevent instabilities in extreme breaking gravity waves

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG 
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jc = i_startidx, i_endidx
              ptr_delp_mc_new(jc,jk,jb) =                        &
            &       MAX(0.1_wp*p_delp_mc_new(jc,jk,jb),      &
            &                  p_delp_mc_new(jc,jk,jb)       &
            &                + pdtime_mod                                         &
            &                * ( p_mflx_contra_v(jc,jk+1,jb) &
            &                - p_mflx_contra_v(jc,jk,jb) )   )
            ENDDO
          ENDDO
!$ACC END PARALLEL
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( ptr_delp_mc_new ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

      END IF

    ELSE  ! if lvadv_tracer=.FALSE.

      ptr_delp_mc_new  => p_delp_mc_new

    ENDIF


    !*************************************!
    ! Horizontal advection                !
    ! compute horizontal tracer flux      !
    !*************************************!
    !
    ! calculate horizontal tracer flux with upwind scheme
    !
    i_rlend        = min_rledge_int-1
    !
    CALL hor_upwind_flux( ptr_current_tracer,                                &! in
      &                  ptr_delp_mc_now,                                    &! in
      &                  p_mflx_contra_h, p_vn_contra_traj, p_dtime, p_patch,&! in
      &                  p_int_state, advection_config(jg)%ihadv_tracer,     &! in
      &                  advection_config(jg)%igrad_c_miura,                 &! in
      &                  advection_config(jg)%itype_hlimit,                &! in
      &                  advection_config(jg)%iadv_slev(:),                  &! in
      &                  advection_config(jg)%iord_backtraj,                 &! in
      &                  p_mflx_tracer_h, opt_rlend=i_rlend                  )! inout,in


!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)

    !
    ! update tracer array
    !
    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk  (i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jt,jc,nt,z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                    i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP GANG PRIVATE( z_fluxdiv_c, jt )
      DO nt = 1, trAdvect%len ! Tracer loop

        jt = trAdvect%list(nt)

        IF ( advection_config(jg)%ihadv_tracer(jt) /= 0 ) THEN

!  compute divergence of the upwind fluxes for tracers
         !$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = advection_config(jg)%iadv_slev(jt), nlev
#else
!CDIR UNROLL=6
          DO jk = advection_config(jg)%iadv_slev(jt), nlev
            DO jc = i_startidx, i_endidx
#endif

! TODO: possible GPU optimization: add p_tracer_new calculation here
              z_fluxdiv_c(jc,jk) =  &
                & p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt)*p_int_state%geofac_div(jc,1,jb) + &
                & p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt)*p_int_state%geofac_div(jc,2,jb) + &
                & p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt)*p_int_state%geofac_div(jc,3,jb)

            ENDDO
          ENDDO

        ELSE  ! horizontal advection switched off

!  compute divergence of the upwind fluxes for tracers
          !$ACC LOOP VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = advection_config(jg)%iadv_slev(jt), nlev
#else
!CDIR UNROLL=6
          DO jk = advection_config(jg)%iadv_slev(jt), nlev
            DO jc = i_startidx, i_endidx
#endif

! TODO: possible GPU optimization: add p_tracer_new calculation here
              z_fluxdiv_c(jc,jk) =  ptr_current_tracer(jc,jk,jb,jt) * ( &
                & p_mflx_contra_h(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_int_state%geofac_div(jc,1,jb) + &
                & p_mflx_contra_h(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_int_state%geofac_div(jc,2,jb) + &
                & p_mflx_contra_h(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_int_state%geofac_div(jc,3,jb) )

            ENDDO
          ENDDO

        ENDIF  ! ihadv_tracer(jt) /= 0



        !$ACC LOOP VECTOR COLLAPSE(2)
        DO jk = advection_config(jg)%iadv_slev(jt), nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx

            p_tracer_new(jc,jk,jb,jt) =                                         &
              &   ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb) &
              &    - p_dtime * z_fluxdiv_c(jc,jk) ) / ptr_delp_mc_new(jc,jk,jb)

          ENDDO
        ENDDO

        ! set tracer(nnew) to tracer(nnow) where advection is turned off
        !$ACC LOOP VECTOR COLLAPSE(2)
        DO jk = 1, advection_config(jg)%iadv_slev(jt)-1
          DO jc = i_startidx, i_endidx
            p_tracer_new(jc,jk,jb,jt) = p_tracer_now(jc,jk,jb,jt)
          END DO
        END DO


        ! Store qv advection tendency for convection scheme. 
        ! Store TKE tendency, if TKE advection is turned on
        !
        IF ( is_present_opt_ddt_tracer_adv .AND. (MOD( k_step, 2 ) == 0) &
          &  .AND. iforcing == inwp ) THEN
          IF ( jt == iqv ) THEN
            !$ACC LOOP VECTOR COLLAPSE(2)
            DO jk = advection_config(jg)%iadv_slev(jt), nlev
              DO jc = i_startidx, i_endidx
                opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                  & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime           
              ENDDO
            ENDDO
          ENDIF  ! jt == iqv
          IF ( advection_config(jg)%iadv_tke > 0 .AND. jt == iqtke ) THEN
            !$ACC LOOP VECTOR COLLAPSE(2)
            DO jk = advection_config(jg)%iadv_slev(jt), nlev
              DO jc = i_startidx, i_endidx
                opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                  & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime           
              ENDDO
            ENDDO
          ENDIF  ! jt == iqtke
        ENDIF

      ENDDO  ! Tracer loop
!$ACC END PARALLEL

    ENDDO

!$OMP END DO

    ! Update lateral boundaries of nested domains with interpolated time tendencies
    IF (l_limited_area .OR. p_patch%id > 1) THEN
      i_rlstart  = 1
      i_rlend    = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jt,jc,nt) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

          ! Tracer values are clipped here to avoid generation of negative values
          ! For mass conservation, a correction has to be applied in the
          ! feedback routine anyway
!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
       !$ACC LOOP GANG PRIVATE(jt)
        DO nt = 1, trAdvect%len ! Tracer loop

          jt = trAdvect%list(nt)

          !$ACC LOOP VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_tracer_new(jc,jk,jb,jt) =                            &
                &     MAX(0._wp, p_tracer_now(jc,jk,jb,jt)           &
                &   + p_dtime * p_grf_tend_tracer(jc,jk,jb,jt) )
            ENDDO
          ENDDO  ! Tracer loop
        ENDDO
!$ACC END PARALLEL
      ENDDO

!$OMP END DO NOWAIT

    ENDIF

!$ACC UPDATE HOST( p_tracer_new ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

!$OMP END PARALLEL


    !*********************************!
    ! vertical tracer advection       !
    !*********************************!
    ! subroutine vertical_upwind_flux will be called if either the
    ! current time step is odd or complete Strang-splitting is chosen.
    !
    IF ( advection_config(jg)%lvadv_tracer .AND.        &
      &  ( MOD( k_step, 2 ) == 1 .OR. advection_config(jg)%lstrang) ) THEN

      ptr_current_tracer => p_tracer_new
      ptr_delp_mc_now    => z_delp_mc2
      ptr_delp_mc_new    => p_delp_mc_new


!      IF (lcompute .and. lupdate_pres) THEN
        ! we should think about, whether we should change(update) the pressure values
        ! in the following CALL, in order to be consistent with the intermediate
        ! z_delp_mc2
!      ENDIF

      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int
      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

      CALL vert_upwind_flux( p_patch, ptr_current_tracer,          &! in
        &              p_mflx_contra_v,                            &! inout
        &              p_w_contra_traj,                            &! in
        &              advection_config(jg)%cSTR*p_dtime,          &! in
        &              p_cellhgt_mc_now,                           &! in
        &              ptr_delp_mc_now,                            &! in
        &              advection_config(jg)%ivadv_tracer,          &! in
        &              advection_config(jg)%itype_vlimit,          &! in
        &              advection_config(jg)%iubc_adv,              &! in
        &              advection_config(jg)%iadv_slev(:),          &! in
        &              .FALSE.,                                    &! do not print CFL number
        &              p_mflx_tracer_v,                            &! out
        &              opt_topflx_tra=opt_topflx_tra,              &! in
        &              opt_q_int=opt_q_int,                        &! out
        &              opt_rlstart=i_rlstart,                      &! in
        &              opt_rlend=i_rlend                           )! in


      ! calculate vertical flux divergence
      !
      ! update tracer array
      !

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jt,nt,i_startidx,i_endidx,ikp1) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG PRIVATE(jt)
        DO nt = 1, trAdvect%len ! Tracer loop

          jt = trAdvect%list(nt)

          !$ACC LOOP VECTOR, COLLAPSE(2)
          DO jk = advection_config(jg)%iadv_slev(jt), nlev

!DIR$ IVDEP
            DO jc = i_startidx, i_endidx
              ! index of top half level
              ikp1 = jk + 1   ! WS: put in here to ensure loops can be collapsed


              p_tracer_new(jc,jk,jb,jt) =                                         &
                &  ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb)  &
                &  - pdtime_mod * ( p_mflx_tracer_v(jc,ikp1,jb,jt)                &
                &               -   p_mflx_tracer_v(jc,jk  ,jb,jt) ) )            &
                &  / ptr_delp_mc_new(jc,jk,jb)

            END DO
          END DO

          ! Store qv advection tendency for convection scheme. 
          ! Store TKE tendency, if TKE advection is turned on
          !
          IF ( is_present_opt_ddt_tracer_adv .AND. (iforcing == inwp) ) THEN
            IF ( jt == iqv ) THEN
              !$ACC LOOP VECTOR, COLLAPSE(2)
              DO jk = advection_config(jg)%iadv_slev(jt), nlev
                DO jc = i_startidx, i_endidx
                  opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                    & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime           
                ENDDO
              ENDDO
            ENDIF  ! jt == iqv
            IF ( advection_config(jg)%iadv_tke > 0 .AND. jt == iqtke ) THEN
              !$ACC LOOP VECTOR, COLLAPSE(2)
              DO jk = advection_config(jg)%iadv_slev(jt), nlev
                DO jc = i_startidx, i_endidx
                  opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                    & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime           
                ENDDO
              ENDDO
            ENDIF  ! jt == iqtke
          ENDIF

        END DO  ! Tracer loop
!$ACC END PARALLEL
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC UPDATE HOST( p_tracer_new ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
      IF ( is_present_opt_ddt_tracer_adv .AND. (MOD( k_step, 2 ) == 0) ) THEN
!$ACC UPDATE HOST( opt_ddt_tracer_adv ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )
      ENDIF

    END IF


    ! For tracer fields which are not advected, we just 
    ! perform a copy from time level now to new.
    !
    IF ( trNotAdvect%len > 0 ) THEN

      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jt,nt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG PRIVATE(jt)
        DO nt = 1, trNotAdvect%len ! Tracer loop

          jt = trNotAdvect%list(nt)

          !$ACC LOOP VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_tracer_new(jc,jk,jb,jt) = p_tracer_now(jc,jk,jb,jt)
            ENDDO  !jc
          ENDDO  !jk
 
        ENDDO  !nt
!$ACC END PARALLEL
       
      ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF  ! trNotAdvect%len > 0


    ! Synchronize tracer array after update. This is only necessary, if 
    ! NWP physics has NOT been selected. Otherwise, the SYNC-operation will 
    ! follow AFTER the call of NWP physics.
    ! For efficiency, the synchronization is applied for all tracers at once

    IF (iforcing /= inwp) THEN
      CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=p_tracer_new)
    ENDIF

    !
    ! compute advective tracer tendency
    !

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)
    IF (is_present_opt_ddt_tracer_adv .AND. iforcing /= inwp) THEN

      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int
      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jt,jc,nt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG PRIVATE(jt)
        DO nt = 1, trAdvect%len

          jt = trAdvect%list(nt)

          !$ACC LOOP VECTOR COLLAPSE(2)
          DO jk = advection_config(jg)%iadv_slev(jt), nlev
            DO jc = i_startidx, i_endidx
              opt_ddt_tracer_adv(jc,jk,jb,jt) =               &
                &           (  p_tracer_new(jc,jk,jb,jt)      &
                &           -  p_tracer_now(jc,jk,jb,jt) )    &
                &            /p_dtime
            ENDDO
          ENDDO
        ENDDO
!$ACC END PARALLEL
      ENDDO
!$OMP END DO

!$ACC UPDATE HOST( opt_ddt_tracer_adv ), IF( acc_validate .AND. i_am_accel_node .AND. acc_on )

    ENDIF


    !
    ! eventually do a clipping of negative values to zero
    !
    IF ( advection_config(jg)%lclip_tracer ) THEN
      CALL negative2zero(p_tracer_new(:,:,:,:))
!$OMP BARRIER
    END IF

!$OMP END PARALLEL

    IF (is_present_opt_ddt_tracer_adv .AND. iforcing /= inwp) THEN
      CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer,  &
                                 & f4din=opt_ddt_tracer_adv )
    ENDIF

!$ACC END DATA

  IF (ltimer) CALL timer_stop(timer_transport)

  END SUBROUTINE step_advection


END MODULE mo_advection_stepping
