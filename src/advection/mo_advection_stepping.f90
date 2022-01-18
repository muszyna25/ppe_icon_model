!>
!! Tracer transport module
!! Contains solver for the tracer mass continuity equation(s)
!!
!! Performs time integration of tracer continuity equations in flux form. 
!! Strang splitting is applied between the horizontal and vertical direction.
!! The air mass continuity equation is re-integrated in the samme splitted manner 
!! in order to achieve consistency with continuity.
!! Nonzero right-hand-sides (i.e. physical tendencies) are not considered here. 
!! For NWP slow-physics tendencies are taken into account in the subroutine 
!! tracer_add_phytend which is called at the beginning of the main physics 
!! interface nwp_nh_interface.
!!
!! @author Daniel Reinert, DWD
!!
!! @Literature:
!! Reinert, D. (2020): A Mass Consistent Finite Volume Approach with Fractional Steps.
!!                     Reports on ICON, Issue 4 
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
!! Modification by Daniel Reinert, DWD (2019-02-18)
!! - removed obsolete 3D field p_w_contra_traj, which was only needed 
!!   by the time-restricted PPM scheme. This scheme is no longer available.
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
  USE mo_exception,           ONLY: finish
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_transport, &
    &                               timer_adv_vert, timer_adv_horz
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ntracer, lvert_nest, ltimer, timers_level, iforcing, iqv, iqtke
  USE mo_advection_hflux,     ONLY: hor_upwind_flux
  USE mo_advection_vflux,     ONLY: vert_upwind_flux
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rlcell, inwp
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array, sync_patch_array_mult
  USE mo_advection_config,    ONLY: advection_config, t_trList
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_initicon_config,     ONLY: is_iau_active, iau_wgt_adv
  USE mo_fortran_tools,       ONLY: negative2zero
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif
  USE mo_dynamics_config,     ONLY: ldeepatmo
  USE mo_upatmo_impl_const,   ONLY: idamtr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: step_advection

#if defined( _OPENACC )
#if defined(__ADVECTION_STEPPING_NOACC)
  LOGICAL, PARAMETER ::  acc_on = .FALSE.
#else
  LOGICAL, PARAMETER ::  acc_on = .TRUE.
#endif
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
    &                        p_cellhgt_mc_now, p_rhodz_new, p_rhodz_now,          &
    &                        p_grf_tend_tracer, p_tracer_new,                     &
    &                        p_mflx_tracer_h, p_mflx_tracer_v, rho_incr,          &
    &                        q_ubc, q_int,                                        &
    &                        opt_ddt_tracer_adv, opt_deepatmo_t1mc,               &
    &                        opt_deepatmo_t2mc                 )
  !
    TYPE(t_patch), TARGET, INTENT(INOUT) ::  &  !< patch on which computation
      &  p_patch                             !< is performed
                                             
    TYPE(t_int_state), INTENT(IN) :: &  !< interpolation state
      &  p_int_state

    REAL(wp), INTENT(IN) ::          &  !< advective time step [s]
      &  p_dtime  

    INTEGER,  INTENT(IN) ::          &  !< time step counter [1]
      &  k_step                         !< necessary for Strang Splitting

    REAL(wp) &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     , INTENT(IN) ::                 &  !< tracer mixing ratios (specific concentrations)
      &  p_tracer_now(:,:,:,:)          !< at current time level n (before transport)
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(IN)  ::         &  !< horizontal mass flux (contravariant)
      &  p_mflx_contra_h(:,:,:)         !< NH: v_n*delta_z*\rho  [kg/m/s]
                                        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(IN)  ::         &  !< horizontal velocity component at n+1/2
      &  p_vn_contra_traj(:,:,:)        !< for calculation of backward trajectories
                                        !< [m/s] 
                                        !< dim: (nproma,nlev,nblks_e)

    REAL(wp), INTENT(INOUT)  ::      &  !< vertical mass flux (contravariant)
      &  p_mflx_contra_v(:,:,:)         !< NH: \rho*w     [kg/m**2/s]
                                        !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::          &  !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)        !< time step n 
                                        !< NH: \Delta z       [m]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::          &  !< NH: density weighted cell height at full levels 
      &  p_rhodz_new(:,:,:)             !< at n+1 [kg/m**2]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN), TARGET ::  &  !< NH: density weighted cell height at full levels
      &  p_rhodz_now(:,:,:)             !< at time step n [kg/m**2]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::          &  !< interpolated tracer time tendencies for
      &  p_grf_tend_tracer(:,:,:,:)     !< updating the lateral boundaries of nested domains
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp) &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
     , INTENT(INOUT) ::          &      !< tracer mixing ratios (specific concentrations)
      &  p_tracer_new(:,:,:,:)          !< at time level n+1 (after transport)
                                        !< [kg/kg]  
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT)  ::  &      !< horizontal tracer mass flux at full level edges
      &  p_mflx_tracer_h(:,:,:,:)       !< NH: [kg/m/s]
                                        !< dim: (nproma,nlev,nblks_e,ntracer)

    REAL(wp), INTENT(INOUT)  ::  &      !< vertical tracer mass flux at half level centers
      &  p_mflx_tracer_v(:,:,:,:)       !< NH: [kg/m**2/s]
                                        !< dim: (nproma,nlevp1,nblks_c,ntracer)

    REAL(vp), INTENT(IN), OPTIONAL:: &  !< density increment due to IAU
      &  rho_incr(:,:,:)                !< [kg/m**3]
                                        ! Note that the OPTIONAL argument is only necessary 
                                        ! due to co-use of step_advection by the hydrostatic 
                                        ! model configuration. Can be removed once the hydrostatic 
                                        ! model configuration is gone.

    REAL(wp), INTENT(INOUT) :: &        !< tracer mass fraction at (nest) upper boundary 
      &  q_ubc(:,:,:)                   !< NH: [kg/kg]
                                        !< dim: (nproma,ntracer,nblks_c)

    REAL(wp), INTENT(OUT)   :: &        !< tracer mass fraction at child nest interface level 
      &  q_int(:,:,:)                   !< NH: [kg/kg]
                                        !< dim: (nproma,ntracer,nblks_c)

    REAL(wp), INTENT(INOUT), OPTIONAL :: & !< advective tendency    [kg/kg/s]
      &  opt_ddt_tracer_adv(:,:,:,:)     !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(IN), OPTIONAL ::  &   !< deep-atmosphere modification factors,
      &  opt_deepatmo_t1mc(:,:),       &   !< type 1: (nlev,nitem1)
      &  opt_deepatmo_t2mc(:,:)            !< type 2: (nitem2,nlev)



    ! Local Variables

    REAL(wp), CONTIGUOUS, POINTER ::  & !< intermediate density times cell thickness [ kg/m**2]
      &  rhodz_ast(:,:,:) => NULL()     !< includes any nonzero source/sink term

    REAL(wp),  TARGET ::  &             !< auxiliary field for rhodz_ast
      &  rhodz_aux(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp) ::  &                      !< intermediate density times cell thickness [ kg/m**2]
      &  rhodz_ast2(nproma,p_patch%nlev,p_patch%nblks_c)
                                        !< compared to rhodz_ast it additionally includes either
                                        !< the horizontal or vertical advective density increment. 

    INTEGER  :: nlev                    !< number of full levels
    INTEGER  :: jb, jk, jt, jc, jg, nt            !< loop indices
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: iadv_slev_jt                      ! Workaround OpenACC limitation

    TYPE(t_trList), POINTER :: trAdvect      !< Pointer to tracer sublist
    TYPE(t_trList), POINTER :: trNotAdvect   !< Pointer to tracer sublist

    REAL(wp) ::   &                     !< Local copies of opt_deepatmo_t1mc and opt_deepatmo_t2mc
      &  deepatmo_vol  (p_patch%nlev), &
      &  deepatmo_divh (p_patch%nlev), &
      &  deepatmo_divzu(p_patch%nlev), &
      &  deepatmo_divzl(p_patch%nlev)

    LOGICAL :: is_present_rho_incr      !< as it says
    LOGICAL :: lstep_even               !< TRUE/FALSE : this is an even/odd timestep

   !-----------------------------------------------------------------------

    IF(ltimer) CALL timer_start(timer_transport)

    ! number of vertical levels
    nlev = p_patch%nlev

    jg = p_patch%id

    ! determine whether this is an even or odd timestep
    lstep_even = MOD( k_step, 2 ) == 0

    is_present_rho_incr = PRESENT( rho_incr )

    ! tracer fields which are advected
    trAdvect => advection_config(jg)%trAdvect       ! 2018-06-05: cray bug do not add to PRESENT list
    !
    ! tracer fields which are not advected
    trNotAdvect => advection_config(jg)%trNotAdvect ! 2018-06-05: cray bug do not add to PRESENT list


    ! deep-atmosphere modification factors: 
    ! direct access of opt_deepatmo_t1mc, opt_deepatmo_t2mc not possible, 
    ! due to co-use of step_advection by hydrostatic model configuration
    IF ( .NOT. (PRESENT(opt_deepatmo_t1mc) .AND. PRESENT(opt_deepatmo_t2mc)) .AND. &
      &  ldeepatmo                                                                 ) THEN
      CALL finish('mo_advection_stepping: step_advection', 'Missing deepatmo argument(s).')
    ELSEIF ( PRESENT(opt_deepatmo_t1mc) .AND. PRESENT(opt_deepatmo_t2mc) .AND.     &
      &      ldeepatmo                                                             ) THEN
      deepatmo_vol(:)   = opt_deepatmo_t1mc(:,idamtr%t1mc%vol)
      deepatmo_divh(:)  = opt_deepatmo_t1mc(:,idamtr%t1mc%divh)
      deepatmo_divzu(:) = opt_deepatmo_t2mc(idamtr%t2mc%divzU,:)
      deepatmo_divzl(:) = opt_deepatmo_t2mc(idamtr%t2mc%divzL,:)
    ELSE
      deepatmo_vol(:)   = 1._wp
      deepatmo_divh(:)  = 1._wp
      deepatmo_divzu(:) = 1._wp
      deepatmo_divzl(:) = 1._wp
    ENDIF


!$ACC DATA  PCOPYIN( p_tracer_now, p_mflx_contra_h, p_mflx_contra_v,    &
!$ACC                p_vn_contra_traj,                                  &
!$ACC                p_cellhgt_mc_now, p_rhodz_now, p_rhodz_new,        &
!$ACC                deepatmo_vol, deepatmo_divh, deepatmo_divzu,       &
!$ACC                deepatmo_divzl),    &
!$ACC       PCOPYOUT( p_tracer_new, p_mflx_tracer_h, p_mflx_tracer_v ), &
!$ACC       CREATE( rhodz_aux, rhodz_ast2 ),                            &
!$ACC       PRESENT( p_int_state, p_grf_tend_tracer, q_ubc),            &
!$ACC       IF( i_am_accel_node .AND. acc_on )

!$ACC DATA  PCOPYIN( rho_incr ),                                        &
!$ACC       IF( i_am_accel_node .AND. acc_on .AND. is_present_rho_incr )
!$ACC DATA  PRESENT( opt_ddt_tracer_adv ),                              &
!$ACC       IF( i_am_accel_node .AND. acc_on .AND. PRESENT(opt_ddt_tracer_adv) )


    ! This vertical mass flux synchronization is necessary, as vertical transport 
    ! includes all halo points (see below)
    IF (lvert_nest .AND. p_patch%nshift > 0) THEN ! vertical nesting

      CALL sync_patch_array_mult(SYNC_C, p_patch, 2, p_mflx_contra_v, q_ubc, &
                                 opt_varname = 'step_advecvtion: p_mflx_contra_v,q_ubc')
    ELSE
      ! note that in cases without vertical nesting q_ubc=0._wp is ensured, as 
      ! * q_int = 0._wp
      ! * parent to child interpolation of q_int is constancy preserving
      !
      CALL sync_patch_array(SYNC_C, p_patch, p_mflx_contra_v, opt_varname='step_advection: p_mflx_contra_v')
    ENDIF


    ! In order to achieve consistency with continuity we follow the method of 
    ! Easter (1993) and (re-)integrate the air mass continuity equation in the 
    ! same split manner as the tracer mass continuity equation.
    !
    ! We start by accounting for any RHS. 
    ! Currently, the only nonzero RHS occurs during the IAU phase.

    ! halo points must be included
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    IF (is_present_rho_incr .AND. is_iau_active) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            rhodz_aux(jc,jk,jb) = p_rhodz_now(jc,jk,jb)                    &
              &                 + iau_wgt_adv*p_cellhgt_mc_now(jc,jk,jb)   &
              &                 * rho_incr(jc,jk,jb) * deepatmo_vol(jk)
          ENDDO  !jc
        ENDDO  !jk
!$ACC END PARALLEL

      ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      rhodz_ast => rhodz_aux
    ELSE
      rhodz_ast => p_rhodz_now
    ENDIF

    ! Integrate tracer continuity equation.
    ! Separate treatment of horizontal and vertical directions via Strang splitting.
    ! Hence the order of the horizontal and vertical transport operator 
    ! is reversed every second time step.
    !
    IF (lstep_even) THEN
      !
      ! Vertical transport precedes horizontal transport
      !

      ! vertical transport includes all halo points in order to avoid 
      ! an additional synchronization step.
      ! nest boundary points are needed as well because the subsequent horizontal transport
      ! accesses part of them.
      !
      i_rlstart  = 2
      i_rlend    = min_rlcell
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend)

        ! compute intermediate density which accounts for the density increment 
        ! due to vertical transport.
!$ACC PARALLEL DEFAULT(NONE) PRESENT(rhodz_ast) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
       !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            rhodz_ast2(jc,jk,jb) = rhodz_ast(jc,jk,jb) + p_dtime                      &
                                 * ( p_mflx_contra_v(jc,jk+1,jb) * deepatmo_divzl(jk) &
                                 -   p_mflx_contra_v(jc,jk,jb)   * deepatmo_divzu(jk) )
          ENDDO  !jc
        ENDDO !jk
!$ACC END PARALLEL
      ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


      ! vertical transport
      !
      CALL vert_adv(p_patch          = p_patch,                 & !in
        &           p_dtime          = p_dtime,                 & !in
        &           k_step           = k_step,                  & !in
        &           p_mflx_contra_v  = p_mflx_contra_v(:,:,:),  & !inout
        &           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,:), & !in
        &           rhodz_now        = rhodz_ast(:,:,:),        & !in
        &           rhodz_new        = rhodz_ast2(:,:,:),       & !in
        &           tracer_now       = p_tracer_now(:,:,:,:),   & !in
        &           tracer_new       = p_tracer_new(:,:,:,:),   & !inout
        &           p_mflx_tracer_v  = p_mflx_tracer_v(:,:,:,:),& !inout
        &           deepatmo_divzl   = deepatmo_divzl(:),       & !in
        &           deepatmo_divzu   = deepatmo_divzu(:),       & !in
        &           i_rlstart        = i_rlstart,               & !in
        &           i_rlend          = i_rlend,                 & !in
        &           q_ubc            = q_ubc(:,:,:),            & !in
        &           q_int            = q_int(:,:,:)             ) !out


      ! horizontal transport
      !
      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int
      !
      CALL hor_adv(p_patch         = p_patch,                  & !in
        &          p_int_state     = p_int_state,              & !in
        &          p_dtime         = p_dtime,                  & !in
        &          p_mflx_contra_h = p_mflx_contra_h(:,:,:),   & !in
        &          p_vn            = p_vn_contra_traj(:,:,:),  & !in 
        &          rhodz_now       = rhodz_ast2(:,:,:),        & !in
        &          rhodz_new       = p_rhodz_new(:,:,:),       & !in
        &          tracer_now      = p_tracer_new(:,:,:,:),    & !in
        &          tracer_new      = p_tracer_new(:,:,:,:),    & !inout
        &          p_mflx_tracer_h = p_mflx_tracer_h(:,:,:,:), & !inout
        &          deepatmo_divh   = deepatmo_divh(:),         & !in
        &          i_rlstart       = i_rlstart,                & !in
        &          i_rlend         = i_rlend                   ) !in

    ELSE  ! odd timestep
      !
      ! Horizontal transport precedes vertical transport
      !
      i_rlstart  = grf_bdywidth_c-1
      i_rlend    = min_rlcell
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend)


        ! compute intermediate density which accounts for the density increment 
        ! due to horizontal transport.
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx
            ! here we apply a security measure in order to ensure that 
            ! the cell is not fully emptied during horizontal transport.
            rhodz_ast2(jc,jk,jb) = MAX(0.1_wp*p_rhodz_new(jc,jk,jb),                    &
                                     p_rhodz_new(jc,jk,jb) - p_dtime                    &
                                 * ( p_mflx_contra_v(jc,jk+1,jb) * deepatmo_divzl(jk)   &
                                 -   p_mflx_contra_v(jc,jk,jb)   * deepatmo_divzu(jk) ) )
          ENDDO  !jc
        ENDDO !jk
!$ACC END PARALLEL
      ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


      ! horizontal transport
      !
      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int
      !
      CALL hor_adv(p_patch         = p_patch,                  & !in
        &          p_int_state     = p_int_state,              & !in
        &          p_dtime         = p_dtime,                  & !in
        &          p_mflx_contra_h = p_mflx_contra_h(:,:,:),   & !in
        &          p_vn            = p_vn_contra_traj(:,:,:),  & !in 
        &          rhodz_now       = rhodz_ast(:,:,:),         & !in
        &          rhodz_new       = rhodz_ast2(:,:,:),        & !in
        &          tracer_now      = p_tracer_now(:,:,:,:),    & !in
        &          tracer_new      = p_tracer_new(:,:,:,:),    & !inout
        &          p_mflx_tracer_h = p_mflx_tracer_h(:,:,:,:), & !inout
        &          deepatmo_divh   = deepatmo_divh(:),         & !in
        &          i_rlstart       = i_rlstart,                & !in
        &          i_rlend         = i_rlend                   ) !in


      ! vertical transport
      !
      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int
      !
      CALL vert_adv(p_patch          = p_patch,                 & !in
        &           p_dtime          = p_dtime,                 & !in
        &           k_step           = k_step,                  & !in
        &           p_mflx_contra_v  = p_mflx_contra_v(:,:,:),  & !inout
        &           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,:), & !in
        &           rhodz_now        = rhodz_ast2(:,:,:),       & !in
        &           rhodz_new        = p_rhodz_new(:,:,:),      & !in
        &           tracer_now       = p_tracer_new(:,:,:,:),   & !in
        &           tracer_new       = p_tracer_new(:,:,:,:),   & !inout
        &           p_mflx_tracer_v  = p_mflx_tracer_v(:,:,:,:),& !inout
        &           deepatmo_divzl   = deepatmo_divzl(:),       & !in
        &           deepatmo_divzu   = deepatmo_divzu(:),       & !in
        &           i_rlstart        = i_rlstart,               & !in
        &           i_rlend          = i_rlend,                 & !in
        &           q_ubc            = q_ubc(:,:,:),            & !in
        &           q_int            = q_int(:,:,:)             ) !out

    ENDIF  ! lstep_even
 


    ! For tracer fields which are not advected (neither horizontally nor vertically), 
    ! we perform a copy from time level now to new.
    !
    IF ( trNotAdvect%len > 0 ) THEN

      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,jt,nt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL DEFAULT(NONE) PRESENT(trNotAdvect) &
!$ACC          ASYNC(1) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP SEQ
        DO nt = 1, trNotAdvect%len ! Tracer loop

          jt = trNotAdvect%list(nt)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
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



    ! Update lateral boundaries of nested domains with interpolated time tendencies
    !
    IF (l_limited_area .OR. p_patch%id > 1) THEN
      i_rlstart  = 1
      i_rlend    = grf_bdywidth_c
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,jt,nt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        ! Tracer values are clipped here to avoid generation of negative values
        ! For mass conservation, a correction has to be applied in the
        ! feedback routine anyway

!$ACC PARALLEL DEFAULT(NONE) PRESENT(trAdvect) &
!$ACC          ASYNC(1) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP SEQ
        DO nt = 1, trAdvect%len ! Tracer loop

          jt = trAdvect%list(nt)

          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              p_tracer_new(jc,jk,jb,jt) =                            &
                &     MAX(0._wp, p_tracer_now(jc,jk,jb,jt)           &
                &   + p_dtime * p_grf_tend_tracer(jc,jk,jb,jt) )
            ENDDO
          ENDDO  !jk
        ENDDO  !Tracer loop
!$ACC END PARALLEL
      ENDDO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF


    ! Synchronize tracer array after update. This is only necessary, if 
    ! the NWP physics package is NOT used. Otherwise, the SYNC-operation will 
    ! follow AFTER the call of NWP physics.
    ! For efficiency, the synchronization is applied for all tracers at once.

    IF (iforcing /= inwp) THEN
      ! This wait is mandatory because of the communication
      !$ACC WAIT
      CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=p_tracer_new, opt_varname='ntracer and p_tracer_new')
    ENDIF


    !
    ! store advective tracer tendencies
    !
    ! If NWP physics are used, advective tendencies are only stored 
    ! for water vapour qv and turbulent kinetic energy TKE.
    ! If any other or no physics package is used, advective tendencies 
    ! are stored for all advected tracers.
    IF ( PRESENT(opt_ddt_tracer_adv) ) THEN

      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int
      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,jt,nt,iadv_slev_jt,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

!$ACC PARALLEL DEFAULT(NONE) PRESENT(trAdvect, advection_config ) &
!$ACC          ASYNC(1) IF( i_am_accel_node .AND. acc_on )
        !$ACC LOOP SEQ
        DO nt = 1, trAdvect%len ! Tracer loop

          jt = trAdvect%list(nt)
          iadv_slev_jt = advection_config(jg)%iadv_slev(jt)

          ! Store qv advection tendency for convection scheme.
          ! Store TKE tendency, if TKE advection is turned on 
          !
          IF ( iforcing == inwp ) THEN

            IF ( jt == iqv ) THEN
              !$ACC LOOP GANG VECTOR COLLAPSE(2)
              DO jk = iadv_slev_jt, nlev
                DO jc = i_startidx, i_endidx
                  opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                    & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime           
                ENDDO
              ENDDO
            ENDIF  ! jt == iqv

            IF ( advection_config(jg)%iadv_tke > 0 .AND. jt == iqtke ) THEN
              !$ACC LOOP GANG VECTOR COLLAPSE(2)
              DO jk = iadv_slev_jt, nlev
                DO jc = i_startidx, i_endidx
                  opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                    & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime           
                ENDDO
              ENDDO
            ENDIF  ! jt == iqtke

          ELSE  ! iforcing /= inwp

            ! store advection tendency for all tracers
            !
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO jk = iadv_slev_jt, nlev
              DO jc = i_startidx, i_endidx
                opt_ddt_tracer_adv(jc,jk,jb,jt) =                               &
                  & (p_tracer_new(jc,jk,jb,jt)-p_tracer_now(jc,jk,jb,jt))/p_dtime   
              ENDDO
            ENDDO
          ENDIF ! iforcing == inwp

        END DO  ! Tracer loop
!$ACC END PARALLEL
      END DO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      IF ( iforcing /= inwp ) THEN
        !$ACC WAIT
        CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer,  f4din=opt_ddt_tracer_adv, &
               &                   opt_varname='ntracer and opt_ddt_tracer_adv' )
      ENDIF

    ENDIF  ! PRESENT(opt_ddt_tracer_adv)


    !
    ! eventually do a clipping of negative values to zero
    !
    IF ( advection_config(jg)%lclip_tracer ) THEN
!$OMP PARALLEL
      CALL negative2zero(p_tracer_new(:,:,:,:), .TRUE.)
!$OMP BARRIER
!$OMP END PARALLEL
    END IF


!$ACC WAIT

!$ACC END DATA
!$ACC END DATA
!$ACC END DATA

    IF (ltimer) CALL timer_stop(timer_transport)

  END SUBROUTINE step_advection




  !>
  !! Vertical Transport
  !!
  !! Computes vertical fluxes based on the prescribed mass flux 
  !! and the current tracer mass fraction (tracer_now). 
  !! Finally the vertical mass flux divergence is computed 
  !! and the tracer mass fraction is updated (tracer_new).
  !!
  !! If vertical transport is switched off, the tracer mass fraction is 
  !! kept constant, i.e. tracer_now is copied to tracer_new.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2020-06-09)
  !!
  SUBROUTINE vert_adv (p_patch, p_dtime, k_step, p_mflx_contra_v,          &
    &                  p_cellhgt_mc_now, rhodz_now, rhodz_new, tracer_now, &
    &                  tracer_new, p_mflx_tracer_v, deepatmo_divzl,        &
    &                  deepatmo_divzu, i_rlstart, i_rlend, q_ubc,          &
    &                  q_int)

    TYPE(t_patch), INTENT(IN   )   ::  & !< compute patch
      &  p_patch

    REAL(wp),      INTENT(IN   )   ::  & !< advective time step
      &  p_dtime                         !< [s]

    INTEGER,       INTENT(IN   )   ::  & !< timestep counter
      &  k_step

    REAL(wp),      INTENT(INOUT)   ::  & !< vertical mass flux
      &  p_mflx_contra_v(:,:,:)          !< [kg/m**2/s]

    REAL(wp),      INTENT(IN   )   ::  & !< cell height defined at full levels
      &  p_cellhgt_mc_now(:,:,:)         !< [m]

    REAL(wp),      INTENT(IN   )   ::  & !< density times cell thickness (current value)
      &  rhodz_now(:,:,:)                !< [kg/m**2]

    REAL(wp),      INTENT(IN   )   ::  & !< density times cell thickness (updated value)
      &  rhodz_new(:,:,:)                !< additionally includes increment due to vertical transport
                                         !< [kg/m**2]

    REAL(wp),      INTENT(IN   )   ::  & !< tracer mass fraction (current value)
      &  tracer_now(:,:,:,:)             !< [kg/kg]

    REAL(wp),      INTENT(INOUT)   ::  & !< tracer mass fraction (updated value)
      &  tracer_new(:,:,:,:)             !< [kg/kg]

    REAL(wp),      INTENT(INOUT)   ::  & !< vertical tracer mass flux at half level centers
      &  p_mflx_tracer_v(:,:,:,:)        !< [kg/m**2/s]

    REAL(wp),      INTENT(IN   )   ::  & !< deep-atmosphere modification factors
      &  deepatmo_divzl(:), deepatmo_divzu(:)

    INTEGER,       INTENT(IN   )   ::  &
      &  i_rlstart, i_rlend

    REAL(wp), INTENT(IN)           ::  & !< tracer mass fraction at (nest) upper boundary 
      &  q_ubc(:,:,:)                    !< NH: [kg/kg]

    REAL(wp), INTENT(OUT)          ::  & !< tracer mass fraction at child nest interface level 
      &  q_int(:,:,:)                    !< NH: [kg/kg]

    ! local variables
    !
    INTEGER :: jc, jk, jb, nt         ! loop indices
    INTEGER :: jt                     ! tracer index
    INTEGER :: jg                     ! domain ID
    INTEGER :: ikp1                   ! jk+1
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: iadv_slev_jt           ! tracer dependent vertical start level

    TYPE(t_trList), POINTER :: trAdvect   !< Pointer to tracer sublist

    LOGICAL :: lprint_cfl             ! TRUE: compute and print vertical CFL number
   !-----------------------------------------------------------------------

    IF (timers_level > 2) CALL timer_start(timer_adv_vert)

    jg = p_patch%id

    ! compute and print vertical CFL number only every second time step, 
    ! as the computational overhead is considerable. 
    lprint_cfl = MOD( k_step, 2 ) == 0

    ! compute vertical tracer flux
    CALL vert_upwind_flux(                                                         &
      &              p_patch             = p_patch,                                & !in
      &              p_cc                = tracer_now(:,:,:,:),                    & !in
      &              p_mflx_contra_v     = p_mflx_contra_v(:,:,:),                 & !inout
      &              p_dtime             = p_dtime,                                & !in
      &              p_cellhgt_mc_now    = p_cellhgt_mc_now(:,:,:),                & !in
      &              p_cellmass_now      = rhodz_now(:,:,:),                       & !in
      &              lprint_cfl          = lprint_cfl,                             & !in
      &              p_upflux            = p_mflx_tracer_v(:,:,:,:),               & !out
      &              q_ubc               = q_ubc(:,:,:),                           & !in
      &              q_int               = q_int(:,:,:),                           & !out
      &              opt_rlstart         = i_rlstart,                              & !in
      &              opt_rlend           = i_rlend                                 ) !in


    ! tracer fields which are advected
    trAdvect => advection_config(jg)%trAdvect       ! 2018-06-05: cray bug do not add to PRESENT list

!$ACC DATA  PRESENT(p_mflx_tracer_v, tracer_now, tracer_new, rhodz_now, rhodz_new, &
!$ACC               deepatmo_divzl, deepatmo_divzu ),                              &
!$ACC       IF( i_am_accel_node .AND. acc_on )

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,jt,nt,iadv_slev_jt,i_startidx,i_endidx,ikp1) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,    &
                     i_startidx, i_endidx, i_rlstart, i_rlend)


      ! compute vertical flux divergences and update tracer mass fractions
      !
!$ACC PARALLEL DEFAULT(NONE) PRESENT( trAdvect, advection_config ) &
!$ACC          ASYNC(1) IF( i_am_accel_node .AND. acc_on )
      !$ACC LOOP SEQ
      TRACERLOOP: DO nt = 1, trAdvect%len ! Tracer loop

        jt = trAdvect%list(nt)
        iadv_slev_jt = advection_config(jg)%iadv_slev(jt)

        IF ( advection_config(jg)%ivadv_tracer(jt) /= 0 ) THEN
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = iadv_slev_jt, p_patch%nlev
            DO jc = i_startidx, i_endidx

              ! index of top half level
              ikp1 = jk + 1

              tracer_new(jc,jk,jb,jt) = ( tracer_now(jc,jk,jb,jt) * rhodz_now(jc,jk,jb)   &
                &    + p_dtime * ( p_mflx_tracer_v(jc,ikp1,jb,jt) * deepatmo_divzl(jk)    &
                &              -   p_mflx_tracer_v(jc,jk  ,jb,jt) * deepatmo_divzu(jk) ) )&
                &    / rhodz_new(jc,jk,jb)

            END DO  !jc
          END DO  !jk

          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          ! set tracer(nnew) to tracer(nnow) at levels where advection is turned off
          DO jk = 1, iadv_slev_jt-1
            DO jc = i_startidx, i_endidx
              tracer_new(jc,jk,jb,jt) = tracer_now(jc,jk,jb,jt)
            END DO
          END DO

        ELSE  ! no vertical transport

          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          ! copy
          DO jk = 1, p_patch%nlev
            DO jc = i_startidx, i_endidx
              tracer_new(jc,jk,jb,jt) = tracer_now(jc,jk,jb,jt)
            ENDDO  !jc
          ENDDO  !jk
        ENDIF

      END DO  TRACERLOOP
!$ACC END PARALLEL

    END DO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

    IF (timers_level > 2) CALL timer_stop(timer_adv_vert)

  END SUBROUTINE vert_adv




  !>
  !! Horizontal Transport
  !!
  !! Computes horizontal tracer mass fluxes based on the prescribed 
  !! air mass flux and the current tracer mass fraction (tracer_now). 
  !! Finally the horizontal mass flux divergence is computed and 
  !! the tracer mass fraction is updated (tracer_new).
  !!
  !! If horizontal transport is switched off, the tracer mass fraction is 
  !! kept constant, i.e. tracer_now is copied to tracer_new.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2020-06-09)
  !!
  SUBROUTINE hor_adv (p_patch, p_int_state, p_dtime, p_mflx_contra_h, p_vn, &
    &                 rhodz_now, rhodz_new, tracer_now, tracer_new,         &
    &                 p_mflx_tracer_h, deepatmo_divh, i_rlstart, i_rlend    )

    TYPE(t_patch), TARGET, INTENT(INOUT)  ::  & !< compute patch
      &  p_patch

    TYPE(t_int_state), INTENT(IN)  ::  & !< interpolation state 
      &  p_int_state

    REAL(wp),      INTENT(IN   )   ::  & !< advective time step
      &  p_dtime                         !< [s]

    REAL(wp),      INTENT(IN   )   ::  & !< horizontal mass flux
      &  p_mflx_contra_h(:,:,:)          !< [kg/m**2/s]

    REAL(wp), INTENT(IN)           ::  & !< horizontal velocity component at n+1/2
      &  p_vn(:,:,:)                     !< for calculation of backward trajectories
                                         !< [m/s] 

    REAL(wp),      INTENT(IN   )   ::  & !< density times cell thickness (current value)
      &  rhodz_now(:,:,:)                !< [kg/m**2]

    REAL(wp),      INTENT(IN   )   ::  & !< density times cell thickness (updated value)
      &  rhodz_new(:,:,:)                !< additionally includes increment due to vertical transport
                                         !< [kg/m**2]

    REAL(wp),      INTENT(IN   )   ::  & !< tracer mass fraction (current value)
      &  tracer_now(:,:,:,:)             !< [kg/kg]

    REAL(wp),      INTENT(INOUT)   ::  & !< tracer mass fraction (updated value)
      &  tracer_new(:,:,:,:)             !< [kg/kg]

    REAL(wp),      INTENT(INOUT)   ::  & !< horizontal tracer mass flux at cell edge midpoints 
      &  p_mflx_tracer_h(:,:,:,:)        !< [kg/m**2/s]

    REAL(wp),      INTENT(IN   )   ::  & !< deep-atmosphere modification factors
      &  deepatmo_divh(:)

    INTEGER,       INTENT(IN   )   ::  &
      &  i_rlstart, i_rlend


    ! local variables
    !
    INTEGER :: jc, jk, jb, nt             ! loop indices
    INTEGER :: jt                         ! tracer index
    INTEGER :: jg                         ! domain ID
    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: iadv_slev_jt               ! tracer dependent vertical start level
    INTEGER  :: nlev                      !< number of full levels   

    REAL(vp) ::  &                        !< flux divergence at cell center
      &  z_fluxdiv_c(nproma,p_patch%nlev) 

    TYPE(t_trList), POINTER :: trAdvect   !< Pointer to tracer sublist

    INTEGER, CONTIGUOUS, POINTER ::   &   !< Pointer to line and block indices (array)
      &  iidx(:,:,:) => NULL(),       &   !< of edges
      &  iblk(:,:,:) => NULL()
   !-----------------------------------------------------------------------

    IF (timers_level > 2) CALL timer_start(timer_adv_horz)

    jg = p_patch%id

    nlev = p_patch%nlev

    ! compute horizontal tracer flux
    CALL hor_upwind_flux(                                            &
      &              p_cc                = tracer_now(:,:,:,:),      & !in
      &              p_rhodz_now         = rhodz_now(:,:,:),         & !in
      &              p_rhodz_new         = rhodz_new(:,:,:),         & !in
      &              p_mass_flx_e        = p_mflx_contra_h(:,:,:),   & !in
      &              p_vn                = p_vn(:,:,:),              & !in
      &              p_dtime             = p_dtime,                  & !in
      &              p_patch             = p_patch,                  & !in
      &              p_int               = p_int_state,              & !in
      &              p_upflux            = p_mflx_tracer_h(:,:,:,:)  ) !out



    ! line and block indices of edges as seen from cells
    iidx => p_patch%cells%edge_idx
    iblk => p_patch%cells%edge_blk

    ! tracer fields which are advected
    trAdvect => advection_config(jg)%trAdvect       ! 2018-06-05: cray bug do not add to PRESENT list


!$ACC DATA  PRESENT(p_int_state, p_mflx_tracer_h, tracer_now, tracer_new,    &
!$ACC               rhodz_now, rhodz_new, deepatmo_divh, iidx, iblk ),       &
!$ACC       CREATE( z_fluxdiv_c ),                                           &
!$ACC       IF( i_am_accel_node .AND. acc_on )


    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block  (i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jk,jb,jt,nt,iadv_slev_jt,i_startidx,i_endidx,z_fluxdiv_c) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                    i_startidx, i_endidx, i_rlstart, i_rlend)

      ! compute horizontal flux divergences and update tracer mass fractions
      !
      TRACERLOOP: DO nt = 1, trAdvect%len

        jt = trAdvect%list(nt)
        iadv_slev_jt = advection_config(jg)%iadv_slev(jt)

        IF ( advection_config(jg)%ihadv_tracer(jt) /= 0 ) THEN

          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG(static:1) VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
          DO jc = i_startidx, i_endidx
            DO jk = iadv_slev_jt, nlev
#else
!NEC$ outerloop_unroll(8)
          DO jk = iadv_slev_jt, nlev
            DO jc = i_startidx, i_endidx
#endif

! TODO: possible GPU optimization: add tracer_new calculation here
              z_fluxdiv_c(jc,jk) =  deepatmo_divh(jk) * (                                              &
                & p_mflx_tracer_h(iidx(jc,jb,1),jk,iblk(jc,jb,1),jt)*p_int_state%geofac_div(jc,1,jb) + &
                & p_mflx_tracer_h(iidx(jc,jb,2),jk,iblk(jc,jb,2),jt)*p_int_state%geofac_div(jc,2,jb) + &
                & p_mflx_tracer_h(iidx(jc,jb,3),jk,iblk(jc,jb,3),jt)*p_int_state%geofac_div(jc,3,jb) )

            ENDDO
          ENDDO
          !$ACC LOOP GANG(static:1) VECTOR COLLAPSE(2)
!NEC$ nofuse
          DO jk = iadv_slev_jt, nlev
!DIR$ IVDEP
            DO jc = i_startidx, i_endidx

              tracer_new(jc,jk,jb,jt) =                                   &
                &   ( tracer_now(jc,jk,jb,jt) * rhodz_now(jc,jk,jb)       &
                &    - p_dtime * z_fluxdiv_c(jc,jk) ) / rhodz_new(jc,jk,jb)

            ENDDO  !jc
          ENDDO  !jk
          !$ACC END PARALLEL

          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          ! set tracer(nnew) to tracer(nnow) at levels where advection is turned off
          DO jk = 1, iadv_slev_jt-1
            DO jc = i_startidx, i_endidx
              tracer_new(jc,jk,jb,jt) = tracer_now(jc,jk,jb,jt)
            END DO
          END DO
          !$ACC END PARALLEL

        ELSE  ! horizontal advection switched off

          !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF( i_am_accel_node .AND. acc_on )
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          ! copy
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              tracer_new(jc,jk,jb,jt) = tracer_now(jc,jk,jb,jt)
            ENDDO  !jc
          ENDDO  !jk
          !$ACC END PARALLEL

        ENDIF  ! ihadv_tracer(jt) /= 0

      ENDDO TRACERLOOP

    END DO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!$ACC END DATA

    IF (timers_level > 2) CALL timer_stop(timer_adv_horz)

  END SUBROUTINE hor_adv

END MODULE mo_advection_stepping
