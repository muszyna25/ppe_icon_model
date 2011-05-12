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
MODULE mo_advection_stepping

  USE mo_kind,                ONLY: wp
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_transport
  USE mo_model_domain,        ONLY: t_patch
  USE mo_math_operators,      ONLY: div
  USE mo_interpolation,       ONLY: t_int_state
  USE mo_run_nml,             ONLY: nproma, ntracer, ltimer, i_cell_type,   &
    &                               iforcing, inwp, iqv
  USE mo_nonhydrostatic_nml,  ONLY: iadv_rcf
  USE mo_advection_hflux,     ONLY: hor_upwind_flux
  USE mo_advection_vflux,     ONLY: vert_upwind_flux
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_impl_constants,      ONLY: min_rlcell_int, min_rledge_int, &
    &                               min_rlcell, min_rledge
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_mpi,                 ONLY: p_pe, p_nprocs
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array_mult
  USE mo_parallel_nml,        ONLY: p_test_pe, p_test_run
  USE mo_advection_nml,       ONLY: ihadv_tracer, ivadv_tracer, lvadv_tracer, &
    &                               lclip_tracer, lstrang, itype_vlimit,      &
    &                               itype_hlimit, iord_backtraj,              &
    &                               igrad_c_miura, iadv_slev, iubc_adv, cSTR, &
    &                               coeff_grid
  USE mo_advection_utils,     ONLY: ptr_delp_mc_now, ptr_delp_mc_new
  USE mo_model_domain_import, ONLY: l_limited_area, lfeedback

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: step_advection


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
  !! 
  !!
  SUBROUTINE step_advection( p_patch, p_int_state, p_dtime, k_step, p_tracer_now, &
    &                        p_mflx_contra_h, p_vn_contra_traj, p_mflx_contra_v,  &
    &                        p_w_contra_traj, p_cellhgt_mc_now, p_delp_mc_new,    &
    &                        p_delp_mc_now, p_pres_mc_now, p_pres_ic_now,         &
    &                        p_grf_tend_tracer, p_tracer_new, p_mflx_tracer_h,    &
    &                        p_mflx_tracer_v, opt_rho_ic, opt_topflx_tra,         &
    &                        opt_q_int, opt_ddt_tracer_adv      )
  !
    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed
                                             

    TYPE(t_int_state), TARGET, INTENT(IN) :: & !< interpolation state
      &  p_int_state                         

    REAL(wp), TARGET, INTENT(IN) ::  &  !< tracer mixing ratios (specific concentrations)
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

    REAL(wp), INTENT(IN)  ::         &  !< NH: full level height [m] 
      &  p_pres_mc_now(:,:,:)           !< HA: full level pressure at time step n [Pa] 
                                        !< dim: (nproma,nlev,nblks_c)
          
    REAL(wp), INTENT(IN)  ::         &  !< NH: half level height [m] 
      &  p_pres_ic_now(:,:,:)           !< HA: half level pressure at time step n [Pa]

    REAL(wp), INTENT(INOUT) ::       &  !< interpolated tracer time tendencies for
      &  p_grf_tend_tracer(:,:,:,:)     !< updating the lateral boundaries of nested domains
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), TARGET, INTENT(INOUT) ::  &  !< tracer mixing ratios (specific concentrations)
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

    REAL(wp), INTENT(IN), OPTIONAL :: & !< NH: half level density at n+1/2 (only necessary
      &  opt_rho_ic(:,:,:)              !< for the NH-core, when muscl_cfl or ppm_cfl
                                        !< is applied in vertical direction)
                                        !< HA: --
                                        !< dim: (nproma,nlevp1,nblks_c)

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

    REAL(wp) ::  &                      !< flux divergence at cell center
      &  z_fluxdiv_c(nproma,p_patch%nlev,p_patch%nblks_c,ntracer) 

    REAL(wp) ::  &                      !< reciprocal of p_cellhgt_mc_now
      &  z_rcellhgt_mc_now(nproma,p_patch%nlev,p_patch%nblks_c)

    REAL(wp), POINTER ::  &
      &  ptr_current_tracer(:,:,:,:) => NULL()  !< pointer to tracer field

    REAL(wp) ::  &                              !< result of intermediate RK2 step
      &  z_estim_c(nproma,p_patch%nlev,p_patch%nblks_c,ntracer)       

    REAL(wp) :: pdtime_mod        !< modified time step
                                  !< (multiplied by cSTR * coeff_grid)

    INTEGER  :: nlev, nlevp1      !< number of full and half levels
    INTEGER  :: jb, jk, jt, jc, jg                !< loop indices
    INTEGER  :: ikp1                              !< vertical level + 1
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

    LOGICAL  :: l_parallel

    INTEGER, POINTER :: i_itype_hlimit(:) => NULL()
    INTEGER, TARGET  :: itype_hlimit_0(ntracer)
     

!AL   REAL(wp) :: z_min, z_max
!DR   CHARACTER(LEN=MAX_CHAR_LENGTH) :: ctrdiag

   !-----------------------------------------------------------------------

    IF(ltimer) CALL timer_start(timer_transport)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg  = p_patch%id

    ! precompute modified timestep (timestep multiplied by Strang-splitting
    ! coefficient and a second coefficient to account for either a
    ! vertical z- or p-system)
    pdtime_mod = cSTR * coeff_grid * p_dtime


    IF (p_nprocs == 1 .OR. p_pe == p_test_pe) THEN
      l_parallel = .FALSE.
    ELSE
      l_parallel = .TRUE.
    ENDIF


    !---------------------------------------------------!
    !                                                   !
    !  time integration of tracer continuity-equation   !
    !                                                   !
    !---------------------------------------------------!

    !
    ! Godunov-type advection with two-time level scheme
    !

    !
    ! calculate field of reciprocal cell height at time step n
    ! (HDC:\delta p; NHDC: \delta z)
    !
    i_rlstart  = 1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,jt,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! Only rdelp is updated here
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          z_rcellhgt_mc_now(jc,jk,jb) = 1._wp/p_cellhgt_mc_now(jc,jk,jb)

        END DO

        DO jt=1,ntracer
          DO jc = i_startidx, i_endidx
            ! IF slev >1, then z_fluxdiv_c must be initialized with 0.
            ! For p_test_run=.TRUE. this is necessary, anyway
            z_fluxdiv_c(jc,jk,jb,jt) = 0._wp
          ENDDO
        ENDDO

      END DO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL



    ptr_current_tracer => p_tracer_now
    ptr_delp_mc_now    => p_delp_mc_now



    !*********************************!
    ! Vertical advection              !
    !*********************************!
    ! calculate vertical tracer flux
    ! subroutine vertical_upwind_flux will be called if either the
    ! current time step is even or complete Strang-splitting is chosen.
    !
    IF ( lvadv_tracer ) THEN

      IF ( MOD( k_step, 2 ) == 0 .OR. lstrang ) THEN


        i_rlstart  = grf_bdywidth_c-1
        i_rlend    = min_rlcell_int
        i_startblk = p_patch%cells%start_blk(i_rlstart,1)
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

        CALL vert_upwind_flux( p_patch, ptr_current_tracer,          &! in
          &              p_mflx_contra_v, p_w_contra_traj,           &! inout,in
          &              cSTR*p_dtime, p_pres_ic_now, p_pres_mc_now, &! in
          &              p_cellhgt_mc_now, z_rcellhgt_mc_now,        &! in
          &              ptr_delp_mc_now, ivadv_tracer,              &! in
          &              itype_vlimit, iubc_adv(jg), iadv_slev(jg,:),&! in
          &              p_mflx_tracer_v,                            &! out
          &              opt_topflx_tra=opt_topflx_tra,              &! in
          &              opt_q_int=opt_q_int,                        &! out
          &              opt_rho_ic=opt_rho_ic,                      &! in
          &              opt_rlstart=i_rlstart,                      &! in
          &              opt_rlend=i_rlend                           )! in


        ! calculation of intermediate layer thickness (density)
        ! necessary for tracer-mass consistency.
        ptr_delp_mc_new  => z_delp_mc1

        ! integration of tracer continuity equation in vertical
        ! direction. This is independent of the tracer and thus must be
        ! computed only once.

!$OMP PARALLEL PRIVATE(i_rlend, i_endblk)

        ! Note that intermediate densities are needed by the horizontal 
        ! limiter as well. That is, why we have to extend the following 
        ! computation to halo points.
        i_rlend    = min_rlcell
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )

          ptr_delp_mc_new(i_startidx:i_endidx,1:nlev,jb) =                      &
            &              ptr_delp_mc_now(i_startidx:i_endidx,1:nlev,jb)       &
            &              - pdtime_mod                                         &
            &              * ( p_mflx_contra_v(i_startidx:i_endidx,2:nlevp1,jb) &
            &              - p_mflx_contra_v(i_startidx:i_endidx,1:nlev,jb) )
        ENDDO
!$OMP END DO

        !
        ! compute vertical flux divergence
        !
        i_rlend    = min_rlcell_int
        i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

        ! Note that we need to start the calculation within the
        ! nest boundary, since the following horizontal flux calculation
        ! uses a non-local stencil.

!$OMP DO PRIVATE(jb,jk,jt,jc,ikp1,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jt = 1, ntracer ! Tracer loop

            DO jk = iadv_slev(jg,jt), nlev

              ! index of top half level
              ikp1 = jk + 1

              DO jc = i_startidx, i_endidx

                p_tracer_new(jc,jk,jb,jt) =                                         &
                  &  ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb)  &
                  &  - pdtime_mod * ( p_mflx_tracer_v(jc,ikp1,jb,jt)                &
                  &               -   p_mflx_tracer_v(jc,jk  ,jb,jt) ) )            &
                  &  / ptr_delp_mc_new(jc,jk,jb)

              END DO
            END DO
          END DO  ! Tracer loop
        END DO
!$OMP END DO
!$OMP END PARALLEL


        ptr_current_tracer => p_tracer_new
        ptr_delp_mc_now    => z_delp_mc1

        IF ( lstrang ) THEN

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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
          DO jb = i_startblk, i_endblk

            CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
              &                 i_startidx, i_endidx, i_rlstart, i_rlend )

            ptr_delp_mc_new(i_startidx:i_endidx,1:nlev,jb) =                   &
              &          p_delp_mc_new(i_startidx:i_endidx,1:nlev,jb)          &
              &          + pdtime_mod                                          &
              &          * ( p_mflx_contra_v(i_startidx:i_endidx,2:nlevp1,jb)  &
              &          - p_mflx_contra_v(i_startidx:i_endidx,1:nlev,jb) )
          ENDDO
!$OMP END DO
!$OMP END PARALLEL

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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx)
        DO jb = i_startblk, i_endblk

          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
            &                 i_startidx, i_endidx, i_rlstart, i_rlend )
 
          ptr_delp_mc_new(i_startidx:i_endidx,1:nlev,jb) =                        &
            &                  p_delp_mc_new(i_startidx:i_endidx,1:nlev,jb)       &
            &                + pdtime_mod                                         &
            &                * ( p_mflx_contra_v(i_startidx:i_endidx,2:nlevp1,jb) &
            &                - p_mflx_contra_v(i_startidx:i_endidx,1:nlev,jb) )
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END IF

      ! IF vertical advection proceeds horizontal advection, synchronize the 
      ! updated tracer array. For efficiency, the synchronization is applied for all 
      ! tracers at once
      CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=p_tracer_new, &
        &                        lpart4d=.TRUE.)


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
    ! For the hexagonal model, we have to perform (at least) an RK2 time integration
    ! without the limiter in the estimation step
    IF (i_cell_type == 6) THEN
      itype_hlimit_0 = 0
      i_itype_hlimit => itype_hlimit_0
      i_rlend        = min_rledge
    ELSE
      i_itype_hlimit => itype_hlimit
      i_rlend        = min_rledge_int-1
    ENDIF    
    !
    CALL hor_upwind_flux( ptr_current_tracer, ptr_current_tracer, p_mflx_contra_h, &! in
      &                 p_vn_contra_traj, p_dtime, p_patch,             &! in
      &                 p_int_state, ihadv_tracer, igrad_c_miura,       &! in
      &                 i_itype_hlimit, iadv_slev(jg,:), iord_backtraj, &! in
      &                 p_mflx_tracer_h, opt_rlend=i_rlend              )! inout


    IF (i_cell_type == 6) THEN
      i_rlend        = min_rlcell
    ELSE
      i_rlend        = min_rlcell_int
    ENDIF    

    !
    !  compute divergence of the upwind fluxes for tracers
    !  using optimized divergence routine for 4D fields
    CALL div( p_patch, p_int_state, p_mflx_tracer_h,               &! in
      &       z_fluxdiv_c,                                         &! inout
      &       ntracer, opt_slev=iadv_slev(jg,:), opt_rlend=i_rlend )! in


    IF (i_cell_type == 6) THEN
!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk,jb)
      !
      ! update tracer array
      !
      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell
      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
!$OMP DO PRIVATE(jk,jt,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,      &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)
        DO jt = 1, ntracer ! Tracer loop
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx
              z_estim_c(jc,jk,jb,jt) =                                         &
                &   ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb) &
                &    - 0.5_wp*p_dtime * z_fluxdiv_c(jc,jk,jb,jt) )                &
                &    / (0.5_wp*(ptr_delp_mc_now(jc,jk,jb)+ptr_delp_mc_new(jc,jk,jb)))
            ENDDO
          ENDDO
        ENDDO  ! Tracer loop
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=z_estim_c)

      CALL hor_upwind_flux(z_estim_c, ptr_current_tracer, p_mflx_contra_h, &! in
        &                 p_vn_contra_traj, p_dtime, p_patch,          &! in
        &                 p_int_state, ihadv_tracer, igrad_c_miura,    &! in
        &                 itype_hlimit, iadv_slev(jg,:), iord_backtraj,&! in
        &                 p_mflx_tracer_h, opt_rlend=min_rledge        )! inout

      !
      !  compute divergence of the upwind fluxes for tracers
      !
      DO jt = 1, ntracer
        CALL div( p_mflx_tracer_h(:,:,:,jt),         &! in
          &       p_patch, p_int_state,              &! in
          &       z_fluxdiv_c(:,:,:,jt),             &! inout
          &       opt_slev=iadv_slev(jg,jt)          )! in
      ENDDO

    ENDIF ! i_cell_type == 6

!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk,jb,i_startidx,i_endidx)
    !
    ! update tracer array
    !
    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk  (i_rlend,i_nchdom)


!$OMP DO PRIVATE(jk,jt,jc)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                    i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jt = 1, ntracer ! Tracer loop

        DO jk = iadv_slev(jg,jt), nlev

          DO jc = i_startidx, i_endidx

            p_tracer_new(jc,jk,jb,jt) =                                         &
              &   ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb) &
              &    - p_dtime * z_fluxdiv_c(jc,jk,jb,jt) )                       &
              &    / ptr_delp_mc_new(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO  ! Tracer loop

    ENDDO
!$OMP END DO


    ! Update lateral boundaries of nested domains with interpolated time tendencies
    IF (l_limited_area .OR. p_patch%id > 1) THEN
      i_rlstart  = 1
      i_rlend    = grf_bdywidth_c

      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,1)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

! OpenMP parallelization is done over jk loop here because jb loop is too short in practice
!$OMP DO PRIVATE(jk,jt,jc)
        DO jk = 1, nlev
          ! Tracer values are clipped here to avoid generation of negative values
          ! For mass conservation, a correction has to be applied in the
          ! feedback routine anyway
          DO jt = 1, ntracer ! Tracer loop
            DO jc = i_startidx, i_endidx
              p_tracer_new(jc,jk,jb,jt) =                            &
                &     MAX(0._wp, p_tracer_now(jc,jk,jb,jt)           &
                &   + p_dtime * p_grf_tend_tracer(jc,jk,jb,jt) )
            ENDDO
          ENDDO  ! Tracer loop
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
!$OMP END PARALLEL


    !*********************************!
    ! vertical tracer advection       !
    !*********************************!
    ! subroutine vertical_upwind_flux will be called if either the
    ! current time step is odd or complete Strang-splitting is chosen.
    !
    IF ( lvadv_tracer .AND. ( MOD( k_step, 2 ) == 1 .OR. lstrang) ) THEN

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
        &              p_mflx_contra_v, p_w_contra_traj,           &! inout,in
        &              cSTR*p_dtime, p_pres_ic_now, p_pres_mc_now, &! in
        &              p_cellhgt_mc_now, z_rcellhgt_mc_now,        &! in
        &              ptr_delp_mc_now, ivadv_tracer,              &! in
        &              itype_vlimit, iubc_adv(jg), iadv_slev(jg,:),&! in
        &              p_mflx_tracer_v,                            &! out
        &              opt_topflx_tra=opt_topflx_tra,              &! in
        &              opt_q_int=opt_q_int,                        &! out
        &              opt_rho_ic=opt_rho_ic,                      &! in
        &              opt_rlstart=i_rlstart,                      &! in
        &              opt_rlend=i_rlend                           )! in


      ! calculate vertical flux divergence
      !
      ! update tracer array
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikp1)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jt = 1, ntracer ! Tracer loop

          DO jk = iadv_slev(jg,jt), nlev

            ! index of top half level
            ikp1 = jk + 1

            DO jc = i_startidx, i_endidx

              p_tracer_new(jc,jk,jb,jt) =                                         &
                &  ( ptr_current_tracer(jc,jk,jb,jt) * ptr_delp_mc_now(jc,jk,jb)  &
                &  - pdtime_mod * ( p_mflx_tracer_v(jc,ikp1,jb,jt)                &
                &               -   p_mflx_tracer_v(jc,jk  ,jb,jt) ) )            &
                &  / ptr_delp_mc_new(jc,jk,jb)

            END DO
          END DO
        END DO  ! Tracer loop
      END DO
!$OMP END DO
!$OMP END PARALLEL

    END IF

    ! Synchronize tracer array after update. This is only necessary, if 
    ! NWP physics has NOT been selected. Otherwise, the SYNC-operation will 
    ! follow AFTER the call of NWP physics.
    ! For efficiency, the synchronization is applied for all tracers at once
    IF (iforcing /= inwp) THEN
      CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=p_tracer_new, &
        &                        lpart4d=.TRUE.)
    ENDIF

    !
    ! compute advective tracer tendency
    !
!$OMP PARALLEL PRIVATE(i_rlstart,i_rlend,i_startblk,i_endblk)
    IF (PRESENT(opt_ddt_tracer_adv)) THEN

      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int
      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jt,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,     &
                          i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jt = 1, ntracer
          DO jk = iadv_slev(jg,jt), nlev
            DO jc = i_startidx, i_endidx
              opt_ddt_tracer_adv(jc,jk,jb,jt) =               &
                &           (  p_tracer_new(jc,jk,jb,jt)      &
                &           -  p_tracer_now(jc,jk,jb,jt) )    &
                &            /p_dtime
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    ENDIF


    !
    ! eventually do a clipping of negative values to zero
    !
    IF ( lclip_tracer ) THEN
!$OMP WORKSHARE
      WHERE ( p_tracer_new(:,:,:,:) < 0._wp )
        p_tracer_new(:,:,:,:) = 0._wp
      END WHERE
!$OMP END WORKSHARE
    END IF
!$OMP END PARALLEL

    IF (PRESENT(opt_ddt_tracer_adv)) THEN
      IF (iforcing /= inwp) THEN
        CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer,  &
                                 & f4din=opt_ddt_tracer_adv, &
                                 & lpart4d=.TRUE.)
      ENDIF
    ENDIF

!   output some diagnostics
!DR     DO jk = 1, nlev
!DR      WRITE(0,'(2X,I2)') jk
!DR      DO jt = 1, ntracer
!DR        z_min = MINVAL( p_tracer_new(:,jk,:,jt) )
!DR        z_max = MAXVAL( p_tracer_new(:,jk,:,jt) )
!DR        WRITE( 0,'(I1,A2,E16.9,2X,E16.9)' )  &
!DR          &  jt, ": ", z_min, z_max
!DR      END DO
!DR      WRITE(6,'(A)') TRIM(ctrdiag)
!DR     END DO

  IF (ltimer) CALL timer_stop(timer_transport)

  END SUBROUTINE step_advection


END MODULE mo_advection_stepping

