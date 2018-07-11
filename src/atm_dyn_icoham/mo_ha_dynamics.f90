!>
!! This module contains subroutines for evaluating the right-hand side
!! of the primitive equations
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan (MPI-M, 2009-11)
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
MODULE mo_ha_dynamics


  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rcpd, rd
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_math_gradients,     ONLY: grad_fd_norm
  USE mo_math_divrot,        ONLY: div, div_avg
  USE mo_dynamics_config,    ONLY: lshallow_water, idiv_method
  USE mo_ha_dyn_config,      ONLY: ha_dyn_config 
  USE mo_parallel_config,    ONLY: nproma, use_icon_comm, p_test_run

  USE mo_run_config,         ONLY: nlev, nlevp1
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_intp,               ONLY: cell_avg, cells2edges_scalar, &
                                   edges2cells_scalar

  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_ha_dynamics_adv,    ONLY: temp_adv_vertical, temp_adv_horizontal, &
                                     vn_adv_vertical,   vn_adv_horizontal
  USE mo_ha_diag_util,       ONLY: update_diag_state
  USE mo_vertical_coord_table, ONLY: nplev, nplvp1, vct_b, alrrdic, rdlnp0i,    &
    &                                rdt0ral, t0icao

  USE mo_sync,                 ONLY: SYNC_C, SYNC_E, sync_patch_array
  USE mo_timer,              ONLY: ltimer, timer_start, timer_stop,&
    & timer_dyn_temp
   
   USE mo_icon_comm_lib,     ONLY: new_icon_comm_variable, &
     & icon_comm_sync, icon_comm_sync_all,  is_ready, &
     & until_sync

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dyn_temp
  PUBLIC :: continuity
  PUBLIC :: energy_conversion_terms

CONTAINS
  !>
  !!
  !! @par Revision History
  !!
  !! First version by Hui Wan (MPI-M, 2009-11)
  !!
  SUBROUTINE dyn_temp( pt_patch, pt_int_state, pt_prog,    & ! input
                       pt_ext_data, pt_diag, pt_tend_dyn   ) ! input, output

  !! Arguments

  TYPE(t_patch), TARGET, INTENT(INOUT) :: pt_patch
  TYPE(t_int_state),TARGET,INTENT(IN)    :: pt_int_state
  TYPE(t_hydro_atm_prog),INTENT(IN)    :: pt_prog
  TYPE(t_external_data),   INTENT(INOUT)    :: pt_ext_data   !< external data

  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag
  TYPE(t_hydro_atm_prog),INTENT(INOUT) :: pt_tend_dyn

  !! Local variables

  REAL(wp) :: z_mdiv_int (nproma,nlevp1,pt_patch%nblks_c)
  REAL(wp) :: z_mdiv     (nproma,nlev  ,pt_patch%nblks_c)

  INTEGER :: nblks_c, nblks_e
  INTEGER :: jb, jbs, is, ie
  INTEGER :: temp_comm, vn_comm

  

  IF (ltimer) CALL timer_start(timer_dyn_temp)
! Dimension parameters related to refinement and MPI parallelisation

   nblks_c = pt_patch%nblks_c
   nblks_e = pt_patch%nblks_e

! Update the diagnostic state vector. This includes the calculation of
! pressure and related quantities, vorticity and divergence, u- and v-wind,
! virtual temperature and geopotential.

   ! includes sync_patch_array(SYNC_V, pt_patch, pt_diag%rel_vort)
   ! essentially async
   CALL update_diag_state( pt_prog, pt_patch, pt_int_state, pt_ext_data, &
     &                     pt_diag )

! Diagnose the mass flux, eta-coordinate vertical velocity,
! and calculate surface pressure tendency
   
   ! includes sync_patch_array(SYNC_E, pt_patch, pt_diag%mass_flux_e)
   ! sync required
   CALL continuity( pt_prog%vn, pt_diag%delp_e,         &! in
                    pt_patch, pt_int_state, .TRUE.,     &! in
                    z_mdiv, z_mdiv_int,                 &! inout
                    pt_diag%mass_flux_e,                &! inout
                    pt_tend_dyn%pres_sfc,               &! inout
                    pt_diag%weta                    )    ! inout

! Initialize velocity and temperature tendencies in the interior of each patch

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      pt_tend_dyn%vn(is:ie,:,jb) = 0._wp
   ENDDO
!$OMP END DO

   IF (.NOT.lshallow_water) THEN
      jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
         CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
         pt_tend_dyn%temp(is:ie,:,jb) = 0._wp
      ENDDO
!$OMP END DO NOWAIT
   ENDIF
!$OMP END PARALLEL

! From now on, after calling each subroutine the individual contribution
! will be added to the tendency state.

   IF (.NOT.lshallow_water) THEN

! Vertical advection of momentum

     CALL vn_adv_vertical( pt_prog%vn, pt_diag%weta, pt_diag%delp_e, &! in
                           pt_patch,   pt_int_state,                 &! in
                           pt_tend_dyn%vn                          )  ! inout

! Vertical and horizontal advection of temperature

     CALL temp_adv_vertical( pt_prog%temp, pt_diag%weta,    &! in
                             pt_diag%rdelp_c, pt_patch,     &! in
                             pt_tend_dyn%temp            )   ! inout

     CALL temp_adv_horizontal( pt_prog%temp, pt_diag%mass_flux_e,  &! in
                               pt_diag%rdelp_c, z_mdiv,            &! in
                               pt_patch, pt_int_state,             &! in
                               pt_tend_dyn%temp                  )  ! inout
   ENDIF

! Horizontal advection and Coriolis force

   ! LL: this is already synced by the continuity
!    CALL sync_patch_array( SYNC_E, pt_patch, pt_diag%mass_flux_e )
   CALL vn_adv_horizontal( pt_prog%vn,                          &! in
                           pt_diag%rel_vort,                    &! in
                           pt_diag%delp_c,                      &! in
                           pt_patch, pt_int_state,              &! in
                           pt_tend_dyn%vn,                      &! inout
                           pt_diag%e_kin, pt_diag%vt,           &! inout
                           pt_diag%delp_v,                      &! inout
                           opt_rlstart=grf_bdywidth_e+1 )        ! for nesting

! Pressure gradient force and adiabatic heating

   CALL energy_conversion_terms( pt_diag, pt_patch, pt_int_state, &! inout,in,in
                                 z_mdiv, z_mdiv_int,              &! in
                                 pt_tend_dyn%temp,                &! inout
                                 pt_tend_dyn%vn                  ) ! inout

! Synchronize tendencies

   ! The following should be aggregated
    IF (use_icon_comm) THEN
      ! halo surface proessure has also been calculated in the continuity
!       pres_sfc_comm = new_icon_comm_variable(pt_tend_dyn%pres_sfc, cells_not_in_domain, pt_patch, &
!         & status=is_ready, scope=until_sync, name="dyn_temp pres_sfc")
      temp_comm = new_icon_comm_variable(pt_tend_dyn%temp, pt_patch%sync_cells_not_in_domain, &
        & status=is_ready, scope=until_sync, name="dyn_temp temp")
      vn_comm = new_icon_comm_variable(pt_tend_dyn%vn, pt_patch%sync_edges_not_owned, &
        & status=is_ready, scope=until_sync, name="dyn_temp vn")
      CALL icon_comm_sync_all()
    ELSE
      ! halo surface proessure has also been calculated in the continuity
!       CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_dyn%pres_sfc )
      CALL sync_patch_array( SYNC_C, pt_patch, pt_tend_dyn%temp     )
      CALL sync_patch_array( SYNC_E, pt_patch, pt_tend_dyn%vn       )
    ENDIF
  
   IF (ltimer) CALL timer_stop(timer_dyn_temp)

  END SUBROUTINE dyn_temp
  !----------------------

  !>
  !! Spatially discretized continuity equation
  !!
  !! Purpose:
  !! Vertically integrate mass divergence to compute surface pressure tendency;
  !! Compute the vertical velocity eta-dot at layer interfaces.)
  !!
  !! @par Revision History
  !! Separated from subroutine dyn by Hui Wan (MPI-M, 2009-11-17)
  !!
  SUBROUTINE continuity( p_vn, p_delp_e, pt_patch, pt_int_state, ldiag_weta, &
                         p_mdiv, p_mdiv_int, p_mflux, p_ddt_psfc, p_weta )

  IMPLICIT NONE

  !! Arguments

  REAL(wp), INTENT(in)  :: p_vn    (:,:,:) !< normal velocity
  REAL(wp), INTENT(in)  :: p_delp_e(:,:,:) !< layer thickness at edges
  TYPE(t_patch), TARGET, INTENT(inout) :: pt_patch   !< grid information
  TYPE(t_int_state),INTENT(in) :: pt_int_state    !< interpolation coefficients

  LOGICAL,INTENT(in) :: ldiag_weta              !< if .true., diagnoise weta

  REAL(wp),INTENT(inout) :: p_mdiv    (:,:,:)   !< mass divergence
  REAL(wp),INTENT(inout) :: p_mdiv_int(:,:,:)   !< mass divergence
                                                !< vertically integrated

!   REAL(wp), POINTER, INTENT(inout) :: p_mflux(:,:,:)  !< mass flux at edges
  REAL(wp), POINTER :: p_mflux(:,:,:)  !< mass flux at edges, out

  REAL(wp),INTENT(inout) :: p_ddt_psfc(:,:) !< tendency of surface pressure

  REAL(wp),INTENT(inout),OPTIONAL :: p_weta (SIZE(p_mdiv_int,1), &!< vertical
                                             SIZE(p_mdiv_int,2), &!< velocity
                                             SIZE(p_mdiv_int,3) ) !< rho*eta-dot

  !! Local variables
!   REAL(wp), POINTER :: p_3d(:,:,:)  

  INTEGER  :: nblks_e, nblks_c
  INTEGER  :: jb, jbs, is,ie, jk,jkp, jc

! Dimension parameters

  nblks_e  = pt_patch%nblks_e
  nblks_c  = pt_patch%nblks_c

! Divergence of mass flux at full levels
  IF (p_test_run) p_mflux(:,:,:) = 0.0_wp
!  CALL sync_patch_array(SYNC_E, pt_patch, p_delp_e)
!  CALL sync_patch_array(SYNC_E, pt_patch, p_vn)


!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 2)
      p_mflux(is:ie,:,jb) = p_delp_e(is:ie,:,jb)*p_vn(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   ! LL The assumpiton is that inner edges are the ones
   !    adjacent to the owned cells. Thus no communication
   !    is required here.
   IF (use_icon_comm) THEN
!      p_3d => p_mflux
     CALL icon_comm_sync(p_mflux, pt_patch%sync_edges_not_owned)
   ELSE
     CALL sync_patch_array(SYNC_E, pt_patch, p_mflux)
   ENDIF

   SELECT CASE(idiv_method)
   CASE(1)

       CALL div( p_mflux, pt_patch, pt_int_state, p_mdiv, opt_rlstart=2 )

   CASE(2)

       CALL div_avg( p_mflux, pt_patch, pt_int_state, &
                     pt_int_state%c_bln_avg, p_mdiv, opt_rlstart=2 )

   END SELECT

! Vertically integrated mass divergence at layer interfaces

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%cells%start_blk(2,1)
!$OMP DO PRIVATE(jb,is,ie,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)

      p_mdiv_int(is:ie,1,jb) = 0.0_wp
      DO jk = 1,nlev
         jkp = jk + 1
         p_mdiv_int(is:ie,jkp,jb) = p_mdiv_int(is:ie,jk,jb)+p_mdiv(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO

! Vertical velocity at the interfaces (half-levels)

   IF ((.NOT.lshallow_water).AND.ldiag_weta) THEN
!$OMP DO PRIVATE(jb,jc,is,ie,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)

!DR Workaround for gcc 4.5.0 internal compiler error
        DO jc=is,ie
          p_weta(jc,1,jb)      = 0._wp   ! upper boundary
          p_weta(jc,nlevp1,jb) = 0._wp   ! lower boundary
        ENDDO

        DO jk = 2,nlev
          DO jc=is,ie
            p_weta(jc,jk,jb)  = -p_mdiv_int(jc,jk,jb) + &
                                    p_mdiv_int(jc,nlevp1,jb)*vct_b(jk)
          ENDDO
        ENDDO

!DR        p_weta(is:ie,1,jb)      = 0._wp   ! upper boundary
!DR        p_weta(is:ie,nlevp1,jb) = 0._wp   ! lower boundary

!DR        DO jk = 2,nlev
!DR           p_weta(is:ie,jk,jb)  = -p_mdiv_int(is:ie,jk,jb) + &
!DR                                   p_mdiv_int(is:ie,nlevp1,jb)*vct_b(jk)
!DR        ENDDO
      ENDDO
!$OMP END DO
   ENDIF

! Tendency of surface pressure
! (For nested domains, tendencies are interpolated from the parent domain
! on a boundary zone with a width of grf_bdywidth_c for cells and
! grf_bdywidth_e for edges, respectively. These tendencies must not be
! overwritten.)

   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
      p_ddt_psfc(is:ie,jb) = -p_mdiv_int(is:ie,nlevp1,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE continuity


  !>
  !! energy_conversion_terms
  !!
  !! Purpose:
  !!
  !! @par Revision History
  !! Separated from subroutine dyn and rewritten by Hui Wan (MPI-M, 2009-11-18)
  !!
  SUBROUTINE energy_conversion_terms( pt_diag, pt_patch, pt_int_state, &
                                      p_mdiv, p_mdiv_int,              &
                                      p_ddt_temp, p_ddt_vn,            &
                                      opt_lseparate, opt_ddt_temp_fast )

  !! Arguments

  TYPE(t_hydro_atm_diag),INTENT(INOUT) :: pt_diag
  TYPE(t_patch),TARGET,    INTENT(IN)    :: pt_patch
  TYPE(t_int_state),       INTENT(IN)    :: pt_int_state

  REAL(wp),INTENT(IN)    :: p_mdiv     (:,:,:)
  REAL(wp),INTENT(IN)    :: p_mdiv_int (:,:,:)
  REAL(wp),INTENT(INOUT) :: p_ddt_temp (:,:,:)
  REAL(wp),INTENT(INOUT) :: p_ddt_vn   (:,:,:)

  LOGICAL, INTENT(IN),   OPTIONAL :: opt_lseparate
  REAL(wp),INTENT(INOUT),OPTIONAL :: opt_ddt_temp_fast(:,:,:)

  !! Local variables

  LOGICAL  :: lseparate
  INTEGER  :: jk,jkp
  INTEGER  :: jb,jbs,is,ie
  INTEGER  :: nblks_e,nblks_c
  REAL(wp) :: z2d (nproma,nlev)

  REAL(wp) :: z_geo_mc( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_tv_c  ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_tv_e  ( nproma,nlev,pt_patch%nblks_e )
  REAL(wp) :: z_tvp_c ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_tvp_e ( nproma,nlev,pt_patch%nblks_e )
  REAL(wp) :: z_tmp_e ( nproma,nlev,pt_patch%nblks_e )
  REAL(wp) :: z_tmp_c ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_rlp_c ( nproma,nlev,pt_patch%nblks_c )
  REAL(wp) :: z_fast  ( nproma,nlev,pt_patch%nblks_c )

! Check optional input

  IF (PRESENT(opt_lseparate).AND.PRESENT(opt_ddt_temp_fast)) THEN
     lseparate = opt_lseparate
  ELSE
     lseparate = .FALSE.
  ENDIF

! Dimension parameters

  nblks_c = pt_patch%nblks_c
  nblks_e = pt_patch%nblks_e

!=====================================================================
! Horizontal gradient of geopotential
!=====================================================================
IF (.NOT.lshallow_water) THEN
! For the hydrostatic model,
! if the use of a reference state is desired (to reduce the numerical
! error near steep topography), we need to first construct the reference
! state and calculate the temperature and geopotential perturbation.
! The reference state we use here is the ICAO(1964) standard atmosphere.
! See MPI-M Report 349, p20.

  IF (ha_dyn_config%lref_temp) THEN

    jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,z2d) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_c
       CALL get_indices_c( pt_patch, jb,jbs,nblks_c, is,ie, 2)

       z2d(is:ie,:) = rd*LOG(pt_diag%pres_mc(is:ie,:,jb)) - rdlnp0i
       z2d(is:ie,:) = EXP( alrrdic*z2d(is:ie,:) )

       z_tv_c  (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb)
       z_tvp_c (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb) -  t0icao*z2d(is:ie,:)
       z_geo_mc(is:ie,:,jb) = pt_diag%geo_mc(is:ie,:,jb) + rdt0ral*z2d(is:ie,:)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ELSE !Do not use reference state

    jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = jbs,nblks_c
       CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)
       z_tv_c  (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb)
       z_tvp_c (is:ie,:,jb) = pt_diag%tempv (is:ie,:,jb)
       z_geo_mc(is:ie,:,jb) = pt_diag%geo_mc(is:ie,:,jb)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

ELSE !shallow water model
!$OMP PARALLEL
!$OMP WORKSHARE
      z_geo_mc(:,:,:) = pt_diag%geo_mc(:,:,:)
!$OMP END WORKSHARE
!$OMP END PARALLEL
ENDIF

!--------------------------------------------------------------------------
! Calculate the gradient of geopotential, and accumulate velocity tendency
!--------------------------------------------------------------------------

  CALL grad_fd_norm( z_geo_mc, pt_patch, z_tmp_e, opt_rlstart=4 )

  jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)
      p_ddt_vn(is:ie,:,jb) = p_ddt_vn(is:ie,:,jb) - z_tmp_e(is:ie,:,jb)
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

!=====================================================================
! The other part of the pressure gradient force and the thermodynamic
! equation only exist in the hydrostatic model
!=====================================================================
IF (.NOT.lshallow_water) THEN

! Average virtual temperature perturbation from cells to edges
! Note: by giving the vertical start index optionally, the computation
! is done only for non-pressure levels

   CALL cells2edges_scalar( z_tvp_c, pt_patch, pt_int_state%c_lin_e, &
                            z_tvp_e, nplvp1, opt_rlstart=4)

! The counterpart to be used in the thermodynamic equation:
! Even if the reference state is involved, the adiabatic heating
! is still calculated using the original variables.

   CALL cells2edges_scalar( z_tv_c, pt_patch, pt_int_state%c_lin_e, &
                            z_tv_e, nplvp1, opt_rlstart=4)

! Pressure gradient at edges

   jbs = pt_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,jk,jkp) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs, nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 2)
      DO jk = nplvp1, nlev
         jkp = jk+1
         z_rlp_c(is:ie,jk,jb) = rd*pt_diag%lnp_ic(is:ie,jkp,jb) &
                              - pt_diag%rdalpha_c(is:ie,jk ,jb)
      ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   CALL grad_fd_norm( z_rlp_c, pt_patch, z_tmp_e, nplvp1, opt_rlstart=4 )

!--------------------------------------------------------------------
! Accumulate velocity tendency
!--------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(jbs)
   jbs = pt_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP DO PRIVATE(jb,is,ie,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, grf_bdywidth_e+1)

      DO jk = nplvp1, nlev
         p_ddt_vn(is:ie,jk,jb) = p_ddt_vn(is:ie,jk,jb) &
                                -z_tvp_e(is:ie,jk,jb)*z_tmp_e(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO

!--------------------------------------------------------------------
! Adiabatic heating terms in the thermodynamic equation
!--------------------------------------------------------------------
! Part I: vn*[Rd*T/p*grad(p)]. Use the grad(lnp) calculated above.

   jbs = pt_patch%edges%start_blk(4,1)
!$OMP DO PRIVATE(jb,is,ie,jk) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_e
      CALL get_indices_e(pt_patch, jb,jbs,nblks_e, is,ie, 4)

      DO jk = 1, nplev
         z_tmp_e(is:ie,jk,jb) = 0._wp  !Pressure gradient vanishes on p-levels
      ENDDO

      DO jk = nplvp1, nlev
         z_tmp_e(is:ie,jk,jb) =  pt_diag%mass_flux_e(is:ie,jk,jb) &
                                *z_tv_e (is:ie,jk,jb)             &
                                *z_tmp_e(is:ie,jk,jb)
      ENDDO
   ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! Interpolation from edges to cell centers

   CALL edges2cells_scalar( z_tmp_e, pt_patch, pt_int_state%e_inn_c, &
                            z_tmp_c, opt_rlstart=3 )

! Add smoothing if desired. (Only when using the original gradient operator.
! Otherwise the model will become unstable.

   SELECT CASE (idiv_method)

   CASE(2)

!$OMP PARALLEL
!$OMP WORKSHARE
       z_rlp_c = z_tmp_c
!$OMP END WORKSHARE
!$OMP END PARALLEL
       CALL cell_avg( z_rlp_c, pt_patch, pt_int_state%c_bln_avg, &
                      z_tmp_c, opt_rlstart=4)

   END SELECT

! Divide by the pseudo-density

!$OMP PARALLEL  PRIVATE(jbs)
   jbs = pt_patch%cells%start_blk(3,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 3)
      z_tmp_c(is:ie,:,jb) = z_tmp_c(is:ie,:,jb) * pt_diag%rdelp_c(is:ie,:,jb)
   ENDDO
!$OMP END DO

!-------------------------------------------------------------------
! Part II: Rd*T/p *[ p-tendency + vertical-adv ]

   jbs = pt_patch%cells%start_blk(3,1)
!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
      CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, 3)

      z_fast(is:ie,:,jb) = -pt_diag%rdelp_c(is:ie,:,jb)*z_tv_c(is:ie,:,jb)      &
                 *( pt_diag%rdlnpr_c (is:ie,:,jb) * p_mdiv_int(is:ie,1:nlev,jb) &
                   +pt_diag%rdalpha_c(is:ie,:,jb) * p_mdiv    (is:ie,:,jb)     )
   ENDDO
!$OMP END DO

!-------------------------------------------------------------------
! Accumulate temperature tendency

   jbs = pt_patch%cells%start_blk(grf_bdywidth_c+1,1)
   IF (lseparate) THEN !Provide the fast and slow components separately

!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,nblks_c
        CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
        p_ddt_temp(is:ie,:,jb) = p_ddt_temp(is:ie,:,jb) + z_tmp_c(is:ie,:,jb)*rcpd
        opt_ddt_temp_fast(is:ie,:,jb) = z_fast(is:ie,:,jb)*rcpd
      ENDDO
!$OMP END DO NOWAIT

  ELSE !Add both fast and slow components to the tendency state

!$OMP DO PRIVATE(jb,is,ie) ICON_OMP_DEFAULT_SCHEDULE
   DO jb = jbs,nblks_c
     CALL get_indices_c(pt_patch, jb,jbs,nblks_c, is,ie, grf_bdywidth_c+1)
     p_ddt_temp(is:ie,:,jb) = p_ddt_temp(is:ie,:,jb) &
                            + ( z_tmp_c(is:ie,:,jb) + z_fast(is:ie,:,jb) )*rcpd
   ENDDO
!$OMP END DO NOWAIT

  ENDIF
!$OMP END PARALLEL

ENDIF !shallow_water vs hydrostatic
!------------------------------------------
  END SUBROUTINE energy_conversion_terms

END MODULE mo_ha_dynamics
