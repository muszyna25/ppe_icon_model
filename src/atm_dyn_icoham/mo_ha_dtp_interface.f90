!>
!! This module contains subroutines that are used at the
!! dynamics-tracer-physics interface in the hydrostatic dynamical core.
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-02-02)
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
MODULE mo_ha_dtp_interface

  USE mo_kind,               ONLY: wp
  USE mo_dynamics_nml,       ONLY: itime_scheme
  USE mo_parallel_configuration,  ONLY: nproma
  USE mo_run_nml,            ONLY: nlev, nlevp1
  USE mo_model_domain,       ONLY: t_patch
  USE mo_ext_data,           ONLY: t_external_data
  USE mo_interpolation,      ONLY: t_int_state, verts2edges_scalar,       &
                                   edges2cells_scalar, verts2cells_scalar
  USE mo_icoham_dyn_types,   ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_ha_dynamics,        ONLY: continuity
  USE mo_ha_diag_util,       ONLY: update_omega,       &
                                   update_pres_delp_c, &
                                   update_delp_e,      &
                                   update_pres,        &
                                   update_tempv_geopot
  USE mo_math_operators,     ONLY: rot_vertex
  USE mo_physical_constants, ONLY: rd_o_cpd, p0ref
  USE mo_sync,               ONLY: SYNC_E, SYNC_V, sync_patch_array
  USE mo_impl_constants,     ONLY: TWO_TL_SI

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: prepare_tracer
  PUBLIC :: prepare_tracer_RK
  PUBLIC :: prepare_tracer_leapfrog
  PUBLIC :: prepare_physics, prepare_echam_phy

CONTAINS

  !>
  !! Diagnose some pressure- and velocity-related quantities
  !! for the tracer transport scheme, under the assumption that
  !! a two time level semi-implicit time stepping scheme is
  !! used by the dynamical core.
  !!
  !! Note that after calling this subroutine, many of the pressure
  !! related components of the diagnostic state will be updated
  !! to time t+dt, while other components, e.g., temperature and
  !! geopotential, are still at time t.
  !!
  !! @par Revision History
  !! First version by Hui Wan, MPI-M (2010-02-28)
  !!
  SUBROUTINE prepare_tracer_RK( p_patch, p_int_state,          &! in
                                p_now, p_new,                  &! in
                                p_diag,                        &! inout
                                p_vn_traj,                     &! out
                                p_delp_mc_now,                 &! out
                                p_pres_mc_now, p_pres_ic_now   )! out

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch
    TYPE(t_int_state),INTENT(in) :: p_int_state

    TYPE(t_hydro_atm_prog),INTENT(in)    :: p_now, p_new
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

    REAL(wp),INTENT(out) :: p_vn_traj     (nproma,nlev,  p_patch%nblks_e)
    REAL(wp),INTENT(out) :: p_delp_mc_now (nproma,nlev  ,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_pres_mc_now (nproma,nlev  ,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_pres_ic_now (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp) :: z_delp_me_now (nproma,nlev,  p_patch%nblks_e)
    INTEGER  :: jb, jbs, nblks_e, is, ie

    ! Diagnose pressure-related quantities of time step n, including
    ! full- and half-level pressure values and layer thickess, all
    ! computed at cell centers (i.e., mass points).
    ! Note that the computed values are given to dummy arguments
    ! rather than components of the diagnostic state vector.

    CALL update_pres_delp_c( p_patch, p_now%pres_sfc,       &! in
                           & p_pres_mc_now, p_pres_ic_now,  &! out
                           & p_delp_mc_now                  )! out

    ! Horizontal velocitiy = mass flux / delp_e
    ! The mass "delp_e" is at time step n, derived from "p_delp_mc_now"
    ! computed just above;
    ! The mass flux used here is the component "%mass_flux_e" of the
    ! diagnostic state vector, which contains the mass flux temporally
    ! integrated over the RK substeps. The integration is done in
    ! "step_RungeKutta", before this subroutine is called. The use of
    ! integrated mass flux is essential for achieving mass-tracer consistency.

    CALL update_delp_e( p_patch, p_int_state, p_delp_mc_now, z_delp_me_now)

    nblks_e = p_patch%nblks_int_e
    jbs     = p_patch%edges%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_e
      CALL get_indices_e(p_patch, jb,jbs,nblks_e, is,ie, 2)
      p_vn_traj(is:ie,:,jb) =  p_diag%mass_flux_e(is:ie,:,jb) &
                            & /z_delp_me_now(is:ie,:,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Pressure values and layer thickness of time step n+1.
    ! The thickness will be used in tracer transport to ensure
    ! conservation of tracer mass.

    CALL update_pres_delp_c( p_patch, p_new%pres_sfc,        &! in
                           & p_diag%pres_mc, p_diag%pres_ic, &! out
                           & p_diag%delp_c                   )! out

  END SUBROUTINE prepare_tracer_RK
  !--------

  !>
  !! Diagnose some pressure- and velocity-related quantities
  !! for the tracer transport scheme, under the assumption that
  !! a two time level semi-implicit time stepping scheme is
  !! used by the dynamical core.
  !!
  !! Note that after calling this subroutine, many of the pressure
  !! related components of the diagnostic state will be updated
  !! to time t+dt, while other components, e.g., temperature and
  !! geopotential, are still at time t.
  !!
  !! @par Revision History
  !! First version by Hui Wan, MPI-M (2010-02-02)
  !!
  SUBROUTINE prepare_tracer( p_patch, p_int_state,          &! in
                             p_now, p_new, palpha,          &! in
                             p_diag,                        &! inout
                             p_mflux_me, p_vn_traj,         &! out
                             p_mflux_ic, p_weta_traj,       &! out
                             p_delp_mc_now,                 &! out
                             p_pres_mc_now, p_pres_ic_now   )! out

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch
    TYPE(t_int_state),INTENT(in) :: p_int_state

    TYPE(t_hydro_atm_prog),INTENT(in)    :: p_now, p_new
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

    REAL(wp),INTENT(in)  :: palpha

    REAL(wp),INTENT(out) :: p_mflux_me   (nproma,nlev,  p_patch%nblks_e)
    REAL(wp),INTENT(out) :: p_vn_traj    (nproma,nlev,  p_patch%nblks_e)
    REAL(wp),INTENT(out) :: p_mflux_ic   (nproma,nlevp1,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_weta_traj  (nproma,nlevp1,p_patch%nblks_c)

    REAL(wp),INTENT(out) :: p_delp_mc_now(nproma,nlev  ,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_pres_mc_now(nproma,nlev  ,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_pres_ic_now(nproma,nlevp1,p_patch%nblks_c)

    ! Local variables

    REAL(wp) :: z_mdiv_int (nproma,nlevp1,p_patch%nblks_c)
    REAL(wp) :: z_mdiv     (nproma,nlev  ,p_patch%nblks_c)
    REAL(wp) :: z_aux_me   (nproma,nlev  ,p_patch%nblks_e)
    REAL(wp) :: z_ddt_psfc (nproma,       p_patch%nblks_c)

    REAL(wp) :: z_delp_me_now (nproma,nlev,p_patch%nblks_e)

    INTEGER  :: jb, jbs, is, ie, nblks_e, ircs
    REAL(wp) :: z1ma                       !< 1 - palpha

    !---

    nblks_e  = p_patch%nblks_int_e
    z1ma     = 1._wp - palpha

    ! Diagnose pressure-related quantities of time step n, including
    ! full- and half-level pressure values, and layer thickess, all
    ! computed at cell centers (mass points).

    CALL update_pres_delp_c( p_patch, p_now%pres_sfc,       &! in
                           & p_pres_mc_now, p_pres_ic_now,  &! out
                           & p_delp_mc_now                  )! out

    ! Layer thickness associated to edges. Needed for calculating mass flux

    CALL update_delp_e( p_patch, p_int_state, p_delp_mc_now, z_delp_me_now)

    !---------------------------------------------------------
    ! Mass fluxes
    !---------------------------------------------------------
    ! Note that here the variable "p_vn_traj" is used as a temporary array.
    ! The actual values passed out of the routine are calculated later.

    ircs = 2
    jbs  = p_patch%edges%start_blk(ircs,1)
!$OMP PARALLEL

    IF (itime_scheme == TWO_TL_SI) THEN

      ! For the two time level semi-implicit scheme, this means using the
      ! normal wind of time step n+alpha, and layer thickness of step n

!$OMP DO PRIVATE(jb,is,ie)
      DO jb = jbs, nblks_e
        CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, ircs )
        p_vn_traj(is:ie,:,jb) =    z1ma*p_now%vn(is:ie,:,jb) &
        &                      + palpha*p_new%vn(is:ie,:,jb)
      ENDDO
!$OMP END DO

    ELSE

!$OMP DO PRIVATE(jb,is,ie)
      DO jb = jbs, nblks_e
        CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, ircs )
        p_vn_traj(is:ie,:,jb) =  0.5_wp * ( p_now%vn(is:ie,:,jb)   &
        &                                 + p_new%vn(is:ie,:,jb) )
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL

    ! Get the horizontal and vertical mass fluxes that are consistent
    ! with the continuity equation. Here the vertical mass flux is
    ! eta-dot*(partial-p/partial-eta).

    CALL continuity( p_vn_traj, z_delp_me_now,      &! in
    &                p_patch, p_int_state, .TRUE.,  &! in
    &                z_mdiv, z_mdiv_int,            &! out, but not used here
    &                p_mflux_me,                    &! out
    &                z_ddt_psfc,                    &! out, but not used here
    &                p_mflux_ic                  )   ! out

    !---------------------------------------------------------
    ! Velocities for trajectory calculation
    !---------------------------------------------------------
    ! Normal wind at time step n+1/2

    jbs = p_patch%edges%start_blk(ircs,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
    DO jb = jbs,nblks_e
      CALL get_indices_e( p_patch, jb,jbs,nblks_e, is,ie, ircs )
      p_vn_traj(is:ie,:,jb) = 0.5_wp*( p_now%vn(is:ie,:,jb)  &
      &                               +p_new%vn(is:ie,:,jb) )

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Vertical velocity at time step n+1/2

    CALL continuity( p_vn_traj, z_delp_me_now,      &! in
    &                p_patch, p_int_state, .TRUE.,  &! in
    &                z_mdiv, z_mdiv_int,            &! out, but not used here
    &                z_aux_me,                      &! out, but not used here
    &                z_ddt_psfc,                    &! out, but not used here
    &                p_weta_traj                 )   ! out

    !---------------------------------------------------------
    ! Pressure values and layer thickness of time step n+1
    !---------------------------------------------------------
    ! The thickness will be used in tracer transport to ensure
    ! conservation of tracer mass.

    CALL update_pres_delp_c( p_patch, p_new%pres_sfc,        &! in
                           & p_diag%pres_mc, p_diag%pres_ic, &! out
                           & p_diag%delp_c                   )! out

  END SUBROUTINE prepare_tracer
  !--------
  !>
  !! Diagnose some pressure- and velocity-related quantities
  !! for the tracer transport scheme, under the assumption that
  !! the leapfrog time stepping scheme is used by the dynamical core.
  !!
  !! Note that after calling this subroutine, many of the pressure
  !! related components of the diagnostic state will be updated
  !! to time t+dt, while other components, e.g., temperature and
  !! geopotential, are still at time t.
  !!
  !! @par Revision History
  !! First version by Hui Wan, MPI-M (2010-07-27)
  !!
  SUBROUTINE prepare_tracer_leapfrog( p_patch,                       &! in
                                      p_old, p_now, p_new,           &! in
                                      p_diag,                        &! inout
                                      p_vn_traj,                     &! out
                                      p_delp_mc_old,                 &! out
                                      p_pres_mc_old, p_pres_ic_old   )! out

    TYPE(t_patch),TARGET,INTENT(in) :: p_patch

    TYPE(t_hydro_atm_prog),INTENT(in)    :: p_old, p_now, p_new
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

    REAL(wp),INTENT(out) :: p_vn_traj     (nproma,nlev,  p_patch%nblks_e)
    REAL(wp),INTENT(out) :: p_delp_mc_old (nproma,nlev  ,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_pres_mc_old (nproma,nlev  ,p_patch%nblks_c)
    REAL(wp),INTENT(out) :: p_pres_ic_old (nproma,nlevp1,p_patch%nblks_c)

    ! Diagnose pressure-related quantities of time step n-1, including
    ! full- and half-level pressure values and layer thickess, all
    ! computed at cell centers (i.e., mass points).
    ! Note that the computed values are given to dummy arguments
    ! rather than components of the diagnostic state vector.

    CALL update_pres_delp_c( p_patch, p_old%pres_sfc,       &! in
                           & p_pres_mc_old, p_pres_ic_old,  &! out
                           & p_delp_mc_old                  )! out

    ! Horizontal velocitiy = velocity at time step n
!$OMP PARALLEL WORKSHARE
    p_vn_traj(:,:,:) =  p_now%vn(:,:,:)
!$OMP END PARALLEL WORKSHARE

    ! Pressure values and layer thickness of time step n+1.
    ! The thickness will be used in tracer transport to ensure
    ! conservation of tracer mass.

    CALL update_pres_delp_c( p_patch, p_new%pres_sfc,        &! in
                           & p_diag%pres_mc, p_diag%pres_ic, &! out
                           & p_diag%delp_c                   )! out

  END SUBROUTINE prepare_tracer_leapfrog
  !--------
  !>
  !! Calculate virtual temperature, geopotential and vertical velocity
  !! using the newest values of temperature and horizontal wind.
  !!
  !! @par Revision History
  !! First version by Hui Wan, MPI-M (2010-02-04)
  !! slight modification regarding the theta2T conversion
  !! by Kristina Froehlich, DWD (2010-03-05)
  !!
  SUBROUTINE prepare_physics( ltheta_dyn, p_patch, p_int_state, &
  &                           p_ext_data, p_prog,               &! in
  &                           p_diag                          )  ! inout

    LOGICAL,INTENT(IN) :: ltheta_dyn

    TYPE(t_patch),TARGET,  INTENT(in) :: p_patch
    TYPE(t_int_state),     INTENT(in) :: p_int_state
    TYPE(t_external_data), INTENT(IN) :: p_ext_data !< external data

    TYPE(t_hydro_atm_prog),INTENT(inout) :: p_prog
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

! Local array bounds:

    INTEGER :: nblks_c       ! number of blocks for cells / edges

! Local scalars:

    INTEGER  :: jb          ! loop indices
    INTEGER  :: is,ie,jbs

    nblks_c   = p_patch%nblks_int_c

    ! Pressure at half- and full-levels;
    ! auxiliary variables %delp_c, %rdelp_c, %lnp_ic, %rdalpha_c

    CALL update_pres( p_prog, p_patch, p_int_state, p_diag )

    jbs = p_patch%cells%start_blk(2,1)
    ! If using theta as prognostic variable, convert it to temperauture
    ! KF suggestion: in order to avoid multiple pressure calculations
    ! put conversion equation right here
    IF (ltheta_dyn) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
      DO jb = jbs,nblks_c
         CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, 2)

         p_prog%temp(is:ie,:,jb) =   p_prog%   theta(is:ie,:,jb) &
&                                   *p_diag% rdelp_c(is:ie,:,jb) &
&                  *EXP(rd_o_cpd*LOG(p_diag% pres_mc(is:ie,:,jb)/p0ref))
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
     ENDIF !theta dynamics

    ! Virtual temperature and geopotential
    CALL update_tempv_geopot( p_prog, p_patch, p_ext_data, p_diag )

    ! Vertical velocity
    CALL update_omega( p_prog%vn, p_diag%delp_e, p_diag%pres_mc,  &! in
                     & p_patch, p_int_state,                      &! in
                     & p_diag%wpres_mc                          )  ! inout

  END SUBROUTINE prepare_physics
  !--------
  !>
  !! Calculate virtual temperature, geopotential and vertical velocity
  !! using the newest values of temperature and horizontal wind.
  !!
  !! @par Revision History
  !! First version by Hui Wan, MPI-M (2010-02-04)
  !! slight modification regarding the theta2T conversion
  !! by Kristina Froehlich, DWD (2010-03-05)
  !!
  SUBROUTINE prepare_echam_phy( ltheta_dyn, p_patch, p_int_state, &! in
  &                             p_ext_data, p_prog,               &! in
  &                             p_diag                          )  ! inout

    LOGICAL,INTENT(IN) :: ltheta_dyn

    TYPE(t_patch),TARGET,  INTENT(in) :: p_patch
    TYPE(t_int_state),     INTENT(in) :: p_int_state
    TYPE(t_external_data), INTENT(IN) :: p_ext_data !< external data

    TYPE(t_hydro_atm_prog),INTENT(inout) :: p_prog
    TYPE(t_hydro_atm_diag),INTENT(inout) :: p_diag

    INTEGER :: nblks_c       ! number of blocks for cells / edges
    INTEGER  :: jb           ! loop indices
    INTEGER  :: is,ie,jbs

    !------------
    nblks_c = p_patch%nblks_int_c

    ! Pressure at half- and full-levels;
    ! auxiliary variables %delp_c, %rdelp_c, %lnp_ic, %rdalpha_c

    CALL update_pres( p_prog, p_patch, p_int_state, p_diag )

    ! If using theta as prognostic variable, convert it to temperauture
    ! KF suggestion: in order to avoid multiple pressure calculations
    ! put conversion equation right here
    IF (ltheta_dyn) THEN
      jbs = p_patch%cells%start_blk(2,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie)
      DO jb = jbs,nblks_c
         CALL get_indices_c( p_patch, jb,jbs,nblks_c, is,ie, 2)

         p_prog%temp(is:ie,:,jb) =   p_prog%   theta(is:ie,:,jb) &
&                                   *p_diag% rdelp_c(is:ie,:,jb) &
&                  *EXP(rd_o_cpd*LOG(p_diag% pres_mc(is:ie,:,jb)/p0ref))
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF !theta dynamics

    ! Virtual temperature and geopotential
    CALL update_tempv_geopot( p_prog, p_patch, p_ext_data, p_diag, &
                            & opt_lgeop_wrt_sfc=.TRUE. )

    ! Vertical velocity
    CALL update_omega( p_prog%vn, p_diag%delp_e, p_diag%pres_mc,  &! in
                     & p_patch, p_int_state,                      &! in
                     & p_diag%wpres_mc                          )  ! inout

    ! Relative vorticity at cell centers
    SELECT CASE (p_patch%cell_type)
    CASE (3)
      CALL rot_vertex( p_prog%vn, p_patch, p_int_state, p_diag%rel_vort )
      CALL sync_patch_array( SYNC_V, p_patch, p_diag%rel_vort )
      CALL verts2cells_scalar( p_diag%rel_vort, p_patch,   &
                             & p_int_state%verts_aw_cells, &
                             & p_diag%rel_vort_c          )
    CASE (6)
      CALL rot_vertex( p_prog%vn, p_patch, p_int_state, p_diag%rel_vort )
      CALL sync_patch_array( SYNC_V, p_patch, p_diag%rel_vort )
      CALL verts2edges_scalar( p_diag%rel_vort, p_patch, &
                             & p_int_state%tria_aw_rhom, &
                             & p_diag%rel_vort_e )
      CALL sync_patch_array( SYNC_E, p_patch, p_diag%rel_vort_e )
      CALL edges2cells_scalar( p_diag%rel_vort_e, p_patch, &
                             & p_int_state%r_aw_c,         &
                             & p_diag%rel_vort_c   )
    END SELECT

  END SUBROUTINE prepare_echam_phy
  !--------

END MODULE mo_ha_dtp_interface
