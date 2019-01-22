!>
!! Contains the interface needed to call AWI FEM sea ice model
!! as well as advection and interpolation routines.
!!
!! @par Revision History
!! Developed  by Einar Olason (2013)
!! Restructured by Vladimir Lapin (2015)
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

MODULE mo_ice_fem_interface
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_momentum, timer_ice_interp, timer_ice_advection

! USE mo_grid_config,         ONLY: l_limited_area, n_dom   ! for now sea-ice works on global domain-only
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary

  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_dynamics_config,     ONLY: nold

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_sea_ice_nml,         ONLY: i_ice_advec
  USE mo_ice_fem_advection,   ONLY: fct_ice_solve, ice_TG_rhs
  USE mo_math_types,          ONLY: t_cartesian_coordinates
  USE mo_math_utilities,      ONLY: gvec2cvec, cvec2gvec
  USE mo_scalar_product,      ONLY: map_cell2edges_3D, map_edges2cell_3D
  USE mo_icon_interpolation_scalar, ONLY: verts2cells_scalar
  USE mo_ice_fem_interpolation, ONLY: map_edges2verts, map_verts2edges,                 &
                                      gvec2cvec_c_2d, cvec2gvec_c_2d,                   &
                                      rotate_cvec_v, gvec2cvec_v_fem, cvec2gvec_v_fem,  &
                                      cells2verts_scalar_seaice

  IMPLICIT NONE

  PUBLIC  :: ice_fem_interface

  PUBLIC  :: ice_fem_init_vel_restart
  PUBLIC  :: ice_fem_update_vel_restart

  PRIVATE :: map_icon2fem_vec
  PRIVATE :: map_fem2icon_vec
  PRIVATE :: map_icon2fem_scalar
  PRIVATE :: map_fem2icon_scalar

  CHARACTER(len=12)           :: str_module    = 'IceFemUtils'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1         ! Level of detail for 1 line debug

CONTAINS

!------------------------------------------------------------------------
!
!>  Wrapper for the call to the AWI FEM ice model
!!  We first remap the neccesary inputs, then call the momentum solver (EVPdynamics)
!!  and map the resulting velocity onto edges and cell centres.
!!
!! @par Revision History
!! Developed by Einar Olason, MPI-M (2013-06-05)
!! Modified   by Vladimir Lapin (2015)
!
  SUBROUTINE ice_fem_interface( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff)

    USE mo_ice_fem_types,     ONLY: sigma11, sigma12, sigma22
    USE mo_ice_fem_evp,       ONLY: EVPdynamics
!    USE mo_ice_fem_evp_old,  ONLY: EVPdynamics_old ! non-optimized, original version of the solver

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    TYPE(t_patch), POINTER :: p_patch
    REAL(wp),      POINTER :: ssh(:,:) ! sea surface height (input only)         [m]

!--------------------------------------------------------------------------------------------------
    p_patch => p_patch_3D%p_patch_2D(1)
    ssh     => p_os%p_prog(nold(1))%h(:,:)

    IF (ltimer) CALL timer_start(timer_ice_momentum)

! Initialization of  u_ice, v_ice in case of a restart was moved to mo_hydro_ocean_run

!--------------------------------------------------------------------------------------------------
! Interpolate and/or copy ICON variables to FEM variables
!--------------------------------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Map scalars to vertices. Obtain: m_ice, m_snow, a_ice, elevation
    CALL map_icon2fem_scalar(p_patch, p_ice, ssh)
    ! Map vectors to vertices on the rotated-pole grid. Obtain: stress_atmice_x, stress_atmice_y, u_w, v_w
    CALL map_icon2fem_vec(p_patch_3D, p_os, atmos_fluxes, p_op_coeff)

    IF (ltimer) CALL timer_stop(timer_ice_interp)

!--------------------------------------------------------------------------------------------------
! Call FEM EVP solver
!--------------------------------------------------------------------------------------------------
    sigma11=0._wp; sigma12=0._wp; sigma22=0._wp
!    CALL EVPdynamics_old
    CALL EVPdynamics

!--------------------------------------------------------------------------------------------------
! FCT advection on FEM grid. Advection on ICON grid is done in ice_slow_interface
!--------------------------------------------------------------------------------------------------
    IF (i_ice_advec == 1) THEN
        IF (ltimer) CALL timer_start(timer_ice_advection)

        call ice_TG_rhs
        call fct_ice_solve

        IF (ltimer) CALL timer_stop(timer_ice_advection)
    ENDIF

!--------------------------------------------------------------------------------------------------
! Interpolate and/or copy FEM variables back to ICON variables
!--------------------------------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Rotate and interpolate ice velocities back to ICON grid
    CALL map_fem2icon_vec( p_patch_3D, p_ice, p_op_coeff )
    ! If advection is on FEM grid, interp ice scalars back to ICON grid
    IF (i_ice_advec == 1) THEN
        CALL map_fem2icon_scalar( p_patch, p_ice )
    ENDIF

    IF (ltimer) CALL timer_stop(timer_ice_interp)

! Check to make sure that u_ice, v_ice are written to the restart file was moved to mo_hydro_ocean_run

    IF (ltimer) CALL timer_stop(timer_ice_momentum)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('femIWrap: ice_u' , p_ice%u, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('femIWrap: ice_v' , p_ice%v, str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_fem_interface

  !-------------------------------------------------------------------------
  !
  !> Initialize u_ice, v_ice in case of a restart
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_init_vel_restart(p_patch, p_ice)
    USE mo_ice_fem_types,       ONLY: u_ice, v_ice

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),       INTENT(IN)    :: p_ice

    ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: jk, jb, jv

  !-------------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
!            ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files
!            ! now, so this does not need to be done every timestep, only after restart file is
!            ! read.
         u_ice(jk) = p_ice%u_prog(jv,jb)
         v_ice(jk) = p_ice%v_prog(jv,jb)
      END DO
    END DO

  END SUBROUTINE ice_fem_init_vel_restart

  !-------------------------------------------------------------------------
  !
  !> Update p_ice%u_prog with last u_ice, v_ice values before writing a restart file
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_update_vel_restart(p_patch, p_ice)
    USE mo_ice_fem_types,       ONLY: u_ice, v_ice

    TYPE(t_patch), TARGET, INTENT(in)       :: p_patch
    TYPE(t_sea_ice),       INTENT(inout)    :: p_ice

    ! Local variables
    ! Patch and ranges
    TYPE(t_subset_range), POINTER :: all_verts
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: jk, jb, jv

  !-------------------------------------------------------------------------

    all_verts => p_patch%verts%all

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
!        ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files now,
!        ! so this does not need to be done every timestep, only before restart file is written.
        p_ice%u_prog(jv,jb) = u_ice(jk)
        p_ice%v_prog(jv,jb) = v_ice(jk)

      END DO
    END DO

  END SUBROUTINE ice_fem_update_vel_restart

  !-------------------------------------------------------------------------
  ! #vla, development note for:
  ! ---------- map_icon2fem_vec and map_fem2icon_vec ----------
  ! Old versoin contained a polar singularity because a cartesian vector was
  ! converted to geographical coords form form and then rotated.
  ! FIX: first rotate the cartesian vector to the rotated pole coordinates
  ! and only then convert to geographic form.
  ! #vla old code was removed during clean-up, 05-2017
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> 1) Interpolate wind stress from cell centers to vertices
  !> 2) Rotate ocean velocities (available on the dual grid, i.e. vertices)
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-07-17)
  !
  SUBROUTINE map_icon2fem_vec( p_patch_3D, p_os, atmos_fluxes, p_op_coeff )

    USE mo_ice_fem_icon_init, ONLY: rot_mat_3D
    USE mo_ice_fem_types,     ONLY: u_w, v_w, stress_atmice_x, stress_atmice_y!, u_ice, v_ice !, a_ice

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    TYPE(t_patch), POINTER :: p_patch

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: p_tau_n_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: tau_n(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
!    TYPE(t_cartesian_coordinates) :: p_vn_dual    (nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)

!--------------------------------------------------------------------------------------------------
    p_patch => p_patch_3D%p_patch_2D(1)

    !**************************************************************
    ! (1) Convert lat-lon wind stress to cartesian coordinates
    !**************************************************************
    CALL gvec2cvec_c_2d(p_patch_3D, atmos_fluxes%stress_x, atmos_fluxes%stress_y, p_tau_n_c)

    !**************************************************************
    ! (2) Interpolate 3D wind stress from cell centers to edges
    !**************************************************************
#ifdef NAGFOR
    tau_n = 0.0_wp
#endif
    CALL map_cell2edges_3D(p_patch_3D, p_tau_n_c, tau_n ,p_op_coeff, 1)
    CALL sync_patch_array(SYNC_E, p_patch, tau_n)

    !**************************************************************
    ! (3) Interpolate 3D wind stress from edges to vertices
    !**************************************************************
#ifdef NAGFOR
    p_tau_n_dual(:,:)%x(1) = 0.0_wp
    p_tau_n_dual(:,:)%x(2) = 0.0_wp
    p_tau_n_dual(:,:)%x(3) = 0.0_wp
#endif
    CALL map_edges2verts(p_patch, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual)
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(1))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(2))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(3))

    !**************************************************************
    ! (4) Rotate the vectors onto the rotated grid
    !     + convert back to geographic coordinates
    !**************************************************************
    ! atmospheric stress
    CALL rotate_cvec_v(p_patch, p_tau_n_dual, rot_mat_3D, p_tau_n_dual_fem)
    CALL cvec2gvec_v_fem(p_patch, p_tau_n_dual_fem, stress_atmice_x, stress_atmice_y)
    ! ocean velocities
    CALL rotate_cvec_v(p_patch, p_os%p_diag%p_vn_dual(:,1,:), rot_mat_3D, p_vn_dual_fem)
    CALL cvec2gvec_v_fem(p_patch, p_vn_dual_fem, u_w, v_w)

  END SUBROUTINE map_icon2fem_vec

  !-------------------------------------------------------------------------
  !> 1) Rotate back ice velocities (available on the FEM grid, i.e. vertices)
  !> 2) Interpolate to cell centers and convert back to geographic coordinates
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-07-17)
  !
  SUBROUTINE map_fem2icon_vec( p_patch_3D, p_ice, p_op_coeff )

    USE mo_ice_fem_icon_init, ONLY: rot_mat_3D!, pollon, pollat
    USE mo_ice_fem_types,          ONLY: u_ice, v_ice

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    TYPE(t_patch), POINTER :: p_patch

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: p_vn_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_c_3D(nproma,1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

!--------------------------------------------------------------------------------------------------
    p_patch => p_patch_3D%p_patch_2D(1)

#ifdef NAGFOR
    p_vn_c_3D(:,:,:)%x(1) = 0.0_wp
    p_vn_c_3D(:,:,:)%x(2) = 0.0_wp
    p_vn_c_3D(:,:,:)%x(3) = 0.0_wp
    p_vn_dual_fem(:,:)%x(1) = 0.0_wp
    p_vn_dual_fem(:,:)%x(2) = 0.0_wp
    p_vn_dual_fem(:,:)%x(3) = 0.0_wp
#endif

    !**************************************************************
    ! (1) Rotate ice vels to ICON variables + convert to cc
    !**************************************************************
    ! Convert the lat-lon vectors to 3d cartesian
    CALL gvec2cvec_v_fem(p_patch, u_ice, v_ice, p_vn_dual_fem)
    ! Rotate the vectors back onto the ICON grid
    CALL rotate_cvec_v(p_patch, p_vn_dual_fem, TRANSPOSE(rot_mat_3D(:,:)), p_vn_dual)

    !**************************************************************
    ! (2) Interpolate ice velocities to edges for advection
    !**************************************************************
    CALL map_verts2edges(p_patch, p_vn_dual, p_op_coeff%edge2vert_coeff_cc_t, p_ice%vn_e)
    CALL sync_patch_array(SYNC_E, p_patch, p_ice%vn_e)

    !**************************************************************
    ! (3) ... and cells for drag calculation and output
    !**************************************************************
    CALL map_edges2cell_3D( p_patch_3D, RESHAPE(p_ice%vn_e, (/ nproma, 1, p_patch%nblks_e /)), &
      &   p_op_coeff, p_vn_c_3D, 1, 1)

    !**************************************************************
    ! (4) Convert back to geographic coordinates
    !**************************************************************
    CALL cvec2gvec_c_2d(p_patch_3D, p_vn_c_3D(:,1,:), p_ice%u, p_ice%v)
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%u)
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%v)

  END SUBROUTINE map_fem2icon_vec

  !-------------------------------------------------------------------------
  !> 1) Interpolate ice scalars from vertices to cell centers
  !> 2) Resahpe result to get vars on FEM grid
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2017-05-05)
  !
  SUBROUTINE map_icon2fem_scalar(p_patch, p_ice, ssh)

    USE mo_ice_fem_icon_init, ONLY: c2v_wgt
    USE mo_ice_fem_types,     ONLY: m_ice, m_snow, a_ice, elevation

    TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch
    TYPE(t_sea_ice),        INTENT(IN)  :: p_ice
    REAL(wp),DIMENSION(nproma,p_patch%alloc_cell_blocks), INTENT(IN) :: ssh

    ! Temporary variables/buffers
    REAL(wp), DIMENSION (nproma, p_ice%kice, p_patch%nblks_v)   :: buffy_array
    REAL(wp), DIMENSION (nproma*p_patch%nblks_v)                :: buffy

    buffy_array(:,:,:)   = 0.0_wp

    ! Interpolate tracers to vertices
    CALL cells2verts_scalar_seaice( p_ice%hi*MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    m_ice = buffy(1:SIZE(m_ice))

    CALL cells2verts_scalar_seaice( p_ice%conc, p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    a_ice = buffy(1:SIZE(a_ice))

    CALL cells2verts_scalar_seaice( p_ice%hs*MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    m_snow= buffy(1:SIZE(m_snow))

    ! Interpolate SSH to vertices
    CALL cells2verts_scalar_seaice( RESHAPE(ssh,(/ nproma,1,p_patch%alloc_cell_blocks /)), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    elevation = buffy(1:SIZE(elevation) )

  END SUBROUTINE map_icon2fem_scalar


  !-------------------------------------------------------------------------
  !> 1) Transform ice vars on FEM grid back to ICON structures
  !> 2) Interpolate from vertices to cell centers
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2017-05-05)
  !
  SUBROUTINE map_fem2icon_scalar(p_patch, p_ice)

    USE mo_ice_fem_icon_init, ONLY: v2c_wgt
    USE mo_ice_fem_types,          ONLY: m_ice, m_snow, a_ice

    TYPE(t_patch), TARGET,  INTENT(INOUT)   :: p_patch
    TYPE(t_sea_ice),        INTENT(INOUT)   :: p_ice

    ! Local variables
    TYPE(t_subset_range), POINTER :: all_verts
    INTEGER  :: i_startidx_v, i_endidx_v, jb, jk, jv

    ! Temporary variables/buffers
    REAL(wp), DIMENSION (nproma, p_ice%kice, p_patch%nblks_v) :: m_ice_buff, m_snow_buff, a_ice_buff

    all_verts => p_patch%verts%all

    m_ice_buff(:,:,:)   = 0.0_wp
    m_snow_buff(:,:,:)  = 0.0_wp
    a_ice_buff(:,:,:)   = 0.0_wp

    !**************************************************************
    ! (1) Transform ice vars on FEM grid back to ICON structures
    !**************************************************************
    jk=0 ! position index for vars on FEM grid
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v,i_endidx_v

        jk = jk+1
        m_ice_buff(jv,1,jb)   = m_ice(jk)
        m_snow_buff(jv,1,jb)  = m_snow(jk)
        a_ice_buff(jv,1,jb)   = a_ice(jk)

      END DO
    END DO
    ! does the same as above
!    m_ice_buff  = RESHAPE(m_ice, SHAPE(m_ice_buff))
!    m_snow_buff = RESHAPE(m_snow, SHAPE(m_snow_buff))
!    a_ice_buff  = RESHAPE(a_ice, SHAPE(a_ice_buff))

    !**************************************************************
    ! (2) Interpolate FEM ice variables to cells
    !**************************************************************
    CALL verts2cells_scalar ( m_ice_buff, p_patch, v2c_wgt, p_ice%vol ) ! multiplied by cell-area below
    CALL verts2cells_scalar ( m_snow_buff, p_patch, v2c_wgt,p_ice%vols) ! multiplied by cell-area below
    CALL verts2cells_scalar ( a_ice_buff, p_patch, v2c_wgt, p_ice%conc )

    CALL sync_patch_array   ( SYNC_C, p_patch, p_ice%vol )
    CALL sync_patch_array   ( SYNC_C, p_patch, p_ice%vols )
    CALL sync_patch_array   ( SYNC_C, p_patch, p_ice%conc )

    !**************************************************************
    ! (3) Calculate ICON ice-variables
    !**************************************************************
!ICON_OMP_WORKSHARE
        WHERE ( p_ice%conc(:,1,:) > 0._wp )
          ! New ice and snow thickness
          p_ice%hi  (:,1,:) = p_ice%vol (:,1,:) / p_ice%conc(:,1,:)
          p_ice%hs  (:,1,:) = p_ice%vols(:,1,:) / p_ice%conc(:,1,:)

          ! multiply vol and vols by cell-area, as by definition
          p_ice%vol (:,1,:) = p_ice%vol (:,1,:) * p_patch%cells%area(:,:)
          p_ice%vols(:,1,:) = p_ice%vols(:,1,:) * p_patch%cells%area(:,:)
        ELSEWHERE
          p_ice%hi  (:,1,:) = 0._wp
          p_ice%hs  (:,1,:) = 0._wp
          p_ice%vol (:,1,:) = 0._wp
          p_ice%vols(:,1,:) = 0._wp
          p_ice%conc(:,1,:) = 0._wp
        ENDWHERE
!ICON_OMP_END_WORKSHARE

  END SUBROUTINE map_fem2icon_scalar

END MODULE mo_ice_fem_interface
