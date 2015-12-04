!>
!! Contains the implementation of interpolation needed for the FEM ice model.
!!
!! @par Revision History
!! Developed  by Einar Olason (2013).
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
!-----------------------------------------------------------------------
!
!  ! Wrapper around the AWI FEM ice model as well as advection and drag calculations
!  ! Also averaging and interpolation routines and routines needed to compute the coefficients
!  ! therein
!
!-----------------------------------------------------------------------
!
MODULE mo_ice_fem_utils
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
!  USE mo_exception,           ONLY: message
!  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_momentum,                  &
    &                               timer_ice_interp, timer_ice_advection

  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_io_config,           ONLY: n_checkpoints

!  USE mo_grid_config,         ONLY: l_limited_area, n_dom   ! restrict sea-ice model to the global domain for the time being
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary!, sea

  USE mo_advection_utils,     ONLY: laxfr_upflux
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_ocean_math_operators,  ONLY: div_oce_3D
  USE mo_dynamics_config,     ONLY: nold

  USE mo_scalar_product,      ONLY: map_cell2edges_3D, map_edges2cell_3D
  USE mo_math_constants,      ONLY: rad2deg, deg2rad

  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_patch_array_mult

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes!, t_sfc_flx, t_atmos_for_ocean
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec, cvec2gvec !, &
!    &                               rotate_latlon, rotate_latlon_vec!, disp_new_vect
!  USE mo_icon_interpolation_scalar, ONLY: cells2verts_scalar
  USE mo_icon_to_fem_interpolation, ONLY: map_edges2verts, map_verts2edges,                 &
                                          gvec2cvec_c_2d, cvec2gvec_c_2d,                   &
                                          rotate_cvec_v, gvec2cvec_v_fem, cvec2gvec_v_fem,  &
                                          map_verts2edges_einar, map_edges2verts_einar,     &
                                          cells2verts_scalar_seaice
  USE mo_ice_fem_init,        ONLY: ice_fem_update_vel_for_restart, ice_fem_update_vel_for_restart


  IMPLICIT NONE

  PUBLIC  :: fem_ice_wrap
  PUBLIC  :: ice_advection, ice_advection_vla
  PUBLIC  :: ice_ocean_stress

  PRIVATE :: upwind_hflux_ice
  PRIVATE :: intrp_to_fem_grid_vec
  PRIVATE :: intrp_from_fem_grid_vec
!  PRIVATE :: intrp_to_fem_grid_vec_old
!  PRIVATE :: intrp_from_fem_grid_vec_old

  CHARACTER(len=12)           :: str_module    = 'IceFem'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1         ! Level of detail for 1 line debug

CONTAINS

!------------------------------------------------------------------------
!
!
!>
!!  Wrapper for the call to the AWI FEM ice model
!!  We first remap the neccesary inputs, then call the momentum solver (EVPdynamics) and map the
!!  resulting velocity onto edges and cell centres.
!!
!!
!! @par Revision History
!! Developed by Einar Olason, MPI-M (2013-06-05)
!!

  SUBROUTINE fem_ice_wrap( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff, jstep)

    USE mo_ice_fem_init,        ONLY: c2v_wgt
    USE mo_ice,      ONLY: u_ice, v_ice, m_ice, a_ice, m_snow, u_w, v_w, &
      &   elevation, sigma11, sigma12, sigma22
    USE mo_ice_evp,  ONLY: EVPdynamics
    USE mo_ice_evp_omp,  ONLY: EVPdynamics_omp

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff
    INTEGER,                  INTENT(IN)     :: jstep

    ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
!    TYPE(t_subset_range), POINTER :: all_cells, all_verts

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jv, jb, jk
    INTEGER  :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers
    REAL(wp) :: buffy(nproma*p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp) :: buffy_array(nproma,p_ice%kice,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp) :: tmp2(2)!, tmp_x, tmp_y, delu, u_change

    REAL(wp), DIMENSION (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: &
			& atm_u, atm_v, oce_u, oce_v

    ! TODO: There are too many calls to cvec2gvec and gvec2cvec ... but premature optimisation is the
    ! root of all evil ...

!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
! Set up patch and ranges
!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
!    all_cells => p_patch%cells%all
!    all_verts => p_patch%verts%all

    ! Initialize u_ice, v_ice in case of a restart
!    CALL ice_fem_init_vel(p_patch, p_ice)
    ! No need to do this every timestep; only after a restart file is read in.
    ! So, was moved to after <prepare_ho_stepping> in <ocean_model>

    IF (ltimer) CALL timer_start(timer_ice_momentum)

!--------------------------------------------------------------------------------------------------
! Interpolate and/or copy ICON variables to FEM variables
!--------------------------------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Interpolate tracers to vertices
    buffy_array = 0._wp
    ! TODO: Replace hi/conc to himean
!    CALL cells2verts_scalar( p_ice%hi/MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL cells2verts_scalar_seaice( p_ice%hi*MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    m_ice = buffy(1:SIZE(m_ice) )

    CALL cells2verts_scalar_seaice( p_ice%conc, p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    a_ice = buffy(1:SIZE(a_ice) )

!    CALL cells2verts_scalar( p_ice%hs/MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL cells2verts_scalar_seaice( p_ice%hs*MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    m_snow= buffy(1:SIZE(m_snow))

    ! Interpolate SSH to vertices
    CALL cells2verts_scalar_seaice( RESHAPE(p_os%p_prog(nold(1))%h(:,:), &
      & (/ nproma, 1, p_patch%alloc_cell_blocks /)), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    elevation = buffy(1:SIZE(elevation) )

!	atm_u=-10._wp
!	atm_v=0._wp
!	atm_u=0._wp
!	atm_v=-10._wp

!    DO jb = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!      DO jc = i_startidx_c, i_endidx_c
!        tmp2(1) = p_patch%cells%center(jc,jb)%lon
!        tmp2(2) = p_patch%cells%center(jc,jb)%lat
!	atmos_fluxes%stress_x(jc,jb) = -SIN(tmp2(1)) / ( SIN(tmp2(1))**2 + (cos(tmp2(1))*sin(tmp2(2)))**2 )**0.5
!	atmos_fluxes%stress_y(jc,jb) =  cos(tmp2(1))*sin(tmp2(2)) / ( SIN(tmp2(1))**2 + (cos(tmp2(1))*sin(tmp2(2)))**2 )**0.5

!    atmos_fluxes%stress_x(jc,jb) = cos(pollat)*sin(tmp2(1)-pollon)
!    atmos_fluxes%stress_y(jc,jb) = cos(tmp2(2))*sin(pollat) - sin(tmp2(2))*cos(pollat)*cos(tmp2(1)-pollon)
!      ENDDO
!    ENDDO

!	p_os%p_diag%p_vn_dual(:,1,:)%x(1)=0._wp
!	p_os%p_diag%p_vn_dual(:,1,:)%x(2)=0._wp
!	p_os%p_diag%p_vn_dual(:,1,:)%x(3)=0._wp
    
    ! Interpolate wind stress to vertices and  onto the rotated-pole lat-lon grid
!    CALL intrp_to_fem_grid_vec_old( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff )
    CALL intrp_to_fem_grid_vec( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff )

    IF (ltimer) CALL timer_stop(timer_ice_interp)

!--------------------------------------------------------------------------------------------------
! Call FEM EVP
!--------------------------------------------------------------------------------------------------

    sigma11=0._wp; sigma12=0._wp; sigma22=0._wp
!    CALL EVPdynamics
    CALL EVPdynamics_omp

!--------------------------------------------------------------------------------------------------
! Post-processing: Copy FEM variables back to ICON variables
!--------------------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Rotate and interpolate ice velocities back to the original ICON grid
!    CALL intrp_from_fem_grid_vec_old( p_patch_3D, p_ice, p_op_coeff )
    CALL intrp_from_fem_grid_vec( p_patch_3D, p_ice, p_op_coeff )

    ! Initialize u_ice, v_ice in case of a restart
    ! This does not need to be done every timestep, only after restart file is read.
    ! TODO: move to after <read_restart_files> in <ocean_model>
    !       ie into <prepare_ho_stepping>
    !!!!!!!!!!!!!!! DONE !!!!!!!!!!!!!!!
!    CALL ice_fem_init_vel(p_patch, p_ice)

      ! write a restart or checkpoint file
      IF (MOD(jstep,n_checkpoints())==0) THEN
        CALL ice_fem_update_vel_for_restart(p_patch, p_ice)
      END IF

    IF (ltimer) CALL timer_stop(timer_ice_interp)
    IF (ltimer) CALL timer_stop(timer_ice_momentum)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('femIWrap: ice_u' , p_ice%u, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('femIWrap: ice_v' , p_ice%v, str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE fem_ice_wrap

!-------------------------------------------------------------------------------
  !>
  !! First order upwind scheme for horizontal tracer advection
  !!
  !! Calculation of horizontal tracer fluxes using the first
  !! order Godunov method.
  !!
  !! @par Revision History
  !! Developed by L.Bonaventura  (2004).
  !! Adapted to new grid structure by L. Bonaventura, MPI-M, August 2005.
  !! Modification by Daniel Reinert, DWD (2010-02-09)
  !! - transferred to separate subroutine
  !! Modification by Stephan Lorenz, MPI (2010-09-06)
  !! - adapted to hydrostatic ocean core
  !! Modification by Einar Olason, MPI (2013-07-30)
  !! - adapted for the FEM ice model
  !!
  SUBROUTINE upwind_hflux_ice( p_patch_3D, pvar_c, pvn_e, pupflux_e, opt_slev, opt_elev )

    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
    REAL(wp), INTENT(IN)              :: pvar_c   (:,:,:) !< advected cell centered variable
    REAL(wp), INTENT(IN)              :: pvn_e    (:,:)   !< normal velocity on edges
    REAL(wp), INTENT(OUT)             :: pupflux_e(:,:,:) !< variable in which the upwind flux is stored
    INTEGER, INTENT(in), OPTIONAL     :: opt_slev    ! optional vertical start level
    INTEGER, INTENT(in), OPTIONAL     :: opt_elev    ! optional vertical end level

    ! local variables
    INTEGER, DIMENSION(:,:,:), POINTER :: iilc,iibc  ! pointer to line and block indices
    INTEGER  :: slev, elev
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: je, jk, jb         !< index of edge, vert level, block
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_patch), POINTER         :: p_patch 

    !-----------------------------------------------------------------------  
    p_patch         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => p_patch%edges%in_domain
    !-----------------------------------------------------------------------

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF
    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = UBOUND(pvar_c,2)
    END IF
    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    ! for no-slip boundary conditions, boundary treatment for tracer (zero at leteral walls)
    !is implicit done via velocity boundary conditions
    !
    ! line and block indices of two neighboring cells
    iilc => p_patch%edges%cell_idx
    iibc => p_patch%edges%cell_blk

    ! loop through all patch edges (and blocks)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block

      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = slev, elev
#else
!CDIR UNROLL=6
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
#endif
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
!          IF ( p_patch_3D%lsm_e(je,1,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  &
            &  laxfr_upflux( pvn_e(je,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
!          ELSE
!            pupflux_e(je,jk,jb) = 0.0_wp
!          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE upwind_hflux_ice

  !-------------------------------------------------------------------------
  !
  !> Advection of sea ice and snow on ice
  !! This uses the upwind_hflux_ice routine and the ocean's div_oce_3D routine to do upwind
  !! advection of the relevant variables.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_advection( p_patch_3D, p_op_coeff, p_ice )
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_operator_coeff),   INTENT(IN)    :: p_op_coeff
    TYPE(t_sea_ice),          INTENT(INOUT) :: p_ice

    ! Local variables
    ! Patch and range
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: cells_in_domain

    ! Indexing
    INTEGER  :: jk, jb, jc
    INTEGER  :: i_startidx_c, i_endidx_c

    ! Temporary variables/buffers
    REAL(wp) :: z_adv_flux_h (nproma,p_ice%kice,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: flux_hi  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_conc(nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_hs  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

!--------------------------------------------------------------------------------------------------

!    IF (ltimer) CALL timer_start(timer_ice_advection)

    p_patch => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

!--------------------------------------------------------------------------------------------------
! Do advection
!--------------------------------------------------------------------------------------------------

    !upwind estimate of tracer flux
    CALL upwind_hflux_ice( p_patch_3D, p_ice%vol,  p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), p_patch, p_op_coeff%div_coeff, flux_hi  (:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    CALL upwind_hflux_ice( p_patch_3D, p_ice%conc, p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), p_patch, p_op_coeff%div_coeff, flux_conc(:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    CALL upwind_hflux_ice( p_patch_3D, p_ice%vols, p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), p_patch, p_op_coeff%div_coeff, flux_hs  (:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    DO jk = 1,p_ice%kice
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
            p_ice%vol (jc,jk,jb)= p_ice%vol (jc,jk,jb)-dtime*flux_hi  (jc,jk,jb)
            p_ice%conc(jc,jk,jb)= p_ice%conc(jc,jk,jb)-dtime*flux_conc(jc,jk,jb)
            p_ice%vols(jc,jk,jb)= p_ice%vols(jc,jk,jb)-dtime*flux_hs  (jc,jk,jb)
          ENDIF
          ! TODO ram - remove p_patch%cells%area(jc,jb) and test
          ! See also thermodyn/mo_sea_ice.f90
          IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN 
            p_ice%hi(jc,jk,jb) = p_ice%vol (jc,jk,jb)   &
              &         /( p_ice%conc(jc,jk,jb)*p_patch%cells%area(jc,jb) )
            p_ice%hs(jc,jk,jb) = p_ice%vols(jc,jk,jb)   &
              &         /( p_ice%conc(jc,jk,jb)*p_patch%cells%area(jc,jb) )
          ENDIF
        END DO
      END DO
    END DO

!--------------------------------------------------------------------------------------------------
! Sync results
!--------------------------------------------------------------------------------------------------

    CALL sync_patch_array(SYNC_C, p_patch, p_ice%vol (:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%vols(:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%conc(:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%hs  (:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%hi  (:,:,:))

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('ice_adv: vol ice'  , p_ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: vol snow' , p_ice%vols, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: hi'       , p_ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: hs'       , p_ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: conc'     , p_ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

!    IF (ltimer) CALL timer_stop(timer_ice_advection)

  END SUBROUTINE ice_advection

  !-------------------------------------------------------------------------
  !
  !> Advection of sea ice and snow on ice.
  !> Modified to advect ice%hi*ice%conc and conc instead of ice%vol.
  !! This uses the upwind_hflux_ice routine and the ocean's div_oce_3D routine to do upwind
  !! advection of the relevant variables.
  !!
  !! @par Revision History
  !! Developed by Vladimir Lapin, MPI-M (2015-06-04)
  !
  SUBROUTINE ice_advection_vla( p_patch_3D, p_op_coeff, p_ice )
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_operator_coeff),   INTENT(IN)    :: p_op_coeff
    TYPE(t_sea_ice),          INTENT(INOUT) :: p_ice

    ! Local variables
    ! Patch and range
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: cells_in_domain

    ! Indexing
    INTEGER  :: jk, jb, jc
    INTEGER  :: i_startidx_c, i_endidx_c

    ! Temporary variables/buffers
    REAL(wp) :: z_adv_flux_h (nproma,p_ice%kice,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: flux_hi  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_conc(nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: flux_hs  (nproma,p_ice%kice, p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

!--------------------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_advection)

    p_patch => p_patch_3D%p_patch_2D(1)
    cells_in_domain => p_patch%cells%in_domain

!--------------------------------------------------------------------------------------------------
! Do advection
!--------------------------------------------------------------------------------------------------

    !upwind estimate of tracer flux
    !CALL upwind_hflux_ice( p_patch_3D, p_ice%vol,  p_ice%vn_e, z_adv_flux_h )
    CALL upwind_hflux_ice( p_patch_3D, p_ice%hi*p_ice%conc,  p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), p_patch, p_op_coeff%div_coeff, flux_hi  (:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    CALL upwind_hflux_ice( p_patch_3D, p_ice%conc, p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), p_patch, p_op_coeff%div_coeff, flux_conc(:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    !CALL upwind_hflux_ice( p_patch_3D, p_ice%vols, p_ice%vn_e, z_adv_flux_h )
    CALL upwind_hflux_ice( p_patch_3D, p_ice%hs*p_ice%conc, p_ice%vn_e, z_adv_flux_h )
    DO jk=1,p_ice%kice
      CALL div_oce_3D( z_adv_flux_h(:,jk,:), p_patch, p_op_coeff%div_coeff, flux_hs  (:,jk,:),&
        & 1, cells_in_domain)
    ENDDO

    DO jk = 1,p_ice%kice
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

            ! will be divided by p_ice%conc(jc,jk,jb) after the -dtime*flux_conc is added
            p_ice%hi(jc,jk,jb)=p_ice%hi(jc,jk,jb)*p_ice%conc(jc,jk,jb)
            p_ice%hs(jc,jk,jb)=p_ice%hs(jc,jk,jb)*p_ice%conc(jc,jk,jb)

            p_ice%conc(jc,jk,jb)= p_ice%conc(jc,jk,jb)-dtime*flux_conc(jc,jk,jb)
            !p_ice%vol (jc,jk,jb)= p_ice%vol (jc,jk,jb)-dtime*flux_hi  (jc,jk,jb)
            !p_ice%vols(jc,jk,jb)= p_ice%vols(jc,jk,jb)-dtime*flux_hs  (jc,jk,jb)

            IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN
              p_ice%hi(jc,jk,jb)= (p_ice%hi(jc,jk,jb)-dtime*flux_hi(jc,jk,jb))/p_ice%conc(jc,jk,jb)
              p_ice%hs(jc,jk,jb)= (p_ice%hs(jc,jk,jb)-dtime*flux_hs(jc,jk,jb))/p_ice%conc(jc,jk,jb)
            ENDIF
          ENDIF
          ! TODO ram - remove p_patch%cells%area(jc,jb) and test
          ! See also thermodyn/mo_sea_ice.f90
          IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN
            p_ice%vol(jc,jk,jb) = p_ice%hi (jc,jk,jb)   &
              &         *( p_ice%conc(jc,jk,jb)*p_patch%cells%area(jc,jb) )
            p_ice%vols(jc,jk,jb) = p_ice%hs(jc,jk,jb)   &
              &         *( p_ice%conc(jc,jk,jb)*p_patch%cells%area(jc,jb) )
          ENDIF
!          IF ( p_ice%conc(jc,jk,jb) > 0.0_wp) THEN
!            p_ice%hi(jc,jk,jb) = p_ice%vol (jc,jk,jb)   &
!              &         /( p_ice%conc(jc,jk,jb)*p_patch%cells%area(jc,jb) )
!            p_ice%hs(jc,jk,jb) = p_ice%vols(jc,jk,jb)   &
!              &         /( p_ice%conc(jc,jk,jb)*p_patch%cells%area(jc,jb) )
!          ENDIF
        END DO
      END DO
    END DO

!--------------------------------------------------------------------------------------------------
! Sync results
!--------------------------------------------------------------------------------------------------

    CALL sync_patch_array(SYNC_C, p_patch, p_ice%vol (:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%vols(:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%conc(:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%hs  (:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%hi  (:,:,:))

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('ice_adv: vol ice'  , p_ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: vol snow' , p_ice%vols, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: hi'       , p_ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: hs'       , p_ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('ice_adv: conc'     , p_ice%conc, str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_ice_advection)

  END SUBROUTINE ice_advection_vla

  !-------------------------------------------------------------------------
  !
  !> Calculate the stress the ocean sees because of the precence of ice
  !! A future version of ths subroutine should also calculate the ice-atmosphere stresses (and have
  !! a different name)
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !<Optimize:inUse>
  SUBROUTINE ice_ocean_stress( p_patch, atmos_fluxes, p_ice, p_os )
    USE mo_physical_constants,  ONLY:  rho_ref, Cd_io
    USE mo_sea_ice_nml,         ONLY: stress_ice_zero

    TYPE(t_patch), TARGET,    INTENT(IN)    :: p_patch
    TYPE (t_atmos_fluxes),    INTENT(INOUT) :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os

    ! Local variables
    ! Ranges
    TYPE(t_subset_range), POINTER :: all_cells

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb

    ! Temporary variables/buffers
    REAL(wp) :: tau, delu, delv

!--------------------------------------------------------------------------------------------------

    all_cells => p_patch%cells%all

!--------------------------------------------------------------------------------------------------
! Modify oceanic stress
!--------------------------------------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc, delu, delv, tau) ICON_OMP_DEFAULT_SCHEDULE

  ! wind-stress is either calculated in bulk-formula or from atmosphere via coupling;
  ! it is stored in stress_xw for open water and stress_x for ice-covered area
  ! ice velocities are calculated using stress_x in ice dynamics
  ! difference of ice and ocean velocities determines ocean stress below sea ice
  ! resulting stress on ocean surface is stored in atmos_fluxes%topBoundCond_windStress_u
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      delu = p_ice%u(jc,jb) - p_os%p_diag%u(jc,1,jb)
      delv = p_ice%v(jc,jb) - p_os%p_diag%v(jc,1,jb)
      ! Ice with concentration lower than 0.01 simply flows with the speed of the ocean and does not alter drag
      ! TODO: The ice-ocean drag coefficient should depend on the depth of the upper most ocean
      ! velocity point: Cd_io = ( kappa/log(z/z0) )**2, with z0 ~= 0.4 cm
      ! Should we multiply with concSum here? 
      !tau = p_ice%concSum(jc,jb)*rho_ref*Cd_io*SQRT( delu**2 + delv**2 )
      ! #slo# - to avoid stress proportional to concSum**2 it is omitted here
      tau = rho_ref*Cd_io*SQRT( delu**2 + delv**2 )
      ! set ocean stress below sea ice to zero wrt concentration for forced runs without ice dynamics;
      ! then ocean gets no stress (no deceleration) below sea ice
      IF (stress_ice_zero) tau = 0.0_wp
      atmos_fluxes%topBoundCond_windStress_u(jc,jb) = atmos_fluxes%stress_xw(jc,jb)*( 1._wp - p_ice%concSum(jc,jb) )   &
        &               + p_ice%concSum(jc,jb)*tau*delu
      atmos_fluxes%topBoundCond_windStress_v(jc,jb) = atmos_fluxes%stress_yw(jc,jb)*( 1._wp - p_ice%concSum(jc,jb) )   &
        &               + p_ice%concSum(jc,jb)*tau*delv
    ENDDO
  ENDDO
!ICON_OMP_END_PARALLEL_DO


!   CALL sync_patch_array_mult(SYNC_C, p_patch, 2, atmos_fluxes%topBoundCond_windStress_u(:,:), &
!     & atmos_fluxes%topBoundCond_windStress_v(:,:))
!   CALL sync_patch_array(SYNC_C, p_patch, atmos_fluxes%topBoundCond_windStress_u(:,:))
!   CALL sync_patch_array(SYNC_C, p_patch, atmos_fluxes%topBoundCond_windStress_v(:,:))

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('IO-Str: windStr-u', atmos_fluxes%topBoundCond_windStress_u, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IO-Str: stress_xw', atmos_fluxes%stress_xw,                 str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IO-Str: ice%u'    , p_ice%u,                                str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IO-Str: ice%concS', p_ice%concSum,                          str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IO-Str: diag%u'   , p_os%p_diag%u,                          str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE ice_ocean_stress

  !-------------------------------------------------------------------------
  !> 1) Interpolate wind stress from cell centers to vertices
  !> 2) Rotate ocean velocities (available on the dual grid, i.e. vertices)
  !-------------------------------------------------------------------------
  ! This is a CORRECTED VERSION OF intrp_to_fem_grid_vec_old
  ! in which a problem occured when a cartesian-coords vector on ICON vertices
  ! was converted onto geogr-coords, i.e. (zonal,meridional) form.
  ! FIX: rotate the cartesian-coords vector to the rotated pole lat-lon coord
  ! and only then convert to geographic form.
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-07-17)
  !
  SUBROUTINE intrp_to_fem_grid_vec( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff ) ! TODO: replace oce_vel by oce_stress in the future.

    USE mo_ice_fem_init, ONLY: rot_mat_3D
    USE mo_ice,          ONLY: u_w, v_w, stress_atmice_x, stress_atmice_y!, u_ice, v_ice !, a_ice
!    USE mo_ice_mesh,     ONLY: coord_nod2D
!    USE mo_physical_constants,    ONLY: Cd_io, rho_ref

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
!    TYPE(t_hydro_ocean_state), target, INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
!    TYPE(t_subset_range), POINTER :: all_cells, all_verts

    ! Indexing
!    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb, jk, jv
!    INTEGER  :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: p_tau_n_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp)                      :: tau_n(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
!    TYPE(t_cartesian_coordinates) :: p_vn_dual    (nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)

!    TYPE(t_cartesian_coordinates), pointer :: tmp(:,:)

!    REAL(wp) :: tmp3(3)!, delu, u_change
!    REAL(wp) :: lat, lon!, lat1, lon1

!--------------------------------------------------------------------------------------------------
! Set up patch and ranges
!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
!    all_cells => p_patch%cells%all
!    all_verts => p_patch%verts%all

    !**************************************************************
    ! (1) Convert lat-lon wind stress to cartesian coordinates
    !**************************************************************
    CALL gvec2cvec_c_2d(p_patch_3D, atmos_fluxes%stress_x, atmos_fluxes%stress_y, p_tau_n_c)
!    DO jb = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!      DO jc = i_startidx_c, i_endidx_c
!        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
!          CALL gvec2cvec(  atmos_fluxes%stress_x(jc,jb), &
!                         & atmos_fluxes%stress_y(jc,jb), &
!                         & p_patch%cells%center(jc,jb)%lon,     &
!                         & p_patch%cells%center(jc,jb)%lat,     &
!                         & p_tau_n_c(jc,jb)%x(1),               &
!                         & p_tau_n_c(jc,jb)%x(2),               &
!                         & p_tau_n_c(jc,jb)%x(3))
!        ELSE
!          p_tau_n_c(jc,jb)%x    = 0.0_wp
!        ENDIF
!      END DO
!    END DO

    !**************************************************************
    ! (2) Interpolate 3D wind stress from cell centers to edges
    !**************************************************************

#ifdef NAGFOR
    ! only for parallel testing with nag
    tau_n = 0.0_wp
#endif
    CALL map_cell2edges_3D( p_patch_3D, p_tau_n_c, tau_n ,p_op_coeff, 1)
    CALL sync_patch_array(SYNC_E, p_patch, tau_n)

    !**************************************************************
    ! (3) Interpolate 3D wind stress from edges to vertices
    !**************************************************************
#ifdef NAGFOR
    ! only for parallel testing with nag
    p_tau_n_dual(:,:)%x(1) = 0.0_wp
    p_tau_n_dual(:,:)%x(2) = 0.0_wp
    p_tau_n_dual(:,:)%x(3) = 0.0_wp
#endif
!    CALL map_edges2verts( p_patch_3D, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual)
    CALL map_edges2verts_einar( p_patch, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual)
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(1))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(2))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(3))

    !**************************************************************
    ! (4) Rotate to the rotated pole grid
    !**************************************************************
    ! (a) atmospheric stress

    ! Rotate the vectors onto the rotated grid
    CALL rotate_cvec_v(p_patch_3D, p_tau_n_dual, rot_mat_3D, p_tau_n_dual_fem)
    ! Convert back to geographic coordinates
    CALL cvec2gvec_v_fem(p_patch_3D, p_tau_n_dual_fem, stress_atmice_x, stress_atmice_y)

    ! (4b) ocean velocities
    ! TODO: Is p_vn_dual updated?
    ! Rotate the vectors onto the rotated grid
    CALL rotate_cvec_v(p_patch_3D, p_os%p_diag%p_vn_dual(:,1,:), rot_mat_3D, p_vn_dual_fem)
!    tmp => p_os%p_diag%p_vn_dual(:,1,:)
!    CALL rotate_cvec_v(p_patch_3D, tmp, rot_mat_3D, p_vn_dual_fem)
    ! Convert back to geographic coordinates
    CALL cvec2gvec_v_fem(p_patch_3D, p_vn_dual_fem, u_w, v_w)

!    jk=0
!    DO jb = all_verts%start_block, all_verts%end_block
!      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
!      DO jv = i_startidx_v, i_endidx_v
!        jk=jk+1
!!        IF(p_patch_3D%surface_vertex_sea_land_mask(jv,jb)<=sea_boundary)THEN
!          !*************************
!          ! (4a) atmospheric stress
!          !*************************
!          ! Rotate the vectors onto the rotated grid
!          tmp3 = MATMUL( rot_mat_3D(:,:), &
!                       & (/ p_tau_n_dual(jv,jb)%x(1),               &
!                       &    p_tau_n_dual(jv,jb)%x(2),               &
!                       &    p_tau_n_dual(jv,jb)%x(3) /) )
!          ! Convert back to geographic coordinates
!          ! NOTE: this is ridiculous, but FEM mesh-coords
!	  ! are first converted to degrees and then back to radians.
!          ! lon = coord_nod2D(1,jk)*deg2rad ! FEM x-coords in degrees -- NOOOOOOOOOO
!          ! lat = coord_nod2D(2,jk)*deg2rad ! FEM x-coords in degrees -- NOOOOOOOOOO
!          lon = coord_nod2D(1,jk)
!          lat = coord_nod2D(2,jk)
!!          lat1 = p_patch%verts%vertex(jv,jb)%lat
!!          lon1 = p_patch%verts%vertex(jv,jb)%lon
!!          CALL rotate_latlon(lat1, lon1, pollat, pollon)
!!	  if ( abs(lon-lon1) + abs(lat-lat1) > 1e-12) then
!!		write(0,*) 'lon, lat: ', lon, lat, 'lon1, lat1: ', lon1, lat1
!!	  endif
!
!          CALL cvec2gvec(tmp3(1), tmp3(2), tmp3(3), lon, lat,       &
!                       & stress_atmice_x(jk), stress_atmice_y(jk))
!
!          !*************************
!          ! (4b) ocean velocities
!          !*************************
!          ! TODO: Is p_vn_dual updated?
!
!          ! Rotate the vectors onto the rotated grid
!          tmp3 = MATMUL( rot_mat_3D(:,:), &
!                       & (/ p_os%p_diag%p_vn_dual(jv,1,jb)%x(1),    &
!                       &    p_os%p_diag%p_vn_dual(jv,1,jb)%x(2),    &
!                       &    p_os%p_diag%p_vn_dual(jv,1,jb)%x(3) /) )
!          ! Convert back to geographic coordinates
!          CALL cvec2gvec(tmp3(1), tmp3(2), tmp3(3), lon, lat,       &
!                       & u_w(jk), v_w(jk))
!
!          !*************************
!          ! (4c) ice free drift
!          !*************************
!!          ! Set the ice speed to free drift speed where concentration is less than 0.01
!!          IF ( a_ice(jk) <= 0.01_wp ) THEN
!!            u_ice(jk) = 0._wp; v_ice(jk) = 0._wp; u_change = 0._wp
!!            ! TODO: Change u_change to speed_change
!!            DO WHILE ( u_change > 1e-6_wp )
!!              u_change = SQRT(u_ice(jk)**2+v_ice(jk)**2)
!!              delu = SQRT( (u_w(jk)-u_ice(jk))**2 + (v_w(jk)-v_ice(jk))**2 )
!!              u_ice(jk) = stress_atmice_x(jk)/( Cd_io*rho_ref*delu )
!!              v_ice(jk) = stress_atmice_y(jk)/( Cd_io*rho_ref*delu )
!!              u_change = ABS(u_change-SQRT(u_ice(jk)**2+v_ice(jk)**2))
!!            ENDDO
!!          ELSE
!!            ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files
!!            ! now, so this does not need to be done every timestep, only after restart file is
!!            ! read.
!            u_ice(jk) = p_ice%u_prog(jv,jb)
!            v_ice(jk) = p_ice%v_prog(jv,jb)
!!          ENDIF
!!        ELSE
!!          stress_atmice_x(jk) = 0._wp
!!          stress_atmice_y(jk) = 0._wp
!!          u_w(jk)             = 0._wp
!!          v_w(jk)             = 0._wp
!!        ENDIF
!      END DO
!    END DO

  END SUBROUTINE intrp_to_fem_grid_vec

    !-------------------------------------------------------------------------
  !> 1) Rotate back ice velocities (available on the FEM grid, i.e. vertices)
  !> 2) Interpolate to cell centers and convert back to geographic coordinates
  !-------------------------------------------------------------------------
  ! This is a CORRECTED VERSION OF intrp_from_fem_grid_vec_old
  ! the problem occurs at a pole verteces when the lat-lon velocities
  ! are rotated back to the ICON grid.
  ! FIX: first convert to the cartesian-coords vector and only then rotate
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Modified by Vladimir Lapin, MPI-M (2015-07-17)
  !
  SUBROUTINE intrp_from_fem_grid_vec( p_patch_3D, p_ice, p_op_coeff ) ! TODO: replace oce_vel by oce_stress in the future.

    USE mo_ice_fem_init, ONLY: rot_mat_3D!, pollon, pollat
    USE mo_ice,          ONLY: u_ice, v_ice
    USE mo_ice_mesh,     ONLY: coord_nod2D

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells, all_verts

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb, jk, jv
    INTEGER  :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: p_vn_dual_fem(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_vn_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: &
      & p_vn_c_3D(nproma,1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
!    REAL(wp) :: tmp3(3)
    REAL(wp) :: lat, lon

!--------------------------------------------------------------------------------------------------
! Set up patch and ranges
!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    all_verts => p_patch%verts%all

#ifdef NAGFOR
    ! only for parallel testing with nag

    p_vn_c_3D(:,:,:)%x(1) = 0.0_wp
    p_vn_c_3D(:,:,:)%x(2) = 0.0_wp
    p_vn_c_3D(:,:,:)%x(3) = 0.0_wp

    p_vn_dual_fem(:,:)%x(1) = 0.0_wp
    p_vn_dual_fem(:,:)%x(2) = 0.0_wp
    p_vn_dual_fem(:,:)%x(3) = 0.0_wp
#endif

    !**************************************************************
    ! (1) Rotate ice velocities to ICON variables and convert to cc
    !**************************************************************
    ! Convert the lat-lon vectors to 3d cartesian
    CALL gvec2cvec_v_fem(p_patch_3D, u_ice, v_ice, p_vn_dual_fem)
    ! Rotate the vectors back onto the ICON grid
    CALL rotate_cvec_v(p_patch_3D, p_vn_dual_fem, TRANSPOSE(rot_mat_3D(:,:)), p_vn_dual)

!    jk=0
!    DO jb = all_verts%start_block, all_verts%end_block
!      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
!      DO jv = i_startidx_v, i_endidx_v
!        jk=jk+1
!        ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files now,
!        ! so this does not need to be done every timestep, only before restart file is written.
!        p_ice%u_prog(jv,jb) = u_ice(jk)
!        p_ice%v_prog(jv,jb) = v_ice(jk)
!        IF(p_patch_3D%surface_vertex_sea_land_mask(jv,jb) <= sea_boundary)THEN
!          ! Rotate the vectors back to geographic grid
!          ! Convert back to geographic coordinates
!          ! NOTE: this is ridiculous, but FEM mesh-coords
!	  ! are first converted to degrees and then back to radians.
!          ! lon = coord_nod2D(1,jk)*deg2rad ! FEM x-coords in degrees -- NOOOOOOOOOO
!          ! lat = coord_nod2D(2,jk)*deg2rad ! FEM x-coords in degrees -- NOOOOOOOOOO
!          lon = coord_nod2D(1,jk)
!          lat = coord_nod2D(2,jk)
!!          lat1 = p_patch%verts%vertex(jv,jb)%lat
!!          lon1 = p_patch%verts%vertex(jv,jb)%lon
!!          CALL rotate_latlon(lat1, lon1, pollat, pollon)
!          CALL gvec2cvec(  u_ice(jk), v_ice(jk),                   &
!                         & lon, lat,                              &
!                         & tmp3(1), tmp3(2), tmp3(3) )
!          ! Rotate the vectors onto the rotated grid
!          p_vn_dual(jv,jb)%x = MATMUL( TRANSPOSE(rot_mat_3D(:,:)), &
!                                    & (/ tmp3(1), tmp3(2), tmp3(3) /) )
!        ELSE
!          p_vn_dual(jv,jb)%x(:) = 0.0_wp
!        ENDIF
!      END DO
!    END DO


    !**************************************************************
    ! (2) Interpolate ice velocities to edges for advection
    !**************************************************************

!    CALL map_verts2edges_einar( p_patch_3D, p_vn_dual, p_op_coeff%edge2cell_coeff_cc_t, p_ice%vn_e )
    CALL map_verts2edges(p_patch_3D, p_vn_dual, p_op_coeff%edge2vert_coeff_cc_t, p_ice%vn_e)
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
!    DO jb = all_cells%start_block, all_cells%end_block
!      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
!      DO jc = i_startidx_c, i_endidx_c
!!        IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN
!          CALL cvec2gvec(p_vn_c_3D(jc,1,jb)%x(1),                &
!                       & p_vn_c_3D(jc,1,jb)%x(2),                &
!                       & p_vn_c_3D(jc,1,jb)%x(3),                &
!                       & p_patch%cells%center(jc,jb)%lon,        &
!                       & p_patch%cells%center(jc,jb)%lat,        &
!                       & p_ice%u(jc,jb),                         &
!                       & p_ice%v(jc,jb))
!!        ELSE
!!          p_ice%u(jc,jb) = 0._wp
!!          p_ice%v(jc,jb) = 0._wp
!!        ENDIF
!      END DO
!    END DO

    CALL sync_patch_array(SYNC_C, p_patch, p_ice%u)
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%v)


  END SUBROUTINE intrp_from_fem_grid_vec

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!!!!!!!!!!!!!!!! Depreciated interpolation routines. !!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------
!------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> 1) Interpolate wind stress from cell centers to vertices
  !> 2) Rotate ocean velocities (available on the dual grid, i.e. vertices)
  !-------------------------------------------------------------------------
  ! IMPORTANT: SHOULD NOT BE USED DUE TO THE POLE SINGULARITY
  ! the problem occurs when a cartesian-coords vector on ICON vertices
  ! is converted onto geogr-coords, i.e. (zonal,meridional) form.
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Moved to a separate subroutine by Vladimir Lapin, MPI-M (2015-07-17)
  !
  SUBROUTINE intrp_to_fem_grid_vec_old( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff ) ! TODO: replace oce_vel by oce_stress in the future.

    USE mo_physical_constants,  ONLY: rho_ref
    USE mo_ice_fem_init, ONLY: rot_mat!, pollon, pollat
    USE mo_ice,         ONLY: a_ice, u_ice, v_ice, u_w, v_w, stress_atmice_x, stress_atmice_y
    USE mo_physical_constants,    ONLY: Cd_io
    USE mo_ice_mesh,    ONLY: coord_nod2D

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells, all_verts

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb, jk, jv
    INTEGER  :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: &
      & p_tau_n_c(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    REAL(wp) :: tau_n(nproma,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: tmp_x, tmp_y, tmp2(2), delu, u_change

!--------------------------------------------------------------------------------------------------
! Set up patch and ranges
!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    all_verts => p_patch%verts%all

    !**************************************************************
    ! (1) Convert lat-lon wind stress to cartesian coordinates
    !**************************************************************

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          CALL gvec2cvec(  atmos_fluxes%stress_x(jc,jb), &
                         & atmos_fluxes%stress_y(jc,jb), &
                         & p_patch%cells%center(jc,jb)%lon,     &
                         & p_patch%cells%center(jc,jb)%lat,     &
                         & p_tau_n_c(jc,jb)%x(1),               &
                         & p_tau_n_c(jc,jb)%x(2),               &
                         & p_tau_n_c(jc,jb)%x(3))
        ELSE
          p_tau_n_c(jc,jb)%x    = 0.0_wp
        ENDIF
      END DO
    END DO

    !**************************************************************
    ! (2) Interpolate 3D wind stress from cell centers to edges
    !**************************************************************

#ifdef NAGFOR
    ! only for parallel testing with nag
    tau_n = 0.0_wp
#endif
    CALL map_cell2edges_3D( p_patch_3D, p_tau_n_c, tau_n ,p_op_coeff, 1)
    CALL sync_patch_array(SYNC_E, p_patch, tau_n)

    !**************************************************************
    ! (3) Interpolate 3D wind stress from edges to vertices
    !**************************************************************

    CALL map_edges2verts_einar( p_patch, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual)
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(1))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(2))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(3))

    !**************************************************************
    ! (4) Convert back to geographic coordinates;
    !     rotate to the rotated grid and copy to FEM model variables.
    ! We also convert ocean forcing to geographic coordinates, rotate and copy to FEM model variables.
    ! We set the ice speed to ocean speed where concentration is low.
    !**************************************************************

    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
        IF(p_patch_3D%surface_vertex_sea_land_mask(jv,jb)<=sea_boundary)THEN
          CALL cvec2gvec(p_tau_n_dual(jv,jb)%x(1),               &
                       & p_tau_n_dual(jv,jb)%x(2),               &
                       & p_tau_n_dual(jv,jb)%x(3),               &
                       & p_patch%verts%vertex(jv,jb)%lon,        &
                       & p_patch%verts%vertex(jv,jb)%lat,        &
                       & tmp_x, tmp_y)
          ! Rotate the vectors onto the rotated grid
          tmp2 = MATMUL( rot_mat(jv,jb,:,:), (/ tmp_x, tmp_y /) )
          stress_atmice_x(jk) = tmp2(1)
          stress_atmice_y(jk) = tmp2(2)

!!!!!!!!!!!!!!!!!! check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          CALL disp_new_vect(tmp_x, tmp_y,          &
!           & p_patch%verts%vertex(jv,jb)%lon,       &
!                        & p_patch%verts%vertex(jv,jb)%lat,       &
!           & pollon, (pollat-90*deg2rad),          &
!               & coord_nod2D(1,jk), coord_nod2D(2,jk), &
!           & tmp2(1), tmp2(2))
!   if (abs(p_patch%verts%vertex(jv,jb)%lat*rad2deg) < 1) then
!     if ( abs(stress_atmice_x(jk)-tmp2(1)) + abs(stress_atmice_y(jk)-tmp2(2)) > 1e-12) then
!       write(0,*) 'lon, lat: ', p_patch%verts%vertex(jv,jb)%lon*rad2deg, p_patch%verts%vertex(jv,jb)%lat*rad2deg
!       write(0,*) 'lon_r, lat_r: ', coord_nod2D(1,jk)*rad2deg, coord_nod2D(2,jk)*rad2deg
!
!       write(0,*) 'strx, stry: ', stress_atmice_x(jk), stress_atmice_y(jk)
!       write(0,*) 'strx_new, stry_new: ', tmp2(1), tmp2(2)
!     endif
!   endif
!!!!!!!!!!!!!!!!!! check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! TODO: Is p_vn_dual_fem updated?
          CALL cvec2gvec(p_os%p_diag%p_vn_dual(jv,1,jb)%x(1),    &
                       & p_os%p_diag%p_vn_dual(jv,1,jb)%x(2),    &
                       & p_os%p_diag%p_vn_dual(jv,1,jb)%x(3),    &
                       & p_patch%verts%vertex(jv,jb)%lon,        &
                       & p_patch%verts%vertex(jv,jb)%lat,        &
                       & tmp_x, tmp_y)
          ! Rotate the vectors onto the rotated grid
          tmp2 = MATMUL( rot_mat(jv,jb,:,:), (/ tmp_x, tmp_y /) )
          u_w(jk) = tmp2(1)
          v_w(jk) = tmp2(2)
          ! Set the ice speed to free drift speed where concentration is less than 0.01
          IF ( a_ice(jk) <= 0.01_wp ) THEN
            u_ice(jk) = 0._wp; v_ice(jk) = 0._wp; u_change = 0._wp
            ! TODO: Change u_change to speed_change
            DO WHILE ( u_change > 1e-6_wp )
              u_change = SQRT(u_ice(jk)**2+v_ice(jk)**2)
              delu = SQRT( (u_w(jk)-u_ice(jk))**2 + (v_w(jk)-v_ice(jk))**2 )
              u_ice(jk) = stress_atmice_x(jk)/( Cd_io*rho_ref*delu )
              v_ice(jk) = stress_atmice_y(jk)/( Cd_io*rho_ref*delu )
              u_change = ABS(u_change-SQRT(u_ice(jk)**2+v_ice(jk)**2))
            ENDDO
          ELSE
            ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files
            ! now, so this does not need to be done every timestep, only after restart file is
            ! read.
            u_ice(jk) = p_ice%u_prog(jv,jb)
            v_ice(jk) = p_ice%v_prog(jv,jb)
          ENDIF
        ELSE
          stress_atmice_x(jk) = 0._wp
          stress_atmice_y(jk) = 0._wp
          u_w(jk)             = 0._wp
          v_w(jk)             = 0._wp
        ENDIF
      END DO
    END DO

  END SUBROUTINE intrp_to_fem_grid_vec_old

  !-------------------------------------------------------------------------
  !> 1) Rotate back ice velocities (available on the FEM grid, i.e. vertices)
  !> 2) Interpolate to cell centers and convert back to geographic coordinates
  !-------------------------------------------------------------------------
  ! IMPORTANT: SHOULD NOT BE USED DUE TO THE POLE SINGULARITY
  ! the problem occurs at a pole verteces when the lat-lon velocities
  ! are rotated back to the ICON grid.
  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !! Moved to a separate subroutine by Vladimir Lapin, MPI-M (2015-07-17)
  !
  SUBROUTINE intrp_from_fem_grid_vec_old( p_patch_3D, p_ice, p_op_coeff ) ! TODO: replace oce_vel by oce_stress in the future.

    USE mo_ice_fem_init, ONLY: rot_mat
    USE mo_ice,          ONLY: u_ice, v_ice

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells, all_verts

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb, jk, jv
    INTEGER  :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: p_vn_dual(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    TYPE(t_cartesian_coordinates) :: &
      & p_vn_c_3D(nproma,1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: tmp2(2)

!--------------------------------------------------------------------------------------------------
! Set up patch and ranges
!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    all_verts => p_patch%verts%all

    !**************************************************************
    ! (1) Rotate, reshape and copy ice velocities to ICON variables
    !     and convert to cartesian coordinates
    !**************************************************************
    jk=0
    DO jb = all_verts%start_block, all_verts%end_block
      CALL get_index_range(all_verts, jb, i_startidx_v, i_endidx_v)
      DO jv = i_startidx_v, i_endidx_v
        jk=jk+1
        ! Strictly speaking p_ice%u_prog and p_ice%v_prog are only used by the restart files now,
        ! so this does not need to be done every timestep, only before restart file is written.
        p_ice%u_prog(jv,jb) = u_ice(jk)
        p_ice%v_prog(jv,jb) = v_ice(jk)
        IF(p_patch_3D%surface_vertex_sea_land_mask(jv,jb) <= sea_boundary)THEN
          ! Rotate the vectors back to geographic grid
          tmp2 = MATMUL( TRANSPOSE(rot_mat(jv,jb,:,:)),(/ u_ice(jk), v_ice(jk) /) )
          CALL gvec2cvec(  tmp2(1),                               &
                         & tmp2(2),                               &
                         & p_patch%verts%vertex(jv,jb)%lon,       &
                         & p_patch%verts%vertex(jv,jb)%lat,       &
                         & p_vn_dual(jv,jb)%x(1),                 &
                         & p_vn_dual(jv,jb)%x(2),                 &
                         & p_vn_dual(jv,jb)%x(3))
        ELSE
          p_vn_dual(jv,jb)%x(:) = 0.0_wp
        ENDIF
      END DO
    END DO


    !**************************************************************
    ! (2) Interpolate ice velocities to edges for advection
    !**************************************************************

    CALL map_verts2edges_einar( p_patch_3D, p_vn_dual, p_op_coeff%edge2cell_coeff_cc_t, p_ice%vn_e )
    CALL sync_patch_array(SYNC_E, p_patch, p_ice%vn_e)

    !**************************************************************
    ! (3) ... and cells for drag calculation and output
    !**************************************************************
    CALL map_edges2cell_3D( p_patch_3D, RESHAPE(p_ice%vn_e, (/ nproma, 1, p_patch%nblks_e /)), &
      &   p_op_coeff, p_vn_c_3D, 1, 1)

    !**************************************************************
    ! (4) Convert back to geographic coordinates
    !**************************************************************

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN
          CALL cvec2gvec(p_vn_c_3D(jc,1,jb)%x(1),                &
                       & p_vn_c_3D(jc,1,jb)%x(2),                &
                       & p_vn_c_3D(jc,1,jb)%x(3),                &
                       & p_patch%cells%center(jc,jb)%lon,        &
                       & p_patch%cells%center(jc,jb)%lat,        &
                       & p_ice%u(jc,jb),                         &
                       & p_ice%v(jc,jb))
        ELSE
          p_ice%u(jc,jb) = 0._wp
          p_ice%v(jc,jb) = 0._wp
        ENDIF
      END DO
    END DO
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%u)
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%v)


  END SUBROUTINE intrp_from_fem_grid_vec_old



END MODULE mo_ice_fem_utils
