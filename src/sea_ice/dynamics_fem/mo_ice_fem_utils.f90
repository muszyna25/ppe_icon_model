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

MODULE mo_ice_fem_utils
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: success, min_rlcell, min_rledge, min_rlvert, min_rlcell_int, &
                                    sea, boundary
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_parallel_config,     ONLY: nproma
!  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_momentum
!    &                               timer_ice_advection, timer_ice_interp
  USE mo_grid_config,         ONLY: n_dom
  USE mo_impl_constants,      ONLY: max_char_length
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, t_atmos_for_ocean
  USE mo_oce_types,           ONLY: t_hydro_ocean_state
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates, gvec2cvec, cvec2gvec, rotate_latlon,&
    &                               rotate_latlon_vec
  USE mo_exception,           ONLY: message
  USE mo_icon_interpolation_scalar, ONLY: cells2verts_scalar
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_impl_constants,      ONLY: sea_boundary, sea
  USE mo_advection_utils,     ONLY: laxfr_upflux
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_oce_math_operators,  ONLY: div_oce_3D
  USE mo_dynamics_config,     ONLY: nold
  USE mo_scalar_product,      ONLY: map_cell2edges_3D, map_edges2cell_3D
  USE mo_math_constants,      ONLY: rad2deg, deg2rad
  USE mo_physical_constants,  ONLY: rhoi, Cd_ia, rho_ref
  USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_util_sort,           ONLY: quicksort

  IMPLICIT NONE

  PUBLIC  :: fem_ice_wrap
  PUBLIC  :: init_fem_wgts
  PUBLIC  :: destruct_fem_wgts
  PUBLIC  :: ice_fem_grid_init
  PUBLIC  :: exchange_nod2D
  PUBLIC  :: exchange_nod2Di
  PUBLIC  :: exchange_elem2D
  PUBLIC  :: ice_fem_grid_post
  PUBLIC  :: ice_advection
  PUBLIC  :: ice_ocean_stress

  PRIVATE :: map_edges2verts
  PRIVATE :: map_verts2edges
  PRIVATE :: upwind_hflux_ice

  REAL(wp), ALLOCATABLE, PRIVATE  :: c2v_wgt(:,:,:), rot_mat(:,:,:,:)
  TYPE(t_patch), POINTER, PRIVATE :: fem_patch
  ! Longitude and latitude of the north pole of the rotated grid
  ! These values put the north pole on the Indonesian/Malasian island Kalimantan and the south pole
  ! in north-west Brazil, near the border to Venezuela and Colombia.
  REAL(wp), PARAMETER, PRIVATE :: pollon = 114._wp*deg2rad, pollat = 0._wp

CONTAINS

!-----------------------------------------------------------------------
!       
!  ! Wrapper around the AWI FEM ice model as well as advection and drag calculations
!  ! Also averaging and interpolation routines and routines needed to compute the coefficients
!  ! therein
!
!-----------------------------------------------------------------------
!
!>

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

  SUBROUTINE fem_ice_wrap( p_patch_3D, p_ice, p_os, atmos_fluxes, p_op_coeff )

    USE mo_ice,      ONLY: u_ice, v_ice, m_ice, a_ice, m_snow, u_w, v_w, &
      &   stress_atmice_x, stress_atmice_y, elevation, sigma11, sigma12, sigma22
    USE mo_ice_iceparam, ONLY: C_d_io

    TYPE(t_patch_3D), TARGET, INTENT(IN)     :: p_patch_3D
    TYPE(t_sea_ice),          INTENT(INOUT)  :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)     :: p_os
    TYPE (t_atmos_fluxes),    INTENT(IN)     :: atmos_fluxes
    TYPE(t_operator_coeff),   INTENT(IN)     :: p_op_coeff

    ! Local variables
    ! Patch and ranges
    TYPE(t_patch), POINTER :: p_patch
    TYPE(t_subset_range), POINTER :: all_cells, all_verts
    INTEGER  :: i

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jv, jb, jk
    INTEGER  :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers
    TYPE(t_cartesian_coordinates) :: &
      & p_tau_n_c(nproma,p_patch_3D%p_patch_2D(n_dom)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: &
      & p_ws_n_c(nproma,p_patch_3D%p_patch_2D(n_dom)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: &
      & p_vn_c_3D(nproma,1,p_patch_3D%p_patch_2D(n_dom)%alloc_cell_blocks)
    TYPE(t_cartesian_coordinates) :: p_vn_dual(nproma,p_patch_3D%p_patch_2D(n_dom)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_tau_n_dual(nproma,p_patch_3D%p_patch_2D(n_dom)%nblks_v)
    TYPE(t_cartesian_coordinates) :: p_ws_n_dual(nproma,p_patch_3D%p_patch_2D(n_dom)%nblks_v)
    REAL(wp) :: buffy(nproma*p_patch_3D%p_patch_2D(n_dom)%nblks_v)
    REAL(wp) :: buffy_array(nproma,p_ice%kice,p_patch_3D%p_patch_2D(n_dom)%nblks_v) 
    REAL(wp) :: ws_n(nproma,p_patch_3D%p_patch_2D(n_dom)%nblks_e)
    REAL(wp) :: tau_n(nproma,p_patch_3D%p_patch_2D(n_dom)%nblks_e)
    REAL(wp) :: tmp_x, tmp_y, tmp2(2), delu, u_change

    ! TODO: There are too many calls to cvec2gvec and gvec2cvec ... but premature optimisation is the
    ! root of all evil ...

!--------------------------------------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_momentum)

!--------------------------------------------------------------------------------------------------
! Set up patch and ranges
!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(n_dom)
    all_cells => p_patch%cells%all
    all_verts => p_patch%verts%all

!--------------------------------------------------------------------------------------------------
! Interpolate and/or copy ICON variables to FEM variables
!--------------------------------------------------------------------------------------------------

!    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Interpolate tracers to vertices
    buffy_array = 0._wp
    ! TODO: Replace hi/conc to himean
    CALL cells2verts_scalar( p_ice%hi/MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    m_ice = buffy(1:SIZE(m_ice) )

    CALL cells2verts_scalar( p_ice%conc, p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    a_ice = buffy(1:SIZE(a_ice) )

    CALL cells2verts_scalar( p_ice%hs/MAX(TINY(1._wp),p_ice%conc), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    m_snow= buffy(1:SIZE(m_snow))

    ! Interpolate SSH to vertices
    CALL cells2verts_scalar( RESHAPE(p_os%p_prog(nold(1))%h(:,:), &
      & (/ nproma, 1, p_patch%alloc_cell_blocks /)), p_patch, c2v_wgt, buffy_array )
    CALL sync_patch_array(SYNC_V, p_patch, buffy_array )
    buffy = RESHAPE(buffy_array, SHAPE(buffy))
    elevation = buffy(1:SIZE(elevation) )
    
    ! Interpolate wind forcing on vertices (wind speed - ws and wind stress - tau)
    ! First we convert to cartesian coordinates
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary)THEN
          CALL gvec2cvec(  atmos_fluxes%stress_x(jc,jb),                &
                         & atmos_fluxes%stress_y(jc,jb),                &
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

    ! Then we map to edges
#ifdef NAGFOR
    ! only for parallel testing with nag
    tau_n = 0.0_wp
#endif
    CALL map_cell2edges_3D( p_patch_3D, p_tau_n_c, tau_n ,p_op_coeff, 1)
    CALL sync_patch_array(SYNC_E, p_patch, tau_n)

    ! Then we map to verts
    CALL map_edges2verts( p_patch, tau_n, p_op_coeff%edge2vert_coeff_cc, p_tau_n_dual)
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(1))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(2))
    CALL sync_patch_array(SYNC_V, p_patch, p_tau_n_dual%x(3))

    ! Now we convert back to geographic coordinates, rotate to the rotated grid and copy to FEM
    ! model variables.
    ! We also convert ocean forcing to geographic coordinates, rotate and copy to FEM model variables.
    ! We set the ice speed to ocean speed where concentration is low.
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
          ! TODO: Is p_vn_dual updated?
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
              u_ice(jk) = stress_atmice_x(jk)/( C_d_io*rho_ref*delu )
              v_ice(jk) = stress_atmice_y(jk)/( C_d_io*rho_ref*delu )
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

!    IF (ltimer) CALL timer_stop(timer_ice_interp)

!--------------------------------------------------------------------------------------------------
! Call FEM EVP
!--------------------------------------------------------------------------------------------------

    sigma11=0._wp; sigma12=0._wp; sigma22=0._wp
    CALL EVPdynamics

!--------------------------------------------------------------------------------------------------
! Post-processing: Copy FEM variables back to ICON variables
!--------------------------------------------------------------------------------------------------

!    IF (ltimer) CALL timer_start(timer_ice_interp)

    ! Rotate, reshape and copy ice velocities to ICON variables and convert to cartesian
    ! coordinates
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

    ! Interpolate ice velocities to edges for advection
    CALL map_verts2edges( p_patch_3D, p_vn_dual, p_op_coeff%edge2cell_coeff_cc_t, p_ice%vn_e )
    CALL sync_patch_array(SYNC_E, p_patch, p_ice%vn_e)

    ! ... and cells for drag calculation and output
    CALL map_edges2cell_3D( p_patch_3D, RESHAPE(p_ice%vn_e, (/ nproma, 1, p_patch%nblks_e /)), &
      &   p_op_coeff, p_vn_c_3D, 1, 1)

    ! Convert back to geographic coordinates
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

!    IF (ltimer) CALL timer_stop(timer_ice_interp)
    IF (ltimer) CALL timer_stop(timer_ice_momentum)

  END SUBROUTINE fem_ice_wrap

  !-------------------------------------------------------------------------
  !
  !> Constructor of weights for FEM sea-ice model
  !! We calculate only c2v_wgt (PRIVATE), which is used by cells2verts_scalar
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !
  SUBROUTINE init_fem_wgts(p_patch_3D)
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D

    !Local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:init_fem_wgts'

    ! Patch
    TYPE(t_patch), POINTER :: p_patch
    INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

    ! Indexing
    INTEGER :: ist, nblks_v
    INTEGER :: jb, jv
    INTEGER :: i_startblk, i_endblk, i_startidx_v, i_endidx_v

    ! The denominator
    REAL(wp) :: rdeno

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch => p_patch_3D%p_patch_2D(n_dom)
    nblks_v = p_patch%nblks_v

    ALLOCATE(c2v_wgt(nproma,6,nblks_v),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating c2v_wgt failed')
    ENDIF

    ! Weights for the cells2verts_scalar call
    i_startblk = p_patch%verts%all%start_block
    i_endblk   = p_patch%verts%all%end_block
    ! Indexing of neighbours
    iidx => p_patch%verts%cell_idx
    iblk => p_patch%verts%cell_blk
    !
    ! loop through all patch edges
    !
    BLK_LOOP: DO jb = i_startblk, i_endblk
 
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v)
 
      IDX_LOOP: DO jv =  i_startidx_v, i_endidx_v
 
        rdeno =  1._wp/MAX( TINY(rdeno),                                &
          &     p_patch_3D%wet_c  (iidx(jv,jb,1),1,iblk(jv,jb,1)) *     &
          &     p_patch%cells%area(iidx(jv,jb,1),  iblk(jv,jb,1)) +     &
          &     p_patch_3D%wet_c  (iidx(jv,jb,2),1,iblk(jv,jb,2)) *     &
          &     p_patch%cells%area(iidx(jv,jb,2),  iblk(jv,jb,2)) +     &
          &     p_patch_3D%wet_c  (iidx(jv,jb,3),1,iblk(jv,jb,3)) *     &
          &     p_patch%cells%area(iidx(jv,jb,3),  iblk(jv,jb,3)) +     &
          &     p_patch_3D%wet_c  (iidx(jv,jb,4),1,iblk(jv,jb,4)) *     &
          &     p_patch%cells%area(iidx(jv,jb,4),  iblk(jv,jb,4)) +     &
          &     p_patch_3D%wet_c  (iidx(jv,jb,5),1,iblk(jv,jb,5)) *     &
          &     p_patch%cells%area(iidx(jv,jb,5),  iblk(jv,jb,5)) +     &
          &     p_patch_3D%wet_c  (iidx(jv,jb,6),1,iblk(jv,jb,6)) *     &
          &     p_patch%cells%area(iidx(jv,jb,6),  iblk(jv,jb,6)) )

        c2v_wgt(jv,1,jb) =                              &
          &     p_patch_3D%wet_c  (iidx(jv,jb,1),1,iblk(jv,jb,1)) *     &
          &     p_patch%cells%area(iidx(jv,jb,1),  iblk(jv,jb,1)) * rdeno 
        c2v_wgt(jv,2,jb) =                              &
          &     p_patch_3D%wet_c  (iidx(jv,jb,2),1,iblk(jv,jb,2)) *     &
          &     p_patch%cells%area(iidx(jv,jb,2),  iblk(jv,jb,2)) * rdeno 
        c2v_wgt(jv,3,jb) =                              &
          &     p_patch_3D%wet_c  (iidx(jv,jb,3),1,iblk(jv,jb,3)) *     &
          &     p_patch%cells%area(iidx(jv,jb,3),  iblk(jv,jb,3)) * rdeno 
        c2v_wgt(jv,4,jb) =                              &
          &     p_patch_3D%wet_c  (iidx(jv,jb,4),1,iblk(jv,jb,4)) *     &
          &     p_patch%cells%area(iidx(jv,jb,4),  iblk(jv,jb,4)) * rdeno 
        c2v_wgt(jv,5,jb) =                              &
          &     p_patch_3D%wet_c  (iidx(jv,jb,5),1,iblk(jv,jb,5)) *     &
          &     p_patch%cells%area(iidx(jv,jb,5),  iblk(jv,jb,5)) * rdeno 
        c2v_wgt(jv,6,jb) =                              &
          &     p_patch_3D%wet_c  (iidx(jv,jb,6),1,iblk(jv,jb,6)) *     &
          &     p_patch%cells%area(iidx(jv,jb,6),  iblk(jv,jb,6)) * rdeno 

      ENDDO  IDX_LOOP
 
    ENDDO BLK_LOOP

    CALL message (TRIM(routine), 'end')        

  END SUBROUTINE init_fem_wgts

  !-------------------------------------------------------------------------
  !
  !> Destructor of patch for FEM sea-ice model, deallocates c2v_wgt
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !
  SUBROUTINE destruct_fem_wgts

    !Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:destruct_fem_wgts'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    DEALLOCATE(c2v_wgt,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating c2v_wgt failed')
    ENDIF

    DEALLOCATE(rot_mat,STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'deallocating rot_mat failed')
    ENDIF

    CALL message (TRIM(routine), 'end')        

  END SUBROUTINE destruct_fem_wgts
  !-------------------------------------------------------------------------

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
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)

#ifdef __SX__
!CDIR UNROLL=6
#endif
      DO jk = slev, elev
        DO je = i_startidx, i_endidx
          !
          ! compute the first order upwind flux; notice
          ! that multiplication by edge length is avoided to
          ! compute final conservative update using the discrete
          ! div operator
          IF ( p_patch_3D%lsm_e(je,1,jb) <= sea_boundary ) THEN
            pupflux_e(je,jk,jb) =  &
            &  laxfr_upflux( pvn_e(je,jb), pvar_c(iilc(je,jb,1),jk,iibc(je,jb,1)), &
            &                                 pvar_c(iilc(je,jb,2),jk,iibc(je,jb,2)) )
          ELSE
            pupflux_e(je,jk,jb) = 0.0_wp
          ENDIF
        END DO  ! end loop over edges
      END DO  ! end loop over levels
    END DO  ! end loop over blocks

  END SUBROUTINE upwind_hflux_ice

  !-------------------------------------------------------------------------
  !
  !> Map vectors from edges to vertices
  !! Based on map_edges2vert_3d in oce_dyn_icohom/mo_oce/math_operators.f90
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE map_edges2verts(p_patch, vn, edge2vert_coeff_cc, p_vn_dual)
    
    TYPE(t_patch), TARGET, INTENT(in)      :: p_patch
    REAL(wp), INTENT(in)           :: vn(:,:)
    TYPE(t_cartesian_coordinates),INTENT(in)  :: edge2vert_coeff_cc(:,:,:,:)
    TYPE(t_cartesian_coordinates),INTENT(out) :: p_vn_dual(:,:)

    ! Local variables

    ! Sub-set
    TYPE(t_subset_range), POINTER :: verts_in_domain

    ! Indexing
    INTEGER :: jv, jb,jev
    INTEGER :: ile, ibe
    INTEGER :: i_startidx_v, i_endidx_v

    !-----------------------------------------------------------------------

    verts_in_domain => p_patch%verts%in_domain
    ! Set to zero for nag compiler
    p_vn_dual(:,:)%x(1) = 0._wp
    p_vn_dual(:,:)%x(2) = 0._wp
    p_vn_dual(:,:)%x(3) = 0._wp

    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v

          p_vn_dual(jv,jb)%x = 0.0_wp
          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)

            p_vn_dual(jv,jb)%x = p_vn_dual(jv,jb)%x        &
            & +edge2vert_coeff_cc(jv,1,jb,jev)%x &
            & *vn(ile,ibe)
          END DO
        END DO ! jv = i_startidx_v, i_endidx_v
    END DO ! jb = verts_in_domain%start_block, verts_in_domain%end_block

  END SUBROUTINE map_edges2verts

!--------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !> Map vectors from vertices to edges
  !! Based on ideas from rot_vertex_ocean_3d in oce_dyn_icohom/mo_oce_math_operators.f90
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE map_verts2edges(p_patch_3D, p_vn_dual, edge2cell_coeff_cc_t, vn)

    TYPE(t_patch_3D), TARGET, INTENT(in)      :: p_patch_3D
    TYPE(t_cartesian_coordinates),INTENT(IN) :: p_vn_dual(:,:)
    TYPE(t_cartesian_coordinates),INTENT(in):: edge2cell_coeff_cc_t(:,:,:,:)
    REAL(wp), INTENT(inOUT)           :: vn(:,:)

    ! Local variables

    ! Patch and sub-set
    TYPE(t_subset_range), POINTER :: verts_in_domain
    TYPE(t_patch), POINTER        :: p_patch 

    ! Indexing
    INTEGER :: jb, jv, jev
    INTEGER :: ile, ibe
    INTEGER :: il_v1, il_v2,ib_v1, ib_v2
    INTEGER :: i_startidx_v, i_endidx_v

    !-----------------------------------------------------------------------

    p_patch         => p_patch_3D%p_patch_2D(1)
    verts_in_domain => p_patch%verts%in_domain


    ! Copied from oce_dyn_icohom/mo_oce_math_operators.f90:rot_vertex_ocean_3d
    ! Replaced the coefficient and sign to get normal velocity instead of tangental
    DO jb = verts_in_domain%start_block, verts_in_domain%end_block
      CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v

          DO jev = 1, p_patch%verts%num_edges(jv,jb)

            ! get line and block indices of edge jev around vertex jv
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)

            IF(p_patch_3D%lsm_e(ile,1,ibe) <= sea_boundary)THEN
              !calculate normal velocity
              il_v1 = p_patch%edges%vertex_idx(ile,ibe,1)
              ib_v1 = p_patch%edges%vertex_blk(ile,ibe,1)
              il_v2 = p_patch%edges%vertex_idx(ile,ibe,2)
              ib_v2 = p_patch%edges%vertex_blk(ile,ibe,2)

              vn(ile,ibe) = &
                &   DOT_PRODUCT(p_vn_dual(il_v1,ib_v1)%x,edge2cell_coeff_cc_t(ile,1,ibe,1)%x) &
                & + DOT_PRODUCT(p_vn_dual(il_v2,ib_v2)%x,edge2cell_coeff_cc_t(ile,1,ibe,2)%x)
            ELSE
              vn(ile,ibe) = 0._wp
            ENDIF
          END DO
      END DO
    END DO
  
  END SUBROUTINE map_verts2edges

  !-------------------------------------------------------------------------
  !
  !> Initialise the FEM grid
  !! This replaces the routine Sergey used to read the grid from file. We give values to
  !! coord_nod2D, nod2D, index_nod2D, elem2D_nodes and elem2D from the ICON grid.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_grid_init(p_patch_3D)
    ! For the AWI FEM
    USE mo_ice_mesh,           ONLY: coord_nod2D, nod2D, index_nod2D
    USE mo_ice_elements,       ONLY: elem2D_nodes, elem2D
   ! USE mo_mpi

    TYPE(t_patch_3D), TARGET, INTENT(in) :: p_patch_3D

    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ice_fem_utils:ice_fem_grid_init'

    ! Patch
    TYPE(t_patch), POINTER :: p_patch

    ! Indexing
    INTEGER :: k, jb, jc, jv
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: i_startidx_v, i_endidx_v

    ! Temporary variables/buffers and flags
    INTEGER :: verts(nproma,p_patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: buffy(nproma*p_patch_3D%p_patch_2D(n_dom)%nblks_v)
    INTEGER :: ist
    REAL(wp) :: lat, lon, cos_d, sin_d

    ! Masking for the halo-points
    INTEGER :: halo_mask(nproma, p_patch_3D%p_patch_2D(n_dom)%nblks_v)

   ! REAL(wp) :: test(nproma, p_patch_3D%p_patch_2D(n_dom)%alloc_cell_blocks, 3)
   ! INTEGER :: my_mpi_id

!--------------------------------------------------------------------------------------------------

    p_patch => p_patch_3D%p_patch_2D(1)
    fem_patch => p_patch
    ! my_mpi_id = get_my_global_mpi_id()

    ! nod2D is the number of nodes = number of vertices
    nod2D = p_patch%n_patch_verts
    ALLOCATE(coord_nod2D(2,nod2D), index_nod2D(nod2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coord_nod2D failed')
    ENDIF
    ALLOCATE(rot_mat(nproma,p_patch%nblks_v,2,2),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating rot_mat failed')
    ENDIF


    ! Go through all the vertices and assign coordinates and land mask to coord_nod2D and
    ! index_nod2D
    k=0
    DO jb = 1,p_patch%nblks_v
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v) 
      DO jv = i_startidx_v,i_endidx_v
        k=k+1
        verts(jv,jb) = k
        lat = p_patch%verts%vertex(jv,jb)%lat
        lon = p_patch%verts%vertex(jv,jb)%lon

        ! Calculate the rotation matrix once
        CALL rotate_latlon_vec( lon, lat, pollon, pollat, sin_d, cos_d )
        rot_mat(jv,jb,1,1) =  cos_d
        rot_mat(jv,jb,1,2) =  sin_d
        rot_mat(jv,jb,2,1) = -sin_d
        rot_mat(jv,jb,2,2) =  cos_d

        ! Use a rotated grid with the north pole at (pollon,pollat)
        CALL rotate_latlon(lat, lon, pollon, pollat)

        ! x-coords in degrees
        coord_nod2D(1,k) = lon*rad2deg
        ! y-coords in degrees
        coord_nod2D(2,k) = lat*rad2deg
        ! border (and land) points (0 or 1)
        index_nod2D(k)   =  MIN(1, MAX(0, 1+p_patch_3D%surface_vertex_sea_land_mask(jv,jb)))
      ENDDO
    ENDDO

    ! elem2D is the number of elements = number of cells
   ! test(:,:,:) = 0.0_wp
    elem2D = p_patch%n_patch_cells
    ALLOCATE(elem2D_nodes(3,elem2D),STAT=ist)
    IF (ist /= SUCCESS) THEN
      CALL finish (routine,'allocating coord_nod2D failed')
    ENDIF

    ! Establish the connectivity between elements and nodes. The array elem2D_nodes lists the
    ! nodes each element is made up of.
    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        elem2D_nodes(1,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,1),p_patch%cells%vertex_blk(jc,jb,1))
        elem2D_nodes(2,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,2),p_patch%cells%vertex_blk(jc,jb,2))
        elem2D_nodes(3,k) =   &
          &   verts(p_patch%cells%vertex_idx(jc,jb,3),p_patch%cells%vertex_blk(jc,jb,3))

!         test(jc, jb, 1) = &
!           & REAL(p_patch%verts%glb_index(verts(p_patch%cells%vertex_idx(jc,jb,1),p_patch%cells%vertex_blk(jc,jb,1))), wp)
!         test(jc, jb, 2) = &
!           & REAL(p_patch%verts%glb_index(verts(p_patch%cells%vertex_idx(jc,jb,2),p_patch%cells%vertex_blk(jc,jb,2))), wp)
!         test(jc, jb, 3) = &
!           & REAL(p_patch%verts%glb_index(verts(p_patch%cells%vertex_idx(jc,jb,3),p_patch%cells%vertex_blk(jc,jb,3))), wp)
!
!         write(0,*) my_mpi_id, "::", p_patch%cells%glb_index(k), ":",  test(jc, jb, 1), test(jc, jb, 2), test(jc, jb, 3)

      ENDDO
    ENDDO

!    write(0,*) "sync test(jc, jb, 1)..."
!    CALL sync_patch_array(SYNC_C, p_patch, test(:, :, 1))
!    write(0,*) "sync test(jc, jb, 2)..."
!    CALL sync_patch_array(SYNC_C, p_patch, test(:, :, 2))
!    write(0,*) "sync test(jc, jb, 3)..."
!    CALL sync_patch_array(SYNC_C, p_patch, test(:, :, 3))
!    write(0,*) "sync test( done."


    ! Mask out halopoints
    ! This is neccesary so that Sergey's EVP routine doesn't try to calculate velocities on the
    ! halo-points. In his code ponts with index_nod2D=0 are used for calculation and the velocity
    ! at points with index_nod2D=1 (land borders) is set to zero (no-slip). The velocity at points
    ! with index_nod2D>1 will therefore not be calculated and not set to zero.
    halo_mask = 0
    DO jb = p_patch%verts%not_owned%start_block, p_patch%verts%not_owned%end_block
      CALL get_index_range(p_patch%verts%not_owned, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v
          halo_mask(jv,jb) = 2
        ENDDO
    ENDDO

    buffy = RESHAPE(halo_mask(:,:), SHAPE(buffy))
    index_nod2D = index_nod2D+buffy(1:SIZE(index_nod2D))

  END SUBROUTINE ice_fem_grid_init

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM
  !! This replaces the routine Sergey used to synchronize arrays across CPUs. It is just a wrapper
  !! around sync_patch_array using the PRIVATE variable fem_patch as patch and reshaping the
  !! variable before and after syncing.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE exchange_nod2D(u_ice)

    REAL(wp), INTENT(INOUT) :: u_ice(fem_patch%n_patch_verts)

    ! local variables
    ! Temporary variables/buffers and pad
    REAL(wp) :: u(nproma, fem_patch%nblks_v)
    REAL(wp) :: buffy(nproma*fem_patch%nblks_v)
    REAL(wp), ALLOCATABLE :: pad(:)

  !-------------------------------------------------------------------------

    ! Reshape and copy ice velocities to ICON variables
    ALLOCATE(pad(nproma*fem_patch%nblks_v - fem_patch%n_patch_verts))
    pad = -9999._wp
    u=RESHAPE(u_ice, SHAPE(u), pad)
    DEALLOCATE(pad)

    CALL sync_patch_array(SYNC_V, fem_patch, u(:,:))

    ! Reshape and copy ice velocity to FEM model variables
    buffy = RESHAPE(u, SHAPE(buffy))
    u_ice = buffy(1:SIZE(u_ice))

  END SUBROUTINE exchange_nod2D

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM having an integer array as input
  !! For testing purposes only
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE exchange_nod2Di(index)

    INTEGER, INTENT(INOUT) :: index(fem_patch%n_patch_verts)

    ! local variables
    ! Temporary variables/buffers and pad
    INTEGER :: i(nproma, fem_patch%nblks_v)
    INTEGER :: buffy(nproma*fem_patch%nblks_v)
    INTEGER, ALLOCATABLE :: pad(:)

  !-------------------------------------------------------------------------

    ! Reshape and copy ice velocities to ICON variables
    ALLOCATE(pad(nproma*fem_patch%nblks_v - fem_patch%n_patch_verts))
    pad = -9999
    i=RESHAPE(index, SHAPE(i), pad)
    DEALLOCATE(pad)

    CALL sync_patch_array(SYNC_V, fem_patch, i(:,:))

    ! Reshape and copy ice velocity to FEM model variables
    buffy = RESHAPE(i, SHAPE(buffy))
    index = buffy(1:SIZE(index))

  END SUBROUTINE exchange_nod2Di

  !-------------------------------------------------------------------------
  !
  !> Synchronisation routine for the FEM elements
  !! For testing purposes only
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE exchange_elem2D(sigma)

    REAL(wp), INTENT(INOUT) :: sigma(fem_patch%n_patch_cells)

    ! local variables
    ! Temporary variables/buffers and pad
    REAL(wp) :: s(nproma, fem_patch%alloc_cell_blocks)
    REAL(wp) :: buffy(nproma*fem_patch%alloc_cell_blocks)
    REAL(wp), ALLOCATABLE :: pad(:)

  !-------------------------------------------------------------------------

    ! Reshape and copy fem variable to ICON variable
    ALLOCATE(pad(nproma*fem_patch%alloc_cell_blocks - fem_patch%n_patch_cells))
    pad = -9999._wp
    s=RESHAPE(sigma, SHAPE(s), pad)
    DEALLOCATE(pad)

    CALL sync_patch_array(SYNC_C, fem_patch, s(:,:))

    ! Reshape and copy ice velocity to FEM model variables
    buffy = RESHAPE(s, SHAPE(buffy))
    sigma = buffy(1:SIZE(sigma))

  END SUBROUTINE exchange_elem2D

  !-------------------------------------------------------------------------
  !
  !> Post-initialisation work for the FEM grid
  !! In order for the FEM model to pass the p-test we need to sort the lists myList_elem2D and
  !! myList_nod2D according to the global index. This is neccesary so that summation in stress2rhs
  !! (in ice_dyn_fem/ice_evp.f90) is always done in the same order on all CPUs.
  !! In addition to this we also need to set lmass_matrix to the dual area size, since lmass_matrix
  !! does not synchronise across CPUs the way it is calculated here. We also set voltriangle to the
  !! cell area, since the cell area is more accurately calculated by ICON than the FEM code.
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-08-05)
  !
  SUBROUTINE ice_fem_grid_post(p_patch)
    USE mo_ice_parsup,          ONLY: myList_nod2D, myList_elem2D
    USE mo_ice,                 ONLY: lmass_matrix
    USE mo_ice_elements,        ONLY: voltriangle
    USE mo_ice_iceparam,        ONLY: coriolis_nod2D
    USE mo_ice_mesh,            ONLY: cos_elem2D, sin_elem2D

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch

    ! Local variables
    ! Indexing
    INTEGER :: i_startidx_v, i_endidx_v
    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: k, jb, jv, jc

    ! Global lists
    INTEGER :: globList_nod2D(p_patch%n_patch_verts)
    INTEGER :: globList_elem2D(p_patch%n_patch_cells)

    ! Temporary variables/buffers
    REAL(wp) :: buffy_v(nproma*p_patch%nblks_v)
    REAL(wp) :: buffy_c(nproma*p_patch%alloc_cell_blocks)

  !-------------------------------------------------------------------------
    
    ! Sort the list myList_nod2D using the global index
    k=0
    DO jb = 1,p_patch%nblks_v
      CALL get_index_range(p_patch%verts%all, jb, i_startidx_v, i_endidx_v) 
      DO jv = i_startidx_v,i_endidx_v
        k=k+1
        globList_nod2D(k) = p_patch%verts%decomp_info%glb_index(k)
      ENDDO
    ENDDO

    CALL quicksort(globList_nod2D,myList_nod2D)

    ! Sort the list myList_elem2D using the global index
    k=0
    DO jb = p_patch%cells%all%start_block, p_patch%cells%all%end_block
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c) 
      DO jc = i_startidx_c,i_endidx_c
        k=k+1
        globList_elem2D(k) = p_patch%cells%decomp_info%glb_index(k)
      ENDDO
    ENDDO

    CALL quicksort(globList_elem2D,myList_elem2D)

    ! Re-create lmass_matrix as one third of the sum of areas of triangles that share vertex n
    ! Or in ICON-terms as the dual area

    buffy_v = RESHAPE(p_patch%verts%dual_area(:,:), SHAPE(buffy_v))
    lmass_matrix = buffy_v(1:SIZE(lmass_matrix))

    ! voltriangle should also be the cell area
    buffy_c = RESHAPE(p_patch%cells%area(:,:), SHAPE(buffy_c))
    voltriangle = buffy_c(1:SIZE(voltriangle))

    ! coriolis_nod2D, cos_elem2D, sin_elem2D need to be reset because of grid rotation

    buffy_v = RESHAPE(p_patch%verts%f_v(:,:), SHAPE(buffy_v))
    coriolis_nod2D = buffy_v(1:SIZE(coriolis_nod2D))

    buffy_c = RESHAPE(p_patch%cells%center(:,:)%lat, SHAPE(buffy_c))
    cos_elem2D = COS(buffy_c(1:SIZE(cos_elem2D)))

    buffy_c = RESHAPE(p_patch%cells%center(:,:)%lat, SHAPE(buffy_c))
    sin_elem2D = SIN(buffy_c(1:SIZE(sin_elem2D)))

  END SUBROUTINE ice_fem_grid_post

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

    p_patch => p_patch_3D%p_patch_2D(n_dom)
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
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%conc(:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%hs  (:,:,:))
    CALL sync_patch_array(SYNC_C, p_patch, p_ice%hi  (:,:,:))

!    IF (ltimer) CALL timer_stop(timer_ice_advection)

  END SUBROUTINE ice_advection

  !-------------------------------------------------------------------------
  !
  !> Calculate the stress the ocean sees because of the precence of ice
  !! A future version of ths subroutine should also calculate the ice-atmosphere stresses (and have
  !! a different name)
  !!
  !! @par Revision History
  !! Developed by Einar Olason, MPI-M (2013-06-05)
  !
  SUBROUTINE ice_ocean_stress( p_patch, atmos_fluxes, p_ice, p_os )
    USE mo_ice_param,    ONLY: density_0
    USE mo_ice_iceparam, ONLY: C_d_io
    TYPE(t_patch), TARGET,    INTENT(IN)    :: p_patch
    TYPE (t_atmos_fluxes),    INTENT(INOUT) :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(IN)    :: p_ice
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os

    ! Local variables
    ! Ranges
    TYPE(t_subset_range), POINTER :: cells_in_domain

    ! Indexing
    INTEGER  :: i_startidx_c, i_endidx_c, jc, jb

    ! Temporary variables/buffers
    REAL(wp) :: tau, delu, delv

!--------------------------------------------------------------------------------------------------

    cells_in_domain => p_patch%cells%in_domain

!--------------------------------------------------------------------------------------------------
! Modify oceanic stress
!--------------------------------------------------------------------------------------------------

  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
    DO jc = i_startidx_c, i_endidx_c
      ! Ice with concentration lower than 0.01 simply flows with the speed of the ocean and does
      ! not alter drag
      ! TODO: The ice-ocean drag coefficient should depend on the depth of the upper most ocean
      ! velocity point: C_d_io = ( kappa/log(z/z0) )**2, with z0 ~= 0.4 cm
      delu = p_ice%u(jc,jb) - p_os%p_diag%u(jc,1,jb)
      delv = p_ice%v(jc,jb) - p_os%p_diag%v(jc,1,jb)

      ! Should we multiply with concSum here?
      tau = p_ice%concSum(jc,jb)*density_0*C_d_io*SQRT( delu**2 + delv**2 )
      atmos_fluxes%topBoundCond_windStress_u(jc,jb) = atmos_fluxes%stress_xw(jc,jb)*( 1._wp - p_ice%concSum(jc,jb) )   &
        &               + p_ice%concSum(jc,jb)*tau*delu
      atmos_fluxes%topBoundCond_windStress_v(jc,jb) = atmos_fluxes%stress_yw(jc,jb)*( 1._wp - p_ice%concSum(jc,jb) )   &
        &               + p_ice%concSum(jc,jb)*tau*delv
    ENDDO
  ENDDO
  CALL sync_patch_array(SYNC_C, p_patch, atmos_fluxes%topBoundCond_windStress_u(:,:))
  CALL sync_patch_array(SYNC_C, p_patch, atmos_fluxes%topBoundCond_windStress_v(:,:))

  END SUBROUTINE ice_ocean_stress

END MODULE mo_ice_fem_utils
