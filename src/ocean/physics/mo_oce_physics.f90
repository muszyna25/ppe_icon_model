!>
!! Provide an implementation of the ocean physics.
!!
!! Provide an implementation of the physical parameters and characteristics
!! for the hydrostatic ocean model.
!!
!! @author Stephan Lorenz, MPI
!! @author Peter Korn, MPI
!!
!! @par Revision History
!!  Original version by Peter Korn, MPI-M (2009)
!!  Modified by Stephan Lorenz,     MPI-M (2010-07)
!!    adapted to structures discussed in 2010-01.
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_oce_physics
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
!
!
!
USE mo_kind,                ONLY: wp
USE mo_ocean_nml,           ONLY: n_zlev, bottom_drag_coeff, k_veloc_h, k_veloc_v,        &
  &                               k_pot_temp_h, k_pot_temp_v, k_sal_h, k_sal_v, no_tracer,&
  &                               MAX_VERT_DIFF_VELOC, MAX_VERT_DIFF_TRAC,                &
  &                               HORZ_VELOC_DIFF_TYPE, veloc_diffusion_order,            &
  &                               N_POINTS_IN_MUNK_LAYER,                                 &
  &                               biharmonic_diffusion_factor,                            &
  &                               richardson_tracer, richardson_veloc,                    &
  &                               l_constant_mixing, l_smooth_veloc_diffusion,            &
  &                               l_wind_mixing    !, l_convection, l_pp_scheme
USE mo_parallel_config,     ONLY: nproma
USE mo_model_domain,        ONLY: t_patch, t_patch_3D
USE mo_impl_constants,      ONLY: success, max_char_length, MIN_DOLIC, SEA
USE mo_exception,           ONLY: message, message_text, finish
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_oce_state,           ONLY: t_hydro_ocean_state, oce_config!, v_base
USE mo_physical_constants,  ONLY: grav, rho_ref, SItodBar,sal_ref
USE mo_math_constants,      ONLY: dbl_eps
USE mo_dynamics_config,     ONLY: nold!, nnew
USE mo_sea_ice_types,       ONLY: t_sfc_flx
USE mo_run_config,          ONLY: dtime
USE mo_linked_list,         ONLY: t_var_list
USE mo_var_list,            ONLY: add_var,                  &
  &                               new_var_list,             &
  &                               delete_var_list,          &
  &                               default_var_list_settings,&
  &                               add_ref, groups
USE mo_cf_convention
USE mo_grib2
USE mo_cdi_constants,       ONLY: GRID_CELL, GRID_EDGE, GRID_REFERENCE,           &
  &                               GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_CELL, &
  &                               ZA_DEPTH_BELOW_SEA, ZA_DEPTH_BELOW_SEA_HALF,    &
  &                               datatype_pack16, datatype_flt32, filetype_nc2
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
USE mo_sync,                ONLY: SYNC_C, SYNC_E, sync_patch_array, global_max
IMPLICIT NONE


PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

CHARACTER(len=*), PARAMETER :: this_mod_name = 'mo_oce_physics'
CHARACTER(len=12)           :: str_module    = 'ocePhysics  '  ! Output of module for 1 line debug
INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

! Public interface
PUBLIC  :: t_ptr3d, t_ho_params

!PUBLIC :: init_ho_physics
PUBLIC  :: construct_ho_params
PUBLIC  :: destruct_ho_params
PUBLIC  :: init_ho_params
PUBLIC  :: update_ho_params

! variables
Type (t_var_list), PUBLIC :: ocean_params_list

TYPE t_ptr3d
  REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
END TYPE t_ptr3d

! Parameters below appear directly in the ocean model/equation. They are eventually
! dynamically updated by using the "ocean-physics" structure. #slo# - not yet
TYPE t_ho_params

  ! diffusion coefficients for horizontal velocity, temp. and salinity, dim=(nproma,n_zlev,nblks_e)
  REAL(wp),POINTER ::     &
    &  K_veloc_h(:,:,:),  & ! coefficient of horizontal velocity diffusion
    &  K_tracer_h(:,:,:,:)  ! coefficient of horizontal tracer diffusion
  TYPE(t_ptr3d),ALLOCATABLE :: tracer_h_ptr(:)

  ! diffusion coefficients for vertical velocity, temp. and salinity, dim=(nproma,n_zlev+1,nblks_e)
  REAL(wp),POINTER ::     &
    &  A_veloc_v(:,:,:),  & ! coefficient of vertical velocity diffusion
    &  A_tracer_v(:,:,:,:)  ! coefficient of vertical tracer diffusion
  TYPE(t_ptr3d),ALLOCATABLE :: tracer_v_ptr(:)

  !constant background values of coefficients above
  REAL(wp) :: K_veloc_h_back, &! coefficient of horizontal velocity diffusion
           &  A_veloc_v_back   ! coefficient of vertical velocity diffusion

  REAL(wp),ALLOCATABLE ::     &
    &  K_tracer_h_back(:),    & ! coefficient of horizontal tracer diffusion dim=no_tracer
    &  A_tracer_v_back(:)       ! coefficient of vertical tracer diffusion dim=no_tracer

  REAL(wp) :: bottom_drag_coeff

END TYPE t_ho_params


TYPE(t_ho_params),PUBLIC,TARGET :: v_params
CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! Initialisation of ocean physics
  !!
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07)
  !
  !! mpi parallelized, syncs p_phys_param%K_tracer_h, p_phys_param%K_veloc_h
  SUBROUTINE init_ho_params(  p_patch_3D, p_phys_param )
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
    TYPE (t_ho_params)                          :: p_phys_param

    ! Local variables
    INTEGER  :: i, i_no_trac 
    INTEGER  :: je,jb
    INTEGER  :: i_startidx_e, i_endidx_e
    REAL(wp) :: z_lower_bound_diff
    REAL(wp) :: z_largest_edge_length ,z_diff_multfac, z_diff_efdt_ratio
    REAL(wp) :: points_in_munk_layer
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    all_edges => p_patch%edges%all
    !-------------------------------------------------------------------------
    points_in_munk_layer = REAL(N_POINTS_IN_MUNK_LAYER,wp)
    !Init from namelist
    p_phys_param%K_veloc_h_back = k_veloc_h
    p_phys_param%A_veloc_v_back = k_veloc_v
    p_phys_param%K_veloc_h      = k_veloc_h
    p_phys_param%A_veloc_v      = k_veloc_v

    z_largest_edge_length = global_max(MAXVAL(p_patch%edges%primal_edge_length))


    !Distinghuish between harmonic and biharmonic laplacian
    !Harmonic laplacian
    IF(veloc_diffusion_order==1)THEN
      SELECT CASE(HORZ_VELOC_DIFF_TYPE)
      CASE(0)!no friction
        p_phys_param%K_veloc_h(:,:,:) = 0.0_wp

      CASE(1)!use uniform viscosity coefficient from namelist
        CALL calc_lower_bound_veloc_diff(  p_patch, z_lower_bound_diff )
        IF(z_lower_bound_diff>p_phys_param%K_veloc_h_back)THEN
          ! SX9 cannot handle messages of that size -> split
          CALL message ('init_ho_params','WARNING: Specified diffusivity&
                        & does not satisfy Munk criterion.')
          CALL message ('init_ho_params','WARNING: This may lead&
                        & to stability problems for experiments with lateral boundaries')
        ENDIF

        p_phys_param%K_veloc_h(:,:,:) = p_phys_param%K_veloc_h_back
        !write(0,*)'lower bound of diffusivity:',z_lower_bound_diff
        WRITE(message_text,'(a,g25.16)') 'Lower bound of diffusivity:',z_lower_bound_diff
        CALL message ('init_ho_params', message_text)

      CASE(2)!calculate uniform viscosity coefficient, according to Munk criterion

        p_phys_param%K_veloc_h(:,:,:) = 3.82E-12_wp&
        &*(points_in_munk_layer*z_largest_edge_length)**3

      CASE(3)! calculate coefficients for each location based on MUNK layer
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO je = i_startidx_e, i_endidx_e
            !calculate lower bound for diffusivity
            !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29)) 
            p_phys_param%K_veloc_h(je,:,jb) = 3.82E-12_wp&
            &*(points_in_munk_layer*p_patch%edges%primal_edge_length(je,jb))**3
          END DO
        END DO
      END SELECT
      CALL dbg_print('horzVelocDiff:',p_phys_param%K_veloc_h ,str_module,0,in_subset=p_patch%cells%owned)
    !Biharmonic laplacian
    ELSEIF(veloc_diffusion_order==2)THEN

        !The general form follows the hydrostatic atmospheric code.
        !The number that controls all that the "z_diff_efdt_ratio"
        !is different. Higher z_diff_efdt_ratio decreases the final
        !diffusion coefficient 
        z_diff_efdt_ratio = 10000.0_wp * biharmonic_diffusion_factor
        z_diff_multfac = (1._wp/ (z_diff_efdt_ratio*64._wp))/3._wp
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO je = i_startidx_e, i_endidx_e
             p_phys_param%K_veloc_h(je,:,jb) = &
              &p_patch%edges%area_edge(je,jb)*p_patch%edges%area_edge(je,jb)*z_diff_multfac
          END DO
        END DO

!          z_diff_multfac = 0.0045_wp*dtime/3600.0_wp
!         DO jb = all_edges%start_block, all_edges%end_block
!            CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
!            DO je = i_startidx_e, i_endidx_e
!              p_phys_param%K_veloc_h(je,:,jb) = z_diff_multfac*&
!              &maxval(p_patch%edges%primal_edge_length)**4
!            END DO
!          END DO


    ENDIF
    IF ( l_smooth_veloc_diffusion ) CALL smooth_lapl_diff( p_patch, p_patch_3D, p_phys_param%K_veloc_h )


    DO i=1,no_tracer

      IF(i==1)THEN!temperature
        p_phys_param%K_tracer_h_back(i) = k_pot_temp_h
        p_phys_param%A_tracer_v_back(i) = k_pot_temp_v

      ELSEIF(i==2)THEN!salinity
        p_phys_param%K_tracer_h_back(2) = k_sal_h
        p_phys_param%A_tracer_v_back(2) = k_sal_v
      ELSE

        CALL finish ('mo_oce_physics:init_ho_params',  &
        &      'number of tracers exceeds number of background values')
      ENDIF
      p_phys_param%K_tracer_h(:,:,:,i) = p_phys_param%K_tracer_h_back(i)
      p_phys_param%A_tracer_v(:,:,:,i) = p_phys_param%A_tracer_v_back(i)
    END DO

    p_phys_param%bottom_drag_coeff = bottom_drag_coeff

    DO i_no_trac=1, no_tracer
      CALL sync_patch_array(SYNC_C,p_patch,p_phys_param%K_tracer_h(:,:,:,i_no_trac))
    END DO
    CALL sync_patch_array(SYNC_E,p_patch,p_phys_param%K_veloc_h(:,:,:))

  END SUBROUTINE init_ho_params

  !-------------------------------------------------------------------------
  !
  !>
  !! Calculation of a lower bound for horizontal velocity diffusion of laplacian type ! that is
  !! required to have N (default =1) points in Munk layer. The lower bound is calculated ! with
  !! respect to the equator.
  !! The code is based on  Griffies, Fundamentals of ocean climate modeling, sect 18, p. 413.  !
  !! The lower bound is given in units [m^2/s].
  !!
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
  !
  !
  SUBROUTINE calc_lower_bound_veloc_diff(  p_patch, lower_bound_diff )
    TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
    REAL(wp), INTENT(inout)              :: lower_bound_diff

    ! Local variables
    REAL(wp) :: points_in_munk_layer
    REAL(wp)            :: z_largest_edge_length
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-------------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    ! Get the largest edge length globally
    z_largest_edge_length = global_max(MAXVAL(p_patch%edges%primal_edge_length))

    !calculate lower bound for diffusivity: The factor cos(lat) is omitted here, because of
    !equatorial reference (cf. Griffies, eq.  (18.29)) 
    points_in_munk_layer = REAL(N_POINTS_IN_MUNK_LAYER,wp)
    lower_bound_diff = 3.82E-12_wp*(points_in_munk_layer*z_largest_edge_length)**3

  END SUBROUTINE calc_lower_bound_veloc_diff

  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
  !
  !! mpi parallelized, no sync required
  SUBROUTINE smooth_lapl_diff( p_patch,p_patch_3D, K_h )
   TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
   TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: p_patch_3D
   REAL(wp), INTENT(INOUT)    :: K_h(:,:,:)
    ! Local variables
    INTEGER  :: je,jv,jb,jk, jev, ile, ibe, i_edge_ctr
    INTEGER  :: il_v1,ib_v1, il_v2,ib_v2
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: i_startidx_v, i_endidx_v
    REAL(wp) :: z_K_ave_v(nproma,n_zlev,p_patch%nblks_v), z_K_max
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, verts_in_domain
    !-------------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain
    verts_in_domain => p_patch%verts%in_domain

    z_K_ave_v(:,:,:) = 0.0_wp

    DO jk = 1, n_zlev
      DO jb = verts_in_domain%start_block, verts_in_domain%end_block
        CALL get_index_range(verts_in_domain, jb, i_startidx_v, i_endidx_v)
        DO jv = i_startidx_v, i_endidx_v
          i_edge_ctr = 0
          z_K_max    = 0.0_wp
          DO jev = 1, p_patch%verts%num_edges(jv,jb)
            ile = p_patch%verts%edge_idx(jv,jb,jev)
            ibe = p_patch%verts%edge_blk(jv,jb,jev)
!             write(0,*) jv,jb, p_patch%verts%num_edges(jv,jb), ":", ile, ibe
            IF ( p_patch_3D%lsm_e(ile,jk,ibe) == sea) THEN
              z_K_ave_v(jv,jk,jb)= z_K_ave_v(jv,jk,jb) + K_h(ile,jk,ibe)
              i_edge_ctr=i_edge_ctr+1
              IF(K_h(ile,jk,ibe)>z_K_max)THEN
              z_K_max=K_h(ile,jk,ibe)
              ENDIF
            ENDIF
          END DO
          IF(i_edge_ctr/=0)THEN!.and.i_edge_ctr== p_patch%verts%num_edges(jv,jb))THEN
            z_K_ave_v(jv,jk,jb)= z_K_ave_v(jv,jk,jb)/REAL(i_edge_ctr,wp)
          ELSEIF(i_edge_ctr==0)THEN
            z_K_ave_v(jv,jk,jb)=0.0_wp
          ENDIF
          !IF(p_patch%verts%num_edges(jv,jb)== 5)THEN
          !  z_K_ave_v(jv,jk,jb)=80000_wp!°z_K_max
          !ENDIF
        END DO
      ENDDO
    END DO


    DO jk = 1, n_zlev
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e

          il_v1 = p_patch%edges%vertex_idx(je,jb,1)
          ib_v1 = p_patch%edges%vertex_blk(je,jb,1)
          il_v2 = p_patch%edges%vertex_idx(je,jb,2)
          ib_v2 = p_patch%edges%vertex_blk(je,jb,2)

          IF ( p_patch_3D%lsm_e(je,jk,jb) == sea) THEN
            K_h(je,jk,jb)= 0.5_wp*(z_K_ave_v(il_v1,jk,ib_v1) + z_K_ave_v(il_v2,jk,ib_v2))
          ELSE
            K_h(je,jk,jb)=0.0_wp
          ENDIF
!          IF(p_patch%verts%num_edges(il_v1,ib_v1)== 5.OR.p_patch%verts%num_edges(il_v2,ib_v2)==5)THEN
!            K_h(je,jk,jb)=max(z_K_ave_v(il_v1,jk,ib_v1),z_K_ave_v(il_v2,jk,ib_v2))
!          ENDIF
        END DO
      ENDDO
    END DO

    !---------Debug Diagnostics-------------------------------------------
    idt_src=0  ! output print level - 0: print in any case
    CALL dbg_print('smoothed Laplac Diff.'     ,K_h                     ,str_module,idt_src, &
      & in_subset=p_patch%edges%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE smooth_lapl_diff
  !
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07)
  !
  !
  SUBROUTINE construct_ho_params(p_patch, params_oce)

    TYPE(t_patch), INTENT(IN)         :: p_patch
    TYPE (t_ho_params), INTENT(INOUT) :: params_oce

    ! Local variables
    INTEGER   :: ist, i,jtrc
    INTEGER   :: alloc_cell_blocks, nblks_e

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':construct_ho_physics'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct hydro ocean physics')

    CALL new_var_list(ocean_params_list, 'ocean_params_list', patch_id=p_patch%id)
    CALL default_var_list_settings( ocean_params_list,         &
      &                             lrestart=.TRUE.,           &
      &                             restart_type=FILETYPE_NC2, &
      &                             model_type='oce' )

    ! determine size of arrays
    alloc_cell_blocks = p_patch%alloc_cell_blocks
    nblks_e = p_patch%nblks_e

    CALL add_var(ocean_params_list, 'K_veloc_h', params_oce%K_veloc_h , GRID_UNSTRUCTURED_EDGE,&
    &            ZA_DEPTH_BELOW_SEA, &
    &            t_cf_var('K_veloc_h', 'kg/kg', 'horizontal velocity diffusion', DATATYPE_FLT32),&
    &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))

    CALL add_var(ocean_params_list, 'A_veloc_v', params_oce%A_veloc_v , GRID_UNSTRUCTURED_EDGE,&
    &            ZA_DEPTH_BELOW_SEA_HALF, &
    &            t_cf_var('A_veloc_v', 'kg/kg', 'vertical velocity diffusion', DATATYPE_FLT32),&
    &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev+1,nblks_e/),in_group=groups("oce_physics"))


    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL add_var(ocean_params_list, 'K_tracer_h', params_oce%K_tracer_h , &
      &            GRID_UNSTRUCTURED_EDGE, ZA_DEPTH_BELOW_SEA, &
      &            t_cf_var('K_tracer_h', '', '1:temperature 2:salinity', DATATYPE_FLT32),&
      &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
      &            ldims=(/nproma,n_zlev,nblks_e,no_tracer/), &
      &            lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(ocean_params_list, 'A_tracer_v', params_oce%A_tracer_v , &
      &            GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA_HALF, &
      &            t_cf_var('A_tracer_v', '', '1:temperature 2:salinity', DATATYPE_FLT32),&
      &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
      &            ldims=(/nproma,n_zlev+1,alloc_cell_blocks,no_tracer/), &
      &            lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! Reference to individual tracer, for I/O

      ALLOCATE(params_oce%tracer_h_ptr(no_tracer))
      ALLOCATE(params_oce%tracer_v_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_params_list, 'K_tracer_h',&
                    & 'K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc)),     &
                    & params_oce%tracer_h_ptr(jtrc)%p,                             &
                    & GRID_UNSTRUCTURED_EDGE, ZA_DEPTH_BELOW_SEA,               &
                    & t_cf_var('K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc)), &
                    &          'kg/kg', &
                    &          TRIM(oce_config%tracer_longnames(jtrc))//'(K_tracer_h_)', &
                    &          DATATYPE_FLT32), &
                    & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
                    & ldims=(/nproma,n_zlev,nblks_e/),in_group=groups("oce_physics"))
        CALL add_ref( ocean_params_list, 'A_tracer_v',&
                    & 'A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)),     &
                    & params_oce%tracer_h_ptr(jtrc)%p,                             &
                    & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA_HALF,            &
                    & t_cf_var('A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)), &
                    &          'kg/kg', &
                    &          TRIM(oce_config%tracer_longnames(jtrc))//'(A_tracer_v)', &
                    &          DATATYPE_FLT32), &
                    & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
                    & ldims=(/nproma,n_zlev+1,alloc_cell_blocks/),in_group=groups("oce_physics"))

      END DO
!TODO     use the following code, if add_var support 1d arrays:
!TODO     CALL add_var(ocean_params_list, 'K_tracer_h_back', params_oce%K_tracer_h_back , &
!TODO     &            GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, &
!TODO     &            t_cf_var('K_tracer_h_back', '', '1:temperature 2:salinity', DATATYPE_FLT32),&
!TODO     &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
!TODO     &            ldims=(/ no_tracer /))
!TODO     CALL add_var(ocean_params_list, 'A_tracer_v_back', params_oce%A_tracer_v_back , &
!TODO     &            GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
!TODO     &            t_cf_var('A_tracer_v_back', '', '1:temperature 2:salinity', DATATYPE_FLT32),&
!TODO     &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
!TODO     &            ldims=(/no_tracer/))
    ENDIF ! no_tracer > 0


    ALLOCATE(params_oce%K_tracer_h_back(no_tracer), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for horizontal background tracer diffusion failed')
     END IF

    ALLOCATE(params_oce%A_tracer_v_back(no_tracer), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for vertical tracer background diffusion failed')
     END IF


    DO i=1,no_tracer
      params_oce%K_tracer_h_back(i)  = 0.0_wp
      params_oce%A_tracer_v_back(i)  = 0.0_wp
    END DO
  END SUBROUTINE construct_ho_params
  !-------------------------------------------------------------------------
  !
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07)
  !
  !
  SUBROUTINE destruct_ho_params(params_oce)

    TYPE (t_ho_params), INTENT(INOUT) :: params_oce

    ! Local variables
    INTEGER :: ist
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':destruct_ho_physics'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'destruct hydro ocean physics')

    CALL delete_var_list(ocean_params_list)

    DEALLOCATE(params_oce%K_tracer_h_back, STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for horizontal tracer &
         &background iffusion failed')
     END IF

    DEALLOCATE(params_oce%A_tracer_v_back, STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for vertical background &
         &temperaure diffusion failed')
     END IF
  END SUBROUTINE destruct_ho_params 
  !-------------------------------------------------------------------------
  !
  !>
  !! Update of parameters
  !!
  !! Update of ocean physics: This routine is used used only if time-dependent
  !! changes of physical parametrizations.
  !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
  !! Dependent on the local Richardson number the diffusivity are calculated
  !! (Large & Gent JPO 29, (1999), 449-464).
  !! The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean
  !! Modelling 5, 2003).
  !! The notational convention is also taken from this paper( cf. eqs (14) and (19)).
  !! What is missing is the fractional ice cover (see eqs. (15-16)).
  !! Eq. (18) is the Redi part that is not implemented, yet
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-02)
  !
  !! mpi parallelized, sync required:
  !!                   params_oce%A_tracer_v, params_oce%A_veloc_v
  SUBROUTINE update_ho_params(p_patch_3D, p_os, p_sfc_flx, params_oce, calc_density_func)

    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: p_patch_3D
    TYPE(t_hydro_ocean_state), TARGET           :: p_os
    TYPE(t_sfc_flx),   INTENT(IN)               :: p_sfc_flx
    TYPE(t_ho_params), INTENT(INOUT)            :: params_oce
    INTERFACE !This contains the function version of the actual EOS as chosen in namelist
      FUNCTION calc_density_func(tpot, sal, press) RESULT(rho)
        USE mo_kind, ONLY: wp
        REAL(wp), INTENT(IN) :: tpot
        REAL(wp), INTENT(IN) :: sal
        REAL(wp), INTENT(IN) :: press
        REAL(wp) :: rho
     ENDFUNCTION calc_density_func
    END INTERFACE

    ! Local variables
    INTEGER  :: jc, jb, je,jk, itracer
   !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER  :: ilc1, ibc1, ilc2,ibc2
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: i_startidx_e, i_endidx_e
    INTEGER  :: z_dolic

    REAL(wp) :: z_rho_up, z_rho_down, z_stabio, z_shear_c, z_av0, z_dv0
    REAL(wp) :: z_vert_density_grad_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_vert_density_grad_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_Ri_c               (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_Ri_e               (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_c                  (nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    REAL(wp) :: dz_inv, z_lambda_frac

    !Below is a set of variables and parameters for tracer and velocity
    !REAL(wp), PARAMETER :: z_beta            = 0.6_wp
    !REAL(wp), PARAMETER :: z_one_minus_beta  = 0.4_wp
    REAL(wp), PARAMETER :: z_lambda          = 0.03_wp   !  wind mixing stability parameter (16)
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_T            = 5.0_wp    !  PP diffusivity tuning constant
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp    !  PP viscosity tuning constant
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: z_grav_rho, z_inv_rho_ref, z_press, press, A_T_tmp, z_s1, z_s2
    REAL(wp) :: density_grad_e, mean_z_r
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_cells!, cells_in_domain
    TYPE(t_patch), POINTER        :: p_patch

    !-------------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => p_patch%edges%in_domain
    !cells_in_domain => p_patch%cells%in_domain
    all_cells       => p_patch%cells%all

    !-------------------------------------------------------------------------
    z_av0 = richardson_veloc
    z_dv0 = richardson_tracer

    !-------------------------------------------------------------------------
    IF (l_constant_mixing) THEN
      !nothing to do!In sbr init_ho_params (see above)
      !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, itracer) is already
      !initialzed with params_oce%A_tracer_v_back(itracer)
      !and velocity diffusion coefficient

    ! prepare independent logicals for PP and convection parametrizations - not yet activated
    ! IF (.NOT. (l_convection .AND. l_pp_scheme)) THEN

    ELSE !(.NOT.l_constant_mixing)
      ! Attention: with l_constant_mixing=.true. there is no application of
      ! convective mixing parameters in case of instability
      ! max_vert_diff_veloc / max_vert_diff_trac
      ! control of convective and constant mixing should be independent

      !Density gradient for 1st level not calulated yet, but required below. Following
      !parxis in MPI-OM we set first layer equal to second layer (cf MPI-OM mo_ocean_vertical_mixing)
      !z_vert_density_grad_c(:,1,:) =  z_vert_density_grad_c(:,2,:)
      z_vert_density_grad_c(:,:,:) = 0.0_wp
      z_vert_density_grad_e(:,:,:) = 0.0_wp
      z_Ri_c(:,:,:)                = 0.0_wp
      z_Ri_e(:,:,:)                = 0.0_wp
      z_stabio                     = 0.0_wp
      z_shear_c                    = 0.0_wp
      z_s1                         = sal_ref
      z_s2                         = sal_ref
      z_grav_rho                   = grav/rho_ref
      z_inv_rho_ref                = 1.0_wp/rho_ref

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

        DO jc = i_startidx_c, i_endidx_c
          DO jk = 2, p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

             !dz_inv = p_os%p_diag%inv_prism_center_dist_c(jc,jk,jb) !1.0_wp/v_base%del_zlev_i(jk)
             !dz_inv = p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,jk,jb)
             !This calculates the localshear at cells
             ! - add small epsilon to avoid division by zero
          z_shear_c = dbl_eps + sum((p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)**2)

          z_press = p_patch_3D%p_patch_1D(1)%zlev_i(jk)*rho_ref*SItodBar
            press = p_patch_3D%p_patch_1D(1)%zlev_m(jk)*rho_ref*SItodBar

          !salinity at upper and lower cell
          IF(no_tracer >= 2) THEN
            z_s1 = p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,2)
            z_s2 = p_os%p_prog(nold(1))%tracer(jc,jk,jb,2) !TODO: this is the current cell, not the lower one!
          ENDIF
          !density of upper and lower cell w.r.t.to pressure at intermediate level
          z_rho_up           = calc_density_func(p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,1), z_s1, z_press)
          z_rho_down         = calc_density_func(p_os%p_prog(nold(1))%tracer(jc,jk,jb,1), z_s2, z_press)

            ! comments from MPIOM
            !! calculate vertical density stabio gradient between upper and lower box
            !! vertical density gradient 1/delta_z * (rho(k)-rho(k-1))

            !! stabio > 0 stable   stratification  (lower layer is havier)  => vert_density_grad  > 0
            !! stabio < 0 instable stratification  (lower layer is lighter) => vert_density_grad  < 0
          z_stabio  = z_rho_down-z_rho_up

             ! z_stabio  = dz_inv*0.5_wp*(z_rho_up(jc,jk,jb)-z_rho_down(jc,jk,jb))
             ! #slo# - think once more about 0.5, and this line in mo_convection of MPIOM:
             ! rhoo(:, j, k-1) = 0.5_wp * (rhoo(:, j, k-1) + rhuppo(:))

          z_vert_density_grad_c(jc,jk,jb) = dbl_eps + z_stabio

             ! taken from: G.R. Stuhne, W.R. Peltier / Journal of Computational Physics 213 (2006), p. 719
             ! Richardson number is positive for stable stratification (rho_down > rho_up)
             ! Richardson number is zero for unstable strat., see switch z_frac below
             ! The expression z_grav_rho*z_vert_density_grad_c/z_shear_c is the
             ! Buoyancy frequency:
             ! z_Ri_c(jc,jk,jb)=MAX(z_grav_rho*z_vert_density_grad_c(jc,jk,jb)/z_shear_c(jc,jk,jb)&
             !                     &,0.0_wp)
             !buoyance_frequence = z_grav_rho*z_vert_density_grad_c/z_shear_c
 !   buoyance_frequence            = z_grav_rho*z_vert_density_grad_c(jc,jk,jb)/z_shear_c(jc,jk,jb)
 !z_Ri_c(jc,jk,jb)=v_base%zlev_i(jk)*buoyance_frequence !TODO this created a difference in results!!
          z_Ri_c(jc,jk,jb)=p_patch_3D%p_patch_1D(1)%zlev_i(jk)*z_grav_rho*z_vert_density_grad_c(jc,jk,jb)/z_shear_c
          END DO
        END DO
      END DO

      !The tracer mixing coefficient at cell centers
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
          IF ( z_dolic >= MIN_DOLIC ) THEN
            DO jk = 2, z_dolic

              dz_inv        = p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,jk,jb)
              z_lambda_frac = z_lambda * dz_inv
              ! p_os%p_diag%inv_prism_center_dist_c(jc,jk,jb) !1.0_wp/v_base%del_zlev_i(jk)

              !calculate vertical tracer mixing based on local Richardson number
              DO itracer = 1, no_tracer

                !Store old diffusivity
                !z_A_tracer_v_old = params_oce%A_tracer_v(jc,jk,jb,itracer)

                !! vert_density_grad == 0 or below threshold ('semi-stable'): use background value
                IF ( ABS(z_vert_density_grad_c(jc,jk,jb)) < z_threshold ) THEN
                  params_oce%A_tracer_v(jc,jk,jb, itracer) = params_oce%A_tracer_v_back(itracer)

                !! vert_density_grad  < 0 or smaler than threshold ('unstable'): use convective mixing parameter
                ELSE IF (z_vert_density_grad_c(jc,jk,jb) < -z_threshold ) THEN
                ! IF (l_convection) &  ! else take calculated or background values
                  params_oce%A_tracer_v(jc,jk,jb, itracer) = MAX_VERT_DIFF_TRAC

                !! vert_density_grad  >= 0 or greater-o-equal then threshold ('stable'): use calculated value
                !! again taken from: G.R. Stuhne, W.R. Peltier / Journal of Computational Physics 213 (2006), p. 719
                ELSE

                  ! IF (l_pp_scheme) THEN

                  z_Ri_c(jc,jk,jb) = MAX(z_Ri_c(jc,jk,jb),0.0_wp)

                  ! This follows (19) in Marsland et al. but with lambda_D=1.0
                  ! but Marsland et al. use lambda_D=0.6 (z_beta, see below), which is not implemented here
                  ! The small wind mixing term D_w=5.0e-4 is missing
                  A_T_tmp = params_oce%A_tracer_v_back(itracer) + &
                    & z_dv0/((1.0_wp + z_c1_T * z_Ri_c(jc,jk,jb))**3)

                  params_oce%A_tracer_v(jc,jk,jb, itracer) = A_T_tmp

                  IF (l_wind_mixing) THEN
             
                    ! This is (16) in Marsland et al. and identical to treatment of velocity
                    ! but it allows to use different parameters

                    !z_A_W_T(jc,jk,jb) = z_A_W_T (jc,jk-1,jb)                      &
                    !&*(z_lambda_frac*exp(-v_base%del_zlev_i(jk)/z_0_T))&
                    !&/(z_lambda_frac+z_vert_density_grad_c)

                    ! For positive Richardson number set vertical mixing coefficient to maximal number 
                    IF(z_Ri_c(jc,jk,jb) <= 0.0_wp)THEN
                      params_oce%A_tracer_v(jc,jk,jb,itracer)                        &
                      & = params_oce%A_tracer_v_back(itracer)!z_one_minus_beta* z_A_tracer_v_old                            &
                      !& + z_beta*(params_oce%A_tracer_v_back(itracer)                 &
                      !& + params_oce%A_tracer_v_back(itracer)/(1.0_wp+z_c1_T*z_Ri_c)**3 &
                      !& + z_A_W_T(jc,jk,jb))
                      !write(123,*)'neg T-Ri number',jc,jk,jb,params_oce%A_tracer_v(jc,jk,jb, itracer)
                    ELSEIF(z_Ri_c(jc,jk,jb) > 0.0_wp)THEN
                      !write(123,*)'pos T-Ri number',jc,jk,jb,z_max_diff_T
                      params_oce%A_tracer_v(jc,jk,jb, itracer) = params_oce%A_tracer_v_back(itracer)
                                                                   !z_max_diff_T!params_oce%A_tracer_v_back(itracer)
                    END IF
                  ENDIF  !  l_wind_mixing
                  
                  ! ELSE - nothing to calculate - background values
                  ! ENDIF  !  pp_scheme

                END IF
              ENDDO
            END DO
          ENDIF
        END DO
      END DO
      !--------------------------------------------

      !--------------------------------------------
      ! Calculate params_oce%A_veloc_v:
      ! use mean values between the two cells; change to min, max if required
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx_e, i_endidx_e)
        DO je = i_startidx_e, i_endidx_e

          ilc1 = p_patch%edges%cell_idx(je,jb,1)
          ibc1 = p_patch%edges%cell_blk(je,jb,1)
          ilc2 = p_patch%edges%cell_idx(je,jb,2)
          ibc2 = p_patch%edges%cell_blk(je,jb,2)

          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_e(je, jb)
          DO jk = 2, z_dolic
            ! Set to zero for land + boundary locations edges
            IF (p_patch_3D%lsm_e(je,jk,jb) > SEA) THEN
               params_oce%A_veloc_v(je,jk,jb) = 0.0_wp
            ELSE
              ! TODO: the following expect equally sized cells
              ! compute density gradient at edges
              density_grad_e = 0.5_wp * &
                (z_vert_density_grad_c(ilc1,jk,ibc1) + z_vert_density_grad_c(ilc2,jk,ibc2))

              !! density gradient smaller then threshold ('semi-stable'): use background value
              IF ( ABS(density_grad_e) < z_threshold ) THEN
                params_oce%A_veloc_v(je,jk,jb) = params_oce%A_veloc_v_back

              !! vert_density_grad below that neg. theshold ('unstable stratification'): use max value
              ELSE IF (density_grad_e < -z_threshold ) THEN
                params_oce%A_veloc_v(je,jk,jb) = MAX_VERT_DIFF_VELOC

              !! vert_density_grad  > 0 stable stratification: use calculated value
              ELSE IF (density_grad_e > z_threshold ) THEN
                ! TODO: the following expect equally sized cells
                 mean_z_r = MAX(0.5_wp * (z_Ri_c(ilc1,jk,ibc1) + z_Ri_c(ilc2,jk,ibc2)),0.0_wp)
                 params_oce%A_veloc_v(je,jk,jb) = &
                   & params_oce%A_veloc_v_back +  &
                   & z_av0 /                      &
                   & ((1.0_wp + z_c1_v * mean_z_r)**2)

              ENDIF
            ENDIF
          END DO ! jk = 2, z_dolic
        ENDDO ! je = i_startidx_e, i_endidx_e
      ENDDO ! jb = edges_in_domain%start_block, edges_in_domain%end_block
    ENDIF !l_constant_mixing

    ! Sync the results, the A_tracer_v is only for checking
    DO itracer = 1, no_tracer
      CALL sync_patch_array(SYNC_C,p_patch,params_oce%A_tracer_v(:,:,:,itracer))
    END DO
    CALL sync_patch_array(SYNC_E,p_patch,params_oce%A_veloc_v(:,:,:))

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    idt_src=4  ! output print level (1-5, fix)
    CALL dbg_print('UpdPar: p_vn%x(1)'         ,p_os%p_diag%p_vn%x(1)    ,str_module,idt_src, &
      & in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdPar: p_vn%x(2)'         ,p_os%p_diag%p_vn%x(2)    ,str_module,idt_src, &
      & in_subset=p_patch%cells%owned)
    idt_src=2  ! output print level (1-5, fix)
    DO itracer = 1, no_tracer
      z_c(:,:,:)=params_oce%A_tracer_v(:,:,:,itracer)
      CALL dbg_print('UpdPar FinalTracerMixing'  ,z_c                    ,str_module,idt_src, &
        & in_subset=p_patch%cells%owned)
    ENDDO
    CALL dbg_print('UpdPar FinalVelocMixing'   ,params_oce%A_veloc_v     ,str_module,idt_src, &
        & in_subset=p_patch%edges%owned)
    !---------------------------------------------------------------------
  END SUBROUTINE update_ho_params
END MODULE mo_oce_physics
