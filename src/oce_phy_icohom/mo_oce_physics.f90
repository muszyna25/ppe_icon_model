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
  &                               biharmonic_diffusion_factor, &
  &                               richardson_factor_tracer, richardson_factor_veloc,      &
  &                               l_constant_mixing, l_smooth_veloc_diffusion
USE mo_parallel_config,     ONLY: nproma
USE mo_model_domain,        ONLY: t_patch, t_patch_3D_oce
USE mo_impl_constants,      ONLY: success, max_char_length, MIN_DOLIC, SEA
USE mo_exception,           ONLY: message, finish
USE mo_util_dbg_prnt,       ONLY: dbg_print
USE mo_oce_state,           ONLY: t_hydro_ocean_state, oce_config!, v_base
USE mo_physical_constants,  ONLY: grav, rho_ref, SItodBar
USE mo_math_constants,      ONLY: dbl_eps
USE mo_dynamics_config,     ONLY: nold!, nnew
USE mo_sea_ice_types,       ONLY: t_sfc_flx
USE mo_run_config,          ONLY: dtime
USE mo_linked_list,         ONLY: t_var_list
USE mo_var_list,            ONLY: add_var,                  &
  &                               new_var_list,             &
  &                               delete_var_list,          &
  &                               default_var_list_settings,&
  &                               add_ref
USE mo_cf_convention
USE mo_grib2
USE mo_cdi_constants
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
    TYPE(t_patch_3D_oce ),TARGET, INTENT(IN) :: p_patch_3D
    TYPE (t_ho_params)                          :: p_phys_param

    ! Local variables
    INTEGER  :: i, i_no_trac 
    INTEGER  :: je,jb
    INTEGER  :: i_startidx_e, i_endidx_e
    REAL(wp) :: z_lower_bound_diff
    REAL(wp) :: z_largest_edge_length ,z_diff_multfac, z_diff_efdt_ratio
    REAL(wp), PARAMETER :: N_POINTS_IN_MUNK_LAYER = 1.0_wp
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_patch), POINTER        :: p_patch 
    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    all_edges => p_patch%edges%all
    !-------------------------------------------------------------------------
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
        write(0,*)'lower bound of diffusivity:',z_lower_bound_diff

      CASE(2)!calculate uniform viscosity coefficient, according to Munk criterion

        p_phys_param%K_veloc_h(:,:,:) = 3.82E-12_wp&
        &*(N_POINTS_IN_MUNK_LAYER*z_largest_edge_length)**3

      CASE(3)!
        DO jb = all_edges%start_block, all_edges%end_block
          CALL get_index_range(all_edges, jb, i_startidx_e, i_endidx_e)
          DO je = i_startidx_e, i_endidx_e
            !calculate lower bound for diffusivity
            !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29)) 
            p_phys_param%K_veloc_h(je,:,jb) = 3.82E-12_wp&
            &*(N_POINTS_IN_MUNK_LAYER*p_patch%edges%primal_edge_length(je,jb))**3
          END DO
        END DO
      END SELECT
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

write(*,*)'max-min coeff',z_diff_multfac, maxval(p_phys_param%K_veloc_h(:,1,:)),&
&minval(p_phys_param%K_veloc_h(:,1,:))


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
    REAL(wp), INTENT(OUT)              :: lower_bound_diff

    ! Local variables
    REAL(wp), PARAMETER :: N_POINTS_IN_MUNK_LAYER = 1.0_wp
    REAL(wp)            :: z_largest_edge_length
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain
    !-------------------------------------------------------------------------
    edges_in_domain => p_patch%edges%in_domain

    ! Get the largest edge length globally
    z_largest_edge_length = global_max(MAXVAL(p_patch%edges%primal_edge_length))

    !calculate lower bound for diffusivity: The factor cos(lat) is omitted here, because of
    !equatorial reference (cf. Griffies, eq.  (18.29)) 
    lower_bound_diff = 3.82E-12_wp*(N_POINTS_IN_MUNK_LAYER*z_largest_edge_length)**3

  END SUBROUTINE calc_lower_bound_veloc_diff

  !-------------------------------------------------------------------------
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
  !
  !! mpi parallelized, no sync required
  SUBROUTINE smooth_lapl_diff( p_patch,p_patch_3D, K_h )
   TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
   TYPE(t_patch_3D_oce ),TARGET, INTENT(IN)   :: p_patch_3D
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
    CALL dbg_print('smoothed Laplac Diff.'     ,k_h                     ,str_module,idt_src)
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
    INTEGER   :: nblks_c, nblks_e

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
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    CALL add_var(ocean_params_list, 'K_veloc_h', params_oce%K_veloc_h , GRID_UNSTRUCTURED_EDGE,&
    &            ZA_DEPTH_BELOW_SEA, &
    &            t_cf_var('K_veloc_h', 'kg/kg', 'horizontal velocity diffusion', DATATYPE_FLT32),&
    &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    CALL add_var(ocean_params_list, 'A_veloc_v', params_oce%A_veloc_v , GRID_UNSTRUCTURED_EDGE,&
    &            ZA_DEPTH_BELOW_SEA_HALF, &
    &            t_cf_var('A_veloc_v', 'kg/kg', 'vertical velocity diffusion', DATATYPE_FLT32),&
    &            t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev+1,nblks_e/))


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
      &            ldims=(/nproma,n_zlev+1,nblks_c,no_tracer/), &
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
                    & ldims=(/nproma,n_zlev,nblks_e/))
        CALL add_ref( ocean_params_list, 'A_tracer_v',&
                    & 'A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)),     &
                    & params_oce%tracer_h_ptr(jtrc)%p,                             &
                    & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA_HALF,            &
                    & t_cf_var('A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)), &
                    &          'kg/kg', &
                    &          TRIM(oce_config%tracer_longnames(jtrc))//'(A_tracer_v)', &
                    &          DATATYPE_FLT32), &
                    & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),&
                    & ldims=(/nproma,n_zlev+1,nblks_c/))

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
  SUBROUTINE update_ho_params(p_patch_3D, p_os, p_sfc_flx, params_oce, calc_density)

    TYPE(t_patch_3D_oce ),TARGET, INTENT(IN) :: p_patch_3D
    TYPE(t_hydro_ocean_state), TARGET           :: p_os
    TYPE(t_sfc_flx),   INTENT(IN)               :: p_sfc_flx
    TYPE(t_ho_params), INTENT(INOUT)            :: params_oce
    INTERFACE !This contains the function version of the actual EOS as chosen in namelist
      FUNCTION calc_density(tpot, sal, press) RESULT(rho)
        USE mo_kind, ONLY: wp
        REAL(wp), INTENT(IN) :: tpot
        REAL(wp), INTENT(IN) :: sal
        REAL(wp), INTENT(IN) :: press
        REAL(wp) :: rho
     ENDFUNCTION calc_density
    END INTERFACE

    ! Local variables
    INTEGER  :: jc, jb, je,jk, itracer
   !INTEGER  :: ile1, ibe1,ile2, ibe2,ile3, ibe3
    INTEGER  :: ilc1, ibc1, ilc2,ibc2!, jj, ible,idxe
    INTEGER  :: i_startidx_c, i_endidx_c
    INTEGER  :: i_startidx_e, i_endidx_e
    REAL(wp) :: z_vert_density_grad_c(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_vert_density_grad_e(nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_stabio             (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_shear_c            (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_rho_up             (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_rho_down           (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_Ri_c               (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_c)
    REAL(wp) :: z_Ri_e               (nproma,n_zlev,p_patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp) :: z_c                  (nproma,n_zlev+1,p_patch_3D%p_patch_2D(1)%nblks_c)

!    REAL(wp) :: dz_inv
!    REAL(wp) :: z_rho_up_c1, z_rho_down_c1,z_rho_up_c2, z_rho_down_c2
!   REAL(wp) :: z_lambda_frac 
!   REAL(wp) :: z_A_tracer_v_old!, z_A_veloc_v_old
    INTEGER  :: z_dolic

    !Below is a set of variables and parameters for tracer and velocity
    !REAL(wp), PARAMETER :: z_beta            = 0.6_wp
    !REAL(wp), PARAMETER :: z_one_minus_beta  = 0.4_wp
    REAL(wp), PARAMETER :: z_lambda          = 0.05_wp
    REAL(wp), PARAMETER :: z_0               = 40.0_wp
    REAL(wp), PARAMETER :: z_c1_T            = 5.0_wp
    REAL(wp), PARAMETER :: z_c1_v            = 5.0_wp
    REAL(wp)            :: z_av0             = 0.5E-2_wp ! later set via nml richardson_factor_veloc
    REAL(wp)            :: z_dv0             = 0.5E-2_wp ! later set via nml richardson_factor_tracer
    REAL(wp), PARAMETER :: z_threshold       = 5.0E-8_wp
    REAL(wp) :: z_grav_rho, z_inv_rho_ref!, z_stabio
    REAL(wp) :: z_press!, z_frac
    REAL(wp) :: A_T_tmp!, A_v_tmp
    REAL(wp) :: z_s1, z_s2, density_grad_e, mean_z_r
    ! REAL(wp) :: tmp_communicate_c(nproma,p_patch%nblks_c)
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: edges_in_domain,cells_in_domain,all_cells
    TYPE(t_patch), POINTER        :: p_patch 
    !-------------------------------------------------------------------------
    z_av0 = richardson_factor_veloc
    z_dv0 = richardson_factor_tracer
    !-------------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    edges_in_domain => p_patch%edges%in_domain
    cells_in_domain => p_patch%cells%in_domain
    all_cells       => p_patch%cells%all

    IF (l_constant_mixing) THEN
      !nothing to do!In sbr init_ho_params (see above)
      !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, itracer) is already
      !initialzed with params_oce%A_tracer_v_back(itracer)
      !and velocity diffusion coefficient

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
      z_stabio (:,:,:)             = 0.0_wp
      z_shear_c(:,:,:)             = 0.0_wp
      z_s1                         = 0.0_wp
      z_s2                         = 0.0_wp
      z_grav_rho                   = grav/rho_ref
      z_inv_rho_ref                = 1.0_wp/rho_ref

    ! #slo# this part is taken from revision r5966
    ! !Following MPI-OM (cf. vertical mixing sbr)
    ! z_w_T = CWT/6.0_wp**3
    ! z_w_v = CWA/6.0_wp**3
    ! 
    ! !The wind part
    ! DO jb = i_startblk_e, i_endblk_e
    !   CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
    !   &                   rl_start_e, rl_end_e)
    !   DO je = i_startidx_e, i_endidx_e 
    !     IF ( v_base%lsm_e(je,1,jb) <= sea_boundary ) THEN
    ! 
    !      ilc1 = p_patch%edges%cell_idx(je,jb,1)
    !      ibc1 = p_patch%edges%cell_blk(je,jb,1)
    !      ilc2 = p_patch%edges%cell_idx(je,jb,2)
    !      ibc2 = p_patch%edges%cell_blk(je,jb,2)
    ! 
    !       !This is (15) in Marsland et al. 
    !       z_10m_wind_e(je,1,jb)= SQRT(&
    !       &0.5_wp*(DOT_PRODUCT(p_sfc_flx%forc_wind_cc(ilc1,ibc1)%x,     &
    !       &                    p_sfc_flx%forc_wind_cc(ilc1,ibc1)%x)     &
    !       &       +DOT_PRODUCT(p_sfc_flx%forc_wind_cc(ilc2,ibc2)%x,     &
    !       &                    p_sfc_flx%forc_wind_cc(ilc2,ibc2)%x)))**3
    !       z_A_W_v (je,1,jb) = z_w_v*z_10m_wind_e(je,1,jb)
    !     ENDIF
    !   END DO
    ! END DO
    !
    ! DO jb = i_startblk_c, i_endblk_c
    !   CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    !   &                   rl_start_c, rl_end_c)
    !   DO jc = i_startidx_c, i_endidx_c 
    !     IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN
    !       !This is (15) in Marsland et al. 
    !       z_10m_wind_c(jc,1,jb)= SQRT(DOT_PRODUCT(p_sfc_flx%forc_wind_cc(jc,jb)%x,&
    !                                              &p_sfc_flx%forc_wind_cc(jc,jb)%x))**3
    !       z_A_W_T (jc,1,jb) = z_w_T*z_10m_wind_c(jc,1,jb)
    !     ENDIF
    !   END DO
    ! END DO

      !Calculate Richardson number and vertical density gradient
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

        DO jc = i_startidx_c, i_endidx_c
          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
          IF ( z_dolic >= MIN_DOLIC ) THEN        
            DO jk = 2, z_dolic 

              !dz_inv = p_os%p_diag%inv_prism_center_dist_c(jc,jk,jb) !1.0_wp/v_base%del_zlev_i(jk)
              !dz_inv = p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,jk,jb)
              !This calculates the localshear at cells
              ! - add small epsilon to avoid division by zero
              z_shear_c(jc,jk,jb) = dbl_eps + &
                & sum((p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)**2)

              z_press = p_patch_3D%p_patch_1D(1)%zlev_i(jk)*rho_ref*SItodBar

              !salinity at upper and lower cell
              IF(no_tracer >= 2) THEN
                z_s1 = p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,2)
                z_s2 = p_os%p_prog(nold(1))%tracer(jc,jk,jb,2) !TODO: this is the current cell, not
                                                               !the lower one!
              ENDIF
              !density of upper and lower cell w.r.t.to pressure at intermediate level
              !z_rho_up(jc,jk,jb)   = p_os%p_diag%rho(jc,jk-1,jb)
              z_rho_up(jc,jk,jb)   = calc_density &
               & (p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,1), z_s1, z_press)

              !z_rho_down(jc,jk,jb)   = p_os%p_diag%rho(jc,jk,jb)
              z_rho_down(jc,jk,jb) = calc_density &
                & (p_os%p_prog(nold(1))%tracer(jc,jk,jb,1), z_s2, z_press) !TODO: local/current
                                                                           !density

              ! comments from MPIOM
              !! calculate vertical density stabio gradient between upper and lower box
              !! vertical density gradient 1/delta_z * (rho(k)-rho(k-1))

              !! stabio > 0 stable   stratification  (lower layer is havier)  => vert_density_grad  > 0
              !! stabio < 0 instable stratification  (lower layer is lighter) => vert_density_grad  < 0
              z_stabio(jc,jk,jb) = (z_rho_down(jc,jk,jb)-z_rho_up(jc,jk,jb))

              ! z_stabio  = dz_inv*0.5_wp*(z_rho_up(jc,jk,jb)-z_rho_down(jc,jk,jb))
              ! #slo# - think once more about 0.5, and this line in mo_convection of MPIOM:
              ! rhoo(:, j, k-1) = 0.5_wp * (rhoo(:, j, k-1) + rhuppo(:))

              z_vert_density_grad_c(jc,jk,jb) = dbl_eps+z_stabio(jc,jk,jb)

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
 z_Ri_c(jc,jk,jb)=p_patch_3D%p_patch_1D(1)%zlev_i(jk)*z_grav_rho*z_vert_density_grad_c(jc,jk,jb)/z_shear_c(jc,jk,jb)
            END DO
          ENDIF
        END DO
      END DO

      !The tracer mixing coefficient at cell centers
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c

          z_dolic = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
          IF ( z_dolic >= MIN_DOLIC ) THEN
            DO jk = 2, z_dolic

             !dz_inv = p_patch_3D%p_patch_1D(1)%inv_prism_center_dist_c(jc,jk,jb)
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
                  params_oce%A_tracer_v(jc,jk,jb, itracer) = MAX_VERT_DIFF_TRAC

                !! vert_density_grad  >= 0 or greater-o-equal then threshold ('stable'): use calculated value
                !! again taken from: G.R. Stuhne, W.R. Peltier / Journal of Computational Physics 213 (2006), p. 719
                ELSE
                  z_Ri_c(jc,jk,jb) = MAX(z_Ri_c(jc,jk,jb),0.0_wp)

                  ! This follows (19) in Marsland et al. but with lambda_D=1.0
                  ! but Marsland et al. use lambda_D=0.6 (z_beta, see below), which is not implemented here
                  ! The small wind mixing term D_w=5.0e-4 is missing
                  A_T_tmp = params_oce%A_tracer_v_back(itracer) + &
                    & z_dv0/((1.0_wp + z_c1_T * z_Ri_c(jc,jk,jb))**3)

                  params_oce%A_tracer_v(jc,jk,jb, itracer) = A_T_tmp

!               !This is (16) in Marsland et al. and identical to treatment of velocity
!               !but it allows to use different parameters
!               z_lambda_frac     = z_lambda_T/v_base%del_zlev_i(jk)
!               z_A_W_T(jc,jk,jb) = z_A_W_T (jc,jk-1,jb)                      &
!               &*(z_lambda_frac*exp(-v_base%del_zlev_i(jk)/z_0_T))&
!               &/(z_lambda_frac+z_vert_density_grad_c)
!             !For positive Richardson number set vertical mixing coefficient to maximal number 
!               IF(z_Ri_c <= 0.0_wp)THEN
!                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac)                        &
!                 & = params_oce%A_tracer_v_back(i_no_trac)!z_one_minus_beta* z_A_tracer_v_old                            &
!                 !& + z_beta*(params_oce%A_tracer_v_back(i_no_trac)                 &
!                 !& + params_oce%A_tracer_v_back(i_no_trac)/(1.0_wp+z_c1_T*z_Ri_c)**3 &
!                 !& + z_A_W_T(jc,jk,jb))
!              !write(123,*)'neg T-Ri number',jc,jk,jb,params_oce%A_tracer_v(jc,jk,jb, i_no_trac)
!               ELSEIF(z_Ri_c > 0.0_wp)THEN
! !                 write(123,*)'pos T-Ri number',jc,jk,jb,z_max_diff_T
!                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = params_oce%A_tracer_v_back(i_no_trac)!z_max_diff_T!params_oce%A_tracer_v_back(i_no_trac)
!               ENDIF

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
    CALL dbg_print('UpdPar: p_vn%x(1)'         ,p_os%p_diag%p_vn%x(1)    ,str_module,idt_src)
    CALL dbg_print('UpdPar: p_vn%x(2)'         ,p_os%p_diag%p_vn%x(1)    ,str_module,idt_src)
    idt_src=3  ! output print level (1-5, fix)
    CALL dbg_print('UpdPar: z_shear_c'         ,z_shear_c                ,str_module,idt_src)
    idt_src=2  ! output print level (1-5, fix)
    DO itracer = 1, no_tracer
      z_c(:,:,:)=params_oce%A_tracer_v(:,:,:,itracer)
      CALL dbg_print('UpdPar FinalTracerMixing'  ,z_c                    ,str_module,idt_src)
    ENDDO
    CALL dbg_print('UpdPar FinalVelocMixing'   ,params_oce%A_veloc_v     ,str_module,idt_src)
    !---------------------------------------------------------------------
  END SUBROUTINE update_ho_params
END MODULE mo_oce_physics
