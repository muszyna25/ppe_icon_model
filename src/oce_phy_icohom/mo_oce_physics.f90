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
  &                               CWA, CWT, HORZ_VELOC_DIFF_TYPE
USE mo_parallel_config,     ONLY: nproma
USE mo_model_domain,        ONLY: t_patch
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, min_rledge,&
  &                               sea_boundary, MIN_DOLIC
USE mo_exception,           ONLY: message, finish
USE mo_oce_index,           ONLY: print_mxmn, jkc, jkdim, ipl_src
USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base, oce_config
USE mo_physical_constants,  ONLY: grav, rho_ref, SItodBar
USE mo_loopindices,         ONLY: get_indices_c,get_indices_e
USE mo_math_constants,      ONLY: dbl_eps
USE mo_dynamics_config,     ONLY: nold, nnew
! USE mo_oce_forcing,         ONLY: t_sfc_flx
USE mo_sea_ice,             ONLY: t_sfc_flx

USE mo_linked_list,         ONLY: t_var_list
USE mo_var_list,            ONLY: add_var,                  &
  &                               new_var_list,             &
  &                               delete_var_list,          &
  &                               default_var_list_settings,&
  &                               add_ref
USE mo_cf_convention
USE mo_grib2
USE mo_cdi_constants

IMPLICIT NONE


PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

CHARACTER(len=*), PARAMETER :: this_mod_name = 'mo_oce_physics'

! Public interface
PUBLIC :: t_ptr3d, t_ho_params

!PUBLIC :: init_ho_physics
PUBLIC :: construct_ho_params
PUBLIC :: destruct_ho_params
PUBLIC :: init_ho_params
PUBLIC :: update_ho_params
PRIVATE :: calc_munk_based_lapl_diff

! variables
Type (t_var_list), PUBLIC :: ocean_params_list

TYPE t_ptr3d
  REAL(wp),POINTER :: p(:,:,:)  ! pointer to 3D (spatial) array
END TYPE t_ptr3d

! Parameters below appear directly in the ocean model/equation. They are eventually
! dynamically updated by using the "ocean-physics" structure. #slo# - not yet
TYPE t_ho_params

  ! diffusion coefficients for horizontal velocity, temp. and salinity, dim=(nproma,n_zlev,nblks_e)
  REAL(wp),POINTER ::    &
    &  K_veloc_h(:,:,:),  & ! coefficient of horizontal velocity diffusion
    &  K_tracer_h(:,:,:,:)     ! coefficient of horizontal tracer diffusion
  TYPE(t_ptr3d),ALLOCATABLE :: tracer_h_ptr(:)

  ! diffusion coefficients for vertical velocity, temp. and salinity, dim=(nproma,n_zlev+1,nblks_e)
  REAL(wp),POINTER ::    &
    &  A_veloc_v(:,:,:),  & ! coefficient of vertical velocity diffusion
    &  A_tracer_v(:,:,:,:)     ! coefficient of vertical tracer diffusion
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
  !
  SUBROUTINE init_ho_params(  ppatch, p_phys_param )
   TYPE(t_patch), INTENT(IN)  :: ppatch
   TYPE (t_ho_params)         :: p_phys_param 
    ! Local variables
    INTEGER :: i
    REAL(wp) :: z_lower_bound_diff
    !-------------------------------------------------------------------------
    !Init from namelist
    p_phys_param%K_veloc_h_back = k_veloc_h
    p_phys_param%A_veloc_v_back = k_veloc_v
    p_phys_param%K_veloc_h = k_veloc_h
    p_phys_param%A_veloc_v = k_veloc_v
    SELECT CASE(HORZ_VELOC_DIFF_TYPE)


    CASE(0)
      p_phys_param%K_veloc_h(:,:,:) = 0.0_wp

    CASE(1)
      CALL calc_lower_bound_veloc_diff(  ppatch, z_lower_bound_diff )
      IF(z_lower_bound_diff>p_phys_param%K_veloc_h_back)THEN

        CALL message ('init_ho_params','WARNING: Specified diffusivity&
                      & does not satisfy Munk criterion. This may lead&
                      &to stability problems for experiments with lateral boundaries')
      ENDIF
      p_phys_param%K_veloc_h(:,:,:) = p_phys_param%K_veloc_h_back
      write(*,*)'lower bound',z_lower_bound_diff

    CASE(2)

      CALL calc_munk_based_lapl_diff(ppatch,p_phys_param%K_veloc_h)

    CASE(3)

      CALL calc_munk_based_lapl_diff(ppatch,p_phys_param%K_veloc_h) 

    END SELECT

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

  END SUBROUTINE init_ho_params
  !-------------------------------------------------------------------------
  !
  !>
  !! Calculation of a lower bound for horizontal velocity diffusion of laplacian type
  !! that is required to have N (default =1) points in Munk layer. The lower bound is calculated
  !! with respect to the equator. 
  !!The code is based on  Griffies, Fundamentals of ocean climate modeling, sect 18, p. 413. 
  !!The lower bound is given in units [m^2/s].
  !!
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
  !
  !
  SUBROUTINE calc_lower_bound_veloc_diff(  p_patch, lower_bound_diff )
   TYPE(t_patch), INTENT(IN)  :: p_patch
   REAL(wp), INTENT(OUT)      :: lower_bound_diff
    ! Local variables
    REAL(wp), PARAMETER  :: N_POINTS_IN_MUNK_LAYER = 1.0_wp
    INTEGER  :: je,jb 
    INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start_e, rl_end_e
    REAL(wp) :: z_largest_edge_length
    !-------------------------------------------------------------------------
    rl_start_e   = 1
    rl_end_e     = min_rledge
    i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

    z_largest_edge_length = 0.0_wp
    DO jb = i_startblk_e, i_endblk_e
      CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
      &                   rl_start_e, rl_end_e)
      DO je = i_startidx_e, i_endidx_e 

        !calculate largest edges length within grid 
        IF(p_patch%edges%dual_edge_length(je,jb)>z_largest_edge_length)THEN
          z_largest_edge_length=p_patch%edges%primal_edge_length(je,jb)
        ENDIF
      END DO
    END DO

   !write(*,*)'largest edge length',z_largest_edge_length

    DO jb = i_startblk_e, i_endblk_e
      CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
      &                   rl_start_e, rl_end_e)
      DO je = i_startidx_e, i_endidx_e 

        !calculate lower bound for diffusivity
        !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29)) 
        lower_bound_diff = 3.82E-12_wp&
        &*(N_POINTS_IN_MUNK_LAYER*z_largest_edge_length)**3&
          &*cos( p_patch%edges%center(je,jb)%lat)!*deg2rad
      END DO
    END DO

  END SUBROUTINE calc_lower_bound_veloc_diff
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Initialisation 
  !!
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-08)
  !
  !
  SUBROUTINE calc_munk_based_lapl_diff( p_patch, K_h )
   TYPE(t_patch), INTENT(IN)  :: p_patch
   REAL(wp), INTENT(OUT)      :: K_h(:,:,:) 
    ! Local variables
    REAL(wp), PARAMETER  :: N_POINTS_IN_MUNK_LAYER = 1.0_wp
    INTEGER  :: je,jb 
    INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start_e, rl_end_e 
    REAL(wp) :: z_largest_edge_length
    !-------------------------------------------------------------------------
    rl_start_e   = 1
    rl_end_e     = min_rledge
    i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)


    SELECT CASE(HORZ_VELOC_DIFF_TYPE)
    CASE(2) 
      z_largest_edge_length = 0.0_wp
      DO jb = i_startblk_e, i_endblk_e
        CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
        &                   rl_start_e, rl_end_e)
        DO je = i_startidx_e, i_endidx_e 

          !calculate largest edges length within grid 
          IF(p_patch%edges%dual_edge_length(je,jb)>z_largest_edge_length)THEN
            z_largest_edge_length=p_patch%edges%primal_edge_length(je,jb)
          ENDIF
        END DO
      END DO
      !write(*,*)'largest edge length',z_largest_edge_length
      DO jb = i_startblk_e, i_endblk_e
        CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
        &                   rl_start_e, rl_end_e)
        DO je = i_startidx_e, i_endidx_e 

          !calculate lower bound for diffusivity
          !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29)) 
          K_h(je,:,jb) = 3.82E-12_wp&
          &*(N_POINTS_IN_MUNK_LAYER*z_largest_edge_length)**3&
          &*cos( p_patch%edges%center(je,jb)%lat)!*deg2rad
        END DO
      END DO

    CASE(3)

      DO jb = i_startblk_e, i_endblk_e
        CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
        &                   rl_start_e, rl_end_e)
        DO je = i_startidx_e, i_endidx_e 

          !calculate lower bound for diffusivity
          !The factor cos(lat) is omitted here, because of equatorial reference (cf. Griffies, eq. (18.29)) 
          K_h(je,:,jb) = 3.82E-12_wp&
          &*(N_POINTS_IN_MUNK_LAYER*p_patch%edges%primal_edge_length(je,jb))**3&
          &*cos( p_patch%edges%center(je,jb)%lat)!*deg2rad
! write(123,*)'visc',K_h(je,1,jb),&
! &3.82E-12_wp&
! &*(N_POINTS_IN_MUNK_LAYER*p_patch%edges%primal_edge_length(je,jb))**3,&
! &cos( p_patch%edges%center(je,jb)%lat)
        END DO
      END DO
    END SELECT

    ipl_src=1  ! output print level (1-5, fix)
    CALL print_mxmn('PHY diffusivity',1,K_h(:,:,:),n_zlev,p_patch%nblks_c,'per',ipl_src)

  END SUBROUTINE calc_munk_based_lapl_diff
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
  SUBROUTINE construct_ho_params(ppatch, params_oce)

    TYPE(t_patch), INTENT(IN)         :: ppatch
    TYPE (t_ho_params), INTENT(INOUT) :: params_oce

    ! Local variables
    INTEGER   :: ist, i,jtrc
    INTEGER   :: nblks_c, nblks_e

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':construct_ho_physics'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct hydro ocean physics')

    CALL new_var_list(ocean_params_list, 'ocean_params_list', patch_id=ppatch%id)
    CALL default_var_list_settings( ocean_params_list,         &
      &                             lrestart=.TRUE.,           &
      &                             restart_type=FILETYPE_NC2, &
      &                             model_type='oce' )

    ! determine size of arrays
    nblks_c = ppatch%nblks_c
    nblks_e = ppatch%nblks_e
    
    CALL add_var(ocean_params_list, 'K_veloc_h', params_oce%K_veloc_h , GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('K_veloc_h', 'kg/kg', 'horizontal velocity diffusion'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    CALL add_var(ocean_params_list, 'A_veloc_v', params_oce%A_veloc_v , GRID_UNSTRUCTURED_EDGE,&
    &            ZAXIS_DEPTH_BELOW_SEA, &
    &            t_cf_var('A_veloc_v', 'kg/kg', 'vertical velocity diffusion'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
    &            ldims=(/nproma,n_zlev+1,nblks_e/))


    !! Tracers
    IF ( no_tracer > 0 ) THEN
      CALL add_var(ocean_params_list, 'K_tracer_h', params_oce%K_tracer_h , &
      &            GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA, &
      &            t_cf_var('K_tracer_h', '', '1:temperature 2:salinity'),&
      &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
      &            ldims=(/nproma,n_zlev,nblks_e,no_tracer/), &
      &            lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      CALL add_var(ocean_params_list, 'A_tracer_v', params_oce%A_tracer_v , &
      &            GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA, &
      &            t_cf_var('A_tracer_v', '', '1:temperature 2:salinity'),&
      &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
      &            ldims=(/nproma,n_zlev+1,nblks_c,no_tracer/), &
      &            lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ! Reference to individual tracer, for I/O

      ALLOCATE(params_oce%tracer_h_ptr(no_tracer))
      ALLOCATE(params_oce%tracer_v_ptr(no_tracer))
      DO jtrc = 1,no_tracer
        CALL add_ref( ocean_params_list, 'K_tracer_h',&
                    & 'K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc)),     &
                    & params_oce%tracer_h_ptr(jtrc)%p,                             &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_DEPTH_BELOW_SEA,            &
                    & t_cf_var('K_tracer_h_'//TRIM(oce_config%tracer_names(jtrc)), &
                    &          'kg/kg', &
                    &          TRIM(oce_config%tracer_longnames(jtrc))//'(K_tracer_h_)'), &
                    & t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
                    & ldims=(/nproma,n_zlev,nblks_e/))
        CALL add_ref( ocean_params_list, 'A_tracer_v',&
                    & 'A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)),     &
                    & params_oce%tracer_h_ptr(jtrc)%p,                             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_DEPTH_BELOW_SEA,            &
                    & t_cf_var('A_tracer_v_'//TRIM(oce_config%tracer_names(jtrc)), &
                    &          'kg/kg', &
                    &          TRIM(oce_config%tracer_longnames(jtrc))//'(A_tracer_v)'), &
                    & t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
                    & ldims=(/nproma,n_zlev+1,nblks_c/))

      END DO
!TODO     use the following code, if add_var support 1d arrays:
!TODO     CALL add_var(ocean_params_list, 'K_tracer_h_back', params_oce%K_tracer_h_back , &
!TODO     &            GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, &
!TODO     &            t_cf_var('K_tracer_h_back', '', '1:temperature 2:salinity'),&
!TODO     &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE),&
!TODO     &            ldims=(/ no_tracer /))
!TODO     CALL add_var(ocean_params_list, 'A_tracer_v_back', params_oce%A_tracer_v_back , &
!TODO     &            GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, &
!TODO     &            t_cf_var('A_tracer_v_back', '', '1:temperature 2:salinity'),&
!TODO     &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL),&
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
  !!(Large & Gent JPO 29, (1999), 449-464).
  !!The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean Modelling, 2002).
  !!The notationnal convection is also taken from this paper( cf. eqs (14)-(19)).
  !! What is missing is the fractional ice cover (see eqs. (15-16)).
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-02)
  !
  !
 SUBROUTINE update_ho_params(p_patch, p_os, p_sfc_flx, params_oce, calc_density)
   TYPE(t_patch),     INTENT(IN)     :: p_patch
   TYPE(t_hydro_ocean_state), TARGET :: p_os
   TYPE(t_sfc_flx),   INTENT(INOUT)  :: p_sfc_flx
   TYPE(t_ho_params), INTENT(INOUT)  :: params_oce  
INTERFACE !This contains the function version of the actual EOS as chosen in namelist
  FUNCTION calc_density(tpot, sal, press) RESULT(rho) 
    USE mo_kind, ONLY: wp
    REAL(wp), INTENT(IN) :: tpot
    REAL(wp), INTENT(IN) :: sal
    REAL(wp), INTENT(IN) :: press
    REAL(wp) :: rho
 ENDFUNCTION calc_density
END INTERFACE

!   ! Local variables
   INTEGER :: jc, jb, je,jk, i_no_trac
  !INTEGER :: ile1, ibe1,ile2, ibe2,ile3, ibe3
   INTEGER :: ilc1,ibc1,ilc2,ibc2,jj, ible,idxe
   INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
   INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start_e, rl_end_e
   REAL(wp) :: z_vert_density_grad_c(nproma,n_zlev,p_patch%nblks_c)!(1:n_zlev)
   REAL(wp) :: z_vert_density_grad_e(nproma,n_zlev,p_patch%nblks_e)!(1:n_zlev)
   REAL(wp) :: z_stabio(nproma,n_zlev,p_patch%nblks_c)!(1:n_zlev)
   REAL(wp) :: z_shear_e
   REAL(wp) :: z_shear_c(nproma,n_zlev,p_patch%nblks_c)
   REAL(wp) :: z_shear2_c(nproma,n_zlev,p_patch%nblks_c)
   REAL(wp) :: z_rho_up(nproma,n_zlev,p_patch%nblks_c)
   REAL(wp) :: z_rho_down(nproma,n_zlev,p_patch%nblks_c)
   REAL(wp) :: z_Ri_c(nproma,n_zlev,p_patch%nblks_c)
   REAL(wp) :: z_Ri_e(nproma,n_zlev,p_patch%nblks_e)
   REAL(wp) :: dz_inv
   REAL(wp) :: z_rho_up_c1, z_rho_down_c1,z_rho_up_c2, z_rho_down_c2
   REAL(wp) :: z_lambda_frac 
   REAL(wp) :: z_A_veloc_v_old, z_A_tracer_v_old
   INTEGER ::  z_dolic
   INTEGER :: idx_c1,ibk_c1, idx_c2,ibk_c2, idx_c3,ibk_c3
   INTEGER :: idx_e1,ibk_e1, idx_e2,ibk_e2, idx_e3,ibk_e3,idx_e4,ibk_e4

   !Below is a set of variables and parameters for tracer and velocity
   REAL(wp), PARAMETER :: z_beta          =0.6_wp
   REAL(wp), PARAMETER :: z_one_minus_beta=0.4_wp
   REAL(wp) :: z_A_W_T(nproma,n_zlev,p_patch%nblks_c)
   REAL(wp) :: z_A_W_v(nproma,n_zlev,p_patch%nblks_e)
   REAL(wp) :: z_10m_wind_c(nproma,1,p_patch%nblks_c)
   REAL(wp) :: z_10m_wind_e(nproma,1,p_patch%nblks_e)
   REAL(wp) :: z_grav_rho, z_inv_rho_ref!, z_stabio
   REAL(wp) :: z_press!, z_frac
   REAL(wp) :: A_v_tmp, A_T_tmp
   REAL(wp) :: z_w_T
   REAL(wp) :: z_w_v
   REAL(wp) :: z_s1, z_s2
   REAL(wp) :: z_c(nproma,n_zlev+1,p_patch%nblks_c)
   REAL(wp), PARAMETER :: z_lambda = 0.05_wp
   REAL(wp), PARAMETER :: z_0      = 40.0_wp
   REAL(wp), PARAMETER :: z_c1_T   = 5.0_wp
   REAL(wp), PARAMETER :: z_c1_v   = 5.0_wp
   REAL(wp), PARAMETER :: z_av0    = 0.5E-2_wp 
   REAL(wp), PARAMETER :: z_dv0    = 0.5E-2_wp
   REAL(wp), PARAMETER :: z_treshold= 5.0E-8_wp
   LOGICAL,  PARAMETER :: l_constant_mixing = .FALSE.
!   !-------------------------------------------------------------------------
! DO, jk=1,n_zlev
! CALL &
! & print_mxmn('(uhp) params_oce%K_veloc_h',jk,params_oce%K_veloc_h,n_zlev,p_patch%nblks_e,'phy',3)
! CALL &
! & print_mxmn('(uhp) params_oce%A_veloc_v',jk,params_oce%A_veloc_v,n_zlev,p_patch%nblks_e,'phy',3)
! ENDDO
!write(0,*)'K_veloc_h_back:',params_oce%K_veloc_h_back
!write(0,*)'A_veloc_v_back:',params_oce%A_veloc_v_back
!DO jk=1,no_tracer
!write(0,*)'K_tracer_h_back(',jk,'):',params_oce%K_tracer_h_back(jk)
!write(0,*)'A_tracer_v_back(',jk,'):',params_oce%A_tracer_v_back(jk)
!ENDDO
  IF(l_constant_mixing)THEN

    !nothing to do!In sbr init_ho_params (see above)
    !tracer mixing coefficient params_oce%A_tracer_v(:,:,:, i_no_trac) is already
    !initialzed with params_oce%A_tracer_v_back(i_no_trac)
    !and velocity diffusion coefficient
    ! params_oce%A_veloc_v(je,jk,jb) is initialzed with params_oce%A_veloc_v_back

  ELSEIF(.NOT.l_constant_mixing)THEN
    rl_start_c   = 1
    rl_end_c     = min_rlcell
    i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
    i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)

    rl_start_e   = 1
    rl_end_e     = min_rledge
    i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
    i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)

    z_A_W_T (:,:,:)              = 0.0_wp
    z_A_W_v (:,:,:)              = 0.0_wp
    z_10m_wind_e(:,:,:)          = 0.0_wp
    z_10m_wind_c(:,:,:)          = 0.0_wp
    z_vert_density_grad_c(:,:,:) = 0.0_wp
    z_vert_density_grad_e(:,:,:) = 0.0_wp
    z_Ri_c(:,:,:)                = 0.0_wp
    z_Ri_e(:,:,:)                = 0.0_wp
    z_stabio (:,:,:)             = 0.0_wp
    z_shear_c(:,:,:)             = 0.0_wp
    z_s1  = 0.0_wp
    z_s2  = 0.0_wp
    z_grav_rho    = grav/rho_ref
    z_inv_rho_ref = 1.0_wp/rho_ref

    !Following MPI-OM (cf. vertical mixing sbr)
    z_w_T = CWT/6.0_wp**3
    z_w_v = CWA/6.0_wp**3

    !The wind part
    DO jb = i_startblk_e, i_endblk_e
      CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
      &                   rl_start_e, rl_end_e)
      DO je = i_startidx_e, i_endidx_e 
        IF ( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN

         ilc1 = p_patch%edges%cell_idx(je,jb,1)
         ibc1 = p_patch%edges%cell_blk(je,jb,1)
         ilc2 = p_patch%edges%cell_idx(je,jb,2)
         ibc2 = p_patch%edges%cell_blk(je,jb,2)

          !This is (15) in Marsland et al. 
          z_10m_wind_e(je,1,jb)= SQRT(&
          &0.5_wp*(DOT_PRODUCT(p_sfc_flx%forc_wind_cc(ilc1,ibc1)%x,     &
          &                    p_sfc_flx%forc_wind_cc(ilc1,ibc1)%x)     &
          &       +DOT_PRODUCT(p_sfc_flx%forc_wind_cc(ilc2,ibc2)%x,     &
          &                    p_sfc_flx%forc_wind_cc(ilc2,ibc2)%x)))**3
          z_A_W_v (je,1,jb) = z_w_v*z_10m_wind_e(je,1,jb)
        ENDIF
      END DO
    END DO

    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                   rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c 
        IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
          !This is (15) in Marsland et al. 
          z_10m_wind_c(jc,1,jb)= SQRT(DOT_PRODUCT(p_sfc_flx%forc_wind_cc(jc,jb)%x,&
                                                 &p_sfc_flx%forc_wind_cc(jc,jb)%x))**3
          z_A_W_T (jc,1,jb) = z_w_T*z_10m_wind_c(jc,1,jb)
        ENDIF
      END DO
    END DO

    !Calculate Richardson number and vertical density gradient
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                   rl_start_c, rl_end_c)

      DO jc = i_startidx_c, i_endidx_c

        z_dolic = v_base%dolic_c(jc,jb)
        IF ( z_dolic>=MIN_DOLIC ) THEN        
          DO jk = 2, z_dolic 
           dz_inv = 1.0_wp/v_base%del_zlev_i(jk)

            !This calculates the localshear at cells
            ! - add small epsilon to avoid division by zero
            z_shear_c(jc,jk,jb) = dbl_eps +&!+ dz_inv*dz_inv*                              &
!             & DOT_PRODUCT(p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x,&
!             &              p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)
            & sum((p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)**2)


!------------------------------
!Just a test
!                 idx_c1=p_patch%cells%edge_idx(jc,jb,1)
!                 ibk_c1=p_patch%cells%edge_blk(jc,jb,1)
!                 idx_c2=p_patch%cells%edge_idx(jc,jb,2)
!                 ibk_c2=p_patch%cells%edge_blk(jc,jb,2)
!                 idx_c3=p_patch%cells%edge_idx(jc,jb,3)
!                 ibk_c3=p_patch%cells%edge_blk(jc,jb,3)
! z_shear2_c(jc,jk,jb) = &
! & dz_inv*(&
! & p_os%p_prog(nold(1))%vn(idx_c1,jk-1,ibk_c1)-p_os%p_prog(nold(1))%vn(idx_c1,jk,ibk_c1)&
! &+p_os%p_prog(nold(1))%vn(idx_c2,jk-1,ibk_c2)-p_os%p_prog(nold(1))%vn(idx_c2,jk,ibk_c2)&
! &+p_os%p_prog(nold(1))%vn(idx_c3,jk-1,ibk_c3)-p_os%p_prog(nold(1))%vn(idx_c3,jk,ibk_c3))/3.0_wp
! 
! z_shear2_c(jc,jk,jb) = dbl_eps+z_shear2_c(jc,jk,jb)*z_shear2_c(jc,jk,jb)
!------------------------------
            z_press = v_base%zlev_i(jk)*rho_ref*SItodBar !*grav!z_press = v_base%zlev_i(jk)*rho_ref*grav

            !salinity at upper and lower cell
            IF(no_tracer >= 2) THEN
              z_s1 = p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,2)
              z_s2 = p_os%p_prog(nold(1))%tracer(jc,jk,jb,2)
            ENDIF
            !density of upper and lower cell w.r.t.to pressure at intermediate level
            z_rho_up(jc,jk,jb) = calc_density &
             & (p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,1), z_s1, z_press)

            z_rho_down(jc,jk,jb) = calc_density &
              & (p_os%p_prog(nold(1))%tracer(jc,jk,jb,1), z_s2, z_press)

            ! comments from MPIOM
            !! calculate vertical density stabio gradient between upper and lower box
            !! vertical density gradient 1/delta_z * (rho(k)-rho(k-1))

            ! #slo# 2011-09-02 correction
            ! rho_up: rho(k-1); rho_down: rho(k)
            !  i.e.: dz_inv*(rho_down-rho_up)
            z_stabio(jc,jk,jb)  =(z_rho_down(jc,jk,jb)-z_rho_up(jc,jk,jb))! *dz_inv
            ! z_stabio  = dz_inv*0.5_wp*(z_rho_up(jc,jk,jb)-z_rho_down(jc,jk,jb))
            ! #slo# - think once more about 0.5, and this line in mo_convection of MPIOM:
            ! rhoo(:, j, k-1) = 0.5_wp * (rhoo(:, j, k-1) + rhuppo(:))

            !! stabio > 0 stable stratification    (lower layer is havier) => vert_density_grad  > 0
            !! stabio < 0 instable stratification  (lower layer is lighter)=>vert_density_grad  < 0

            !! set negative values to zero for switch below
            !z_vert_density_grad_c(jc,jk,jb) = MAX(z_stabio(jc,jk,jb), 0.0_wp)

            z_vert_density_grad_c(jc,jk,jb) = dbl_eps+z_stabio(jc,jk,jb)

            ! Richardson number is positive for stable stratification (rho_down>rho_up)
            ! Richardson number is zero for unstable strat., see switch z_frac below
            ! The expression z_grav_rho*z_vert_density_grad_c/z_shear_c is the
            ! Buoyancy frequency
!             z_Ri_c(jc,jk,jb)=MAX(z_grav_rho*z_vert_density_grad_c(jc,jk,jb)/z_shear_c(jc,jk,jb)&
!                                 &,0.0_wp)
            z_Ri_c(jc,jk,jb)=v_base%zlev_i(jk)*z_grav_rho*z_vert_density_grad_c(jc,jk,jb)&
                            &/z_shear_c(jc,jk,jb)
          END DO
        ENDIF
      END DO
    END DO    
!Do jk=1,n_zlev
!write(*,*)'max-min Ri-Nr:densgrad:shear:',&
!&maxval( z_Ri_c(:,jk,:)),&!minval( z_Ri_c(:,jk,:)),&
!&maxval(z_vert_density_grad_c(:,jk,:)),&!minval(z_vert_density_grad_c(:,jk,:)),&
!&maxval( z_shear_c(:,jk,:))!,minval( z_shear_c(:,jk,:))
!END DO

   !Density gradient for 1st level not calulated yet, but required below. Following
   !parxis in MPI-OM we set first layer equal to second layer (cf MPI-OM mo_ocean_vertical_mixing)
   !z_vert_density_grad_c(:,1,:) =  z_vert_density_grad_c(:,2,:)

   !The tracer mixing coefficient at cell centers
    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
      &                   rl_start_c, rl_end_c)

      DO jc = i_startidx_c, i_endidx_c

        z_dolic = v_base%dolic_c(jc,jb)
        IF ( z_dolic>=MIN_DOLIC ) THEN        
          DO jk = 2, z_dolic
 
           dz_inv = 1.0_wp/v_base%del_zlev_i(jk)           
  
            !calculate vertical tracer mixing based on local Richardson number
            DO i_no_trac=1, no_tracer

              !Store old diffusivity
              z_A_tracer_v_old = params_oce%A_tracer_v(jc,jk,jb,i_no_trac)

              !! vert_density_grad  = 0 'semi-stable': use background value
              IF (     z_vert_density_grad_c(jc,jk,jb) > -z_treshold&
                &.AND. z_vert_density_grad_c(jc,jk,jb) < z_treshold) THEN
                params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = params_oce%A_tracer_v_back(i_no_trac)
                DO jj=1,3
                  idxe = p_patch%cells%edge_idx(jc,jb,jj)
                  ible = p_patch%cells%edge_blk(jc,jb,jj)
                  params_oce%A_veloc_v(idxe,jk,ible) = params_oce%A_veloc_v_back
                ENDDO 
              !! vert_density_grad  < 0 instable stratification: use convective mixing parameter
              ELSE IF (z_vert_density_grad_c(jc,jk,jb) < -z_treshold ) THEN
                params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = MAX_VERT_DIFF_TRAC
!  write(246,*)'instab',&
!  & jc,jk,jb,z_vert_density_grad_c(jc,jk,jb),&
!  & z_rho_down(jc,jk,jb),z_rho_up(jc,jk,jb),&
!  & z_rho_down(jc,jk,jb)-z_rho_up(jc,jk,jb),&
!  &p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,1),&
!  &p_os%p_prog(nold(1))%tracer(jc,jk,jb,1),&
! &z_Ri_c(jc,jk,jb)
                DO jj=1,3
                  idxe = p_patch%cells%edge_idx(jc,jb,jj)
                  ible = p_patch%cells%edge_blk(jc,jb,jj)
                  params_oce%A_veloc_v(idxe,jk,ible) = MAX_VERT_DIFF_VELOC
                ENDDO 
!                 idx_c1=p_patch%cells%neighbor_idx(jc,jb,1)
!                 ibk_c1=p_patch%cells%neighbor_blk(jc,jb,1)
!                 idx_c2=p_patch%cells%neighbor_idx(jc,jb,2)
!                 ibk_c2=p_patch%cells%neighbor_blk(jc,jb,2)
!                 idx_c3=p_patch%cells%neighbor_idx(jc,jb,3)
!                 ibk_c3=p_patch%cells%neighbor_blk(jc,jb,3)
! 
!                 params_oce%A_tracer_v(idx_c1,jk,ibk_c1,i_no_trac)=MAX_VERT_DIFF_TRAC
!                 params_oce%A_tracer_v(idx_c2,jk,ibk_c2,i_no_trac)=MAX_VERT_DIFF_TRAC
!                 params_oce%A_tracer_v(idx_c3,jk,ibk_c3,i_no_trac)=MAX_VERT_DIFF_TRAC

              !! vert_density_grad  > 0 stable stratification: use calculated value
              ELSE IF (z_vert_density_grad_c(jc,jk,jb) > z_treshold ) THEN

                !This is (16) in Marsland et al. and identical to treatment of velocity
                !but it allows to use different parameters
!                 z_lambda_frac     = z_lambda/v_base%del_zlev_i(jk)
!                 z_A_W_T(jc,jk,jb) = z_A_W_T (jc,jk-1,jb)*z_lambda_frac      &
!                   &*(exp(-v_base%del_zlev_i(jk)/z_0))&
!                   &/(z_lambda_frac+0.5_wp*(z_vert_density_grad_c(jc,jk,jb)&
!                                        &+z_vert_density_grad_c(jc,jk-1,jb)))
!                 ! This is (19) in Marsland et al. valid for stable stratification, with
!                 !   with: z_c1_T=CRD=5.0, z_av0=DVO=0.005, A_tracer_v=Db=1.0e-5, Dw=0.0
!                 A_T_tmp = &
!                 &z_one_minus_beta&
!                 &*MIN(z_A_tracer_v_old, z_dv0+params_oce%A_tracer_v_back(i_no_trac))&
!                 & +z_beta*(z_A_W_T(jc,jk,jb)+z_dv0/((1.0_wp+z_c1_T*z_Ri_c(jc,jk,jb))**3)&
!                 &                  +params_oce%A_tracer_v_back(i_no_trac))
!                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = MIN(MAX_VERT_DIFF_TRAC, A_T_tmp)

                z_Ri_c(jc,jk,jb) = max(z_Ri_c(jc,jk,jb),0.0_wp)

                A_T_tmp = params_oce%A_tracer_v_back(i_no_trac)&
                &+ z_dv0/((1.0_wp+z_c1_T*z_Ri_c(jc,jk,jb))**3)
                
                params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = A_T_tmp

              
                DO jj=1,3
                  idxe = p_patch%cells%edge_idx(jc,jb,jj)
                  ible = p_patch%cells%edge_blk(jc,jb,jj)

                  params_oce%A_veloc_v(idxe,jk,ible) = params_oce%A_veloc_v_back&
                                    &+ z_av0/((1.0_wp+z_c1_v*z_Ri_c(jc,jk,jb))**2)
               END DO

              END IF
            ENDDO
          END DO
        ENDIF
      END DO
    END DO
    !END DO
!    DO jb = i_startblk_c, i_endblk_c
!       CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
!       &                   rl_start_c, rl_end_c)
!         DO jc = i_startidx_c, i_endidx_c
!  write(*,*)'Vert Trac diff:',&!A_T_tmp,&
!    &params_oce%A_tracer_v(jc,:,jb, i_no_trac)
!        END DO
!     END DO 


!     !Viscosity at edges
!     !This calculates vertical density gradient and local Richardson number at edges
!     !as average of neighbouring cell values
!     !This calulation is done separately to allow loop-start from 1st layer: the  loop
!     !for calculation of vertical diffusivity below required edge-density gradient at first layer.  
!     DO jb = i_startblk_e, i_endblk_e
!       CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
!       &                   rl_start_e, rl_end_e)
!       DO jk = 1, n_zlev
!         DO je = i_startidx_e, i_endidx_e
!           z_dolic = v_base%dolic_e(je,jb)
! 
!           IF ( z_dolic >=MIN_DOLIC ) THEN
!             !This calculates the local Richardson number at edges
!             !and  shear at edges as average of shear at centers
!             !indices of neighboring edges 
!             ilc1 = p_patch%edges%cell_idx(je,jb,1)
!             ibc1 = p_patch%edges%cell_blk(je,jb,1)
!             ilc2 = p_patch%edges%cell_idx(je,jb,2)
!             ibc2 = p_patch%edges%cell_blk(je,jb,2)
! 
!             z_vert_density_grad_e(je,jk,jb) = 0.5_wp*&
!             &( z_vert_density_grad_c(ilc1,jk,ibc1)&
!             & +z_vert_density_grad_c(ilc2,jk,ibc2))
! 
!             z_Ri_e(je,jk,jb) = 0.5_wp*(z_Ri_c(ilc1,jk,ibc1)+z_Ri_c(ilc2,jk,ibc2))
!           ENDIF
!         END DO
!       ENDDO
!     END DO
!     !calculate vertical diffusivity
!     DO jb = i_startblk_e, i_endblk_e
!       CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
!       &                   rl_start_e, rl_end_e)
!       DO jk = 2, n_zlev
!         DO je = i_startidx_e, i_endidx_e
!           z_dolic = v_base%dolic_e(je,jb)
! 
!           IF ( z_dolic >=MIN_DOLIC ) THEN
!             !dz_inv  = 1.0_wp/v_base%del_zlev_i(jk)
! 
!             !Store old viscosity
!             z_A_veloc_v_old = params_oce%A_veloc_v(je,jk,jb)
! 
!             !! vert_density_grad  = 0 'semi-stable': use background value
!             IF (z_vert_density_grad_e(je,jk,jb) == 0.0_wp ) THEN
!               params_oce%A_veloc_v(je,jk,jb) = params_oce%A_veloc_v_back
!             !! vert_density_grad  < 0 instable stratification: use convective mixing parameter
!             ELSE IF (z_vert_density_grad_e(je,jk,jb) < 0.0_wp ) THEN
!               params_oce%A_veloc_v(je,jk,jb) = MAX_VERT_DIFF_VELOC
! !                 !Possible averaging of mixing coefficient
! !                 idx_e1=p_patch%edges%quad_idx(je,jb,1)
! !                 ibk_e1=p_patch%edges%quad_blk(je,jb,1)
! !                 idx_e2=p_patch%edges%quad_idx(je,jb,2)
! !                 ibk_e2=p_patch%edges%quad_blk(je,jb,2)
! !                 idx_e3=p_patch%edges%quad_idx(je,jb,3)
! !                 ibk_e3=p_patch%edges%quad_blk(je,jb,3)
! !                 idx_e4=p_patch%edges%quad_idx(je,jb,4)
! !                 ibk_e4=p_patch%edges%quad_blk(je,jb,4)
! !                 params_oce%A_veloc_v(idx_e1,jk,ibk_e1) = MAX_VERT_DIFF_VELOC
! !                 params_oce%A_veloc_v(idx_e2,jk,ibk_e3) = MAX_VERT_DIFF_VELOC
! !                 params_oce%A_veloc_v(idx_e3,jk,ibk_e3) = MAX_VERT_DIFF_VELOC
! !                 params_oce%A_veloc_v(idx_e4,jk,ibk_e4) = MAX_VERT_DIFF_VELOC
!             !! vert_density_grad  > 0 stable stratification: use calculated value
!             ELSE IF (z_vert_density_grad_e(je,jk,jb) > 0.0_wp ) THEN
! 
!               !This is (16) in Marsland et al.
! !               z_lambda_frac = z_lambda/v_base%del_zlev_i(jk)
! !               z_A_W_v (je,jk,jb) = z_A_W_v (je,jk-1,jb)*z_lambda_frac     &
! !               &*(z_lambda_frac*exp(-v_base%del_zlev_i(jk)/z_0))&
! !               &/(z_lambda_frac+0.5_wp*(z_vert_density_grad_e(je,jk,jb)&
! !                                      &+z_vert_density_grad_e(je,jk-1,jb)))
! !               A_v_tmp = z_one_minus_beta*MIN(z_A_veloc_v_old, z_av0+params_oce%A_veloc_v_back)&
! !                       & +z_beta*(z_A_W_v (je,jk,jb)+z_av0/((1.0_wp+z_c1_v*z_Ri_e(je,jk,jb))**2)         &
! !                       &         +params_oce%A_veloc_v_back)
! 
!               z_Ri_e(je,jk,jb) = max(z_Ri_e(je,jk,jb),0.0_wp)
!               params_oce%A_veloc_v(je,jk,jb) = params_oce%A_veloc_v_back&
!                                             &+ z_av0/((1.0_wp+z_c1_v*z_Ri_e(je,jk,jb))**2)
! 
!             END IF
! 
!           !ENDIF
!           ENDIF
!         END DO
!       END DO
!     END DO

    ! set to background value
    !params_oce%A_veloc_v(:,1,:) = params_oce%A_veloc_v_back!params_oce%A_veloc_v(:,2,:)
  ENDIF!l_constant_mixing


DO i_no_trac=1, no_tracer
 z_c(:,:,:)=params_oce%A_tracer_v(:,:,:,i_no_trac)
 DO jk=1,n_zlev
  ipl_src=3  ! output print level (1-5, fix)
  CALL print_mxmn('PHY trac mixing',jk,z_c(:,:,:),n_zlev+1,p_patch%nblks_c,'phy',ipl_src)
  CALL print_mxmn('z_A_W_v',jk,z_A_W_v,n_zlev,p_patch%nblks_e,'phy',ipl_src)
  CALL print_mxmn('z_A_W_T',jk,z_A_W_T,n_zlev,p_patch%nblks_c,'phy',ipl_src)
  CALL print_mxmn('p_vn%x(1)',jk,p_os%p_diag%p_vn%x(1),n_zlev,p_patch%nblks_c,'phy',ipl_src)
  CALL print_mxmn('p_vn%x(2)',jk,p_os%p_diag%p_vn%x(2),n_zlev,p_patch%nblks_c,'phy',ipl_src)
  CALL print_mxmn('z_shear_c',jk,z_shear_c,n_zlev,p_patch%nblks_c,'phy',ipl_src)
  !write(*,*)'max/min trac mixing',jk,maxval(params_oce%A_tracer_v(:,jk,:,i_no_trac)),&
  !&minval(params_oce%A_tracer_v(:,jk,:,i_no_trac))
  !write(*,*)'max/min veloc mixing',jk,maxval(params_oce%A_veloc_v(:,jk,:)),&
  !minval(params_oce%A_veloc_v(:,jk,:))

  !write(123,*)'max/min trac mixing',jk,maxval(params_oce%A_tracer_v(:,jk,:,i_no_trac)),&
  !&minval(params_oce%A_tracer_v(:,jk,:,i_no_trac))
 END DO
END DO
 DO jk=1,n_zlev
   ipl_src=3  ! output print level (1-5, fix)
   CALL print_mxmn('PHY veloc mixing',jk,params_oce%A_veloc_v(:,:,:),n_zlev+1, &
     & p_patch%nblks_e,'phy',ipl_src)
! write(*,*)'max/min veloc mixing',jk,maxval(params_oce%A_veloc_v(:,jk,:)),&
! &minval(params_oce%A_veloc_v(:,jk,:))
 END DO

 END SUBROUTINE update_ho_params
! !   !-------------------------------------------------------------------------
! !   !
! !   !>
! !   !! Update of parameters
! !   !!
! !   !! Update of ocean physics: This routine is used used only if time-dependent
! !   !! changes of physical parametrizations.
! !   !! Currently vertical mixing coefficients for tracers and vertical diffusivity are updated.
! !   !! Dependent on the local Richardson number the diffusivity are calculated
! !   !!(Large & Gent JPO 29, (1999), 449-464).
! !   !!The formulation follows the MPI-OM implementation as described in Marsland et al. (Ocean Modelling, 2002).
! !   !!The notationnal convection is also taken from this paper( cf. eqs (14)-(19)).
! !   !! What is missing is the fractional ice cover (see eqs. (15-16)).
! !   !!
! !   !! @par Revision History
! !   !! Initial release by Peter Korn, MPI-M (2011-02)
! !   !
! !   !
! !  SUBROUTINE update_ho_params_old(p_patch, p_os, p_sfc_flx, params_oce, calc_density)
! !    TYPE(t_patch),     INTENT(IN)     :: p_patch
! !    TYPE(t_hydro_ocean_state), TARGET :: p_os
! !    TYPE(t_sfc_flx),   INTENT(INOUT)  :: p_sfc_flx
! !    TYPE(t_ho_params), INTENT(INOUT)  :: params_oce  
! ! INTERFACE !This contains the function version of the actual EOS as chosen in namelist
! !   FUNCTION calc_density(tpot, sal, press) RESULT(rho) 
! !     USE mo_kind, ONLY: wp
! !     REAL(wp), INTENT(IN) :: tpot
! !     REAL(wp), INTENT(IN) :: sal
! !     REAL(wp), INTENT(IN) :: press
! !     REAL(wp) :: rho
! !  ENDFUNCTION calc_density
! ! END INTERFACE
! ! 
! ! !   ! Local variables
! !    INTEGER :: jc, jb, je,jk, i_no_trac
! !   !INTEGER :: ile1, ibe1,ile2, ibe2,ile3, ibe3
! !    INTEGER :: ilc1,ibc1,ilc2,ibc2
! !    INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
! !    INTEGER  :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, rl_start_e, rl_end_e
! !    REAL(wp) :: z_vert_density_grad_c(1:n_zlev)
! !    REAL(wp) :: z_vert_density_grad_e(1:n_zlev)
! !    REAL(wp) :: z_shear_e
! !    REAL(wp) :: z_shear_c(nproma,n_zlev,p_patch%nblks_c)
! !    REAL(wp) :: z_rho_up(nproma,n_zlev,p_patch%nblks_c)
! !    REAL(wp) :: z_rho_down(nproma,n_zlev,p_patch%nblks_c)
! !    REAL(wp) :: dz_inv
! !    REAL(wp) :: z_Ri_e, z_Ri_c
! !    REAL(wp) :: z_rho_up_c1, z_rho_down_c1,z_rho_up_c2, z_rho_down_c2
! !    REAL(wp) :: z_lambda_frac 
! !    REAL(wp) :: z_A_veloc_v_old, z_A_tracer_v_old
! !    INTEGER :: z_dolic
! ! 
! !    !Below is a set of variables and parameters for tracer and velocity
! !    REAL(wp), PARAMETER :: z_beta          =0.6_wp
! !    REAL(wp), PARAMETER :: z_one_minus_beta=0.4_wp
! !    REAL(wp) :: z_A_W_T(nproma,n_zlev,p_patch%nblks_c)
! !    REAL(wp) :: z_A_W_v(nproma,n_zlev,p_patch%nblks_e)
! !    REAL(wp) :: z_10m_wind_c(nproma,1,p_patch%nblks_c)
! !    REAL(wp) :: z_10m_wind_e(nproma,1,p_patch%nblks_e)
! !    REAL(wp) :: z_grav_rho, z_inv_rho_ref, z_stabio
! !    REAL(wp) :: z_press, z_frac
! !    REAL(wp) :: A_v_tmp, A_T_tmp
! !    REAL(wp) :: z_w_T
! !    REAL(wp) :: z_w_v
! !    REAL(wp) :: z_s1
! !    REAL(wp) :: z_c(nproma,n_zlev+1,p_patch%nblks_c)
! !    REAL(wp),PARAMETER  :: z_lambda = 0.05_wp
! !    REAL(wp), PARAMETER :: z_0      = 40.0_wp
! !    REAL(wp), PARAMETER :: z_c1_T     = 5.0_wp
! !    REAL(wp), PARAMETER :: z_c1_v     = 5.0_wp
! !    REAL(wp), PARAMETER :: z_av0      = 0.5E-2_wp
! ! !   !-------------------------------------------------------------------------
! ! !    IF (no_tracer == 0) RETURN
! !     rl_start_c   = 1
! !     rl_end_c     = min_rlcell
! !     i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
! !     i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
! ! 
! !     rl_start_e   = 1
! !     rl_end_e     = min_rledge
! !     i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! !     i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)
! ! 
! !     z_A_W_T (:,:,:)         = 0.0_wp
! !     z_A_W_v (:,:,:)         = 0.0_wp
! !     z_10m_wind_e(:,:,:)     = 0.0_wp
! !     z_10m_wind_c(:,:,:)     = 0.0_wp
! !     z_vert_density_grad_c(:)= 0.0_wp
! !     z_vert_density_grad_e(:)= 0.0_wp
! !     z_shear_c(:,:,:)        = 0.0_wp
! !     z_s1                    = 0.0_wp
! ! 
! !     z_grav_rho    = grav/rho_ref
! !     z_inv_rho_ref = 1.0_wp/rho_ref
! ! 
! !     !Following MPI-OM (cf. vertical mixing sbr)
! !     z_w_T = CWT/6.0_wp**3
! !     z_w_v = CWA/6.0_wp**3
! ! 
! !     !The wind part
! !     DO jb = i_startblk_e, i_endblk_e
! !       CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
! !       &                   rl_start_e, rl_end_e)
! !       DO je = i_startidx_e, i_endidx_e 
! !         IF ( v_base%lsm_oce_e(je,1,jb) <= sea_boundary ) THEN
! ! 
! !          ilc1 = p_patch%edges%cell_idx(je,jb,1)
! !          ibc1 = p_patch%edges%cell_blk(je,jb,1)
! !          ilc2 = p_patch%edges%cell_idx(je,jb,2)
! !          ibc2 = p_patch%edges%cell_blk(je,jb,2)
! ! 
! !           !This is (15) in Marsland et al. 
! !           z_10m_wind_e(je,1,jb)= SQRT(&
! !           &0.5_wp*(DOT_PRODUCT(p_sfc_flx%forc_wind_cc(ilc1,ibc1)%x,     &
! !           &                    p_sfc_flx%forc_wind_cc(ilc1,ibc1)%x)     &
! !           &       +DOT_PRODUCT(p_sfc_flx%forc_wind_cc(ilc2,ibc2)%x,     &
! !           &                    p_sfc_flx%forc_wind_cc(ilc2,ibc2)%x)))**3
! !           z_A_W_v (je,1,jb) = z_w_v*z_10m_wind_e(je,1,jb)
! !         ENDIF
! !       END DO
! !     END DO
! ! 
! !     DO jb = i_startblk_c, i_endblk_c
! !       CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !       &                   rl_start_c, rl_end_c)
! !       DO jc = i_startidx_c, i_endidx_c 
! !         IF ( v_base%lsm_oce_c(jc,1,jb) <= sea_boundary ) THEN
! !           !This is (15) in Marsland et al. 
! !           z_10m_wind_c(jc,1,jb)= SQRT(DOT_PRODUCT(p_sfc_flx%forc_wind_cc(jc,jb)%x,&
! !                                                  &p_sfc_flx%forc_wind_cc(jc,jb)%x))**3
! !           z_A_W_T (jc,1,jb) = z_w_T*z_10m_wind_c(jc,1,jb)
! !         ENDIF
! !       END DO
! !     END DO
! ! 
! !     !The tracer mixing coefficient at cell centers
! !     DO jb = i_startblk_c, i_endblk_c
! !       CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! !       &                   rl_start_c, rl_end_c)
! !       DO jk = 2, n_zlev
! !  
! !        dz_inv = 1.0_wp/v_base%del_zlev_i(jk)           !(v_base%zlev_m(jk)-v_base%zlev_m(jk-1)) 
! ! 
! !         DO jc = i_startidx_c, i_endidx_c
! !           z_dolic = v_base%dolic_c(jc,jb)
! !           IF ( z_dolic>=MIN_DOLIC ) THEN        
! ! 
! !             !This calculates the local Richardson number at cells
! !             ! - add small epsilon to avoid division by zero
! !             z_shear_c(jc,jk,jb) = dbl_eps + dz_inv*dz_inv                              &
! !             &* DOT_PRODUCT(p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x,&
! !             &              p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)
! ! 
! !             z_press = v_base%zlev_i(jk)*rho_ref*SItodBar !*grav
! !             !z_press = v_base%zlev_i(jk)*rho_ref*grav
! !             !density of upper cell w.r.t.to pressure at intermediate level
! !             ! #slo# 2011-08-10: tracer(:,:,:,2) not defined if no_tracer=1
! !             IF (no_tracer >= 2) z_s1 = p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,2)
! ! 
! !             z_rho_up(jc,jk,jb) = calc_density &
! !              & (p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,1), z_s1, z_press)
! !           ! z_rho_up(jc,jk,jb) = calc_density(&
! !           ! & p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,1),&
! !           ! & p_os%p_prog(nold(1))%tracer(jc,jk-1,jb,2),&
! !           ! & z_press)
! !             !density of lower cell w.r.t.to pressure at intermediate level
! !             IF (no_tracer >= 2) z_s1 = p_os%p_prog(nold(1))%tracer(jc,jk,jb,2)
! ! 
! !             z_rho_down(jc,jk,jb) = calc_density &
! !               & (p_os%p_prog(nold(1))%tracer(jc,jk,jb,1), z_s1, z_press)
! !            ! z_rho_down(jc,jk,jb) = calc_density(&
! !            ! & p_os%p_prog(nold(1))%tracer(jc,jk,jb,1),&
! !            ! & p_os%p_prog(nold(1))%tracer(jc,jk,jb,2),&
! !            ! & z_press)
! !             ! comments from MPIOM
! !             !! calculate vertical density stabio gradient between upper and lower box
! !             !! vertical density gradient 1/delta_z * (rho(k)-rho(k-1))
! !             ! #slo# 2011-09-02 correction
! !             ! rho_up: rho(k-1); rho_down: rho(k)
! !             !  i.e.: dz_inv*(rho_down-rho_up)
! !             z_stabio  = dz_inv*(z_rho_down(jc,jk,jb)-z_rho_up(jc,jk,jb))
! !             ! z_stabio  = dz_inv*0.5_wp*(z_rho_up(jc,jk,jb)-z_rho_down(jc,jk,jb))
! !             ! #slo# - think once more about 0.5, and this line in mo_convection of MPIOM:
! !             ! rhoo(:, j, k-1) = 0.5_wp * (rhoo(:, j, k-1) + rhuppo(:))
! !             !! stabio > 0 stable stratification    (lower layer is havier)
! !             !! stabio < 0 instable stratification  (lower layer is lighter)
! !             !! set negative values to zero for switch below
! ! 
! !             z_vert_density_grad_c(jk) = MAX(z_stabio, 0.0_wp)
! ! 
! !             ! #slo# 2011-09-02 correction
! !             ! Richardson number is positive for stable stratification (rho_down>rho_up)
! !             ! Richardson number is zero for unstable strat., see switch z_frac below
! !             z_Ri_c = MAX(z_grav_rho*z_vert_density_grad_c(jk)/z_shear_c(jc,jk,jb),0.0_wp)
! !        !    write(*,*)'rup/do,stabio,drho,Ri,shear',jk,jc,z_rho_up(jc,jk,jb),&
! !        !    &z_rho_down(jc,jk,jb), z_stabio,z_vert_density_grad_c(jk), z_Ri_c, z_shear_c(jc,jk,jb)
! ! 
! !             !calculate vertical tracer mixing based on local Richardson number
! !             DO i_no_trac=1, no_tracer
! ! 
! !               !Store old diffusivity
! !               z_A_tracer_v_old = params_oce%A_tracer_v(jc,jk,jb,i_no_trac)
! !  
! !               !This is (16) in Marsland et al. and identical to treatment of velocity
! !               !but it allows to use different parameters
! !               z_lambda_frac     = z_lambda/v_base%del_zlev_i(jk)
! ! 
! !               z_A_W_T(jc,jk,jb) = z_A_W_T (jc,jk-1,jb)*z_lambda_frac      &
! !                 &*(exp(-v_base%del_zlev_i(jk)/z_0))&
! !                 &/(z_lambda_frac+0.5_wp*(z_vert_density_grad_c(jk)+z_vert_density_grad_c(jk-1)))
! ! 
! ! 
! !               ! This is (19) in Marsland et al. valid for stable stratification, with
! !               !   with: z_c1_T=CRD=5.0, z_av0=DVO=0.005, A_tracer_v=Db=1.0e-5, Dw=0.0
! !               A_T_tmp = &
! !               &z_one_minus_beta*MIN(z_A_tracer_v_old, z_av0+params_oce%A_tracer_v_back(i_no_trac))&
! !               & +z_beta*(z_A_W_T(jc,jk,jb)+z_av0/((1.0_wp+z_c1_T*z_Ri_c)**3)                      &
! !               &         +params_oce%A_tracer_v_back(i_no_trac))
! !   ! write(*,*)'AWT,tmp1,tmp2,Ttmp',jk,jc,&
! !   ! &z_one_minus_beta*MIN(z_A_tracer_v_old, z_av0+params_oce%A_tracer_v_back(i_no_trac)),&
! !   ! & +z_beta*(z_A_W_T(jc,jk,jb)+z_av0/((1.0_wp+z_c1_T*z_Ri_c)**3)                      &
! !   ! & +params_oce%A_tracer_v_back(i_no_trac)),A_T_tmp
! ! 
! !               ! #slo# In the following code z_frac is used as a switch
! !               ! vert_density_grad  > 0 stable stratification   -> z_frac=-1.0
! !               ! vert_density_grad  = 0 unstable stratification -> z_frac=+1.0
! !               z_frac=(1.0E-11_wp-z_vert_density_grad_c(jk))&
! !               &/(1.0E-11_wp+ABS(z_vert_density_grad_c(jk)))
! ! 
! !  
! !               params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = MAX(MAX_VERT_DIFF_TRAC*z_frac, A_T_tmp)!params_oce%A_tracer_v_back(i_no_trac)!
! ! 
! !               ! z_frac is used to avoid an if-condition
! !               !  IF (z_vert_density_grad_c(jk) == 0.0_wp ) THEN
! !               !    params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = MAX(MAX_VERT_DIFF_TRAC, A_T_tmp)
! !               !  ELSE IF (z_vert_density_grad_c(jk) > 0.0_wp ) THEN
! !               !    params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = A_T_tmp
! !               !  ELSE
! !               !    CALL FINISH(...)
! !               !  END IF
! !     ! write(*,*)'Ri number',jk,jc,jb,z_Ri_c, z_vert_density_grad_c(jk),z_rho_up(jc,jk,jb),&
! !   ! &z_rho_down(jc,jk,jb),&
! !   ! &z_press,z_frac,&
! !   ! & params_oce%A_tracer_v(jc,jk,jb, i_no_trac)
! ! !               !This is (16) in Marsland et al. and identical to treatment of velocity
! ! !               !but it allows to use different parameters
! ! !               z_lambda_frac     = z_lambda_T/v_base%del_zlev_i(jk)
! ! !               z_A_W_T(jc,jk,jb) = z_A_W_T (jc,jk-1,jb)                      &
! ! !               &*(z_lambda_frac*exp(-v_base%del_zlev_i(jk)/z_0_T))&
! ! !               &/(z_lambda_frac+z_vert_density_grad_c)
! !               !For positive Richardson number set vertical mixing coefficient to maximal number 
! ! !               IF(z_Ri_c <= 0.0_wp)THEN
! ! !                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac)                        &
! ! !                 & = params_oce%A_tracer_v_back(i_no_trac)!z_one_minus_beta* z_A_tracer_v_old                            &
! ! !                 !& + z_beta*(params_oce%A_tracer_v_back(i_no_trac)                 &
! ! !                 !& + params_oce%A_tracer_v_back(i_no_trac)/(1.0_wp+z_c1_T*z_Ri_c)**3 &
! ! !                 !& + z_A_W_T(jc,jk,jb))
! ! !              !write(123,*)'neg T-Ri number',jc,jk,jb,params_oce%A_tracer_v(jc,jk,jb, i_no_trac)
! ! !               ELSEIF(z_Ri_c > 0.0_wp)THEN
! ! ! !                 write(123,*)'pos T-Ri number',jc,jk,jb,z_max_diff_T
! ! !                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = params_oce%A_tracer_v_back(i_no_trac)!z_max_diff_T!params_oce%A_tracer_v_back(i_no_trac)
! ! !               ENDIF
! !             ENDDO
! !           !ENDIF
! !           ENDIF
! !         END DO
! !       END DO
! !     END DO
! !     !END DO
! ! 
! ! !    DO jb = i_startblk_c, i_endblk_c
! ! !       CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
! ! !       &                   rl_start_c, rl_end_c)
! ! !         DO jc = i_startidx_c, i_endidx_c
! ! ! 
! ! !  write(*,*)'Vert Trac diff:',&!A_T_tmp,&
! ! !    &params_oce%A_tracer_v(jc,:,jb, i_no_trac)
! ! !        END DO
! ! !     END DO 
! ! 
! !     DO i_no_trac=1, no_tracer
! !       params_oce%A_tracer_v(:,1,:, i_no_trac) = params_oce%A_tracer_v_back(1) !params_oce%A_tracer_v(:,2,:,i_no_trac)
! !     END DO
! !     !Viscosity at edges
! !     DO jb = i_startblk_e, i_endblk_e
! !       CALL get_indices_e( p_patch, jb, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e, &
! !       &                   rl_start_e, rl_end_e)
! !       DO jk = 2, n_zlev
! !         DO je = i_startidx_e, i_endidx_e
! !           z_dolic = v_base%dolic_e(je,jb)
! !          ! write(*,*)'z_dolic e',z_dolic,v_base%lsm_oce_e(je,jk,jb)
! !           IF ( z_dolic >=MIN_DOLIC ) THEN
! ! 
! !             dz_inv = 1.0_wp/v_base%del_zlev_i(jk)
! ! 
! !             !indices of neighboring edges 
! !             ilc1 = p_patch%edges%cell_idx(je,jb,1)
! !             ibc1 = p_patch%edges%cell_blk(je,jb,1)
! !             ilc2 = p_patch%edges%cell_idx(je,jb,2)
! !             ibc2 = p_patch%edges%cell_blk(je,jb,2)
! ! 
! !             !This calculates the local Richardson number at edges
! !             !Shear at edges as average of shear at centers
! !             z_shear_e = 0.5_wp*(z_shear_c(ilc1,jk,ibc1)+z_shear_c(ilc2,jk,ibc2))
! !  
! !             z_press = v_base%zlev_i(jk)*rho_ref*grav
! ! 
! !             !density of two adjacent upper cell centers
! !             z_rho_up_c1 = z_rho_up(ilc1,jk,ibc1)
! !             z_rho_up_c2 = z_rho_up(ilc2,jk,ibc2)
! !             !density of two adjacent lower cell centers
! !             z_rho_down_c1 = z_rho_down(ilc1,jk,ibc1)
! !             z_rho_down_c2 = z_rho_down(ilc2,jk,ibc2)
! ! 
! !             ! #slo# correction, see above:
! !             ! z_stabio  = dz_inv&
! !             ! &*0.5_wp*(z_rho_up_c1 + z_rho_up_c2-z_rho_down_c1-z_rho_down_c2)
! !             z_stabio  = dz_inv*0.5_wp*(z_rho_down_c1 + z_rho_down_c2-z_rho_up_c1-z_rho_up_c2)
! ! 
! !             z_vert_density_grad_e(jk) = MAX(z_stabio, 0.0_wp)
! ! 
! !             z_Ri_e = MAX(z_grav_rho*z_vert_density_grad_e(jk)/z_shear_e,0.0_wp)
! !             !z_vert_density_grad_e/(z_shear_e+0.0005_wp)! 
! !             !write(*,*)'Ri number',jk,jc,jb,z_Ri, z_vert_density_grad,z_shear
! ! 
! !             !calculate vertical viscosity based on local Richardson number
! !             !Store old viscosity
! !             z_A_veloc_v_old = params_oce%A_veloc_v(je,jk,jb)
! ! 
! !             !This is (16) in Marsland et al.
! !             z_lambda_frac = z_lambda/v_base%del_zlev_i(jk)
! ! 
! !             z_A_W_v (je,jk,jb) = z_A_W_v (je,jk-1,jb)*z_lambda_frac     &
! !             &*(z_lambda_frac*exp(-v_base%del_zlev_i(jk)/z_0))&
! !             &/(z_lambda_frac+0.5_wp*(z_vert_density_grad_e(jk)+z_vert_density_grad_e(jk-1)))
! ! 
! !             A_v_tmp = z_one_minus_beta*MIN(z_A_veloc_v_old, z_av0+params_oce%A_veloc_v_back)&
! !                     & +z_beta*(z_A_W_v (je,jk,jb)+z_av0/((1.0_wp+z_c1_v*z_Ri_e)**2)         &
! !                     &         +params_oce%A_veloc_v_back)
! ! 
! !             z_frac=(1.0E-11_wp-z_vert_density_grad_e(jk))&
! !             &/(1.0E-11_wp+ABS(z_vert_density_grad_e(jk)))
! ! 
! !             params_oce%A_veloc_v(je,jk,jb) = MAX(MAX_VERT_DIFF_VELOC*z_frac, A_v_tmp)
! !             !params_oce%A_veloc_v_back
! ! !  write(*,*)'Ri number',jk,jc,jb,z_Ri_e, z_vert_density_grad_e(jk),&
! ! !  &z_frac,&
! ! !  & params_oce%A_veloc_v(je,jk,jb)
! !           !ENDIF
! !           ENDIF
! !         END DO
! !       END DO
! !     END DO
! !     params_oce%A_veloc_v(:,1,:) = params_oce%A_veloc_v_back!params_oce%A_veloc_v(:,2,:)
! ! 
! ! DO i_no_trac=1, no_tracer
! !  z_c(:,:,:)=params_oce%A_tracer_v(:,:,:,i_no_trac)
! !  DO jk=1,n_zlev
! !   ipl_src=3  ! output print level (1-5, fix)
! !   CALL print_mxmn('PHY trac mixing',jk,z_c(:,:,:),n_zlev+1,p_patch%nblks_c,'phy',ipl_src)
! !   !write(*,*)'max/min trac mixing',jk,maxval(params_oce%A_tracer_v(:,jk,:,i_no_trac)),&
! !   !&minval(params_oce%A_tracer_v(:,jk,:,i_no_trac))
! !   write(123,*)'max/min trac mixing',jk,maxval(params_oce%A_tracer_v(:,jk,:,i_no_trac)),&
! !   &minval(params_oce%A_tracer_v(:,jk,:,i_no_trac))
! !  END DO
! ! END DO
! !  DO jk=1,n_zlev
! !   ipl_src=3  ! output print level (1-5, fix)
! !   CALL print_mxmn('PHY veloc mixing',jk,params_oce%A_veloc_v(:,:,:),n_zlev, &
! !     & p_patch%nblks_e,'phy',ipl_src)
! !   write(123,*)'max/min veloc mixing',jk,maxval(params_oce%A_veloc_v(:,jk,:)),&
! !   &minval(params_oce%A_veloc_v(:,jk,:))
! !  END DO
! ! 
! !  END SUBROUTINE update_ho_params_old



 
END MODULE mo_oce_physics
