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
USE mo_ocean_nml,           ONLY: n_zlev, bottom_drag_coeff, k_veloc_h, k_veloc_v,&
                            &  k_pot_temp_h, k_pot_temp_v, k_sal_h, k_sal_v, no_tracer,&
                            &  expl_vertical_velocity_diff
USE mo_parallel_config,  ONLY: nproma
USE mo_model_domain,        ONLY: t_patch
USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, &!! min_rledge, min_rlvert,      &
 &                               land, sea_boundary    
USE mo_exception,           ONLY: message, finish
USE mo_oce_state,           ONLY: t_hydro_ocean_state!, t_hydro_ocean_diag
USE mo_physical_constants,  ONLY: grav
USE mo_loopindices,         ONLY: get_indices_c
USE mo_math_constants,      ONLY: pi 
IMPLICIT NONE


PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

CHARACTER(len=*), PARAMETER :: this_mod_name = 'mo_oce_physics'


! Public interface

! public subroutines
PUBLIC :: construct_ho_physics
PUBLIC :: destruct_ho_physics
PUBLIC :: init_ho_physics
PUBLIC :: construct_ho_params
PUBLIC :: destruct_ho_params
PUBLIC :: init_ho_params
PUBLIC :: update_ho_params

! Parameters below appear directly in the ocean model/equation. They are eventually
! dynamically updated by using the "ocean-physics" structure. #slo# - not yet
TYPE t_ho_params

  REAL(wp),ALLOCATABLE ::  &
  ! diffusion coefficients for horizontal velocity, temp. and salinity, dim=(nproma, n_zlev,nblks_c)
    &  K_veloc_h   (:,:,:),  & ! coefficient of horizontal velocity diffusion
    &  K_tracer_h(:,:,:,:),  & ! coefficient of horizontal tracer diffusion
  ! diffusion coefficients for vertical velocity, temp. and salinity, dim=(n_zlev)
    &  A_veloc_v   (:,:,:),    & ! coefficient of vertical velocity diffusion
    &  A_tracer_v(:,:,:,:)       ! coefficient of vertical tracer diffusion

  !constant background values of coefficients above
  REAL(wp) :: K_veloc_h_back, &! coefficient of horizontal velocity diffusion
           &  A_veloc_v_back ! coefficient of vertical velocity diffusion


  REAL(wp),ALLOCATABLE ::  &
    &  K_tracer_h_back(:), &! coefficient of horizontal tracer i diffusion
    &  A_tracer_v_back(:)   ! coefficient of vertical temperature diffusion
  ! coefficients in linear EOS
  REAL(wp) :: a_T = 2.55E-04_wp,     & ! thermal expansion coefficient (kg/m3/K)
    &         b_S = 7.64E-01_wp         ! haline contraction coefficient (kg/m3/psu)

  ! density reference values, to be constant in Boussinesq ocean models
  REAL(wp) :: rho_ref = 1025.022_wp         ! reference density [kg/m^3]
  REAL(wp) :: rho_inv = 0.0009755881663_wp  ! inverse reference density [m^3/kg]

  REAL(wp) :: bottom_drag_coeff

END TYPE t_ho_params


!Data used for parametrization of unresolved processes. The list below is uncomplete.
TYPE t_ho_physics

! INTEGER,  ALLOCATABLE ::  &
!   &  Jwtype     (:,:)  &  ! Jerlov water type (clarity)
!   &  ksbl       (:,:)  &
!   &  kbbl       (:,:)  &

  REAL(wp), ALLOCATABLE ::   &
    &  Akv        (:,:,:),   &  ! vertical mixing coefficient (m2/s) for momentum
    &  Akt        (:,:,:,:)! &  ! vertical mixing coefficient (m2/s) for tracers
!   &  alpha       (:,:),    &  ! surface thermal expansion coefficient (1/Celsius)
!   &  beta       (:,:),     &  ! surface saline contraction coefficient (1/PSU)
!   &  bvf        (:,:,:),   &  ! Brunt-Vaisala frequency squared (1/s2)
!   &  shear2     (:,:),     &  ! velocity shear squared
!   &  Rig        (:,:,:),   &  ! Richardson gradient number
!   &  alphaobeta (:,:,:),   &  ! ratio of thermal expansion and saline contraction
!                            &  ! coefficients (Celsius/PSU) used in double diffusion
!   &  hsbl       (:,:),     &  ! depth of surface oceanic boundary layer (m)
!   &  hbbl       (:,:),     &  ! depth of bottom oceanic boundary layer (m)
!   &  ghats      (:,:,:,:), &  ! boundary layer nonlocal transport (T units/m)
!   &  Kv         (:,:,:),   &  ! vertical viscosity scratch array
!   &  Kt         (:,:,:),   &  ! vertical diffusion scratch array, temperature
!   &  Ks         (:,:,:),   &  ! vertical diffusion scratch array, salinity
!   &  aux_e      (:,:),     &  ! aux array with edge dimensions
!   &  aux_c      (:,:,:),   &  ! aux_array with cell dimensions
!   &  sl_dpth    (:,:),     &  ! depth of surface mixed layer (work array)
!   &  bl_dpth    (:,:),     &  !
!   &  Bo         (:,:),     &  ! surface buoyancy forcing
!   &  Bosol      (:,:),     &  ! surface radiative buoyancy forcing
!   &  swdk       (:,:),     &  ! shortwave (radiation) fractional decay
!   &  Bflux      (:,:,:),   &  ! total buoyancy flux
!   &  dR         (:,:,:),   &  ! delta Rho
!   &  d_U        (:,:,:),   &  ! delta vn (at cells)
!   &  d_V        (:,:,:),   &  ! delta vt (at cells)
!   &  FC         (:,:,:),   &  ! critical function
!   &  wm         (:,:),     &  ! turbulent velocity scale (m/s) for momentum
!   &  ws         (:,:),     &  ! turbulent velocity scale (m/s) for tracer
!   &  zgrid      (:,:),     &  ! aux. depth array
!   &  Bfsfc      (:,:),     &  ! total buoyancy flux at surface boundary layer depth
!   &  Bfbot      (:,:),     &  !
!   &  f1         (:,:),     &  !
!   &  Gm1        (:,:),     &  ! shape function for momentum
!   &  Gt1        (:,:),     &  ! shape function for temperature
!   &  Gs1        (:,:),     &  ! shape function for salinity
!   &  dGm1dS     (:,:),     &  ! derivative of shape function for momentum
!   &  dGt1dS     (:,:),     &  ! derivative of shape function for temperature
!   &  dGs1dS     (:,:)      &  ! derivative of shape function for salinity

END TYPE t_ho_physics


PUBLIC :: t_ho_physics
PUBLIC :: t_ho_params

CONTAINS
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Initialisation of ocean physics
  !!
  !! Initialisation of ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE init_ho_physics (  physics_oce )
   TYPE(t_ho_physics) :: physics_oce
    ! Local variables
    !-------------------------------------------------------------------------
    physics_oce%akv     (:,:,:) = 0.0_wp
  END SUBROUTINE init_ho_physics

  !-------------------------------------------------------------------------
  !
  !>
  !! Initialisation of ocean physics
  !!
  !! Initialisation of ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE init_ho_params(  p_phys_param )
   TYPE (t_ho_params) :: p_phys_param 
    ! Local variables
    INTEGER :: i

    !-------------------------------------------------------------------------

    !Init from namelist
    p_phys_param%K_veloc_h_back = k_veloc_h
    p_phys_param%A_veloc_v_back = k_veloc_v

    p_phys_param%K_veloc_h(:,:,:) = p_phys_param%K_veloc_h_back
    p_phys_param%A_veloc_v(:,:,:) = p_phys_param%A_veloc_v_back

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

    p_phys_param%rho_inv        = 1.0_wp / p_phys_param%rho_ref


    p_phys_param%bottom_drag_coeff = bottom_drag_coeff

  END SUBROUTINE init_ho_params
  !
  !-------------------------------------------------------------------------
  !
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE construct_ho_physics (ppatch, physics_oce)

    TYPE(t_patch), INTENT(IN) :: ppatch
    TYPE(t_ho_physics) :: physics_oce
    !TYPE (t_ho_params),        INTENT(INOUT)  :: p_phys_param 
    ! Local variables

    INTEGER   :: ist
    INTEGER   :: nblks_c

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':construct_ho_physics'

    !-------------------------------------------------------------------------

    CALL message(TRIM(routine), 'construct hydro ocean physics')

    ! determine size of arrays
    nblks_c = ppatch%nblks_c

    ! preliminary Akv
    ALLOCATE(physics_oce%akv(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vertical mixing coefficient failed')
    END IF

    ! preliminary Kh
    !ALLOCATE(params_oce%k_veloc_h(nproma,nblks_c), STAT=ist)
!     IF (ist/=SUCCESS) THEN
!       CALL finish(TRIM(routine), 'allocation for horizontal velocity diffusion failed')
!     END IF

  END SUBROUTINE construct_ho_physics
  !-------------------------------------------------------------------------
  !
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE construct_ho_params(ppatch, params_oce)

    TYPE(t_patch), INTENT(IN)         :: ppatch
    TYPE (t_ho_params), INTENT(INOUT) :: params_oce

    ! Local variables
    INTEGER   :: ist, i
    INTEGER   :: nblks_c, nblks_e

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':construct_ho_physics'

    !-------------------------------------------------------------------------

    CALL message(TRIM(routine), 'construct hydro ocean physics')

    ! determine size of arrays
    nblks_c = ppatch%nblks_c
    nblks_e = ppatch%nblks_e

    ! preliminary Kh
    ALLOCATE(params_oce%K_veloc_h(nproma,n_zlev,nblks_e), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for horizontal velocity diffusion failed')
     END IF

    ALLOCATE(params_oce%K_tracer_h(nproma,n_zlev,nblks_e, no_tracer), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for horizontal tracer diffusion failed')
     END IF

    ALLOCATE(params_oce%K_tracer_h_back(no_tracer), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for horizontal background tracer diffusion failed')
     END IF

    IF(expl_vertical_velocity_diff==0)THEN !explicit
      ALLOCATE(params_oce%A_veloc_v(nproma,n_zlev,nblks_c), STAT=ist)
       IF (ist/=SUCCESS) THEN
         CALL finish(TRIM(routine), 'allocation for vertical velocity diffusion failed')
       END IF
    ELSEIF(expl_vertical_velocity_diff==1)THEN !implicit
      ALLOCATE(params_oce%A_veloc_v(nproma,n_zlev,nblks_e), STAT=ist)
       IF (ist/=SUCCESS) THEN
         CALL finish(TRIM(routine), 'allocation for vertical velocity diffusion failed')
       END IF
    ENDIF


    ALLOCATE(params_oce%A_tracer_v(nproma,n_zlev,nblks_c,no_tracer), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for vertical tracer diffusion failed')
     END IF

    ALLOCATE(params_oce%A_tracer_v_back(no_tracer), STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'allocation for vertical tracer background diffusion failed')
     END IF

    params_oce%K_veloc_h    = 0.0_wp
    params_oce%A_veloc_v    = 0.0_wp

    DO i=1,no_tracer
      params_oce%K_tracer_h(:,:,:,i) = 0.0_wp
      params_oce%A_tracer_v(:,:,:,i) = 0.0_wp
     params_oce%K_tracer_h_back(i)   = 0.0_wp
      params_oce%A_tracer_v_back(i)  = 0.0_wp
    END DO
  END SUBROUTINE construct_ho_params
  !-------------------------------------------------------------------------
  !
  !>
  !! Destruction of arrays for ocean physics
  !!
  !! Destruction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE destruct_ho_physics(physics_oce)
    TYPE(t_ho_physics) :: physics_oce


    ! Local variables
    INTEGER   :: ist
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':destruct_ho_physics'

    !-------------------------------------------------------------------------

    CALL message(TRIM(routine), 'destruct hydro ocean physics')

    ! preliminary Akv
    DEALLOCATE(physics_oce%akv, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for vertical mixing coefficient failed')
    END IF
 
  END SUBROUTINE destruct_ho_physics
  !-------------------------------------------------------------------------
  !
  !>
  !! Construction of arrays for ocean physics
  !!
  !! Construction of arrays for ocean physics ...
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-07)
  !
  !
  SUBROUTINE destruct_ho_params(params_oce)

    TYPE (t_ho_params), INTENT(INOUT) :: params_oce

    ! Local variables
    INTEGER :: ist
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = this_mod_name//':destruct_ho_physics'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'construct hydro ocean physics')

    ! preliminary Kh
    DEALLOCATE(params_oce%K_veloc_h, STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for horizontal velocity diffusion failed')
     END IF

    DEALLOCATE(params_oce%K_tracer_h, STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for horizontal tracer diffusion failed')
     END IF

    DEALLOCATE(params_oce%A_veloc_v, STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for vertical velocity diffusion failed')
     END IF

    DEALLOCATE(params_oce%A_tracer_v, STAT=ist)
     IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for vertical temperaure diffusion failed')
     END IF


  END SUBROUTINE destruct_ho_params 

 !
  !-------------------------------------------------------------------------
  !
  !>
  !! Update of parameters
  !!
  !! Update of ocean physics: This routine is used used only if time-dependent
  !! changes of physical parametrizations.
  !! Currently vertical mixing coefficients for tracers are updated.
  !! Dependent on the Richardson number the diffusivity are calculated
  !!(Large & Gent JPO 29, (1999), 449-464).
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2011-02)
  !
  !
 SUBROUTINE update_ho_params(p_patch, p_os, params_oce)
   TYPE(t_patch), INTENT(IN)         :: p_patch
   TYPE(t_hydro_ocean_state), TARGET :: p_os
   TYPE (t_ho_params), INTENT(INOUT) :: params_oce
!   ! Local variables
   INTEGER  :: ist, jc, jb, jk, i_no_trac
   INTEGER  :: il_e1, ib_e1,il_e2, ib_e2,il_e3, ib_e3
   INTEGER  :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, rl_start_c, rl_end_c
   REAL(wp) :: z_bv, z_a, dz_inv, z_shear
   REAL(wp) :: z_Ri
   REAL(wp) :: z_edge_ave 
!   !-------------------------------------------------------------------------
    rl_start_c   = 1
    rl_end_c     = min_rlcell
    i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
    i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)


    z_bv   = 0.0_wp
    z_a    = 0.0_wp
    dz_inv = 0.0_wp
    z_shear= 0.0_wp
    DO i_no_trac=1, no_tracer
      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c( p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
        &                   rl_start_c, rl_end_c)

        DO jk = 2, n_zlev
          dz_inv = 1.0_wp/(p_patch%patch_oce%zlev_m(jk-1)-p_patch%patch_oce%zlev_m(jk))

          DO jc = i_startidx_c, i_endidx_c
            IF ( p_patch%patch_oce%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN

               z_bv  = dz_inv&
               &*(p_os%p_diag%rho(jc,jk-1,jb)-p_os%p_diag%rho(jc,jk,jb))&
               &/params_oce%rho_ref
               z_shear =&
               & DOT_PRODUCT(p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x,&
                              &p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)

               z_Ri = z_bv/(z_shear + 0.05_wp)
               IF(z_Ri>0.0_wp)THEN
                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac) &
                 &= params_oce%A_tracer_v_back(i_no_trac)&
                 & + ( params_oce%A_veloc_v(jc,jk,jb)/(1.0_wp+5.0_wp*z_Ri))

!                write(*,*)'Mixing coeffs',jc,jk,jb,&
!                &params_oce%A_tracer_v(jc,jk,jb,1),&
!                &params_oce%A_tracer_v_back(i_no_trac), z_Ri 

               ELSE
                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = params_oce%A_tracer_v_back(i_no_trac)
               ENDIF

! il_e1 = p_patch%cells%edge_idx(jc,jb,1)
! ib_e1 = p_patch%cells%edge_blk(jc,jb,1)
! il_e2 = p_patch%cells%edge_idx(jc,jb,2)
! ib_e2 = p_patch%cells%edge_blk(jc,jb,2)
! il_e3 = p_patch%cells%edge_idx(jc,jb,3)
! ib_e3 = p_patch%cells%edge_blk(jc,jb,3)
! 
! z_edge_ave = (p_patch%edges%primal_edge_length(il_e1,ib_e1)&
!             &+p_patch%edges%primal_edge_length(il_e2,ib_e2)&
!             &+p_patch%edges%primal_edge_length(il_e3,ib_e3))/3.0_wp
! params_oce%A_veloc_v(jc,jk,jb)= 4.0_wp*pi*cos(p_patch%cells%center(jc,jb)%lat&
! &* (sqrt(3.0_wp)*z_edge_ave/pi)*(sqrt(3.0_wp)*z_edge_ave/pi)*(sqrt(3.0_wp)*z_edge_ave/pi))
! write(*,*)' A_v',jc,jk,jb,params_oce%A_veloc_v(jc,jk,jb)

!               !------------------------------------
!               z_bv  = -grav*dz_inv&
!               &*(p_os%p_diag%rho(jc,jk-1,jb)-p_os%p_diag%rho(jc,jk,jb))&
!               &/params_oce%rho_ref
! 
!               IF (z_bv<0.0_wp)THEN
!                 z_a=1.0_wp 
!               ELSE
!                 z_shear =&
!                 & DOT_PRODUCT(p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x,&
!                              &p_os%p_diag%p_vn(jc,jk-1,jb)%x-p_os%p_diag%p_vn(jc,jk,jb)%x)       
!sum((Uelem(:,nz-1,nelem) - Uelem(:, nz, nelem))**2)
! !write(*,*)'shear',z_shear
!                 z_shear = z_shear*dz_inv*dz_inv
! 
!                 IF(z_shear==0.0_wp)THEN 
!                   z_a=1.0_wp
!                 ELSE
!                   z_a     = z_shear/(z_shear+10.0_wp*z_bv)
!                 ENDIF
!                 params_oce%A_tracer_v(jc,jk,jb, i_no_trac) = params_oce%A_tracer_v_back(i_no_trac)&
!                                                        & +0.05_wp*z_a*z_a*z_a
!                 write(*,*)'Mixing coeffs',jc,jk,jb,&
!                 &params_oce%A_tracer_v(jc,jk,jb,1),&
!                 &params_oce%A_tracer_v_back(i_no_trac), z_a 
               !------------------------------------------
               !ENDIF
            ENDIF
          END DO
        END DO
      END DO
    END DO

 END SUBROUTINE update_ho_params

END MODULE mo_oce_physics
