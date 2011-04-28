!>
!!  Contains the data structures for the hydrostatic ocean model.
!!
!!  Contains the data structures to store the hydrostatic & boussinesq ocean model state.
!!  Implementation is based on ICON-Shallow-Water model
!!  to store the shallow water model state and other auxiliary variables.
!!  Constructors and destructors for these data structures are also defined here.
!!
!! @par Revision History
!!  Initial version by Peter Korn (MPI-M), (2006).
!!  Big recoding by P. Korn (MPI-M), (2009/2010)
!!  Modification by Stephan Lorenz, MPI-M, (2010-03-19):
!!   - renaming and adjustment to ocean domain and patch_oce
!
!
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
!!
MODULE mo_oce_state
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,                ONLY: wp
  USE mo_run_nml,             ONLY: nproma
  USE mo_impl_constants,      ONLY: success, max_char_length, min_rledge, min_rlcell, boundary, sea
  USE mo_ocean_nml,           ONLY: n_zlev, ntrac_oce, t_ref, s_ref
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch
  USE mo_model_domain_import, ONLY: n_dom
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c

  USE mo_linked_list,  ONLY: t_var_list
  USE mo_var_list,     ONLY: default_var_list_settings, &
                           & add_var,                   &
                           & new_var_list,              &
                           & delete_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants

  IMPLICIT NONE
  PRIVATE

  ! subroutines
  PUBLIC :: construct_hydro_ocean_state
  PUBLIC :: destruct_hydro_ocean_state
  PUBLIC :: set_lateral_boundary_values
  !
  !variables
  PUBLIC :: t_hydro_ocean_state
  PUBLIC :: t_hydro_ocean_prog
  PUBLIC :: t_hydro_ocean_aux
  PUBLIC :: t_hydro_ocean_diag

  !constructors
  PRIVATE :: construct_hydro_ocean_diag
  PRIVATE :: construct_hydro_ocean_prog
  PRIVATE :: construct_hydro_ocean_aux
  !destructors
  PRIVATE :: destruct_hydro_ocean_diag
  PRIVATE :: destruct_hydro_ocean_prog
  PRIVATE :: destruct_hydro_ocean_aux


!
!! data structure defining model states
!

!
!! prognostic variables
!
  TYPE t_hydro_ocean_prog

    REAL(wp), POINTER ::    &
      &  h(:,:)                ,& ! height of the free surface. Unit: [m]
                                  ! dimension:(nproma, nblks_c)
      &  vn(:,:,:)             ,& ! velocity component normal to cell edge. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  tracer(:,:,:,:)          ! tracer concentration.
                                  ! dimension: (nproma, n_zlev, nblks_c, ntrac_oce)
                                  ! Ordering of tracers:
                                  !   1) pot_temp:= potential temperature, Unit: [deg C]
                                  !   2) salinity:= salinity, Unit [psu]

  END TYPE t_hydro_ocean_prog

!
!! diagnostic variables
!
  TYPE t_hydro_ocean_diag

    REAL(wp), POINTER ::        &
      &  vt(:,:,:)             ,& ! tangential velocity component at edges. Unit [m/s].
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  rho(:,:,:)            ,& ! density. Unit: [kg/m^3]
                                  ! dimension: (nproma,n_zlev, nblks_c)
      &  h_e(:,:)              ,& ! surface height at cell edges. Unit [m].
                                  ! dimension: (nproma, nblks_e)
      &  thick_c(:,:)          ,& ! individual fluid column thickness at cells. Unit [m].
                                  ! dimension: (nproma, nblks_c)
      &  thick_e(:,:)          ,& ! individual fluid column thickness at edges. Unit [m].
                                  ! dimension: (nproma, nblks_e)
      &  w(:,:,:)              ,& ! vertical velocity. Unit [m/s].
                                  ! dimension: (nproma, n_zlev+1, nblks_c)
      &  w_old(:,:,:)          ,& ! vertical velocity from previous timestep. Unit [m/s].
                                  ! dimension: (nproma, n_zlev+1, nblks_c)
      &  w_e(:,:,:)            ,& ! vertical velocity at edges. Unit [m/s]
                                  ! dimension: (nproma, n_zlev+1, nblks_e)
      &  w_prev(:,:,:)         ,& ! vertical velocity at cells, from previous timestep. Unit [m/s]
                                  ! dimension: (nproma, n_zlev+1, nblks_c)
      &  u(:,:,:)              ,& ! reconstructed zonal velocity component. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  v(:,:,:)              ,& ! reconstructed meridional velocity component. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  ptp_vn(:,:,:)         ,& ! normal velocity after mapping P^T P
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_pred(:,:,:)        ,& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vn_impl_vert_diff(:,:,:),& ! predicted normal velocity vector at edges.
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  vort(:,:,:)           ,& ! vorticity at triangle vertices. Unit [1/s]
                                  ! dimension: (nproma, n_zlev, nblks_v)
      &  vort_e(:,:,:)         ,& ! vorticity interpolated to triangle edges. Unit [1/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  kin(:,:,:)            ,& ! kinetic energy. Unit [m/s].
                                  ! (nproma, n_zlev, nblks_c)
      &  veloc_adv_horz(:,:,:) ,& ! horizontal velocity advection
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  veloc_adv_vert(:,:,:) ,& ! vertical velocity advection
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  laplacian_horz(:,:,:) ,& ! horizontal diffusion of horizontal velocity
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  laplacian_vert(:,:,:) ,& ! vertical diffusion of horizontal velocity
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  grad(:,:,:)           ,& ! gradient of kinetic energy. Unit [m/s]
                                  ! dimension: (nproma,n_zlev, nblks_e)
      &  div(:,:,:)            ,& ! divergence. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  press_hyd(:,:,:)      ,& ! hydrostatic pressure. Unit [m]
                                  ! dimension: (nproma, n_zlev, nblks_c)
      &  press_grad(:,:,:)       ! hydrostatic pressure gradient term. Unit [m/s]
                                  ! dimension: (nproma, n_zlev, nblks_e)

 TYPE(t_cartesian_coordinates), POINTER :: p_vn(:,:,:) ! reconstructed velocity at cell center in cartesian coordinates
                                                           ! dimension: (nproma, n_zlev, nblks_c)

  END TYPE t_hydro_ocean_diag

!
!! auxiliary data
!
  TYPE t_hydro_ocean_aux

    REAL(wp), POINTER ::    &
      &  g_n(:,:,:)            ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                  ! at timelevel n
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  g_nm1(:,:,:)          ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                  ! at timelevel n-1
                                  ! dimension: (nproma, n_zlev, nblks_e)
      &  g_nimd(:,:,:)         ,& ! explicit velocity term in Adams-Bashford time marching routines,
                                  ! located at intermediate timelevel
      &  g_n_c_h(:,:,:,:)      ,& ! explicit tracer term in Adams-Bashford time marching routines,
                                  ! at timelevel n for each tracer, horizontal
                                  ! dimension: (nproma, n_zlev, nblks_c, ntrac_oce )
      &  g_nm1_c_h(:,:,:,:)    ,& ! explicit tracer term in Adams-Bashford time marching routines,
                                  ! at timelevel n-1 for each tracer, horizontal
                                  ! dimension: (nproma, n_zlev, nblks_c, ntrac_oce)
      &  g_nimd_c_h(:,:,:,:)   ,& ! explicit tracer term in Adams-Bashford time marching routines,
                                  ! located at intermediate timelevel for each tracer, horizontal
                                  ! dimension: (nproma, n_zlev, nblks_c,ntrac_oce )
      &  g_n_c_v(:,:,:,:)      ,& ! explicit tracer term in Adams-Bashford time marching routines,
                                  ! at timelevel n for each tracer, vertical
                                  ! dimension: (nproma, n_zlev, nblks_c, ntrac_oce )
      &  g_nm1_c_v(:,:,:,:)    ,& ! explicit tracer term in Adams-Bashford time marching routines,
                                  ! at timelevel n-1 for each tracer, vertical
                                  ! dimension: (nproma, n_zlev, nblks_c, ntrac_oce)
      &  g_nimd_c_v(:,:,:,:)   ,& ! explicit tracer term in Adams-Bashford time marching routines,
                                  ! located at intermediate timelevel for each tracer, vertical
                                  ! dimension: (nproma, n_zlev, nblks_c,ntrac_oce )
      &  bc_top_vn(:,:)         ,& ! normal velocity boundary condition at surface
                                  ! dimension: (nproma,nblks_e)
      &  bc_bot_vn(:,:)         ,& ! normal velocity boundary condition at bottom
                                  ! dimension: (nproma,nblks_c)
      &  bc_top_u(:,:)         ,& ! zonal velocity boundary condition at surface
                                  ! dimension: (nproma,nblks_c)
      &  bc_top_v(:,:)         ,& ! meridional velocity boundary condition at surface
                                  ! dimension: (nproma,nblks_c)
      &  bc_bot_u(:,:)         ,& ! zonal velocity boundary condition at bottom
                                  ! dimension: (nproma,nblks_c)
      &  bc_bot_v(:,:)         ,& ! meridional velocity boundary condition at bottom
                                  ! dimension: (nproma,nblks_c)
      &  bc_top_w(:,:)         ,& ! vertical velocity boundary condition at surface
                                  ! dimension: (nproma,nblks_c)
      &  bc_bot_w(:,:)         ,& ! vertical velocity boundary condition at bottom
      &  bc_top_tracer(:,:,:)  ,& ! vertical velocity boundary condition at surface
                                  ! dimension: (nproma,nblks_c)
      &  bc_bot_tracer(:,:,:)  ,& ! vertical velocity boundary condition at bottom
      &  p_rhs_sfc_eq(:,:)!,      & ! right hand side of surface equation
                                         ! dimension: (nproma,nblks_c)
     TYPE(t_cartesian_coordinates), POINTER :: bc_top_veloc_cc(:,:), &
                                  &                bc_bot_veloc_cc(:,:)
  END TYPE t_hydro_ocean_aux

!
!! array of states
!
  TYPE t_hydro_ocean_state

    TYPE(t_hydro_ocean_prog), POINTER :: p_prog(:)    ! time array of prognostic states at
                                                        ! different time levels
    TYPE(t_hydro_ocean_diag) :: p_diag
    TYPE(t_hydro_ocean_aux)  :: p_aux

  END TYPE t_hydro_ocean_state

!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
!

!>
!! Constructor for hydrostatic ocean state + diagnostic and auxiliary  states.
!!
!! Constructor for hydrostatic ocean state
!! It calls  constructors to single time level
!! auxiliary and diagnostic states. Then it constructs state array,
!! whose components (representing multiple time levels).
!! Initialization of all components with zero.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2007).
!!  Modification by Stephan Lorenz, MPI-M, (2010-06-01) - no temporary memory array
!

  SUBROUTINE construct_hydro_ocean_state( p_patch, p_os )

    TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)
    TYPE(t_hydro_ocean_state), TARGET :: p_os(n_dom)

    ! local variables
    INTEGER           :: jg
    INTEGER           :: i_status, jp, prlength ! local prognostic array length

    ! variabels for dynamic variable list construction
    TYPE(t_var_list),POINTER :: list

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_state'

    CALL message(TRIM(routine), 'start to construct hydro_ocean state' )

    ! Using Adams-Bashforth semi-implicit timestepping with 3 prognostic time levels:
    prlength = 3

    CALL new_var_list(list, 'ocean_var_list')
    DO jg = 1, n_dom

      ALLOCATE(p_os(jg)%p_prog(1:prlength), STAT=i_status)
      IF (i_status/=SUCCESS) THEN
         CALL finish(TRIM(routine), 'allocation of progn. state array failed')
      END IF

      ! construction loop: create components of state array
      DO jp = 1, prlength
         !CALL construct_hydro_ocean_prog(p_patch(jg), p_os(jg)%p_prog(jp))
         CALL add_hydro_ocean_prog_vars(list, p_patch(jg), p_os(jg)%p_prog(jp))
      END DO

     !CALL construct_hydro_ocean_prog(p_patch(jg), p_os(jg)%p_prog)

      CALL construct_hydro_ocean_diag(p_patch(jg), p_os(jg)%p_diag)
      CALL add_hydro_ocean_diag_vars( list, p_patch(jg), p_os(jg)%p_diag)

      CALL construct_hydro_ocean_aux(p_patch(jg), p_os(jg)%p_aux)

      CALL message(TRIM(routine),'construction of hydrostatic ocean state finished')

    END DO

  END SUBROUTINE construct_hydro_ocean_state

!-------------------------------------------------------------------------
!>
!!               Destructor for hydrostatic ocean state.
!
!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2006).
!!
  SUBROUTINE destruct_hydro_ocean_state(p_os)
    TYPE(t_hydro_ocean_state), TARGET,INTENT(inout)   :: p_os(n_dom)

    ! local variables

    INTEGER                                   :: jg, prlength, ist, jp

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_state'


!-------------------------------------------------------------------------

    CALL message(TRIM(routine), 'start to destruct hydro ocean state ')

    prlength = SIZE(p_os(1)%p_prog)

    IF (prlength==0) THEN
      CALL finish(TRIM(routine),'prog array has length zero')
    END IF

    DO jg=1, n_dom

      !destruction loop for array components
      DO jp = 1, prlength
         CALL destruct_hydro_ocean_prog(p_os(jg)%p_prog(jp))
      END DO

!     CALL destruct_hydro_ocean_prog(p_os(jg)%p_prog)
      CALL destruct_hydro_ocean_diag(p_os(jg)%p_diag)
      CALL destruct_hydro_ocean_aux (p_os(jg)%p_aux)

      ! destruct state array
      ist = 1
      DEALLOCATE(p_os(jg)%p_prog, STAT=ist)
      IF (ist/=SUCCESS) THEN
         CALL finish(TRIM(routine),'deallocation of state array failed')
      END IF

    END DO

    CALL message(TRIM(routine),'destruction of hydrostatic ocean state finished')

  END SUBROUTINE destruct_hydro_ocean_state

!-------------------------------------------------------------------------
!>
!!               Allocation of components of hydrostatic ocean prognostic state.
!!               Initialization of components with zero.
!
!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2006).
!!
  SUBROUTINE construct_hydro_ocean_prog(p_patch, p_os_prog)

!
    TYPE(t_patch), INTENT(in), TARGET         :: p_patch
    TYPE(t_hydro_ocean_prog), INTENT(inout)   :: p_os_prog


    ! local variables

    INTEGER  :: n
    INTEGER  :: nblks_c, nblks_e !, nblks_v
    INTEGER  :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_prog'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine), 'start to construct hydro_ocean prognostic state')

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! allocate prognostic part of state vector
    !
    ! height
!   ALLOCATE(p_os_prog%h(nproma,nblks_c), STAT=ist)
!   IF (ist/=SUCCESS) THEN
!      CALL finish(TRIM(routine),'allocation for h on cells failed')
!   END IF

    ! normal velocity component
!   ALLOCATE(p_os_prog%vn(nproma,n_zlev,nblks_e), STAT=ist)
!   IF (ist/=SUCCESS) THEN
!      CALL finish(TRIM(routine),'allocation for normal velocity on edges failed')
!   END IF

    ! Tracer
    !ALLOCATE(p_os_prog%tracer(nproma,n_zlev,nblks_c, ntrac_oce), STAT=ist)
    !IF (ist/=SUCCESS) THEN
       !CALL finish(TRIM(routine),'allocation for tracers temp. and sal. on cells failed')
    !END IF


    ! initialize all components with zero
    !
!    p_os_prog%h(:,:)    = 0.0_wp
!    p_os_prog%vn(:,:,:) = 0.0_wp
    !DO n=1,ntrac_oce
      !p_os_prog%tracer(:,:,:,n) = 0.0_wp
    !END DO

    ! initialize tracers temperature and salinity with reference values from namelist
    !
    !p_os_prog%tracer(:,:,:,1) = t_ref  !  temperature
    !p_os_prog%tracer(:,:,:,2) = s_ref  !  salinity

    !CALL message(TRIM(routine), 'construction of hydrostatic ocean prognostic state finished')

  END SUBROUTINE construct_hydro_ocean_prog

  SUBROUTINE add_hydro_ocean_prog_vars(list, p_patch, p_os_prog)

    TYPE(t_var_list), POINTER :: list
    TYPE(t_patch), INTENT(in), TARGET         :: p_patch
    TYPE(t_hydro_ocean_prog), INTENT(inout)   :: p_os_prog


    INTEGER  :: n
    INTEGER  :: nblks_c, nblks_e !, nblks_v
    INTEGER  :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_prog'

    !-------------------------------------------------------------------------
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! height
    CALL add_var(list, 'h', p_os_prog%h , GRID_UNSTRUCTURED, &
    &            t_cf_var('h', 'm', 'surface elevation at cell center'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL, ZAXIS_HYBRID),&
    &            ldims=(/nproma,nblks_c/))

    !! normal velocity component
    CALL add_var(list, 'vn', p_os_prog%vn , GRID_UNSTRUCTURED, &
    &            t_cf_var('vn', 'm/s', 'normale velocity on edge,m'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_EDGE, ZAXIS_HYBRID),&
    &            ldims=(/nproma,n_zlev,nblks_e/))

    !! Tracers
    CALL add_var(list, 'tracers', p_os_prog%tracer , GRID_UNSTRUCTURED, &
    &            t_cf_var('tracers', '', '1:temperature 2:salinity'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL, ZAXIS_HYBRID),&
    &            ldims=(/nproma,n_zlev,nblks_c,ntrac_oce/) )

  END SUBROUTINE

  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean diagnostic state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!

  SUBROUTINE construct_hydro_ocean_diag(p_patch,p_os_diag)

    TYPE(t_patch), TARGET, INTENT(IN)          :: p_patch
    TYPE(t_hydro_ocean_diag), INTENT(INOUT)    :: p_os_diag

    ! local variables

    INTEGER :: ist
    INTEGER :: nblks_c, nblks_e, nblks_v
    INTEGER ::  jc,jb,jk, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_diag'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine), 'start to construct diagnostic hydro ocean state')

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! density
!   ALLOCATE(p_os_diag%rho(nproma,n_zlev,nblks_c), STAT=ist)
!   IF (ist/=SUCCESS) THEN
!   CALL finish(TRIM(routine), 'allocation for density failed')
!   END IF

    ! tangential velocity component
    ALLOCATE(p_os_diag%vt(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for tangential velocity at edges failed')
    END IF

    ALLOCATE(p_os_diag%h_e(nproma,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for surface height at edges failed')
    END IF

    ALLOCATE(p_os_diag%thick_c(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for fluid column thickness at cells failed')
    END IF
    ALLOCATE(p_os_diag%thick_e(nproma,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for fluid column thickness at edges failed')
    END IF
    ALLOCATE(p_os_diag%w(nproma,n_zlev+1,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vertical velocity at cells failed')
    END IF
    ALLOCATE(p_os_diag%w_old(nproma,n_zlev+1,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vertical velocity at cells failed')
    END IF

    ALLOCATE(p_os_diag%w_e(nproma,n_zlev+1,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vertical velocity at edges failed')
    END IF

    ALLOCATE(p_os_diag%w_prev(nproma,n_zlev+1,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vert. veloc. at cells at prev. timestep failed')
    END IF

    ! reconstructed u velocity component
    ALLOCATE(p_os_diag%u(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for u at cells failed')
    END IF

    ! reconstructed v velocity component
    ALLOCATE(p_os_diag%v(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for v at cells failed')
    END IF

    !reconstrcuted velocity in cartesian coordinates
    ALLOCATE(p_os_diag%p_vn(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for p_vn at cells failed')
    END IF

    !reconstrcuted velocity in cartesian coordinates
    ALLOCATE(p_os_diag%ptp_vn(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for ptp_vn at edges failed')
    END IF

    ! predicted vn normal velocity component
    ALLOCATE(p_os_diag%vn_pred(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vn_pred at cells failed')
    END IF
    ! predicted vn normal velocity component
    ALLOCATE(p_os_diag%vn_impl_vert_diff(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vn_impl_vert_diff failed')
    END IF

    ! vorticity
    write(*,*) 'allocation of vort in shape ',(/nproma,n_zlev,nblks_v/)
    ALLOCATE(p_os_diag%vort(nproma,n_zlev,nblks_v), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for vorticity at vertices failed')
    END IF

    ALLOCATE(p_os_diag%vort_e(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for vorticity at edges failed')
    END IF

    ! kinetic energy component
    ALLOCATE(p_os_diag%kin(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for kinetic energy at cells failed')
    END IF

    ! gradient term
    ALLOCATE(p_os_diag%grad(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for gradient term at edges failed')
    END IF

    ! divergence component
    ALLOCATE(p_os_diag%div(nproma,n_zlev,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for divergence at cells failed')
    END IF

    ! hydrostatic pressure
    ALLOCATE(p_os_diag%press_hyd(nproma,n_zlev,nblks_c) , STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for hydrostatic pressure at cells failed')
    END IF

    ! pressure gradient
    ALLOCATE(p_os_diag%press_grad(nproma,n_zlev,nblks_e) , STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'allocation for pressure gradient at edges failed')
    END IF

!     ! mass flux
!     ALLOCATE(p_os_diag%flux_mass(nproma,n_zlev,nblks_e) , STAT=ist)
!     IF (ist/=SUCCESS) THEN
!       CALL finish(TRIM(routine),'allocation for mass flux at edges failed')
!     END IF
!     ! tracer flux
!     ALLOCATE(p_os_diag%flux_tracer(nproma,n_zlev,nblks_e), STAT=ist)
!     IF (ist/=SUCCESS) THEN
!       CALL finish(TRIM(routine),'allocation for tracer flux at edges failed')
!     END IF

    ! horizontal velocity advection
    ALLOCATE(p_os_diag%veloc_adv_horz(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation for horizontal velocity advection at edges failed')
    ENDIF

    ! vertical velocity advection
    ALLOCATE(p_os_diag%veloc_adv_vert(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation for vertical velocity advection at edges failed')
    ENDIF

    ! horizontal diffusion
    ALLOCATE(p_os_diag%laplacian_horz(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation for horizontal diff. of hor. veloc. at edges failed')
    ENDIF

    ALLOCATE(p_os_diag%laplacian_vert(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation for vertical diff. of hor. veloc. at edges failed')
    ENDIF

    ! initialize all components with zero (this is preliminary)
    p_os_diag%vt        = 0.0_wp
!    p_os_diag%rho       = 0.0_wp
    p_os_diag%h_e       = 0.0_wp
    p_os_diag%thick_c   = 0.0_wp
    p_os_diag%thick_e   = 0.0_wp
    p_os_diag%w         = 0.0_wp
    p_os_diag%w_old     = 0.0_wp
    p_os_diag%w_e       = 0.0_wp
    p_os_diag%w_prev    = 0.0_wp
    p_os_diag%u         = 0.0_wp
    p_os_diag%v         = 0.0_wp
    p_os_diag%ptp_vn    = 0.0_wp
    p_os_diag%vn_pred   = 0.0_wp
    p_os_diag%vort      = 0.0_wp
    p_os_diag%vort_e    = 0.0_wp
    p_os_diag%kin       = 0.0_wp
    p_os_diag%grad      = 0.0_wp
    p_os_diag%div       = 0.0_wp
    p_os_diag%press_hyd = 0.0_wp
    p_os_diag%press_grad= 0.0_wp
    p_os_diag%vn_impl_vert_diff= 0.0_wp

    p_os_diag%veloc_adv_horz=0.0_wp
    p_os_diag%veloc_adv_vert=0.0_wp
    p_os_diag%laplacian_horz=0.0_wp
    p_os_diag%laplacian_vert=0.0_wp

   rl_start = 1
   rl_end = min_rlcell

   i_startblk = p_patch%cells%start_blk(rl_start,1)
   i_endblk   = p_patch%cells%end_blk(rl_end,1)
   DO jk=1,n_zlev
     DO jb = i_startblk, i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
                        & i_startidx, i_endidx, rl_start, rl_end)
       DO jc = i_startidx, i_endidx
         p_os_diag%p_vn(jc,jk,jb)%x=0.0_wp
       END DO
    END DO
  END DO


    !CALL message(TRIM(routine), 'construction of hydrostatic ocean diagnostic state finished')

  END SUBROUTINE construct_hydro_ocean_diag

  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean auxiliary state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!

  SUBROUTINE add_hydro_ocean_diag_vars(list,p_patch,p_os_diag)

    TYPE(t_var_list), POINTER :: list
    TYPE(t_patch), TARGET, INTENT(IN)          :: p_patch
    TYPE(t_hydro_ocean_diag), INTENT(INOUT)    :: p_os_diag

    ! local variables

    INTEGER :: ist
    INTEGER :: nblks_c, nblks_e, nblks_v
    INTEGER ::  jc,jb,jk, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_diag'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine), 'start to construct diagnostic hydro ocean state')

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! density
    CALL add_var(list, 'rho', p_os_diag%rho , GRID_UNSTRUCTURED, &
    &            t_cf_var('rho', 'kg/m^3', 'density'),&
    &            t_grib2_var(255, 255, 255, 16, GRID_REFERENCE, GRID_CELL, ZAXIS_HYBRID),&
    &            ldims=(/nproma,n_zlev,nblks_c/))

  ! ! tangential velocity component
  ! ALLOCATE(p_os_diag%vt(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for tangential velocity at edges failed')
  ! END IF

  ! ALLOCATE(p_os_diag%h_e(nproma,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for surface height at edges failed')
  ! END IF

  ! ALLOCATE(p_os_diag%thick_c(nproma,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for fluid column thickness at cells failed')
  ! END IF
  ! ALLOCATE(p_os_diag%thick_e(nproma,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for fluid column thickness at edges failed')
  ! END IF
  ! ALLOCATE(p_os_diag%w(nproma,n_zlev+1,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vertical velocity at cells failed')
  ! END IF
  ! ALLOCATE(p_os_diag%w_old(nproma,n_zlev+1,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vertical velocity at cells failed')
  ! END IF

  ! ALLOCATE(p_os_diag%w_e(nproma,n_zlev+1,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vertical velocity at edges failed')
  ! END IF

  ! ALLOCATE(p_os_diag%w_prev(nproma,n_zlev+1,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vert. veloc. at cells at prev. timestep failed')
  ! END IF

  ! ! reconstructed u velocity component
  ! ALLOCATE(p_os_diag%u(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for u at cells failed')
  ! END IF

  ! ! reconstructed v velocity component
  ! ALLOCATE(p_os_diag%v(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for v at cells failed')
  ! END IF

  ! !reconstrcuted velocity in cartesian coordinates
  ! ALLOCATE(p_os_diag%p_vn(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for p_vn at cells failed')
  ! END IF

  ! !reconstrcuted velocity in cartesian coordinates
  ! ALLOCATE(p_os_diag%ptp_vn(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for ptp_vn at edges failed')
  ! END IF

  ! ! predicted vn normal velocity component
  ! ALLOCATE(p_os_diag%vn_pred(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vn_pred at cells failed')
  ! END IF
  ! ! predicted vn normal velocity component
  ! ALLOCATE(p_os_diag%vn_impl_vert_diff(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vn_impl_vert_diff failed')
  ! END IF

  ! ! vorticity
  ! ALLOCATE(p_os_diag%vort(nproma,n_zlev,nblks_v), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine),'allocation for vorticity at vertices failed')
  ! END IF

  ! ALLOCATE(p_os_diag%vort_e(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for vorticity at edges failed')
  ! END IF

  ! ! kinetic energy component
  ! ALLOCATE(p_os_diag%kin(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine),'allocation for kinetic energy at cells failed')
  ! END IF

  ! ! gradient term
  ! ALLOCATE(p_os_diag%grad(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine),'allocation for gradient term at edges failed')
  ! END IF

  ! ! divergence component
  ! ALLOCATE(p_os_diag%div(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine),'allocation for divergence at cells failed')
  ! END IF

  ! ! hydrostatic pressure
  ! ALLOCATE(p_os_diag%press_hyd(nproma,n_zlev,nblks_c) , STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for hydrostatic pressure at cells failed')
  ! END IF

  ! ! pressure gradient
  ! ALLOCATE(p_os_diag%press_grad(nproma,n_zlev,nblks_e) , STAT=ist)
  ! IF (ist/=SUCCESS) THEN
  !   CALL finish(TRIM(routine), 'allocation for pressure gradient at edges failed')
  ! END IF

! !   ! mass flux
! !   ALLOCATE(p_os_diag%flux_mass(nproma,n_zlev,nblks_e) , STAT=ist)
! !   IF (ist/=SUCCESS) THEN
! !     CALL finish(TRIM(routine),'allocation for mass flux at edges failed')
! !   END IF
! !   ! tracer flux
! !   ALLOCATE(p_os_diag%flux_tracer(nproma,n_zlev,nblks_e), STAT=ist)
! !   IF (ist/=SUCCESS) THEN
! !     CALL finish(TRIM(routine),'allocation for tracer flux at edges failed')
! !   END IF

  ! ! horizontal velocity advection
  ! ALLOCATE(p_os_diag%veloc_adv_horz(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation for horizontal velocity advection at edges failed')
  ! ENDIF

  ! ! vertical velocity advection
  ! ALLOCATE(p_os_diag%veloc_adv_vert(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation for vertical velocity advection at edges failed')
  ! ENDIF

  ! ! horizontal diffusion
  ! ALLOCATE(p_os_diag%laplacian_horz(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation for horizontal diff. of hor. veloc. at edges failed')
  ! ENDIF

  ! ALLOCATE(p_os_diag%laplacian_vert(nproma,n_zlev,nblks_e), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation for vertical diff. of hor. veloc. at edges failed')
  ! ENDIF

  ! ! initialize all components with zero (this is preliminary)
  ! p_os_diag%vt        = 0.0_wp
  ! p_os_diag%rho       = 0.0_wp
  ! p_os_diag%h_e       = 0.0_wp
  ! p_os_diag%thick_c   = 0.0_wp
  ! p_os_diag%thick_e   = 0.0_wp
  ! p_os_diag%w         = 0.0_wp
  ! p_os_diag%w_old     = 0.0_wp
  ! p_os_diag%w_e       = 0.0_wp
  ! p_os_diag%w_prev    = 0.0_wp
  ! p_os_diag%u         = 0.0_wp
  ! p_os_diag%v         = 0.0_wp
  ! p_os_diag%ptp_vn    = 0.0_wp
  ! p_os_diag%vn_pred   = 0.0_wp
  ! p_os_diag%vort      = 0.0_wp
  ! p_os_diag%vort_e    = 0.0_wp
  ! p_os_diag%kin       = 0.0_wp
  ! p_os_diag%grad      = 0.0_wp
  ! p_os_diag%div       = 0.0_wp
  ! p_os_diag%press_hyd = 0.0_wp
  ! p_os_diag%press_grad= 0.0_wp
  ! p_os_diag%vn_impl_vert_diff= 0.0_wp

  ! p_os_diag%veloc_adv_horz=0.0_wp
  ! p_os_diag%veloc_adv_vert=0.0_wp
  ! p_os_diag%laplacian_horz=0.0_wp
  ! p_os_diag%laplacian_vert=0.0_wp

  !rl_start = 1
  !rl_end = min_rlcell

  !i_startblk = p_patch%cells%start_blk(rl_start,1)
  !i_endblk   = p_patch%cells%end_blk(rl_end,1)
  !DO jk=1,n_zlev
  !  DO jb = i_startblk, i_endblk
  !    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
  !                     & i_startidx, i_endidx, rl_start, rl_end)
  !    DO jc = i_startidx, i_endidx
  !      p_os_diag%p_vn(jc,jk,jb)%x=0.0_wp
  !    END DO
  ! END DO
  !END DO


  ! !CALL message(TRIM(routine), 'construction of hydrostatic ocean diagnostic state finished')

  END SUBROUTINE

  !-------------------------------------------------------------------------
  !>
  !!               Allocation of components of hydrostatic ocean auxiliary state.
  !!               Initialization of components with zero.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
  SUBROUTINE construct_hydro_ocean_aux(p_patch, p_os_aux)

    TYPE(t_patch),TARGET, INTENT(IN)                :: p_patch
    TYPE(t_hydro_ocean_aux), TARGET,INTENT(INOUT)   :: p_os_aux
!   INTEGER, INTENT(IN)                           :: k_no_temp_mem

    ! local variables

    INTEGER ::  ist, jc,jb, rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER ::  nblks_c, nblks_e, nblks_v

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:construct_hydro_ocean_aux'
!-------------------------------------------------------------------------
    !CALL message(TRIM(routine), 'start to construct hydro ocean auxiliary state')

    ! determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! allocation for Adam-Bashford time stepping
    ALLOCATE(p_os_aux%g_n(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of expl. term at step n in AB scheme failed')
    ENDIF

    ALLOCATE(p_os_aux%g_nm1(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of expl. term at step n-1 in AB scheme failed')
    ENDIF

    ALLOCATE(p_os_aux%g_nimd(nproma,n_zlev,nblks_e), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of expl. term at intermed. step in AB scheme failed')
    ENDIF

    ALLOCATE(p_os_aux%g_n_c_h(nproma,n_zlev,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of horz expl. tracer term at&
                                & step n in AB scheme failed')
    ENDIF
    ALLOCATE(p_os_aux%g_nm1_c_h(nproma,n_zlev,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of horz expl. tracer term at &
                                   &step n-1 in AB scheme failed')
    ENDIF
    ALLOCATE(p_os_aux%g_nimd_c_h(nproma,n_zlev,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of horz expl. tracer term at&
                              & intermed. step in AB scheme failed')
    ENDIF

    ALLOCATE(p_os_aux%g_n_c_v(nproma,n_zlev,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of vert expl. tracer term at&
                                & step n in AB scheme failed')
    ENDIF
    ALLOCATE(p_os_aux%g_nm1_c_v(nproma,n_zlev,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of vert expl. tracer term at &
                                   &step n-1 in AB scheme failed')
    ENDIF
    ALLOCATE(p_os_aux%g_nimd_c_v(nproma,n_zlev,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of vert expl. tracer term at&
                              & intermed. step in AB scheme failed')
    ENDIF

    ALLOCATE(p_os_aux%p_rhs_sfc_eq(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'allocation of RHS of surface eq. failed')
    ENDIF

    ! allocation for boundary conditions
    ALLOCATE(p_os_aux%bc_top_u(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond u failed')
    END IF
    ALLOCATE(p_os_aux%bc_top_v(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond v failed')
    END IF
    ALLOCATE(p_os_aux%bc_top_veloc_cc(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    ALLOCATE(p_os_aux%bc_top_vn(nproma,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond vn failed')
    END IF

    ALLOCATE(p_os_aux%bc_bot_u(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of bottom boundary cond u failed')
    END IF
    ALLOCATE(p_os_aux%bc_bot_v(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of bottom boundary cond v failed')
    END IF
    ALLOCATE(p_os_aux%bc_bot_veloc_cc(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond cc failed')
    END IF
    ALLOCATE(p_os_aux%bc_bot_vn(nproma,nblks_e), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of bottom boundary cond vn failed')
    END IF

    ALLOCATE(p_os_aux%bc_bot_w(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of bottom boundary cond w failed')
    END IF
    ALLOCATE(p_os_aux%bc_top_w(nproma,nblks_c), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond w failed')
    END IF
    ALLOCATE(p_os_aux%bc_bot_tracer(nproma,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of bottom boundary cond tracer failed')
    END IF
    ALLOCATE(p_os_aux%bc_top_tracer(nproma,nblks_c,ntrac_oce), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation of top boundary cond tracer failed')
    END IF

  ! ! allocation for divergence of fluxes
  ! ALLOCATE(p_os_aux%p_div_flux_horiz_act(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation of div flux horiz act failed')
  ! ENDIF
  !
  ! ALLOCATE(p_os_aux%p_div_flux_vert_act(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation of div flux vert act failed')
  ! ENDIF
  !
  ! ALLOCATE(p_os_aux%p_div_flux_horiz_prev(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation of div flux horiz prev failed')
  ! ENDIF
  !
  ! ALLOCATE(p_os_aux%p_div_flux_vert_prev(nproma,n_zlev,nblks_c), STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'allocation of div flux vert prev failed')
  ! ENDIF

    ! initialize all components with zero (this is preliminary)
    p_os_aux%g_n                   = 0.0_wp
    p_os_aux%g_nm1                 = 0.0_wp
    p_os_aux%g_nimd                = 0.0_wp
    p_os_aux%g_n_c_h               = 0.0_wp
    p_os_aux%g_nm1_c_h             = 0.0_wp
    p_os_aux%g_nimd_c_h            = 0.0_wp
    p_os_aux%g_n_c_v               = 0.0_wp
    p_os_aux%g_nm1_c_v             = 0.0_wp
    p_os_aux%g_nimd_c_v            = 0.0_wp
    p_os_aux%p_rhs_sfc_eq          = 0.0_wp
    p_os_aux%bc_top_u              = 0.0_wp
    p_os_aux%bc_top_v              = 0.0_wp
    p_os_aux%bc_top_vn             = 0.0_wp
    p_os_aux%bc_bot_u              = 0.0_wp
    p_os_aux%bc_bot_v              = 0.0_wp
    p_os_aux%bc_bot_vn             = 0.0_wp
    p_os_aux%bc_bot_w              = 0.0_wp
    p_os_aux%bc_top_w              = 0.0_wp
    p_os_aux%bc_bot_tracer         = 0.0_wp
    p_os_aux%bc_top_tracer         = 0.0_wp
   
    rl_start = 1
    rl_end = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,1)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,&
                       & i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        p_os_aux%bc_top_veloc_cc(jc,jb)%x = 0.0_wp
        p_os_aux%bc_bot_veloc_cc(jc,jb)%x = 0.0_wp
      END DO
   END DO
    !CALL message(TRIM(routine),'construction of hydrostatic oceans auxiliary state finished')

  END SUBROUTINE construct_hydro_ocean_aux

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of hydrostatic ocean prognostic state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2006).
  !!
    SUBROUTINE destruct_hydro_ocean_prog(p_os_prog)

    TYPE(t_hydro_ocean_prog),INTENT(INOUT)      :: p_os_prog

    ! local variables

    INTEGER                                   :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_prog'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine),' start to destruct hydrostatic ocean prognostic state')


    ! deallocate prognostic part of state vector
    !
    ! surface height
    DEALLOCATE(p_os_prog%h, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish('mo_ocean_state:destruct_hydro_ocean_prog', 'deallocation for h failed')
    END IF

    !normal velocity
    DEALLOCATE(p_os_prog%vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine),'deallocation for normal velocity failed')
    END IF

    DEALLOCATE(p_os_prog%tracer, STAT=ist)
    IF (ist/=SUCCESS) THEN
        CALL finish(TRIM(routine), 'deallocation for tracer failed')
    END IF

  END SUBROUTINE destruct_hydro_ocean_prog

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of diagnostic hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
  SUBROUTINE destruct_hydro_ocean_diag(p_os_diag)

    TYPE(t_hydro_ocean_diag), INTENT(INOUT) :: p_os_diag

    ! local variables

    INTEGER :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_diag'

!-------------------------------------------------------------------------

    !CALL message(TRIM(routine),' start to destruct hydrostatic ocean diagnostic state')

    ! density
    DEALLOCATE(p_os_diag%rho, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for density failed')
    END IF

    ! reconstructed velocity
    DEALLOCATE(p_os_diag%u, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine),'deallocation for u failed')
    END IF

    DEALLOCATE(p_os_diag%v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for v failed')
    END IF

    DEALLOCATE(p_os_diag%p_vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for p_vn failed')
    END IF
    DEALLOCATE(p_os_diag%ptp_vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for ptp_vn failed')
    END IF

    ! predicted normal velocity
    DEALLOCATE(p_os_diag%vn_pred, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for predicted vn failed')
    END IF
    DEALLOCATE(p_os_diag%vn_impl_vert_diff, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for vn_impl_vert_diff failed')
    END IF

    DEALLOCATE(p_os_diag%w, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for vertical velocity failed')
    END IF
    DEALLOCATE(p_os_diag%w_old, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for old vertical velocity failed')
    END IF

    DEALLOCATE(p_os_diag%w_e, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for vertical velocity at edges failed')
    END IF

    DEALLOCATE(p_os_diag%w_prev, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for vertical veloc. (cells, prev. timestep) failed')
    END IF

    ! tangential velocity
    DEALLOCATE(p_os_diag%vt, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine),'deallocation for tangential velocity failed')
    END IF

    DEALLOCATE(p_os_diag%h_e, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for surface height at edges failed')
    END IF

    DEALLOCATE(p_os_diag%thick_c, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for fluid column at cells failed')
    END IF

    DEALLOCATE(p_os_diag%thick_e, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine), 'deallocation for fluid column at edges failed')
    END IF

    ! vorticity
    DEALLOCATE(p_os_diag%vort, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for vorticity failed')
    END IF

    ! vorticity interpolated to edges
    DEALLOCATE(p_os_diag%vort_e, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine), 'deallocation for vorticity at edges failed')
    END IF

    ! kinetic energy
    DEALLOCATE(p_os_diag%kin, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for kinetic energy failed')
    END IF
    ! gradient
    DEALLOCATE(p_os_diag%grad, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for gradient term failed')
    END IF
    ! divergence
    DEALLOCATE(p_os_diag%div, STAT=ist)
    IF (ist/=SUCCESS) THEN
        CALL finish(TRIM(routine),'deallocation for divergence failed')
    END IF

    ! hydrostatic pressure
    DEALLOCATE(p_os_diag%press_hyd, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine),'deallocation for press_hyd failed')
    END IF

    DEALLOCATE(p_os_diag%press_grad, STAT=ist)
    IF (ist/=SUCCESS) THEN
       CALL finish(TRIM(routine),'deallocation for press_grad failed')
    END IF

!     DEALLOCATE(p_os_diag%flux_mass, STAT=ist)
!     IF (ist/=SUCCESS) THEN
!        CALL finish(TRIM(routine),'deallocation for mass flux failed')
!     END IF
!     DEALLOCATE(p_os_diag%flux_tracer, STAT=ist)
!     IF (ist/=SUCCESS) THEN
!        CALL finish(TRIM(routine),'deallocation for tracer flux failed')
!     END IF

    DEALLOCATE(p_os_diag%veloc_adv_horz, STAT=ist)
    IF (ist/=SUCCESS)THEN
        CALL finish(TRIM(routine),'deallocation for horizontal veloc advection term failed')
    ENDIF

    DEALLOCATE(p_os_diag%veloc_adv_vert, STAT=ist)
    IF (ist/=SUCCESS)THEN
        CALL finish(TRIM(routine),'deallocation for vertical veloc advection term failed')
    ENDIF

    DEALLOCATE(p_os_diag%laplacian_horz, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation for horizontal diff. of hor. veloc. term failed')
    ENDIF

    DEALLOCATE(p_os_diag%laplacian_vert, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation for vertical diff. of hor. veloc. term failed')
    ENDIF

  END SUBROUTINE destruct_hydro_ocean_diag

  !-------------------------------------------------------------------------
  !>
  !!               Deallocation of auxilliary hydrostatic ocean state.
  !
  !! @par Revision History
  !! Developed  by  Peter Korn, MPI-M (2005).
  !!
  SUBROUTINE destruct_hydro_ocean_aux(p_os_aux)

    TYPE(t_hydro_ocean_aux), INTENT(INOUT)      :: p_os_aux

    ! local variables

    INTEGER                                   :: ist

    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:destruct_hydro_ocean_aux'

!-------------------------------------------------------------------------


    !CALL message(TRIM(routine), 'start to destruct auxiliary hydro ocean state')


!   #slo# - no dynamic allocation
!   CALL deallocate_temp_memory_patch_2d(p_os_aux%p_mem_2d)

    ! deallocation for Adam-Bashford time stepping
    DEALLOCATE(p_os_aux%g_n, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of expl. term at step n in AB scheme failed')
    ENDIF

    DEALLOCATE(p_os_aux%g_nm1, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of expl. term at step n-1 in AB scheme failed')
    ENDIF

    DEALLOCATE(p_os_aux%g_nimd, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'dealloc. of expl. term at intermed. step in AB scheme failed')
    ENDIF

    DEALLOCATE(p_os_aux%g_n_c_h, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of horz expl. tracer term at&
                                & step n in AB scheme failed')
    ENDIF
    DEALLOCATE(p_os_aux%g_nm1_c_h, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of horz expl. tracer term at &
                                   &step n-1 in AB scheme failed')
    ENDIF
    DEALLOCATE(p_os_aux%g_nimd_c_h, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of horz expl. tracer term at&
                              & intermed. step in AB scheme failed')
    ENDIF
    DEALLOCATE(p_os_aux%g_n_c_v, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of vert expl. tracer term at&
                                & step n in AB scheme failed')
    ENDIF
    DEALLOCATE(p_os_aux%g_nm1_c_v, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of vert expl. tracer term at &
                                   &step n-1 in AB scheme failed')
    ENDIF
    DEALLOCATE(p_os_aux%g_nimd_c_v, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of vert expl. tracer term at&
                              & intermed. step in AB scheme failed')
    ENDIF



    DEALLOCATE(p_os_aux%p_rhs_sfc_eq, STAT=ist)
    IF (ist/=SUCCESS)THEN
      CALL finish(TRIM(routine),'deallocation of RHS of surface eq. failed')
    ENDIF

    DEALLOCATE(p_os_aux%bc_top_u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond u failed')
    END IF
    DEALLOCATE(p_os_aux%bc_top_v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond v failed')
    END IF
    DEALLOCATE(p_os_aux%bc_top_vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond vn failed')
    END IF

    DEALLOCATE(p_os_aux%bc_top_veloc_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of top boundary cond cc failed')
    END IF

    DEALLOCATE(p_os_aux%bc_bot_u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of bottom boundary cond u failed')
    END IF
    DEALLOCATE(p_os_aux%bc_bot_v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of bottom boundary cond v failed')
    END IF
    DEALLOCATE(p_os_aux%bc_bot_vn, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of bottom boundary cond vn failed')
    END IF
   DEALLOCATE(p_os_aux%bc_bot_veloc_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation of bot boundary cond cc failed')
    END IF

  ! DEALLOCATE(p_os_aux%p_div_flux_horiz_act, STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !    CALL finish(TRIM(routine),'deallocation of div flux horiz term failed')
  ! ENDIF
  ! DEALLOCATE(p_os_aux%p_div_flux_vert_act, STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'deallocation of div flux vert term failed')
  ! ENDIF
  ! DEALLOCATE(p_os_aux%p_div_flux_horiz_prev, STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !   CALL finish(TRIM(routine),'deallocation of prev div flux horiz term failed')
  ! ENDIF
  ! DEALLOCATE(p_os_aux%p_div_flux_vert_prev, STAT=ist)
  ! IF (ist/=SUCCESS)THEN
  !      CALL finish(TRIM(routine),'deallocation of prev div flux vert term failed')
  ! ENDIF

    !semi-implicit coefficient  at edges
!    IF ( i_oce_stepping == semi_impl_ab ) THEN
!       DEALLOCATE(p_os_aux%gamma_e, STAT=ist)
!       IF (ist/=SUCCESS) THEN
!         CALL finish(TRIM(routine),'deallocation for gamma at edges failed')
!       END IF
!       DEALLOCATE(p_os_aux%sigma_nml_e, STAT=ist)
!       IF (ist/=SUCCESS) THEN
!         CALL finish(TRIM(routine),'deallocation for sigma_nml_e at edges failed')
!       END IF
!       DEALLOCATE(p_os_aux%sigma_tgt_e, STAT=ist)
!       IF (ist/=SUCCESS) THEN
!         CALL finish(TRIM(routine),'deallocation for sigma_tgt_e at edges failed')
!       END IF
!    END IF


  END SUBROUTINE destruct_hydro_ocean_aux
!-------------------------------------------------------------------------
!
!
!>
!! Sbr set boundary values for velocity field.
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!
!
  SUBROUTINE  set_lateral_boundary_values( p_patch, vn )

    TYPE(t_patch), INTENT(in) :: p_patch
    REAL(wp)                  :: vn(:,:,:)

    ! local variables
    INTEGER :: jb,je, jk
    INTEGER :: i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
    INTEGER :: rl_start_e, rl_end_e
    INTEGER :: slev,elev 
    CHARACTER(len=max_char_length), PARAMETER :: &
      &      routine = 'mo_oce_state:set_lateral_boundary_values'
!---------------------------------------------------------------

! blocking
i_startblk_e = p_patch%edges%start_blk(1,1)
i_endblk_e   = p_patch%edges%end_blk(min_rledge,1)
rl_start_e   = 1
rl_end_e     = min_rledge
slev         = 1
elev         = n_zlev

DO jb=i_startblk_e, i_endblk_e
  CALL get_indices_e(p_patch, jb, i_startblk_e, i_endblk_e, &
                     i_startidx_e, i_endidx_e, rl_start_e, rl_end_e)
    DO jk = slev, elev
      DO je= i_startidx_e, i_endidx_e
!       IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) /= sea ) THEN
! #slo# 2011-02-28 - correction (identical)
        IF ( p_patch%patch_oce%lsm_oce_e(je,jk,jb) >= boundary ) THEN
          vn(je,jk,jb) = 0.0_wp
        ENDIF
      END DO
  END DO
END DO



  END SUBROUTINE  set_lateral_boundary_values


END MODULE mo_oce_state
