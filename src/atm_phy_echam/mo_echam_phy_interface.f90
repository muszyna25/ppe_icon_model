!>
!! @brief Interface between ICOHAM dynamics and ECHAM physics
!!
!! @author Kristina Froehlich (DWD)
!! @author Marco Giorgetta (MPI-M)
!! @author Hui Wan (MPI-M)
!!
!! @par Revision History
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
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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
MODULE mo_echam_phy_interface

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlcell
  USE mo_datetime,          ONLY: t_datetime, print_datetime, add_time
  USE mo_math_constants,    ONLY: pi
  USE mo_model_domain,      ONLY: t_patch
  USE mo_master_nml,        ONLY: lrestart
  USE mo_echam_phy_config,  ONLY: phy_config => echam_phy_config
  USE mo_echam_phy_memory,  ONLY: prm_field, prm_tend, mean_charlen
  USE mo_icoham_dyn_types,  ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_interpolation,     ONLY: t_int_state, rbf_vec_interpol_cell, & 
                                & edges2cells_scalar
  USE mo_exception,        ONLY: finish
                                
  USE mo_parallel_config,   ONLY: nproma, use_icon_comm
  
!   USE mo_icon_comm_lib,     ONLY: icon_comm_sync, icon_comm_sync_all, icon_comm_var_is_ready, &
!     & on_cells, comm_dyn_tend_temp, comm_dyn_tend_tracers
   USE mo_icon_comm_lib,     ONLY: new_icon_comm_variable, delete_icon_comm_variable, &
     & icon_comm_var_is_ready, icon_comm_sync, icon_comm_sync_all, on_cells, is_ready, &
     & until_sync
  
  USE mo_run_config,        ONLY: nlev, ltimer, ntracer
  USE mo_radiation_config,  ONLY: dt_rad,izenith
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_e, grf_bdywidth_c
  USE mo_eta_coord_diag,    ONLY: half_level_pressure, full_level_pressure
  USE mo_echam_phy_main,    ONLY: physc
  USE mo_sync,              ONLY: SYNC_C, SYNC_E, sync_patch_array, sync_patch_array_mult
  USE mo_timer,             ONLY: timer_start, timer_stop, timers_level,  &
    & timer_dyn2phy, timer_phy2dyn, timer_echam_phy, timer_coupling, &
    & timer_echam_sync_temp , timer_echam_sync_tracers
                                
  USE mo_coupling_config,    ONLY: is_coupled_run
  USE mo_icon_cpl_exchg,    ONLY: ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_def_field, ONLY:ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids

  USE mo_icoham_sfc_indices,ONLY: iwtr, iice

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_interface

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_interface'

CONTAINS
  !>
  !! SUBROUTINE echam_physics -- the Interface between ICON dynamics and
  !! ECHAM physics
  !!
  !! This subroutine is called in the time loop of the ICOHAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! <li> tendency of the prognostic varibles induced by adiabatic dynamics;
  !! <li> time step;
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients.
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables caused by
  !! the parameterisations.
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE echam_phy_interface( datetime,             &! in
    &                             pdtime, psteplen, jg, &! in
    &                             p_patch,              &! in
    &                             p_int_state,          &! in
    &                             dyn_prog_old,         &! in
    &                             dyn_diag_old,         &! in
    &                             dyn_prog_new,         &! in
    &                             dyn_tend             ) ! inout

    ! Arguments
 
    TYPE(t_datetime),      INTENT(IN)    :: datetime
    REAL(wp),              INTENT(IN)    :: pdtime        !< time step
    REAL(wp),              INTENT(IN)    :: psteplen      !< 2*time step in case of leapfrog
    INTEGER,               INTENT(IN)    :: jg            !< grid level/domain index
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch
    TYPE(t_int_state),TARGET,INTENT(IN)  :: p_int_state

    TYPE(t_datetime)                     :: datetime_radtran

    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: dyn_prog_old
    TYPE(t_hydro_atm_diag),INTENT(IN)    :: dyn_diag_old
    TYPE(t_hydro_atm_prog),INTENT(IN)    :: dyn_prog_new

    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: dyn_tend

    ! Local variables

    LOGICAL  :: l1st_phy_call = .TRUE.
    LOGICAL  :: any_uv_tend, ltrig_rad
    REAL(wp) :: dsec
    REAL(wp) :: ztime_radtran, ztime_radheat  !< time instance (in radians) at which
                                              !< radiative transfer/heating is computed
    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    INTEGER :: jb,jbs   !< block index and its staring value
    INTEGER :: jcs,jce  !< start/end column index within each block
!     INTEGER :: nblks    !< number of blocks for which parameterisations are computed
    INTEGER :: jc, jk   !< column index, vertical level index
    INTEGER :: jcn,jbn  !< column and block indices of a neighbour cell

    INTEGER               :: nbr_fields
    INTEGER               :: nbr_hor_points ! = inner and halo points
    INTEGER               :: nbr_points     ! = nproma * nblks
    INTEGER               :: field_shape(3)
    INTEGER, ALLOCATABLE  :: field_id(:)
    REAL(wp), ALLOCATABLE :: buffer(:,:)

    INTEGER               :: info, ierror !< return values form cpl_put/get calls

    !--------------
    INTEGER:: i_nchdom       !< number of child patches
    INTEGER:: rl_start, rl_end, i_startblk, i_endblk
    !--------------
    INTEGER:: temp_comm, tracers_comm  ! communicators
    INTEGER:: return_status
    CHARACTER(*), PARAMETER :: method_name = "echam_phy_interface"
    
    !-------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_dyn2phy)

    ! Inquire current grid level and the total number of grid cells
    i_nchdom  = MAX(1,p_patch%n_childdom)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
!     nblks = p_patch%nblks_int_c

    !-------------------------------------------------------------------------
    ! Dynamics to physics: remap dynamics variables to physics grid
    !-------------------------------------------------------------------------
    ! Currently this includes
    !  - reconstructing of u- and v-wind and their tendencies at cell centers
    !  - copying scalar fields from the dynamics state and from the tendency
    !    state to the physics state and physics tendencies, repsectively.
    !  - computing pressure values at the "new" time step
    ! Once a physics grid of different resolution is intruduced, conservative
    ! re-mapping will be called here.



    SELECT CASE (p_patch%cell_type)
    CASE (3)
!     LL The physics runs only on the owned cells
!        but the following rbf_vec_interpol_cell may use the halos(?)
!      CALL sync_patch_array( SYNC_E, p_patch, dyn_prog_old%vn )
      CALL rbf_vec_interpol_cell( dyn_prog_old%vn,      &! in
        &                         p_patch, p_int_state, &! in
        &                         prm_field(jg)%u,      &! out
        &                         prm_field(jg)%v,      &! out
        &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in
    CASE (6)
      CALL sync_patch_array( SYNC_E, p_patch, dyn_prog_old%vn )
      CALL edges2cells_scalar(dyn_prog_old%vn,p_patch, &
        &                     p_int_state%hex_east ,prm_field(jg)%u)
      CALL edges2cells_scalar(dyn_prog_old%vn,p_patch, &
        &                     p_int_state%hex_north,prm_field(jg)%v)
    END SELECT

!     LL The physics runs only on the owned cells
!        thus no knowledge of the halo u, v is needed.
!     CALL sync_patch_array( SYNC_C, p_patch, prm_field(jg)%u )
!     CALL sync_patch_array( SYNC_C, p_patch, prm_field(jg)%v )
    

!$OMP PARALLEL WORKSHARE
    prm_field(jg)%         q(:,:,:,:) = dyn_prog_old%     tracer(:,:,:,:)
    prm_field(jg)%      temp(:,:,:)   = dyn_prog_old%       temp(:,:,:)

    prm_field(jg)% presi_old(:,:,:)   = dyn_diag_old%    pres_ic(:,:,:)
    prm_field(jg)% presm_old(:,:,:)   = dyn_diag_old%    pres_mc(:,:,:)

    prm_field(jg)%        qx(:,:,:)   = dyn_diag_old%         qx(:,:,:)
    prm_field(jg)%        tv(:,:,:)   = dyn_diag_old%      tempv(:,:,:)
    prm_field(jg)%      geom(:,:,:)   = dyn_diag_old%     geo_mc(:,:,:)
    prm_field(jg)%      geoi(:,:,:)   = dyn_diag_old%     geo_ic(:,:,:)
    prm_field(jg)%     omega(:,:,:)   = dyn_diag_old%   wpres_mc(:,:,:)
    prm_field(jg)%       vor(:,:,:)   = dyn_diag_old% rel_vort_c(:,:,:)
!$OMP END PARALLEL WORKSHARE

    !---------------------------------
    ! Additional diagnostic variables

!     jbs   = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
!     nblks = p_patch%nblks_int_c
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      ! Pressure at time step "new" (i.e., n+1)

      CALL half_level_pressure( dyn_prog_new%pres_sfc(:,jb),     nproma, jce, &! in
                              & prm_field(jg)%presi_new(:,:,jb)               )! out

      CALL full_level_pressure( prm_field(jg)%presi_new(:,:,jb), nproma, jce, &! in
                              & prm_field(jg)%presm_new(:,:,jb)               )! out
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    !--------------------------------
    ! transfer tendencies


    SELECT CASE (p_patch%cell_type)
    CASE (3)
!     LL The physics runs only on the owned cells
!        but the following rbf_vec_interpol_cell may use the halos(?)
!     CALL sync_patch_array( SYNC_E, p_patch, dyn_tend%vn )
      CALL rbf_vec_interpol_cell( dyn_tend%vn,          &! in
        &                         p_patch, p_int_state, &! in
        &                         prm_tend(jg)%u,       &! out
        &                         prm_tend(jg)%v,       &! out
        &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in
    CASE (6)
      CALL sync_patch_array( SYNC_E, p_patch, dyn_tend%vn )
      CALL edges2cells_scalar(dyn_tend%vn,p_patch, &
        &                     p_int_state%hex_east ,prm_tend(jg)%u)
      CALL edges2cells_scalar(dyn_tend%vn,p_patch, &
        &                     p_int_state%hex_north,prm_tend(jg)%v)
    END SELECT

!   LL The physics runs only on the owned cells
!       thus no knowledge of the halo u, v is needed.
!    CALL sync_patch_array( SYNC_C, p_patch, prm_tend(jg)%u )
!    CALL sync_patch_array( SYNC_C, p_patch, prm_tend(jg)%v )

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      prm_tend(jg)% temp(jcs:jce,:,jb)   = dyn_tend%   temp(jcs:jce,:,jb)
      prm_tend(jg)%    q(jcs:jce,:,jb,:) = dyn_tend% tracer(jcs:jce,:,jb,:)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF (ltimer)  THEN
      CALL timer_stop (timer_dyn2phy)
      CALL timer_start(timer_echam_phy)
    ENDIF

    !-------------------------------------------------------------------------
    ! Prepare some global parameters or parameter arrays
    !-------------------------------------------------------------------------
    ! Check whether it is time to do radiative transfer calculation;
    ! Set the time instance at which radiative transfer and radiative
    ! heating will be computed. The current implementation can not handle
    ! seasonal cycle!

    IF (phy_config%lrad) THEN

      datetime_radtran = datetime                ! copy current date and time
      dsec = 0.5_wp*(dt_rad-pdtime)              ! [s] time increment for zenith angle computation
      CALL add_time(dsec,0,0,0,datetime_radtran) ! add time increment to get date and
      !                                          ! time information for the zenith angle comp.

      ltrig_rad   = ( l1st_phy_call.AND.(.NOT.lrestart)            ).OR. &
                    ( MOD(NINT(datetime%daysec),NINT(dt_rad)) == 0 )
      l1st_phy_call = .FALSE.

      ztime_radheat = 2._wp*pi * datetime%daytim 
      ztime_radtran = 2._wp*pi * datetime_radtran%daytim

!      IF (ltrig_rad) THEN
!        CALL message('mo_echam_phy_interface:physc','Radiative transfer called at:')
!        CALL print_datetime(datetime)
!        CALL message('mo_echam_phy_interface:physc','Radiative transfer computed for:')
!        CALL print_datetime(datetime_radtran)
!      ENDIF
    ELSE
      ltrig_rad = .FALSE.
      ztime_radheat = 0._wp
      ztime_radtran = 0._wp
    ENDIF

    ! Read and interpolate SST for AMIP simulations; prescribed ozone and
    ! aerosol concentrations.
!    WRITE(0,*)'radiation=',ltrig_rad, dt_rad
!    WRITE(0,*)' vor PYHSC rad fluxes sw sfc',  MAXVAL(prm_field(jg)% swflxsfc_avg(:,:))
!    WRITE(0,*)' vor PYHSC rad fluxes lw sfc', MINVAL(prm_field(jg)% lwflxsfc_avg(:,:))

    !-------------------------------------------------------------------------
    ! For each block, call "physc" to compute various parameterised processes
    !-------------------------------------------------------------------------
!     jbs   = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
!     nblks = p_patch%nblks_int_c
!$OMP PARALLEL
#ifdef __SX__
!$OMP DO PRIVATE(jb,jcs,jce), SCHEDULE(guided)
#else 
!$OMP DO PRIVATE(jb,jcs,jce), SCHEDULE(runtime)
#endif
    DO jb = i_startblk,i_endblk

      CALL get_indices_c(p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      ! Prepare some block (slice) specific parameters

      ! CALL xyz(...)

      ! Call parameterisations.
      ! Like in ECHAM, the subroutine *physc* has direct access to the memory
      ! buffers prm_field and prm_tend. In addition it can also directly access
      ! the grid/patch information on which the computations are performed.
      ! Thus the argument list contains only
      ! - the domain index jg;
      ! - the block index jb, and the corresponding staring and ending column
      !   indices jcs and jce;
      ! - a few other (global) constants (nproma, pdtime, psteplen, etc).

      CALL physc( jg,jb,jcs,jce,nproma,pdtime,psteplen,ltrig_rad, &
                & ztime_radtran,ztime_radheat )

    END DO
!$OMP END DO
!$OMP END PARALLEL


!    WRITE(0,*)' nach PYHSC rad fluxes sw sfc', MAXVAL( prm_field(jg)% swflxsfc_avg(:,:))
!    WRITE(0,*)' nach PYHSC rad fluxes lw sfc', MINVAL( prm_field(jg)% lwflxsfc_avg(:,:))


    IF (ltimer)  THEN
      CALL timer_stop (timer_echam_phy)
      CALL timer_start(timer_phy2dyn)
    ENDIF

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information 
    !-------------------------------------------------------------------------

    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_field(jg)% u_stress_tile(:,:,iwtr)  and 
    !    prm_field(jg)% v_stress_tile(:,:,iwtr)  which are the wind stress components;
    !
    ! 2. prm_field(jg)% evap_tile(:,:,iwtr) evaporation rate
    !
    ! 3. prm_field(jg)%rsfl + prm_field(jg)%rsfc + prm_field(jg)%ssfl + prm_field(jg)%ssfc
    !    which gives the precipitation rate;
    !
    ! 4. prm_field(jg)% temp(:,nlev,:)  temperature at the lowest model level, or
    !    prm_field(jg)% temp_2m(:,:)    2-m temperature, not available yet, or
    !    prm_field(jg)% shflx_tile(:,:,iwtr) sensible heat flux
    !
    ! 5  prm_field(jg)% lhflx_tile(:,:,iwtr) latent heat flux
    ! 6. shortwave radiation flux at the surface
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. prm_field(jg)% tsfc_tile(:,:,iwtr)   SST
    ! 2. prm_field(jg)% ocu(:,:) and ocv(:,:) ocean surface current
    ! 
    IF ( is_coupled_run() ) THEN
       IF (ltimer) CALL timer_start(timer_coupling)
    
       nbr_hor_points = p_patch%n_patch_cells
       nbr_points     = nproma * p_patch%nblks_c
       
       ALLOCATE(buffer(nproma*p_patch%nblks_c,4))
       buffer(:,:) = 0.0_wp

       !
       !  see drivers/mo_atmo_model.f90:
       !
       !   field_id(1) represents "TAUX"   wind stress component
       !   field_id(2) represents "TAUY"   wind stress component
       !   field_id(3) represents "SFWFLX" surface fresh water flux
       !   field_id(4) represents "SFTEMP" surface temperature
       !   field_id(5) represents "THFLX"  total heat flux
       !
       !   field_id(6) represents "SST"    sea surface temperature
       !   field_id(7) represents "OCEANU" u component of ocean surface current
       !   field_id(8) represents "OCEANV" v component of ocean surface current
       !
       CALL ICON_cpl_get_nbr_fields ( nbr_fields )
       ALLOCATE(field_id(nbr_fields))
       CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
       !
       !
       field_shape(1) = 1
       field_shape(2) = nbr_hor_points
       field_shape(3) = 1

       !
       ! Send fields away
       ! ----------------
       !
       ! TAUX
       buffer(:,1) = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iwtr), (/ nbr_points /) )
       CALL ICON_cpl_put ( field_id(1), field_shape, buffer(1:nbr_hor_points,1:1), ierror )
       !
       ! TAUY
       buffer(:,1) = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iwtr), (/ nbr_points /) )
       CALL ICON_cpl_put ( field_id(2), field_shape, buffer(1:nbr_hor_points,1:1), ierror )
       !
       ! SFWFLX Note: the evap_tile should be properly updated and added
        buffer(:,1) = RESHAPE ( prm_field(jg)%rsfl(:,:), (/ nbr_points /) ) + &
             &        RESHAPE ( prm_field(jg)%rsfc(:,:), (/ nbr_points /) ) + &
             &        RESHAPE ( prm_field(jg)%ssfl(:,:), (/ nbr_points /) ) + &
             &        RESHAPE ( prm_field(jg)%ssfc(:,:), (/ nbr_points /) )
        buffer(:,2) = RESHAPE ( prm_field(jg)%evap_tile(:,:,iwtr), (/ nbr_points /) )
 
       field_shape(3) = 2
       CALL ICON_cpl_put ( field_id(3), field_shape, buffer(1:nbr_hor_points,1:2), ierror )
       !
       ! SFTEMP
       buffer(:,1) =  RESHAPE ( prm_field(jg)%temp(:,nlev,:), (/ nbr_points /) )
       field_shape(3) = 1
       CALL ICON_cpl_put ( field_id(4), field_shape, buffer(1:nbr_hor_points,1:1), ierror )
       !
       ! THFLX, total heat flux

       buffer(:,1) =  RESHAPE ( prm_field(jg)%swflxsfc  (:,:)     , (/ nbr_points /) ) !net shortwave flux at sfc
       buffer(:,2) =  RESHAPE ( prm_field(jg)%lwflxsfc  (:,:)     , (/ nbr_points /) ) !net longwave flux at sfc
       buffer(:,3) =  RESHAPE ( prm_field(jg)%shflx_tile(:,:,iwtr), (/ nbr_points /) ) !sensible heat flux
       buffer(:,4) =  RESHAPE ( prm_field(jg)%lhflx_tile(:,:,iwtr), (/ nbr_points /) ) !latent heat flux

       field_shape(3) = 4
       CALL ICON_cpl_put ( field_id(5), field_shape, buffer(1:nbr_hor_points,1:4), ierror )
       field_shape(3) = 1

       !
       ! Receive fields, only assign values if something was received ( info > 0 )
       ! -------------------------------------------------------------------------
       !
       ! SST
       CALL ICON_cpl_get ( field_id(6), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
         prm_field(jg)%tsfc_tile(:,:,iwtr) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
         prm_field(jg)%tsfc     (:,:)      = prm_field(jg)%tsfc_tile(:,:,iwtr)
       ENDIF
       !
       ! OCEANU
       CALL ICON_cpl_get ( field_id(7), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
         prm_field(jg)%ocu(:,:) = RESHAPE (buffer(:,1), (/ nproma,  p_patch%nblks_c /) )
       ENDIF
       !
       ! OCEANV
       CALL ICON_cpl_get ( field_id(8), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
         prm_field(jg)%ocv(:,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
       ENDIF

       CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%tsfc_tile(:,:,iwtr))
       CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%tsfc     (:,:))

       CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
       CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))

       DEALLOCATE(buffer)
       DEALLOCATE(field_id)
       IF (ltimer) CALL timer_stop(timer_coupling)

    ENDIF
    
    !-------------------------------------------------------------------------
    ! Physics to dynamics: remap tendencies to the dynamics grid
    !-------------------------------------------------------------------------
    ! Currently this includes a simple copying of the temperature and tracer
    ! tendencies to the dynamics data structure, and convert the u- and v-wind
    ! tendencies to the normal wind tendency.
    ! Once a physics grid of different resolution is intruduced,
    ! conservative re-mapping will be called here.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      dyn_tend%   temp(jcs:jce,:,jb)   = prm_tend(jg)% temp(jcs:jce,:,jb)
      dyn_tend% tracer(jcs:jce,:,jb,:) = prm_tend(jg)%    q(jcs:jce,:,jb,:)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

     IF (use_icon_comm) THEN
       temp_comm = new_icon_comm_variable(dyn_tend%temp, on_cells, p_patch, &
         & status=is_ready, scope=until_sync, name="echam dyn_tend temp")
       tracers_comm = new_icon_comm_variable(dyn_tend%tracer, on_cells, p_patch, &
         & status=is_ready, scope=until_sync, name="echam dyn_tend tracer")
     ELSE
!     IF (timers_level > 5) CALL timer_start(timer_echam_sync_temp)
       CALL sync_patch_array( SYNC_C, p_patch, dyn_tend%temp )
!      IF (timers_level > 5) CALL timer_stop(timer_echam_sync_temp)
!      IF (timers_level > 5) CALL timer_start(timer_echam_sync_tracers)
       CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=dyn_tend% tracer, &
         &                        lpart4d=.TRUE.)
!      IF (timers_level > 5) CALL timer_stop(timer_echam_sync_tracers)
     ENDIF
         
    any_uv_tend = phy_config%lconv.OR.phy_config%lvdiff.OR. &
                 & phy_config%lgw_hines !.OR.phy_config%lssodrag

    IF (any_uv_tend) THEN

       ALLOCATE(zdudt(nproma,nlev,p_patch%nblks_c), &
         & zdvdt (nproma,nlev,p_patch%nblks_c) ,stat=return_status)
       IF (return_status > 0) &
         CALL finish (method_name, 'ALLOCATE(zdudt,zdvdt)')

      ! Accumulate wind tendencies contributed by various parameterized processes.

!       jbs   = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
!       nblks = p_patch%nblks_int_c
!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        zdudt(jcs:jce,:,jb) =   prm_tend(jg)% u_cnv(jcs:jce,:,jb) &
                            & + prm_tend(jg)% u_vdf(jcs:jce,:,jb) &
                            & + prm_tend(jg)% u_gwh(jcs:jce,:,jb)
        zdvdt(jcs:jce,:,jb) =   prm_tend(jg)% v_cnv(jcs:jce,:,jb) &
                            & + prm_tend(jg)% v_vdf(jcs:jce,:,jb) &
                            & + prm_tend(jg)% v_gwh(jcs:jce,:,jb)
      ENDDO
!$OMP END PARALLEL DO

      ! Now derive the physics-induced normal wind tendency, and add it to the
      ! total tendency.
      IF (use_icon_comm) THEN
        CALL icon_comm_sync(zdudt, zdvdt, on_cells, p_patch)
      ELSE
        CALL sync_patch_array_mult(SYNC_C, p_patch, 2, zdudt, zdvdt)
      ENDIF
      
      jbs   = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
!       nblks = p_patch%nblks_int_e
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce,jcn,jbn,zvn1,zvn2)
      DO jb = jbs,p_patch%nblks_int_e
        CALL get_indices_e( p_patch, jb,jbs,p_patch%nblks_int_e, &
          & jcs,jce, grf_bdywidth_e+1)
        DO jk = 1,nlev
          DO jc = jcs,jce

            jcn  =   p_patch%edges%cell_idx(jc,jb,1)
            jbn  =   p_patch%edges%cell_blk(jc,jb,1)
            zvn1 =   zdudt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,1)%v1 &
                 & + zdvdt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,1)%v2

            jcn  =   p_patch%edges%cell_idx(jc,jb,2)
            jbn  =   p_patch%edges%cell_blk(jc,jb,2)
            zvn2 =   zdudt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,2)%v1 &
                 & + zdvdt(jcn,jk,jbn)*p_patch%edges%primal_normal_cell(jc,jb,2)%v2

            dyn_tend%vn(jc,jk,jb) =   dyn_tend%vn(jc,jk,jb)             &
                                  & + p_int_state%c_lin_e(jc,1,jb)*zvn1 &
                                  & + p_int_state%c_lin_e(jc,2,jb)*zvn2
          ENDDO !column loop
        ENDDO !vertical level loop
      ENDDO !block loop
!$OMP END DO
!$OMP END PARALLEL
      DEALLOCATE(zdudt, zdvdt)
 
  END IF !any_uv_tend

   IF (use_icon_comm) THEN
     CALL icon_comm_sync_all()
   ENDIF

  IF (ltimer) CALL timer_stop(timer_phy2dyn)
  !--------------------------
  END SUBROUTINE echam_phy_interface

END MODULE mo_echam_phy_interface
