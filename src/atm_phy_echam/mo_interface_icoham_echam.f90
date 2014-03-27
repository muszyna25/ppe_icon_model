!>
!! @brief Interface between ICOHAM dynamics+transport and ECHAM physics
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_icoham_echam

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_impl_constants,    ONLY: min_rlcell_int, min_rledge_int, min_rlcell, io3_amip
  USE mo_datetime,          ONLY: t_datetime, print_datetime, add_time
  USE mo_math_constants,    ONLY: pi
  USE mo_model_domain,      ONLY: t_patch
  USE mo_master_nml,        ONLY: lrestart
  USE mo_echam_phy_config,  ONLY: phy_config => echam_phy_config
  USE mo_echam_phy_memory,  ONLY: prm_field, prm_tend
  USE mo_icoham_dyn_types,  ONLY: t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_intp_data_strc,    ONLY: t_int_state
  USE mo_intp_rbf,          ONLY: rbf_vec_interpol_cell
  USE mo_intp,              ONLY: edges2cells_scalar
                                
  USE mo_parallel_config,   ONLY: nproma, use_icon_comm, p_test_run
  
  USE mo_icon_comm_lib,     ONLY: new_icon_comm_variable, delete_icon_comm_variable, &
     & icon_comm_var_is_ready, icon_comm_sync, icon_comm_sync_all, is_ready, until_sync
  
  USE mo_run_config,        ONLY: nlev, ltimer, ntracer
  USE mo_radiation_config,  ONLY: ighg, izenith, irad_o3, irad_aero, isolrad, &
                                  tsi, tsi_radt, ssi_radt
  USE mo_loopindices,       ONLY: get_indices_c, get_indices_e
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_e, grf_bdywidth_c
  USE mo_eta_coord_diag,    ONLY: half_level_pressure, full_level_pressure
  USE mo_echam_phy_main,    ONLY: echam_phy_main
  USE mo_sync,              ONLY: SYNC_C, SYNC_E, sync_patch_array, sync_patch_array_mult
  USE mo_timer,             ONLY: timer_start, timer_stop, timers_level,  &
    & timer_dyn2phy, timer_phy2dyn, timer_echam_phy, timer_coupling, &
    & timer_echam_sync_temp , timer_echam_sync_tracers
  USE mo_time_interpolation,ONLY: time_weights_limm
  USE mo_time_interpolation_weights, ONLY: wi_limm, wi_limm_radt      
  USE mo_coupling_config,    ONLY: is_coupled_run
#ifdef YAC_coupling
  USE finterface_description ONLY: yac_fput, yac_fget, yac_fget_nbr_fields, yac_fget_field_ids
#else
  USE mo_icon_cpl_exchg,     ONLY: ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_def_field, ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids
  USE mo_icon_cpl_restart,   ONLY: icon_cpl_write_restart
#endif
  USE mo_icoham_sfc_indices, ONLY: iwtr, iice
  USE mo_o3,                 ONLY: read_amip_o3
  USE mo_aero_kinne,         ONLY: read_aero_kinne
  USE mo_aero_stenchikov,    ONLY: read_aero_stenchikov
  USE mo_amip_bc,            ONLY: read_amip_bc, amip_time_weights, amip_time_interpolation, &
    &                              get_current_amip_bc_year
  USE mo_greenhouse_gases,   ONLY: read_ghg_bc, ghg_time_interpolation, ghg_file_read
  USE mo_solar_irradiance,     ONLY: read_ssi_bc, ssi_time_interpolation

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_icoham_echam

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_interface_icoham_echam'

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

  SUBROUTINE interface_icoham_echam( datetime,             &! in
    &                                pdtime, psteplen, jg, &! in
    &                                p_patch,              &! in
    &                                p_int_state,          &! in
    &                                dyn_prog_old,         &! in
    &                                dyn_diag_old,         &! in
    &                                dyn_prog_new,         &! in
    &                                dyn_tend             ) ! inout

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
    REAL(wp) :: ztime_radtran  !< time instance (in radian) at which radiative transfer is computed
    REAL(wp) :: zvn1, zvn2
    REAL(wp), POINTER :: zdudt(:,:,:), zdvdt(:,:,:)

    INTEGER :: jb,jbs   !< block index and its staring value
    INTEGER :: jcs,jce  !< start/end column index within each block
!     INTEGER :: nblks    !< number of blocks for which parameterisations are computed
    INTEGER :: jc, jk   !< column index, vertical level index
    INTEGER :: jcn,jbn  !< column and block indices of a neighbour cell

    LOGICAL               :: write_coupler_restart
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
    CHARACTER(*), PARAMETER :: method_name = "interface_icoham_echam"

!++jsr
!temporary local variables
!    INTEGER:: inm1, inm2
!    REAL(wp):: wgt1, wgt2
!--jsr

    !-------------------------------------------------------------------------
    IF (ltimer) CALL timer_start(timer_dyn2phy)

    ! Inquire current grid level and the total number of grid cells
    i_nchdom  = MAX(1,p_patch%n_childdom)
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
!     nblks = p_patch%nblks_c

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



    ! LL The physics runs only on the owned cells
    !  but the following rbf_vec_interpol_cell may use the halos(?)
    CALL sync_patch_array( SYNC_E, p_patch, dyn_prog_old%vn )

    CALL rbf_vec_interpol_cell( dyn_prog_old%vn,      &! in
      &                         p_patch, p_int_state, &! in
      &                         prm_field(jg)%u,      &! out
      &                         prm_field(jg)%v,      &! out
      &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in    

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
!     nblks = p_patch%nblks_c
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      ! Pressure at time step "new" (i.e., n+1)

      CALL half_level_pressure( dyn_prog_new%pres_sfc(:,jb),     nproma, jce, &! in
                              & prm_field(jg)%presi_new(:,:,jb)               )! out

      CALL full_level_pressure( prm_field(jg)%presi_new(:,:,jb), nproma, jce, &! in
                              & prm_field(jg)%presm_new(:,:,jb)               )! out
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !--------------------------------
    ! transfer tendencies


    ! LL The physics runs only on the owned cells
    !    but the following rbf_vec_interpol_cell may use the halos(?)
    CALL sync_patch_array( SYNC_E, p_patch, dyn_tend%vn )

    CALL rbf_vec_interpol_cell( dyn_tend%vn,          &! in
      &                         p_patch, p_int_state, &! in
      &                         prm_tend(jg)%u,       &! out
      &                         prm_tend(jg)%v,       &! out
      &   opt_rlstart=rl_start, opt_rlend=rl_end     ) ! in

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      prm_tend(jg)% temp(jcs:jce,:,jb)   = dyn_tend%   temp(jcs:jce,:,jb)
      prm_tend(jg)%    q(jcs:jce,:,jb,:) = dyn_tend% tracer(jcs:jce,:,jb,:)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (ltimer)  THEN
      CALL timer_stop (timer_dyn2phy)
      CALL timer_start(timer_echam_phy)
    END IF

    !-------------------------------------------------------------------------
    ! Prepare some global parameters or parameter arrays
    !-------------------------------------------------------------------------
    ! Check whether it is time to do radiative transfer calculation;
    ! Set the time instance at which radiative transfer and radiative
    ! heating will be computed. The current implementation can not handle
    ! seasonal cycle!

    IF (phy_config%lrad) THEN

      ! different to echam the time for which the radiation has to be calculated has
      ! not to be adjusted by pdtime because the middle of the time interval is hit
      ! perfectly. Only the first one is slightly of centered! Do not correct.

      datetime_radtran = datetime                ! copy current date and time
      dsec = 0.5_wp*phy_config%dt_rad            ! [s] time increment for zenith angle computation
      CALL add_time(dsec,0,0,0,datetime_radtran) ! add time increment to get date and
      !                                          ! time information for the zenith angle comp.

      ltrig_rad   = ( l1st_phy_call.AND.(.NOT.lrestart)            ).OR. &
                    ( MOD(NINT(datetime%daysec),NINT(phy_config%dt_rad)) == 0 )
      ! IF (ltrig_rad) THEN
      !   CALL message('mo_interface_icoham_echam:echam_phy_main','Radiative transfer called at:')
      !   CALL print_datetime(datetime)
      !   CALL message('mo_interface_icoham_echam:echam_phy_main','Radiative transfer computed for:')
      !   CALL print_datetime(datetime_radtran)
      ! END IF

      l1st_phy_call = .FALSE.

      ztime_radtran = 2._wp*pi * datetime_radtran%daytim

      !TODO: Luis Kornblueh
      ! add interpolation of greenhouse gases here, only if radiation is going to be calculated
      IF (ltrig_rad .AND. (ighg > 0)) THEN
        CALL ghg_time_interpolation(datetime_radtran)
      END IF

    ELSE
      ltrig_rad = .FALSE.
      ztime_radtran = 0._wp
    END IF

    ! Read and interpolate SST for AMIP simulations; solar irradiation, prescribed ozone, and
    ! aerosol concentrations.
    !LK:
    !TODO: pass timestep as argument
    !COMMENT: lsmask == slm and is not slf!
    IF (phy_config%lamip) THEN
      IF (datetime%year /= get_current_amip_bc_year()) THEN
        CALL read_amip_bc(datetime%year, p_patch)
      END IF
      CALL amip_time_weights(datetime)
      CALL amip_time_interpolation(prm_field(jg)%seaice(:,:), &
!           &                       prm_field(jg)%tsfc_tile(:,:,:), &
           &                       prm_field(jg)%tsurfw(:,:), &
           &                       prm_field(jg)%siced(:,:), &
           &                       prm_field(jg)%lsmask(:,:))
      prm_field(jg)%tsfc_tile(:,:,iwtr) = prm_field(jg)%tsurfw(:,:)
! The ice model should be able to handle different thickness classes, but for AMIP we only use one
      prm_field(jg)%conc(:,1,:) = prm_field(jg)%seaice(:,:)
      prm_field(jg)%hi(:,1,:)   = prm_field(jg)%siced(:,:)
    END IF
! Calculate interpolation weights for linear interpolation
! of monthly means onto the actual integration time step
    CALL time_weights_limm(datetime, wi_limm)
    IF (ltrig_rad) THEN
! Calculate interpolation weights for linear interpolation
! of monthly means onto the radiation time steps
      CALL time_weights_limm(datetime_radtran,wi_limm_radt)   
    END IF
    IF (isolrad==1) THEN
      CALL read_ssi_bc(datetime%year,.FALSE.)
      CALL ssi_time_interpolation(wi_limm,.FALSE.,tsi)
      IF (ltrig_rad) THEN
         CALL read_ssi_bc(datetime_radtran%year,.TRUE.)
         CALL ssi_time_interpolation(wi_limm_radt,.TRUE.,tsi_radt,ssi_radt)
      END IF
    END IF !isolrad
    IF (ltrig_rad .AND. irad_o3 == io3_amip) THEN
      CALL read_amip_o3(datetime%year, p_patch)
    END IF
    IF (ltrig_rad .AND. irad_aero == 13) THEN
      CALL read_aero_kinne(datetime%year, p_patch)
    END IF
    IF (ltrig_rad .AND. irad_aero == 15) THEN
      CALL read_aero_kinne(datetime%year, p_patch)
      CALL read_aero_stenchikov(datetime%year, p_patch)
    END IF

!    WRITE(0,*)'radiation=',ltrig_rad, dt_rad
!    WRITE(0,*)' vor PYHSC rad fluxes sw sfc',  MAXVAL(prm_field(jg)% swflxsfc_avg(:,:))
!    WRITE(0,*)' vor PYHSC rad fluxes lw sfc', MINVAL(prm_field(jg)% lwflxsfc_avg(:,:))

    !-------------------------------------------------------------------------
    ! For each block, call "echam_phy_main" to compute various parameterised processes
    !-------------------------------------------------------------------------
!     jbs   = p_patch%cells%start_blk(grf_bdywidth_c+1,1)
!     nblks = p_patch%nblks_c
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce),  ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_c(p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      ! Prepare some block (slice) specific parameters

      ! CALL xyz(...)

      ! Call parameterisations.
      ! Like in ECHAM, the subroutine *echam_phy_main* has direct access to the memory
      ! buffers prm_field and prm_tend. In addition it can also directly access
      ! the grid/patch information on which the computations are performed.
      ! Thus the argument list contains only
      ! - the domain index jg;
      ! - the block index jb, and the corresponding staring and ending column
      !   indices jcs and jce;
      ! - a few other (global) constants (nproma, pdtime, psteplen, etc).

      CALL echam_phy_main( jg,jb,jcs,jce,nproma,     &
        &                  pdtime,psteplen,          &
        &                  ltrig_rad,ztime_radtran )

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


!    WRITE(0,*)' nach PYHSC rad fluxes sw sfc', MAXVAL( prm_field(jg)% swflxsfc_avg(:,:))
!    WRITE(0,*)' nach PYHSC rad fluxes lw sfc', MINVAL( prm_field(jg)% lwflxsfc_avg(:,:))


    IF (ltimer)  THEN
      CALL timer_stop (timer_echam_phy)
      CALL timer_start(timer_phy2dyn)
    END IF

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
       
       ALLOCATE(buffer(nproma*p_patch%nblks_c,5))
       buffer(:,:) = 0.0_wp

       !
       !  see drivers/mo_atmo_model.f90:
       !
       !   field_id(1) represents "TAUX"   wind stress component
       !   field_id(2) represents "TAUY"   wind stress component
       !   field_id(3) represents "SFWFLX" surface fresh water flux
       !   field_id(4) represents "SFTEMP" surface temperature
       !   field_id(5) represents "THFLX"  total heat flux
       !   field_id(6) represents "ICEATM" ice temperatures and melt potential
       !
       !   field_id(7) represents "SST"    sea surface temperature
       !   field_id(9) represents "OCEANU" u component of ocean surface current
       !   field_id(9) represents "OCEANV" v component of ocean surface current
       !   field_id(10)represents "ICEOCE" ice thickness, concentration and temperatures
       !
       !
#ifdef YAC_Coupling
       CALL yac_fget_nbr_fields ( nbr_fields )
       ALLOCATE(field_id(nbr_fields))
       CALL yac_fget_field_ids ( nbr_fields, field_id )
#else
       CALL ICON_cpl_get_nbr_fields ( nbr_fields )
       ALLOCATE(field_id(nbr_fields))
       CALL ICON_cpl_get_field_ids ( nbr_fields, field_id )
#endif
       !
       !
       field_shape(1) = 1
       field_shape(2) = nbr_hor_points
       field_shape(3) = 1

       !
       ! Send fields away
       ! ----------------
       !
       write_coupler_restart = .FALSE.
       !
       ! TAUX
       !
       buffer(:,:) = 0.0_wp
       buffer(:,1) = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iwtr), (/ nbr_points /) )
       buffer(:,2) = RESHAPE ( prm_field(jg)%u_stress_tile(:,:,iice), (/ nbr_points /) )

#ifdef YAC_coupling
       CALL yac_fput ( field_id(1), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
       field_shape(3) = 2
       CALL ICON_cpl_put ( field_id(1), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
       IF ( info == 2 ) write_coupler_restart = .TRUE.
       !
       ! TAUY
       !
       buffer(:,1) = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iwtr), (/ nbr_points /) )
       buffer(:,2) = RESHAPE ( prm_field(jg)%v_stress_tile(:,:,iice), (/ nbr_points /) )
#ifdef YAC_coupling
       CALL yac_fput ( field_id(2), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
       CALL ICON_cpl_put ( field_id(2), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
       IF ( info == 2 ) write_coupler_restart = .TRUE.
       !
       ! SFWFLX Note: the evap_tile should be properly updated and added
       !
!       write(0,*)  prm_field(jg)%rsfl(:,:)
!       write(0,*)  prm_field(jg)%rsfc(:,:)
!       write(0,*)  prm_field(jg)%ssfl(:,:)
!       write(0,*)  prm_field(jg)%ssfc(:,:)
!       write(0,*)  prm_field(jg)%evap_tile(:,:,iwtr)

        buffer(:,1) = RESHAPE ( prm_field(jg)%rsfl(:,:), (/ nbr_points /) ) + &
             &        RESHAPE ( prm_field(jg)%rsfc(:,:), (/ nbr_points /) ) + &
             &        RESHAPE ( prm_field(jg)%ssfl(:,:), (/ nbr_points /) ) + &
             &        RESHAPE ( prm_field(jg)%ssfc(:,:), (/ nbr_points /) )
        buffer(:,2) = RESHAPE ( prm_field(jg)%evap_tile(:,:,iwtr), (/ nbr_points /) )
 
#ifdef YAC_coupling
       CALL yac_fput ( field_id(3), nbr_hor_points, 2, 1, 1, buffer, ierror )
#else
       CALL ICON_cpl_put ( field_id(3), field_shape, buffer(1:nbr_hor_points,1:2), info, ierror )
#endif
       IF ( info == 2 ) write_coupler_restart = .TRUE.
       !
       ! SFTEMP
       !
       buffer(:,1) =  RESHAPE ( prm_field(jg)%temp(:,nlev,:), (/ nbr_points /) )
#ifdef YAC_coupling
       CALL yac_fput ( field_id(4), nbr_hor_points, 1, 1, 1, buffer, ierror )
#else
       field_shape(3) = 1
       CALL ICON_cpl_put ( field_id(4), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
       IF ( info == 2 ) write_coupler_restart = .TRUE.
       !
       ! THFLX, total heat flux
       !
       buffer(:,1) =  RESHAPE ( prm_field(jg)%swflxsfc_tile(:,:,iwtr), (/ nbr_points /) ) !net shortwave flux for ocean
       buffer(:,2) =  RESHAPE ( prm_field(jg)%lwflxsfc_tile(:,:,iwtr), (/ nbr_points /) ) !net longwave flux
       buffer(:,3) =  RESHAPE ( prm_field(jg)%shflx_tile(:,:,iwtr),    (/ nbr_points /) ) !sensible heat flux
       buffer(:,4) =  RESHAPE ( prm_field(jg)%lhflx_tile(:,:,iwtr),    (/ nbr_points /) ) !latent heat flux for ocean
#ifdef YAC_coupling
       CALL yac_fput ( field_id(5), nbr_hor_points, 4, 1, 1, buffer, ierror )
#else
       field_shape(3) = 4
       CALL ICON_cpl_put ( field_id(5), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
       !
       ! ICEATM, Ice state determined by atmosphere
       !
       buffer(:,1) =  RESHAPE ( prm_field(jg)%Qtop(:,1,:), (/ nbr_points /) ) !Melt-potential for ice - top
       buffer(:,2) =  RESHAPE ( prm_field(jg)%Qbot(:,1,:), (/ nbr_points /) ) !Melt-potential for ice - bottom
       buffer(:,3) =  RESHAPE ( prm_field(jg)%T1  (:,1,:), (/ nbr_points /) ) !Temperature of upper ice layer
       buffer(:,4) =  RESHAPE ( prm_field(jg)%T2  (:,1,:), (/ nbr_points /) ) !Temperature of lower ice layer
#ifdef YAC_coupling
       CALL yac_fput ( field_id(6), nbr_hor_points, 4, 1, 1, buffer, ierror )
#else
       field_shape(3) = 4
       CALL ICON_cpl_put ( field_id(6), field_shape, buffer(1:nbr_hor_points,1:4), info, ierror )
#endif
       IF ( info == 2 ) write_coupler_restart = .TRUE.
#ifdef YAC_coupling
  TODO
#else
       IF ( write_coupler_restart ) CALL icon_cpl_write_restart ( 6, field_id(1:6), ierror )
#endif
       !
       ! Receive fields, only assign values if something was received ( info > 0 )
       ! -------------------------------------------------------------------------
       !
       ! SST
       !
#ifdef YAC_coupling
       CALL yac_fget ( field_id(7), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
       field_shape(3) = 1
       CALL ICON_cpl_get ( field_id(7), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
         prm_field(jg)%tsfc_tile(:,:,iwtr) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%tsfc_tile(:,:,iwtr))
       END IF
       !
       ! OCEANU
       !
#ifdef YAC_coupling
       CALL yac_fget ( field_id(8), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
       CALL ICON_cpl_get ( field_id(8), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
         prm_field(jg)%ocu(:,:) = RESHAPE (buffer(:,1), (/ nproma,  p_patch%nblks_c /) )
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
       END IF
       !
       ! OCEANV
       !
#ifdef YAC_coupling
       CALL yac_fget ( field_id(9), nbr_hor_points, 1, 1, 1, buffer, info, ierror )
#else
       CALL ICON_cpl_get ( field_id(9), field_shape, buffer(1:nbr_hor_points,1:1), info, ierror )
#endif
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:1) = 0.0_wp
         prm_field(jg)%ocv(:,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
       END IF
       !
       ! ICEOCE
       !
#ifdef YAC_coupling
       CALL yac_fget ( field_id(7), nbr_hor_points, 4, 1, 1, buffer, info, ierror )
#else
       field_shape(3) = 5
       CALL ICON_cpl_get ( field_id(10), field_shape, buffer(1:nbr_hor_points,1:5), info, ierror )
#endif
       IF ( info > 0 ) THEN
         buffer(nbr_hor_points+1:nbr_points,1:4) = 0.0_wp
         prm_field(jg)%hi  (:,1,:) = RESHAPE (buffer(:,1), (/ nproma, p_patch%nblks_c /) )
         prm_field(jg)%hs  (:,1,:) = RESHAPE (buffer(:,2), (/ nproma, p_patch%nblks_c /) )
         prm_field(jg)%conc(:,1,:) = RESHAPE (buffer(:,3), (/ nproma, p_patch%nblks_c /) )
         prm_field(jg)%T1  (:,1,:) = RESHAPE (buffer(:,4), (/ nproma, p_patch%nblks_c /) )
         prm_field(jg)%T2  (:,1,:) = RESHAPE (buffer(:,5), (/ nproma, p_patch%nblks_c /) )
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hi  (:,1,:))
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hs  (:,1,:))
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%seaice(:,:))
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%T1  (:,1,:))
         CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%T2  (:,1,:))
         prm_field(jg)%seaice(:,:) = prm_field(jg)%conc(:,1,:)
       END IF

       DEALLOCATE(buffer)
       DEALLOCATE(field_id)
       IF (ltimer) CALL timer_stop(timer_coupling)

    END IF
    
    !-------------------------------------------------------------------------
    ! Physics to dynamics: remap tendencies to the dynamics grid
    !-------------------------------------------------------------------------
    ! Currently this includes a simple copying of the temperature and tracer
    ! tendencies to the dynamics data structure, and convert the u- and v-wind
    ! tendencies to the normal wind tendency.
    ! Once a physics grid of different resolution is intruduced,
    ! conservative re-mapping will be called here.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c( p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      dyn_tend%   temp(jcs:jce,:,jb)   = prm_tend(jg)% temp(jcs:jce,:,jb)
      dyn_tend% tracer(jcs:jce,:,jb,:) = prm_tend(jg)%    q(jcs:jce,:,jb,:)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

     IF (use_icon_comm) THEN
       temp_comm = new_icon_comm_variable(dyn_tend%temp, p_patch%sync_cells_not_in_domain,  &
         & status=is_ready, scope=until_sync, name="echam dyn_tend temp")
       tracers_comm = new_icon_comm_variable(dyn_tend%tracer, p_patch%sync_cells_not_in_domain, &
         & status=is_ready, scope=until_sync, name="echam dyn_tend tracer")
     ELSE
!     IF (timers_level > 5) CALL timer_start(timer_echam_sync_temp)
       CALL sync_patch_array( SYNC_C, p_patch, dyn_tend%temp )
!      IF (timers_level > 5) CALL timer_stop(timer_echam_sync_temp)
!      IF (timers_level > 5) CALL timer_start(timer_echam_sync_tracers)
       CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer, f4din=dyn_tend% tracer)
!      IF (timers_level > 5) CALL timer_stop(timer_echam_sync_tracers)
     END IF
         
    any_uv_tend = phy_config%lconv.OR.phy_config%lvdiff.OR. &
                 & phy_config%lgw_hines .OR.phy_config%lssodrag

    IF (any_uv_tend) THEN

       ALLOCATE(zdudt(nproma,nlev,p_patch%nblks_c), &
         &      zdvdt(nproma,nlev,p_patch%nblks_c), &
         &      stat=return_status)
       IF (return_status > 0) THEN
         CALL finish (method_name, 'ALLOCATE(zdudt,zdvdt)')
       END IF
       IF (p_test_run) THEN
         zdudt(:,:,:) = 0.0_wp
         zdvdt(:,:,:) = 0.0_wp
       END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jcs,jce) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(p_patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        zdudt(jcs:jce,:,jb) =   prm_tend(jg)% u_phy(jcs:jce,:,jb)
        zdvdt(jcs:jce,:,jb) =   prm_tend(jg)% v_phy(jcs:jce,:,jb)
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      ! Now derive the physics-induced normal wind tendency, and add it to the
      ! total tendency.
      IF (use_icon_comm) THEN
        CALL icon_comm_sync(zdudt, zdvdt, p_patch%sync_cells_not_in_domain)
      ELSE
        CALL sync_patch_array_mult(SYNC_C, p_patch, 2, zdudt, zdvdt)
      END IF
      
      jbs   = p_patch%edges%start_blk(grf_bdywidth_e+1,1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,jcs,jce,jcn,jbn,zvn1,zvn2) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = jbs,p_patch%nblks_e
        CALL get_indices_e( p_patch, jb,jbs,p_patch%nblks_e, &
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
          END DO ! jc
        END DO ! jk
      END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      DEALLOCATE(zdudt, zdvdt)
 
  END IF !any_uv_tend

   IF (use_icon_comm) THEN
     CALL icon_comm_sync_all()
   END IF

  IF (ltimer) CALL timer_stop(timer_phy2dyn)
  !--------------------------
  END SUBROUTINE interface_icoham_echam

END MODULE mo_interface_icoham_echam
