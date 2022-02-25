!>
!! Atmospheric diagnostics around air parcel dynamics
!!
!! This module contains subroutines for the computation of: 
!! * Square of Brunt-Vaisala frequency
!! * Square of general air parcel oscillation frequency
!!
!! @par Revision History
!! Initial revision by Sebastian Borchert, DWD (2020-05-29)
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
!----------------------------
#include "omp_definitions.inc"
!----------------------------
!
MODULE mo_diag_atmo_air_parcel

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_impl_constants,         ONLY: MAX_CHAR_LENGTH, SUCCESS, &
    &                                  start_prog_cells, end_prog_cells
  USE mo_physical_constants,     ONLY: grav, rd, rdv, alv, cpd, alvdcp, rcpd, &
    &                                  earth_angular_velocity
  USE mo_upatmo_impl_const,      ONLY: idamtr
  USE mo_model_domain,           ONLY: t_patch
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_parallel_config,        ONLY: nproma
  USE mo_run_config,             ONLY: iqv, iqc
  USE mo_grid_config,            ONLY: grid_sphere_radius, n_dom
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_nh_diagnose_pres_temp,  ONLY: diag_temp, diag_pres
  USE mo_upatmo_config,          ONLY: upatmo_config
  USE mo_nh_vert_interp_les,     ONLY: intrpl_full2half_inblk
  USE mo_util_phys,              ONLY: exner_from_pres, theta_from_temp_and_exner
  USE mo_fortran_tools,          ONLY: init
  USE mo_util_string,            ONLY: int2string
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_opt_diag_atmo_bvf2, &
    &                                  timer_opt_diag_atmo_parcelfreq2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sqr_of_Brunt_Vaisala_freq
  PUBLIC :: sqr_of_parcel_freq

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_diag_atmo_air_parcel'

  !-----------------------------------------------------
  !               For horizontal loops
  !-----------------------------------------------------

  INTEGER, PARAMETER :: start_bdry_lev1_cells = start_prog_cells - 1
  INTEGER, PARAMETER :: end_halo_lev2_cells   = end_prog_cells   - 2

  !-----------------------------------------------------
  !               Support type for module
  !-----------------------------------------------------

  TYPE :: t_spt4airParcel
    LOGICAL :: initializedSqrOfParcelFreq = .FALSE.
  END TYPE t_spt4airParcel

  !-----------------------------------------------------
  !  Support type for subroutine: sqr_of_parcel_freq
  !-----------------------------------------------------

  TYPE :: t_spt4sqrOfParcelFreq
    INTEGER,  ALLOCATABLE :: trglNbhdBlk(:,:,:)    ! (nproma,nstencil_max,nblks_c) 
                                                   ! blocks of triangles in neighborhood
    INTEGER,  ALLOCATABLE :: trglNbhdIdx(:,:,:)    ! (nproma,nstencil_max,nblks_c) 
                                                   ! indices of triangles in neighborhood
    ! Please see NEC user guide why we need to introduce 
    ! each component separately and cannot work with indices:
    REAL(wp), ALLOCATABLE :: X(:,:,:)              ! (nproma,nlev,nblks_c) Cartesian x-coordinate
                                                   ! of grid cell center                           [m] 
    REAL(wp), ALLOCATABLE :: Y(:,:,:)              ! (nproma,nlev,nblks_c) Cartesian y-coordinate
                                                   ! of grid cell center                           [m] 
    REAL(wp), ALLOCATABLE :: Z(:,:,:)              ! (nproma,nlev,nblks_c) Cartesian z-coordinate
                                                   ! of grid cell center                           [m] 
    !
    REAL(wp), ALLOCATABLE :: GradGradPhi_XX(:,:,:) ! (nproma,nlev,nblks_c) Cartesian xx-component 
                                                   ! of the symmetric Hessian of the geopotential  [s-2]
    REAL(wp), ALLOCATABLE :: GradGradPhi_YY(:,:,:) ! (nproma,nlev,nblks_c) Cartesian yy-component 
                                                   ! of the symmetric Hessian of the geopotential  [s-2]
    REAL(wp), ALLOCATABLE :: GradGradPhi_ZZ(:,:,:) ! (nproma,nlev,nblks_c) Cartesian zz-component 
                                                   ! of the symmetric Hessian of the geopotential  [s-2]
    REAL(wp), ALLOCATABLE :: GradGradPhi_XY(:,:,:) ! (nproma,nlev,nblks_c) Cartesian xy-component 
                                                   ! of the symmetric Hessian of the geopotential  [s-2]
    REAL(wp), ALLOCATABLE :: GradGradPhi_XZ(:,:,:) ! (nproma,nlev,nblks_c) Cartesian xz-component 
                                                   ! of the symmetric Hessian of the geopotential  [s-2]
    REAL(wp), ALLOCATABLE :: GradGradPhi_YZ(:,:,:) ! (nproma,nlev,nblks_c) Cartesian yz-component 
                                                   ! of the symmetric Hessian of the geopotential  [s-2]
    !
    INTEGER :: nstencil_max                       ! maximum number of triangles in neighborhood
    LOGICAL :: initialized = .FALSE.
  END TYPE t_spt4sqrOfParcelFreq

  TYPE(t_spt4airParcel)                            :: spt4airParcel
  TYPE(t_spt4sqrOfParcelFreq), TARGET, ALLOCATABLE :: spt4sqrOfParcelFreq(:) ! (n_dom)

CONTAINS

  !>
  !! Square of Brunt-Vaisala frequency
  !!
  !! Compute square of Brunt-Vaisala frequency after one of the following versions:
  !! * standard version
  !! * hydrostatic version
  !! * version after Durran & Klemp (1982)
  !!
  !! @par Revision History
  !! This is a modified copy of:
  !! * /src/atm_phy_les/mo_les_utilities: brunt_vaisala_freq
  !! by Anurag Dipankar, MPI-M 
  !!
  SUBROUTINE sqr_of_Brunt_Vaisala_freq( p_patch,             & !inout
    &                                   p_metrics,           & !in
    &                                   p_prog,              & !in
    &                                   p_prog_rcf,          & !in
    &                                   p_diag,              & !in
    &                                   bvf2,                & !inout
    &                                   bvf2_mode,           & !in
    &                                   opt_startlev,        & !optin
    &                                   opt_endlev,          & !optin
    &                                   opt_update_temp,     & !optin
    &                                   opt_condensate_list, & !optin
    &                                   opt_kstart_moist,    & !optin
    &                                   opt_timer,           & !optin
    &                                   opt_verbose          ) !optin

    ! in/out variables

    TYPE(t_patch),             INTENT(INOUT) :: p_patch                !< horizontal grid
    TYPE(t_nh_metrics),        INTENT(IN)    :: p_metrics              !< vertical grid
    TYPE(t_nh_prog),           INTENT(IN)    :: p_prog                 !< prognostic variables
    TYPE(t_nh_prog),           INTENT(IN)    :: p_prog_rcf             !< tracers (concentration of water phases)
    TYPE(t_nh_diag),           INTENT(INOUT) :: p_diag                 !< diagnostic variables
    REAL(wp),                  INTENT(INOUT) :: bvf2(:,:,:)            !< square of Brunt-Vaisala frequency 
                                                                       !< (nproma, nlev, nblks_c) [s-2]
    INTEGER,                   INTENT(IN)    :: bvf2_mode              !< computation mode
    INTEGER,         OPTIONAL, INTENT(IN)    :: opt_startlev           !< start level for vertical loops
    INTEGER,         OPTIONAL, INTENT(IN)    :: opt_endlev             !< end level for vertical loops
    LOGICAL,         OPTIONAL, INTENT(IN)    :: opt_update_temp        !< update diagnostic variable temperature
    INTEGER, TARGET, OPTIONAL, INTENT(IN)    :: opt_condensate_list(:) !< list with identifiers 
                                                                       !< of condensed water phases
    INTEGER,         OPTIONAL, INTENT(IN)    :: opt_kstart_moist       !< above this level there are 
                                                                       !< no condensed water phases 
    LOGICAL,         OPTIONAL, INTENT(IN)    :: opt_timer              !< switch for timer monitoring
    LOGICAL,         OPTIONAL, INTENT(IN)    :: opt_verbose            !< switch for message output

    ! local variables

    INTEGER,  POINTER :: condensate_list(:)

    REAL(wp), ALLOCATABLE :: theta_v_ifc(:,:), & 
      &                      theta_ifc(:,:),   &
      &                      temp_v_ifc(:,:),  &
      &                      theta(:,:),       &
      &                      qv_ifc(:,:),      &
      &                      qc_ifc(:,:)

    REAL(wp) :: grav_deepatmo(p_patch%nlev)

    REAL(wp) :: factor
    INTEGER  :: nlev, nlevp1
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: i_startlev, i_endlev
    INTEGER  :: jg, jb, jc, jk
    INTEGER  :: kstart_moist
    INTEGER  :: istat
    LOGICAL  :: mode_standard, mode_hydrostatic, mode_dk1982
    LOGICAL  :: verbose, timer, update_temp
    LOGICAL  :: present_kstart_moist, present_condensate_list
    CHARACTER(LEN=50) :: mode_info_str
    CHARACTER(LEN=2)  :: dom_str

    REAL(wp), PARAMETER :: alv_o_rd          = alv / rd
    REAL(wp), PARAMETER :: rdv_alv2_o_cpd_rd = rdv * alv**2  / ( cpd * rd ) 

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':sqr_of_Brunt_Vaisala_freq'
    
    !------------------------------------------------

    !-------------------------------------------
    !              Description
    !-------------------------------------------

    ! For undersaturated air, the square of the Brunt-Vaisala frequency N**2 [s-2] is computed according to:
    !
    !           g     dtheta_v
    ! N**2 = ------- ---------,   (1)
    !        theta_v     dz
    !
    ! where g is gravity, theta_v is the virtual potential temperature and z is the geometric height.
    !
    ! For saturated air, N**2 is computed according to the approximate formula (36) in Durran & Klemp (1982): 
    ! "On the effects of moisture on the Brunt-Vaisala frequency", JAS, 39, 2152-2158. 
    !
    !            {        1 + ( L * r_s / R * T )             (   1   dtheta      L    dr_s )   dr_w }
    ! N**2 = g * { ---------------------------------------- * ( ----- ------ + ------- ---- ) - ---- } ,   (2)
    !            { 1 + ( eps * L**2 * r_s / c_p * R * T**2)   ( theta   dz     c_p * T  dz  )    dz  }
    ! 
    ! where T = temperature of dry air, theta = potential temperature of dry air, 
    ! R = gas constant of dry air, eps = R / R_v, R_v = gas constant of water vapor, 
    ! L = latent heat of vaporization, c_p = heat capacity of dry air at const. pressure, 
    ! r_s = saturation mixing ratio, r_w = r_s + r_l at (super)saturation, r_l = liquid water mixing ratio,
    ! this formula holds in case of (super)saturated air, effects of ice phase and precipitation 
    ! have not been incorporated.
    !
    ! Some simplifications are applied to formula (2):
    !
    ! * Water phase concentrations are stored in ICON in terms of 
    !   specific concentrations (q_x = mass of water phase x / ( mass of dry air + mass of water phases)). 
    !   From this, water phase mixing ratios (r_x = mass of water phase x / mass of dray air) 
    !   could be computed according to: r_x = q_x / ( 1 - Sum(q_x) ), 
    !   but for efficiency reasons r_x ~ q_x is assumed in the following (which is generally 
    !   regarded as a good approximation).
    !
    ! * The saturation mixing ratio r_s = r_s(T,p) = eps * E / ( p - E ) ~ eps * E / p, 
    !   where E is the saturation vapor pressure, a function of temperature, which can be computed 
    !   using '/src/atm_phy_schemes/mo_satad: sat_pres_water'. 
    !   But this comprises the computationally expensive call of the exponential function, 
    !   so to simplify matters, the following assumptions are made: 
    !   ** q_c (cloud/liquid water) > 0 means that the grid cell volume is (super)saturated.
    !   ** Water phases are always close to thermodynamic equilibrium, 
    !      so if q_c > 0 -> r_s ~ q_s ~ q_v (specific humidity).
    !   ** Ice phase and precipitation are neglected.
    !
    ! Finally, for the square of the hydrostatic Brunt-Vaisala frequency the following formula is used:
    !
    !           g     dtheta_v       {  1  dT_v      R    dp  }                g    { dT_v    g  }
    ! N**2 = ------- --------- = g * { --- ---- - ------- --- }             = --- * { ---- + --- },   (3)
    !        theta_v     dz          { T_v  dz    c_p * p dz  }_hydrostatic   T_v   {  dz    c_p }
    !
    ! where T_v is virtual temperature, p is pressure, and dp/dz = -grav * rho = -grav * p / (R * T_v)
    ! is used in the last step.
    !

    !-------------------------------------------
    !               Computation
    !-------------------------------------------

    !------------
    ! preparation
    !------------

    ! domain index
    jg = p_patch%id

    ! number of vertical grid levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! check presence of optional arguments
    present_kstart_moist    = PRESENT(opt_kstart_moist)
    present_condensate_list = PRESENT(opt_condensate_list) 

    ! type of computation
    mode_standard    = .FALSE.
    mode_hydrostatic = .FALSE.
    mode_dk1982      = .FALSE.
    SELECT CASE(bvf2_mode)
    CASE(1)
      mode_standard    = .TRUE.
      mode_info_str    = "standard mode"
    CASE(2) 
      mode_hydrostatic = .TRUE.
      mode_info_str    = "hydrostatic mode"
    CASE(3)
      mode_dk1982      = .TRUE.
      mode_info_str    = "after Durran & Klemp 1982"
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid computation mode for bvf2')
    END SELECT

    ! please note that we take into account opt_startlev and opt_endlev 
    ! in a poor man's way only in that they are applied where 
    ! the final assignment of the output field bvf2 takes place
    ! (no adjustment of opt_startlev > opt_endlev to opt_startlev <= opt_endlev)
    IF (PRESENT(opt_startlev)) THEN
      i_startlev = MAX(1, opt_startlev)
    ELSE
      i_startlev = 1
    ENDIF
    IF (PRESENT(opt_endlev)) THEN
      i_endlev   = MIN(opt_endlev, SIZE(bvf2, 2), nlev)
    ELSE
      i_endlev   = nlev
    ENDIF

    IF (i_startlev > i_endlev) RETURN

    ! timer monitoring
    IF (PRESENT(opt_timer)) THEN
      timer = opt_timer
    ELSE
      timer = .FALSE.
    ENDIF

    ! transmit messages if required
    IF (PRESENT(opt_verbose)) THEN
      verbose = opt_verbose
    ELSE
      verbose = .FALSE.
    ENDIF

    ! switch on timer monitoring if desired
    IF(timer) CALL timer_start(timer_opt_diag_atmo_bvf2)

    dom_str = TRIM(int2string(jg))

    IF (verbose) CALL message(TRIM(routine), &
      & 'Computation of square of Brunt-Vaisala frequency started ('//TRIM(mode_info_str)//') on domain '//dom_str)

    ! shall we update p_diag%pres and p_diag%tempv?
    IF (PRESENT(opt_update_temp)) THEN
      update_temp = opt_update_temp
    ELSE
      update_temp = .FALSE.
    ENDIF

    ! check for optional index of grid layer above which 
    ! no condensed water phases exist in the model
    IF (present_kstart_moist) THEN
      kstart_moist = MAX(1, MIN(opt_kstart_moist, nlev))
    ELSE
      kstart_moist = 0
    ENDIF

    ! for the temperature update we need the condensate list 
    ! and the start level for the water physics 
    IF (update_temp .AND. .NOT. (present_condensate_list .AND. present_kstart_moist)) &
      & CALL finish(TRIM(routine), 'opt_condensate_list and opt_kstart_moist are required if opt_update_temp = .true.')

    ! account for decrease of gravitational acceleration in upward vertical direction 
    ! if deep-atmosphere mode is switched on: 
    ! grav -> grav * a**2 / r**2, with a = radius of Earth and r = a + z
    ! (please note that this is no more than a poor man's modification for the deep atmosphere, 
    ! if N**2 would be derived in spherical geometry from scratch, 
    ! additional modifications could appear which we neglect)
    DO jk = 1, nlev
      grav_deepatmo(jk) = grav * p_metrics%deepatmo_t1mc(jk,idamtr%t1mc%gradh)**2
    ENDDO

    IF (present_condensate_list) THEN
      condensate_list => opt_condensate_list
    ELSE
      condensate_list => NULL()
    ENDIF

    ! loop over prognostic domain
    rl_start   = start_prog_cells
    rl_end     = end_prog_cells
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !--------------------------------
    ! update temperature if required
    !--------------------------------

    ! only "hydrostatic" and "dk1982" modes make use of the temperature
    IF (update_temp .AND. (mode_hydrostatic .OR. mode_dk1982)) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        CALL diag_temp( pt_prog         = p_prog,          & !in
          &             pt_prog_rcf     = p_prog_rcf,      & !in
          &             condensate_list = condensate_list, & !in
          &             pt_diag         = p_diag,          & !inout
          &             jb              = jb,              & !in
          &             i_startidx      = i_startidx,      & !in
          &             i_endidx        = i_endidx,        & !in
          &             slev            = 1,               & !in
          &             slev_moist      = kstart_moist,    & !in
          &             nlev            = nlev             ) !in   

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL      

    ENDIF  !IF (update_temp .AND. (mode_hydrostatic .OR. mode_dk1982))

    !-------------------------
    ! main computation,
    ! select computation mode
    !-------------------------

    IF (mode_standard) THEN !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      !--------------------------------------------------------
      ! "standard" computation according to Eq. (1) (see above)
      !--------------------------------------------------------

      ALLOCATE( theta_v_ifc(nproma,nlevp1), & 
        &       STAT=istat                  )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of auxiliary fields failed')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, theta_v_ifc) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
                
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      ! interpolate virtual potential temperature from full to half levels
      CALL intrpl_full2half_inblk( i_startidx      = i_startidx,                             & !in
        &                          i_endidx        = i_endidx,                               & !in
        &                          i_startlev_half = 1,                                      & !in
        &                          i_endlev_half   = nlevp1,                                 & !in
        &                          nproma          = nproma,                                 & !in
        &                          nlev            = nlev,                                   & !in
        &                          p_wgtfac        = REAL(p_metrics%wgtfac_c(:,:,jb),   wp), & !in (vp->wp)
        &                          p_wgtfacq       = REAL(p_metrics%wgtfacq_c(:,:,jb),  wp), & !in (vp->wp)
        &                          p_wgtfacq1      = REAL(p_metrics%wgtfacq1_c(:,:,jb), wp), & !in (vp->wp)
        &                          p_var_full      = p_prog%theta_v(:,:,jb),                 & !in
        &                          p_var_half      = theta_v_ifc(:,:)                        ) !inout

      ! compute square of Brunt-Vaisala frequency
      DO jk = i_startlev, i_endlev
        DO jc = i_startidx, i_endidx
          
          ! formula (1) (see above)
          bvf2(jc,jk,jb) = grav_deepatmo(jk) * REAL(p_metrics%inv_ddqz_z_full(jc,jk,jb), wp) * &
            &              (theta_v_ifc(jc,jk) - theta_v_ifc(jc,jk+1)) / p_prog%theta_v(jc,jk,jb)

        ENDDO  !jc
      ENDDO  !jk
      
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      DEALLOCATE( theta_v_ifc, & 
        &         STAT=istat   )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of auxiliary fields failed')

    ELSEIF (mode_hydrostatic) THEN !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      !-----------------------------------------------------
      ! hydrostatic balance according to Eq. (3) (see above)
      !-----------------------------------------------------

      ALLOCATE( temp_v_ifc(nproma,nlevp1), & 
        &       STAT=istat                 )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of auxiliary fields failed')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, temp_v_ifc) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
                
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      ! interpolate virtual temperature from full to half levels
      CALL intrpl_full2half_inblk( i_startidx      = i_startidx,                             & !in
        &                          i_endidx        = i_endidx,                               & !in
        &                          i_startlev_half = 1,                                      & !in
        &                          i_endlev_half   = nlevp1,                                 & !in
        &                          nproma          = nproma,                                 & !in
        &                          nlev            = nlev,                                   & !in
        &                          p_wgtfac        = REAL(p_metrics%wgtfac_c(:,:,jb),   wp), & !in (vp->wp)
        &                          p_wgtfacq       = REAL(p_metrics%wgtfacq_c(:,:,jb),  wp), & !in (vp->wp)
        &                          p_wgtfacq1      = REAL(p_metrics%wgtfacq1_c(:,:,jb), wp), & !in (vp->wp)
        &                          p_var_full      = p_diag%tempv(:,:,jb),                   & !in
        &                          p_var_half      = temp_v_ifc(:,:)                         ) !inout

      ! compute square of Brunt-Vaisala frequency
      DO jk = i_startlev, i_endlev
        DO jc = i_startidx, i_endidx
          
          ! formula (3) (see above)
          bvf2(jc,jk,jb) = grav_deepatmo(jk) * (                             &
            &                REAL(p_metrics%inv_ddqz_z_full(jc,jk,jb), wp) * &
            &                (temp_v_ifc(jc,jk) - temp_v_ifc(jc,jk+1))       &
            &              + rcpd * grav_deepatmo(jk)                        &
            &              ) / p_diag%tempv(jc,jk,jb)

        ENDDO  !jc
      ENDDO  !jk

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      DEALLOCATE( temp_v_ifc, & 
        &         STAT=istat  )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of auxiliary fields failed')

    ELSEIF (mode_dk1982) THEN !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      !------------------------------------------
      ! computation after Durran & Klemp (1982) 
      ! according to Eqs. (1) and (2) (see above)
      !------------------------------------------

      ALLOCATE( theta_v_ifc(nproma,nlevp1), & 
        &       theta_ifc(nproma,nlevp1),   &
        &       qv_ifc(nproma,nlevp1),      &
        &       qc_ifc(nproma,nlevp1),      &
        &       theta(nproma,nlev),         &
        &       STAT=istat                  )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of auxiliary fields failed')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, theta, theta_v_ifc, theta_ifc, &
!$OMP            qv_ifc, qc_ifc, factor) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
                
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      ! compute potential temperature
      DO jk = i_startlev, i_endlev
        theta(i_startidx:i_endidx,jk) = theta_from_temp_and_exner( p_diag%temp(i_startidx:i_endidx,jk,jb), &
          &                                                        p_prog%exner(i_startidx:i_endidx,jk,jb) )
      ENDDO

      ! interpolate virtual potential temperature from full to half levels
      CALL intrpl_full2half_inblk( i_startidx      = i_startidx,                             & !in
        &                          i_endidx        = i_endidx,                               & !in
        &                          i_startlev_half = 1,                                      & !in
        &                          i_endlev_half   = nlevp1,                                 & !in
        &                          nproma          = nproma,                                 & !in
        &                          nlev            = nlev,                                   & !in
        &                          p_wgtfac        = REAL(p_metrics%wgtfac_c(:,:,jb),   wp), & !in (vp->wp)
        &                          p_wgtfacq       = REAL(p_metrics%wgtfacq_c(:,:,jb),  wp), & !in (vp->wp)
        &                          p_wgtfacq1      = REAL(p_metrics%wgtfacq1_c(:,:,jb), wp), & !in (vp->wp)
        &                          p_var_full      = p_prog%theta_v(:,:,jb),                 & !in
        &                          p_var_half      = theta_v_ifc(:,:)                        ) !inout

      ! interpolate potential temperature from full to half levels
      CALL intrpl_full2half_inblk( i_startidx      = i_startidx,                             & !in
        &                          i_endidx        = i_endidx,                               & !in
        &                          i_startlev_half = 1,                                      & !in
        &                          i_endlev_half   = nlevp1,                                 & !in
        &                          nproma          = nproma,                                 & !in
        &                          nlev            = nlev,                                   & !in
        &                          p_wgtfac        = REAL(p_metrics%wgtfac_c(:,:,jb),   wp), & !in (vp->wp)
        &                          p_wgtfacq       = REAL(p_metrics%wgtfacq_c(:,:,jb),  wp), & !in (vp->wp)
        &                          p_wgtfacq1      = REAL(p_metrics%wgtfacq1_c(:,:,jb), wp), & !in (vp->wp)
        &                          p_var_full      = theta(:,:),                             & !in
        &                          p_var_half      = theta_ifc(:,:)                          ) !inout

      ! interpolate specific humidity from full to half levels
      CALL intrpl_full2half_inblk( i_startidx      = i_startidx,                             & !in
        &                          i_endidx        = i_endidx,                               & !in
        &                          i_startlev_half = 1,                                      & !in
        &                          i_endlev_half   = nlevp1,                                 & !in
        &                          nproma          = nproma,                                 & !in
        &                          nlev            = nlev,                                   & !in
        &                          p_wgtfac        = REAL(p_metrics%wgtfac_c(:,:,jb),   wp), & !in (vp->wp)
        &                          p_wgtfacq       = REAL(p_metrics%wgtfacq_c(:,:,jb),  wp), & !in (vp->wp)
        &                          p_wgtfacq1      = REAL(p_metrics%wgtfacq1_c(:,:,jb), wp), & !in (vp->wp)
        &                          p_var_full      = p_prog_rcf%tracer(:,:,jb,iqv),          & !in
        &                          p_var_half      = qv_ifc(:,:)                             ) !inout

      ! interpolate specific cloud water content from full to half levels
      CALL intrpl_full2half_inblk( i_startidx      = i_startidx,                             & !in
        &                          i_endidx        = i_endidx,                               & !in
        &                          i_startlev_half = 1,                                      & !in
        &                          i_endlev_half   = nlevp1,                                 & !in
        &                          nproma          = nproma,                                 & !in
        &                          nlev            = nlev,                                   & !in
        &                          p_wgtfac        = REAL(p_metrics%wgtfac_c(:,:,jb),   wp), & !in (vp->wp)
        &                          p_wgtfacq       = REAL(p_metrics%wgtfacq_c(:,:,jb),  wp), & !in (vp->wp)
        &                          p_wgtfacq1      = REAL(p_metrics%wgtfacq1_c(:,:,jb), wp), & !in (vp->wp)
        &                          p_var_full      = p_prog_rcf%tracer(:,:,jb,iqc),          & !in
        &                          p_var_half      = qc_ifc(:,:)                             ) !inout

      ! compute square of Brunt-Vaisala frequency ...
      ! ... for grid layers without condensed water phases ...
      DO jk = i_startlev, MIN(kstart_moist, i_endlev)
        DO jc = i_startidx, i_endidx      

          ! formula (1) (see above)
          bvf2(jc,jk,jb) = grav_deepatmo(jk) * REAL(p_metrics%inv_ddqz_z_full(jc,jk,jb), wp) * &
            &              (theta_v_ifc(jc,jk) - theta_v_ifc(jc,jk+1)) / p_prog%theta_v(jc,jk,jb)          

        ENDDO  !jc
      ENDDO  !jk
      
      ! ... and for grid layers with condensed water phases
      DO jk = MAX(kstart_moist + 1, i_startlev), i_endlev
        DO jc = i_startidx, i_endidx      

          IF (p_prog_rcf%tracer(jc,jk,jb,iqc) > 0._wp) THEN
            ! probably (super)saturated air:
            ! use formula (2) (see above)
            factor         = ( 1._wp + alv_o_rd * p_prog_rcf%tracer(jc,jk,jb,iqv) / p_diag%temp(jc,jk,jb) ) / &
              &              ( 1._wp + rdv_alv2_o_cpd_rd * p_prog_rcf%tracer(jc,jk,jb,iqv) /                  &
              &                p_diag%temp(jc,jk,jb)**2 )

            bvf2(jc,jk,jb) = grav_deepatmo(jk) * REAL(p_metrics%inv_ddqz_z_full(jc,jk,jb), wp) *    &
              &              (                                                                      &
              &                factor *                                                             &
              &                (                                                                    &
              &                  (theta_ifc(jc,jk) - theta_ifc(jc,jk+1)) / theta(jc,jk)             &
              &                + alvdcp * (qv_ifc(jc,jk) - qv_ifc(jc,jk+1)) / p_diag%temp(jc,jk,jb) &
              &                )                                                                    &
              &              - (qv_ifc(jc,jk) + qc_ifc(jc,jk) - qv_ifc(jc,jk+1) - qc_ifc(jc,jk+1))  &
              &              )
          ELSE
            ! probably undersaturated air:
            ! use formula (1) (see above)
            bvf2(jc,jk,jb) = grav_deepatmo(jk) * REAL(p_metrics%inv_ddqz_z_full(jc,jk,jb), wp) * &
              &              (theta_v_ifc(jc,jk) - theta_v_ifc(jc,jk+1)) / p_prog%theta_v(jc,jk,jb)          
          ENDIF

        ENDDO  !jc
      ENDDO  !jk

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      DEALLOCATE( theta_v_ifc, &
        &         theta_ifc,   & 
        &         qv_ifc,      & 
        &         qc_ifc,      & 
        &         theta,       &
        &         STAT=istat   )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of auxiliary fields failed')

    ENDIF !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    !----------
    ! clean-up
    !----------

    IF (present_condensate_list) condensate_list => NULL()

    IF (verbose) CALL message(TRIM(routine), &
      & 'Computation of square of Brunt-Vaisala frequency finished on domain '//dom_str)

    IF(timer) CALL timer_stop(timer_opt_diag_atmo_bvf2)

  END SUBROUTINE sqr_of_Brunt_Vaisala_freq

  !=====================================================================================================

  !>
  !! Square of general air parcel frequency (real part)
  !!
  !! Compute (real part of) square of general air parcel frequency after:
  !! H. Ertel, J.-J. Jaw, S.-Z. Li (1941) Tensorielle Theorie der Stabilitaet, 
  !! Meteorologische Zeitschrift, 58, 389 - 392.
  !! * standard version
  !! * hydrostatic version
  !!
  !! Please note: if you like to use this subroutine outside the context
  !! of the optional output of diagnostics, you have to implement 
  !! an additional switch in:
  !! /src/shr_horizontal/mo_intp_state: allocate_int_state
  !! that enables the set-up of p_int_state%cell_environ in:
  !! /src/shr_horizontal/mo_intp_state: allocate_int_state.
  !! 
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD (2020-06-08)
  !!
  SUBROUTINE sqr_of_parcel_freq( p_patch,              & !inout
    &                            p_int_state,          & !in
    &                            p_metrics,            & !in
    &                            p_prog,               & !in
    &                            p_diag,               & !inout
    &                            parcelfreq2,          & !inout
    &                            parcelfreq2_mode,     & !in
    &                            opt_startlev,         & !optin
    &                            opt_endlev,           & !optin
    &                            opt_update_pres_temp, & !optin
    &                            opt_p_prog_rcf,       & !optin
    &                            opt_condensate_list,  & !optin
    &                            opt_kstart_moist,     & !optin
    &                            opt_lastcall,         & !optin
    &                            opt_timer,            & !optin
    &                            opt_verbose,          & !optin
    &                            opt_minute            ) !optin

    ! in/out variables

    TYPE(t_patch),      TARGET,   INTENT(INOUT) :: p_patch               !< horizontal grid
    TYPE(t_int_state),  TARGET,   INTENT(IN)    :: p_int_state           !< for horizontal interpolation
    TYPE(t_nh_metrics),           INTENT(IN)    :: p_metrics             !< vertical grid
    TYPE(t_nh_prog),    TARGET,   INTENT(IN)    :: p_prog                !< prognostic variables
    TYPE(t_nh_diag),              INTENT(INOUT) :: p_diag                !< diagnostic variables
    REAL(wp),                     INTENT(INOUT) :: parcelfreq2(:,:,:)    !< square of parcel frequency
                                                                         !< (nproma, nlev, nblks_c) [s-2]
    INTEGER,                      INTENT(IN)    :: parcelfreq2_mode      !< computation mode
    INTEGER,            OPTIONAL, INTENT(IN)    :: opt_startlev          !< start level for vertical loops
    INTEGER,            OPTIONAL, INTENT(IN)    :: opt_endlev            !< end level for vertical loops
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_update_pres_temp  !< update diagnostic variables 
                                                                         !< pressure and temperature
    TYPE(t_nh_prog),    TARGET,   OPTIONAL, INTENT(IN) :: opt_p_prog_rcf !< tracers (concentration of water phases)
    INTEGER,    TARGET, OPTIONAL, INTENT(IN)    :: opt_condensate_list(:)!< list with identifiers 
                                                                         !< for condensed water phases
    INTEGER,            OPTIONAL, INTENT(IN)    :: opt_kstart_moist      !< above this level there are 
                                                                         !< no condensed water phases 
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_lastcall          !< last call?
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_timer             !< switch for timer monitoring
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_verbose           !< switch for message output
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_minute            !< switch for non-essential checks

    ! local variables

    INTEGER,  PARAMETER :: nwork = 10

    TYPE(t_nh_prog), POINTER :: p_prog_rcf

    INTEGER,  POINTER :: condensate_list(:)
    INTEGER,  POINTER :: blk(:,:,:)
    INTEGER,  POINTER :: idx(:,:,:)
    REAL(wp), POINTER :: eR_X(:,:)
    REAL(wp), POINTER :: eR_Y(:,:)
    REAL(wp), POINTER :: eR_Z(:,:)
    REAL(wp), POINTER :: X(:,:,:)
    REAL(wp), POINTER :: Y(:,:,:)
    REAL(wp), POINTER :: Z(:,:,:)
    REAL(wp), POINTER :: GradGradPhi_XX(:,:,:)
    REAL(wp), POINTER :: GradGradPhi_YY(:,:,:)
    REAL(wp), POINTER :: GradGradPhi_ZZ(:,:,:)
    REAL(wp), POINTER :: GradGradPhi_XY(:,:,:)
    REAL(wp), POINTER :: GradGradPhi_XZ(:,:,:)
    REAL(wp), POINTER :: GradGradPhi_YZ(:,:,:)
    REAL(wp), POINTER :: exner(:,:,:)
    REAL(wp), POINTER :: thetav(:,:,:)

    REAL(wp), TARGET, ALLOCATABLE :: exner_hydro(:,:,:)
    REAL(wp), TARGET, ALLOCATABLE :: thetav_hydro(:,:,:)

    REAL(wp) :: GradGradExner_XX(nproma,p_patch%nlev)
    REAL(wp) :: GradGradExner_YY(nproma,p_patch%nlev)
    REAL(wp) :: GradGradExner_ZZ(nproma,p_patch%nlev)
    REAL(wp) :: GradGradExner_XY(nproma,p_patch%nlev)
    REAL(wp) :: GradGradExner_XZ(nproma,p_patch%nlev)
    REAL(wp) :: GradGradExner_YZ(nproma,p_patch%nlev)

    REAL(wp) :: stab_XX, stab_YY, stab_ZZ, stab_XY, stab_XZ, stab_YZ
    REAL(wp) :: cubic_coeff_A, cubic_coeff_B, cubic_coeff_C

    REAL(wp) :: work(nwork)
    REAL(wp) :: comp_mat(3,3,nproma), eval_Re(3,nproma), eval_Im(3,nproma)
    REAL(wp) :: schurvec(1,3)

    INTEGER  :: ierr_hessian(p_patch%nblks_c)
    INTEGER  :: ierr_dhseqr(p_patch%nblks_c)

    INTEGER  :: nlev, nblks_c, nstencil_max
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: i_startlev, i_endlev
    INTEGER  :: jg, jb, jc, jk
    INTEGER  :: kstart_moist
    INTEGER  :: istat
    LOGICAL  :: mode_standard, mode_hydrostatic
    LOGICAL  :: verbose, timer, lastcall
    LOGICAL  :: switch_off_hor_percept, update_pres_temp
    LOGICAL  :: present_kstart_moist, present_condensate_list, present_p_prog_rcf, minute
    CHARACTER(LEN=50) :: mode_info_str
    CHARACTER(LEN=2)  :: dom_str

    INTEGER, PARAMETER :: stencil_type = 13 ! 10 -> edge stencil, 
                                            ! 13 -> vertex stencil
    REAL(wp), PARAMETER :: earth_angular_velocity2_t_4 = 4._wp * earth_angular_velocity**2  ![rad2 s-2]

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':sqr_of_parcel_freq'

    !------------------------------------------------

    !-------------------------------------------
    !              Description
    !-------------------------------------------

    ! 
    ! 1st parcel frequency model:
    !-----------------------------
    ! Here the square of the parcel frequency is the result of a parcel-dynamical stability analysis
    ! described by H. Ertel, J.-J. Jaw, S.-Z. Li (1941) Tensorielle Theorie der Stabilitaet, 
    ! Meteorologische Zeitschrift, 58, 389 - 392. 
    ! Starting point is the momentum budget of an air parcel divided by density rho 
    ! to which a displacement is applied (without specifying the cause of this displacement):
    !
    !     d**2[r]/dt**2 + 2{W} * d[r]/dt + [G]phi + [G]p / rho = 0   | D                          (1)
    !
    ! =>  d**2D[r]/dt**2 + 2{W} * dD[r]/dt + D([G]phi) + D([G]p) / rho - [G]p / rho**2 Drho = 0   (2)
    !
    ! d denotes the Lagrangian (parcel-fixed) derivative, [r] is the position vector, 
    ! {W} = [W] x {1} is the Coriolis tensor, with angular velocity vector [W], identity tensor {1}, 
    ! and cross product x. In addition, * denotes the scalar product, 
    ! [G] is the gradient (Nabla operator), phi is the geopotential, p is the pressure, 
    ! and D stands for the displacement, which should commute with the material derivative d. 
    ! In addition, the displacement is assumed to be so small that 
    !
    !     D(ab) = bDa + aDb                                                                       (3)
    !
    ! is a good approximation. Combined with the quotient rule this has been used, for instance, for:
    !
    !      ( [G]p )   D([G]p)     [G]p
    !     D(------) = ------- - ------- Drho.                                                     (4)
    !      ( rho  )     rho      rho**2
    !
    ! In order to replace the displacements of gradients, we use a spatial Taylor expansion:
    !
    ! D([G]X) = D[r] * [G][G]X = ([G][G]X)**T * D[r] = ([G][G]X) * D[r], with X = phi, p,         (5)
    !
    ! where [G][G]X denotes the Hessian derivative tensor of X. 
    ! This tensor is symmetric, i.e. ([G][G]X)**T = [G][G]X, where (...)**T denotes the transpose of the tensor.
    ! Next, we assume the displacement to be an adiabatic change of state, so that:
    !
    !     Drho / rho = (1 / kappa) / (Dp / p),                                                    (6)
    !
    ! with kappa = cp / cv. Finally, the pressure of the parcel shall adjust to its new environment 
    ! after the displacement, which reads according to its spatial Taylor expansion:
    !
    !     Dp = D[r] * [G]p.                                                                       (7)
    !
    ! The final result for (2) can be written as an operator acting on the displacement D[r]:
    !
    !     ( {1}d**2/dt**2 + 2{W}d/dt + ( [G][G]phi + [G][G]p / rho - ([G]p)([G]p) / (kappa rho p) ) ) * D[r] = 0,  (8)
    !                                  |__________________________________________________________|
    !                                                                |
    !                                                              =: {S}
    !
    ! where Eq. (2) was divided by rho, and ([G]p)([G]p) is the dyadic product of the pressure gradient [G]p with itself.
    ! {S} is called the stability tensor, it is symmetric ({S}**T = {S}).
    ! The operator (...) acting on D[r] in Eq. (8) is a tensor, and for Eq. (8) to hold, the following equation 
    ! has to hold, too:
    ! 
    !     |(...)| D[r] = 0,                                                                       (9) 
    !
    ! where |(...)| denotes the determinant of the tensor operator. After some calculation, we find for the formulation (9) 
    ! of Eq. (8):
    !
    ! ( d**6/dt**6 + ({S}. + 4[W]**2) d**4/dt**4 + (({S}.**2 - {S}..{S}) / 2 + 4[W][W]..{S}) d**2/dt**2 + |{S}| ) D[r] = 0,  (10)
    !                |______________|              |_______________________________________|              |___|
    !                        |                                         |                                    |
    !                      =: a                                      =: b                                 =: c
    !
    ! where {S}. denotes the scalar of {S} (the first invariant of {S}), [W]**2 = [W]*[W], 
    ! {S}..{S} is the double scalar product of {S} with itself (the result of which is a scalar, 
    ! ({S}.**2 - {S}..{S}) / 2 is the second invariant of {S}), 
    ! [W][W]..{S} is the double scalar product of the dyad [W][W] with {S}, 
    ! and |{S}| denotes the determinant of {S} (its third invariant). 
    !
    ! Next, we assume an oscillation for the parcel displacement D[r]:
    !
    !     D[r](t) = <[r]> exp(-i omega t),                                                        (11)
    !
    ! where <[r]> is a time-independent complex vector amplitude, i denotes the imaginary unit, 
    ! and omega is the oscillation frequency. 
    ! With (11) we find for (10):
    !
    !     (omega**2)**3 - a (omega**2)**2 + b (omega**2) - c = 0,                                 (12)
    !
    ! a cubic equation in the square of the oscillation frequency.
    ! In general Eq. (12) has three roots, (omega**2)_1,2,3, of which at least one is real.
    ! Within this model instability would occur, if one of the roots of (omega**2)_1,2,3
    ! has a positive imaginary part. So our goal is to find the one root omega
    ! with the greates imaginary part and output the negative square of that imaginary part.
    !
    ! So much the original formulation. Within ICON, Eq. (1) is replaced by the following reformulation:
    !
    !     d**2[r]/dt**2 + 2{W} * d[r]/dt + [G]phi + cp theta [G]exner = 0   | D                   (13)
    !
    ! =>  ( {1}d**2/dt**2 + 2{W}d/dt + ( [G][G]phi + cp theta [G][G]exner ) ) * D[r] = 0,         (14) 
    !                                  |__________________________________|
    !                                                    |
    !                                                  =: {S}
    !
    ! where we have used that Dtheta = 0 for adiabatic processes, and D([G]exner) = ([G][G]exner) * D[r].
    ! Both formulations, (8) and (14), are equivalent, since:
    !
    !     cp theta [G][G]exner = [ [G][G]p - ([G]p)([G]p) / (kappa p) ] / rho.                    (15)
    !
    ! In principle it would be desirable to account for a hydrostatically balanced background: 
    !
    !     ( [G]phi + cp theta0 [G]exner0 ) * [e_r] = 0                                            (16)
    !
    ! in Eq. (13), where [e_r] denotes the radial unit vector, in order to make the computation more accurate. 
    ! Applying the displacement to the resulting formulation of Eq. (13) is, however, relatively complicated. 
    ! The resulting formulation for {S} would require significantly more computational effort. 
    ! For this reason we use formulation (14).
    !
    ! Measuring the tensorial quantities in the coordinate system of the triangle edges 
    ! would be far too complicated. 
    ! One possible alternative is global Cartesian coordinates with their point of origin in the center of Earth 
    ! and the positive z-axis crossing the North Pole:
    !
    !     x = r cos(lat) cos(lon),  y = r cos(lat) sin(lon),  z = r sin(lat),                     (16)
    !
    ! Measuring {S} in these coordinates and arranging the coefficients in a matrix, 
    ! denoted by ({.})_Cart yields:
    !
    !                  ( S_x,x  S_x,y  S_x,z )   ( S_x,x  S_x,y  S_x,z )
    !     ({S})_Cart = ( S_y,x  S_y,y  S_y,z ) = ( S_x,y  S_y,y  S_y,z ),                         (17)
    !                  ( S_z,x  S_z,y  S_z,z )   ( S_x,z  S_y,z  S_z,z )
    !
    ! In the last step of (17) we used that {S} is symmetric. The six independent coefficients read:
    !
    !             grav r_E**2 ( r**2 - 3 x**2 )            d**2exner 
    !     S_x,x = ----------------------------- + cp theta ---------,                             (18.1)
    !                        r**5                            dx**2   
    !
    !               3 grav r_E**2 x y                      d**2exner
    !     S_x,y = - -----------------           + cp theta ---------,                             (18.2)
    !                      r**5                              dx dy 
    !
    !               3 grav r_E**2 x z                      d**2exner
    !     S_x,z = - -----------------           + cp theta ---------,                             (18.3)
    !                      r**5                              dx dz
    !   
    !             grav r_E**2 ( r**2 - 3 y**2 )            d**2exner
    !     S_y,y = ----------------------------- + cp theta ---------,                             (18.4)
    !                        r**5                            dy**2
    !
    !               3 grav r_E**2 y z                      d**2exner
    !     S_y,z = - -----------------           + cp theta ---------,                             (18.5)
    !                      r**5                              dy dz
    !
    !             grav r_E**2 ( r**2 - 3 z**2 )            d**2exner
    !     S_z,z = ----------------------------- + cp theta ---------,                             (18.6)
    !                        r**5                            dz**2
    !
    ! where grav is the gravitational acceleration and r = sqrt( x**2 + y**2 + z**2 ).
    ! Terms containing the factor grav are the components of ([G][G]phi)_Cart. 
    ! For the shallow-atmosphere approximation they would be absent ([G][G]phi = 0).
    ! In addition, we find for the vectorial Coriolis factor:
    !
    !     (2[W])_Cart = (0, 0, 2 W_z) = (0, 0, 2 Omega) = (0, 0, f_z).                            (19)
    !
    ! where Omega is the magnitude of [W]. 
    ! From (18) and (19) we can compute the coefficients a and b in Eq. (12) according to:
    !
    !     a = S_x,x + S_y,y + S_z,z + 4 Omega**2,                                                 (20)
    !
    !     b = S_x,x S_y,y - S_x,y**2 + S_y,y S_z,z - S_y,z**2 + S_x,x S_z,z - S_x,z**2 + f_z**2 S_z,z.   (21)
    !
    ! The coefficient c = |{S}| can be computed from a Laplace expansion of (17) or some other rule 
    ! for the computation of the determinant of a 3 x 3 matrix.
    ! 
    ! The second alternative is spherical geographic coordinates.  
    ! Measuring {S} in these coordinates and arranging the coefficients in a matrix, 
    ! denoted by ({.})_spheric yields:
    !
    !                     ( S_lon,lon  S_lon,lat  S_lon,r )   ( S_lon,lon  S_lon,lat  S_lon,r )
    !     ({S})_spheric = ( S_lat,lon  S_lat,lat  S_lat,r ) = ( S_lon,lat  S_lat,lat  S_lat,r ),  (22)
    !                     ( S_r,lon    S_r,lat    S_r,r   )   ( S_lon,r    S_lat,r    S_r,r   )
    !
    ! with:
    !
    !                  r_E**2            [        1        d**2exner   tan(lat) dexner    1  dexner ]
    ! S_lon,lon = grav ------ + cp theta [ --------------- --------- - -------- ------ + --- ------ ], (23.1)
    !                   r**3             [ (r cos(lat))**2  dlon**2      r**2    dlat     r    dr   ]
    !
    !                      [       1       d**2exner      tan(lat)   dexner ]
    ! S_lon,lat = cp theta [ ------------- --------- + ------------- ------ ],                    (23.2)
    !                      [ r**2 cos(lat) dlon dlat   r**2 cos(lat)  dlon  ]
    !
    !                    [     1      d**2exner         1       dexner ]
    ! S_lon,r = cp theta [ ---------- --------- - ------------- ------ ],                         (23.3)
    !                    [ r cos(lat)  dlon dr    r**2 cos(lat)  dlon  ]
    !
    !                  r_E**2            [  1   d**2exner    1  dexner ]
    ! S_lat,lat = grav ------ + cp theta [ ---- --------- + --- ------ ],                         (23.4)
    !                   r**3             [ r**2  dlat**2     r    dr   ]
    !
    !                    [  1   d**2exner    1   dexner ]
    ! S_lat,r = cp theta [ --- ---------- - ---- ------ ],                                        (23.5)
    !                    [  r    dlat dr    r**2  dlat  ]
    !
    !                 r_E**2            [ d**2exner ]
    ! S_r,r = -2 grav ------ + cp theta [ --------- ],                                            (23.6)
    !                  r**3             [   dr**2   ]
    !
    ! r_E is the radius of Earth. 
    ! The vectorial Coriolis factor is:
    !
    !     (2[W])_spheric = (2 W_lon, 2 W_lat, 2 W_r) = (0, 2 Omega cos(lat), 2 Omega sin(lat)) = (0, f_lat, f_r),   (24)
    !
    ! where Omega is the magnitude of [W]. 
    ! The coefficients a and b in Eq. (12) are:
    !
    !     a = S_lon,lon + S_lat,lat + S_r,r + 4 Omega**2,                                            (25)
    !
    !     b = S_lon,lon S_lat,lat - S_lon,lat**2 + S_lat,lat S_r,r - S_lat,r**2
    !       + S_lon,lon S_r,r - S_lon,r**2 + f_lat**2 S_lat,lat + f_r**2 S_r,r + 2 f_lat f_r S_lat,  (26)
    !
    ! Both alternatives for the coordinate system have their pros and cons.
    ! In our first try, we used the spherical coordinates. 
    ! This has the advantage that we can account for the atmospheric geometry with the vertical direction 
    ! as the distinct direction of the system. In addition the storage of many time-independent auxiliary fields 
    ! typically requires significantly less memory as compared to the Cartesian coordinates, 
    ! because these fields have to be 3d only up to the height where the terrain imprint 
    ! on the grid layer interfaces vanishes (-> nflatlev). Above vertical 1d profiles are sufficient.
    ! Nevertheless, this first try was not successful. 
    ! Problems occur, for instance, with the computation of the stability tensor at the poles. 
    ! So we came to use the Cartesian coordinates. 
    ! This means a great simplification of the computation which is therefore less error-prone, 
    ! but the memory consumption of the auxiliary fields is considerably larger.
    !
    ! Further technical details follow below.

    !-------------------------------------------
    !               Computation
    !-------------------------------------------

    !------------
    ! preparation
    !------------

    ! domain index
    jg = p_patch%id

    ! number of vertical grid levels
    nlev = p_patch%nlev

    ! number of blocks
    nblks_c = p_patch%nblks_c

    ! check presence of optional arguments
    present_p_prog_rcf      = PRESENT(opt_p_prog_rcf)
    present_kstart_moist    = PRESENT(opt_kstart_moist)
    present_condensate_list = PRESENT(opt_condensate_list)

    ! type of computation
    mode_standard    = .FALSE.
    mode_hydrostatic = .FALSE.
    SELECT CASE(parcelfreq2_mode)
    CASE(11, 12)
      ! we investigate the unmodified 
      ! prognostic state of the model atmosphere
      mode_standard = .TRUE.
      mode_info_str = "standard mode"
    CASE(21, 22)
      ! we investigate the diagnostic 
      ! hydrostatically balanced state of the model atmosphere
      mode_hydrostatic = .TRUE.
      mode_info_str    = "hydrostatic mode"
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid computation mode for parcelfreq2')
    END SELECT

    ! how does the air parcel perceive the atmosphere?
    switch_off_hor_percept = .FALSE.
    SELECT CASE(parcelfreq2_mode)
    CASE(11, 21)
      ! parcel perceives horizontal 
      ! and vertical variation of atmosphere
      switch_off_hor_percept = .FALSE.
    CASE(12, 22)
      ! parcel perceives only vertical variation of atmosphere
      ! (e.g. to compare the result with the square of the Brunt-Vaisala frequency)
      switch_off_hor_percept = .TRUE.
    END SELECT

    ! please note that we take into account opt_startlev and opt_endlev 
    ! in a poor man's way only in that they are applied where 
    ! the final assignment of the output field parcelfreq2 takes place
    ! (no adjustment of opt_startlev > opt_endlev to opt_startlev <= opt_endlev)
    IF (PRESENT(opt_startlev)) THEN
      i_startlev = MAX(1, opt_startlev)
    ELSE
      i_startlev = 1
    ENDIF
    IF (PRESENT(opt_endlev)) THEN
      i_endlev   = MIN(opt_endlev, SIZE(parcelfreq2, 2), nlev)
    ELSE
      i_endlev   = nlev
    ENDIF

    IF (i_startlev > i_endlev) RETURN

    ! timer monitoring
    IF (PRESENT(opt_timer)) THEN
      timer = opt_timer
    ELSE
      timer = .FALSE.
    ENDIF

    ! transmit messages if required
    IF (PRESENT(opt_verbose)) THEN
      verbose = opt_verbose
    ELSE
      verbose = .FALSE.
    ENDIF

    ! switch on timer monitoring if desired
    IF(timer) CALL timer_start(timer_opt_diag_atmo_parcelfreq2)

    dom_str = TRIM(int2string(jg))

    IF (verbose) CALL message(TRIM(routine), &
      & 'Computation of square of parcel frequency ('//TRIM(mode_info_str)//') started on domain '//dom_str)

    ! perform additional checks
    IF (PRESENT(opt_minute)) THEN
      minute = opt_minute
    ELSE
      minute = .FALSE.
    ENDIF

    ! shall we update p_diag%pres and p_diag%tempv?
    IF (PRESENT(opt_update_pres_temp)) THEN
      update_pres_temp = opt_update_pres_temp
    ELSE
      update_pres_temp = .FALSE.
    ENDIF

    ! check for optional index of grid layer above which 
    ! no condensed water phases exist in the model
    IF (present_kstart_moist) THEN
      kstart_moist = MAX(1, MIN(opt_kstart_moist, nlev))
    ELSE
      kstart_moist = 0
    ENDIF

    ! for the pressure and temperature updates we need the concentrations of the water phases,
    ! the condensate list and the start level for the water physics 
    IF ( update_pres_temp .AND.                                                             &
      & .NOT. (present_p_prog_rcf .AND. present_condensate_list .AND. present_kstart_moist) ) &
      & CALL finish(TRIM(routine), 'opt_update_pres_temp = .true. requires presence of '&
      & //'opt_p_prog_rcf, opt_condensate_list and opt_kstart_moist')

    ! is this the last call on this domain?
    IF (PRESENT(opt_lastcall)) THEN
      lastcall = opt_lastcall
    ELSE
      lastcall = .FALSE.
    ENDIF

    ! set up auxiliary fields which have to be stored for the entire run
    ! (no explicit destruction may take place!)
    IF (.NOT. spt4airParcel%initializedSqrOfParcelFreq) THEN
      ! allocate spt4sqrOfParcelFreq and compute auxiliary fields
      ! for the very first call
      CALL setup_sqrOfParcelFreq( p_patch      = p_patch,      & !inout
        &                         p_int_state  = p_int_state,  & !in
        &                         p_metrics    = p_metrics,    & !in
        &                         stencil_type = stencil_type, & !in
        &                         opt_n_dom    = n_dom,        & !optin
        &                         opt_verbose  = opt_verbose,  & !optin
        &                         opt_minute   = minute        ) !optin
    ELSE
      IF (.NOT. spt4sqrOfParcelFreq(jg)%initialized) THEN
        ! compute auxiliary fields for the first calls on the other domains
        CALL setup_sqrOfParcelFreq( p_patch      = p_patch,      & !inout
          &                         p_int_state  = p_int_state,  & !in
          &                         p_metrics    = p_metrics,    & !in
          &                         stencil_type = stencil_type, & !in
          &                         opt_verbose  = opt_verbose,  & !optin
          &                         opt_minute   = minute        ) !optin
      ENDIF
    ENDIF

    !------------------------
    ! start main computation
    !------------------------
    
    ! which triangle neighborhood to use for the computation
    ! of Hessians of scalar fields?
    ! * 10 -> "edge stencil":   current triangle + nearest neighbors + their nearest neighbors
    ! * 13 -> "vertex stencil": current triangle + all surrounding triangles that share 
    !                           one edge or one vertex with the current triangle
    nstencil_max = spt4sqrOfParcelFreq(jg)%nstencil_max

    ! set convenience pointers
    blk            => spt4sqrOfParcelFreq(jg)%trglNbhdBlk
    idx            => spt4sqrOfParcelFreq(jg)%trglNbhdIdx
    eR_X           => p_patch%cells%cartesian_center(:,:)%x(1)
    eR_Y           => p_patch%cells%cartesian_center(:,:)%x(2)
    eR_Z           => p_patch%cells%cartesian_center(:,:)%x(3)
    X              => spt4sqrOfParcelFreq(jg)%X
    Y              => spt4sqrOfParcelFreq(jg)%Y
    Z              => spt4sqrOfParcelFreq(jg)%Z
    GradGradPhi_XX => spt4sqrOfParcelFreq(jg)%GradGradPhi_XX
    GradGradPhi_YY => spt4sqrOfParcelFreq(jg)%GradGradPhi_YY
    GradGradPhi_ZZ => spt4sqrOfParcelFreq(jg)%GradGradPhi_ZZ
    GradGradPhi_XY => spt4sqrOfParcelFreq(jg)%GradGradPhi_XY
    GradGradPhi_XZ => spt4sqrOfParcelFreq(jg)%GradGradPhi_XZ
    GradGradPhi_YZ => spt4sqrOfParcelFreq(jg)%GradGradPhi_YZ

    IF (present_condensate_list) THEN
      condensate_list => opt_condensate_list
    ELSE
      condensate_list => NULL()
    ENDIF

    IF (present_p_prog_rcf) THEN
      p_prog_rcf => opt_p_prog_rcf
    ELSE
      ! this is just a dummy assignment, 
      ! p_prog_rcf will not be used in this case
      p_prog_rcf => p_prog
    ENDIF

    ! loop over prognostic domain
    ! + first halo cell row (i.e. cells that share one edge 
    ! or only one vertex with a prognostic cell )
    ! + first boundary cell row (in case of nests or limited areas)
    rl_start   = start_bdry_lev1_cells
    rl_end     = end_halo_lev2_cells
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !---------------------------------------------
    ! update pressure and temperature if required
    !---------------------------------------------

    IF (update_pres_temp) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        CALL diag_pres( pt_prog        = p_prog,                           & !in
          &             pt_diag        = p_diag,                           & !inout
          &             p_metrics      = p_metrics,                        & !in
          &             jb             = jb,                               & !in
          &             i_startidx     = i_startidx,                       & !in
          &             i_endidx       = i_endidx,                         & !in
          &             slev           = 1,                                & !in
          &             nlev           = nlev,                             & !in
          &             opt_lconstgrav = upatmo_config(jg)%dyn%l_constgrav ) !optin

        CALL diag_temp( pt_prog         = p_prog,          & !in
          &             pt_prog_rcf     = p_prog_rcf,      & !in
          &             condensate_list = condensate_list, & !in
          &             pt_diag         = p_diag,          & !inout
          &             jb              = jb,              & !in
          &             i_startidx      = i_startidx,      & !in
          &             i_endidx        = i_endidx,        & !in
          &             slev            = 1,               & !in
          &             slev_moist      = kstart_moist,    & !in
          &             nlev            = nlev             ) !in

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL      

    ENDIF  !IF (update_pres_temp)

    !-------------------------
    ! select computation mode
    !-------------------------

    IF (mode_standard) THEN

      ! standard mode

      exner  => p_prog%exner
      thetav => p_prog%theta_v

    ELSEIF (mode_hydrostatic) THEN

      ! hydrostatic mode

      ALLOCATE( exner_hydro(nproma,nlev,nblks_c),  &
        &       thetav_hydro(nproma,nlev,nblks_c), &
        &       STAT=istat                         )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of exner_hydro and thetav_hydro failed')

      IF (minute) THEN 
!$OMP PARALLEL
        CALL init(exner_hydro(:,:,:))
        CALL init(thetav_hydro(:,:,:))
!$OMP END PARALLEL
      ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jk =  i_startlev, i_endlev

          ! compute hydrostatic Exner pressure 
          ! from hydrostatic pressure
          exner_hydro(i_startidx:i_endidx,jk,jb) = exner_from_pres( p_diag%pres(i_startidx:i_endidx,jk,jb) )
          
          ! compute hydrostatic virtual potential temperature
          ! from virtual temperature and hydrostatic Exner pressure
          thetav_hydro(i_startidx:i_endidx,jk,jb) = theta_from_temp_and_exner( &
            & p_diag%tempv(i_startidx:i_endidx,jk,jb), &
            & exner_hydro(i_startidx:i_endidx,jk,jb)   )

        ENDDO  !jk
      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      ! set convenience pointers
      exner  => exner_hydro
      thetav => thetav_hydro

    ENDIF  !IF (mode_standard/hydrostatic)

    ! initialize error indicator for computation of Hessian 
    ! and LAPACK routine DHSEQR
    ierr_hessian(:) = SUCCESS
    ierr_dhseqr(:)  = SUCCESS

    ! loop over prognostic domain
    rl_start   = start_prog_cells
    rl_end     = end_prog_cells
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, istat,                                &
!$OMP            GradGradExner_XX, GradGradExner_YY, GradGradExner_ZZ, GradGradExner_XY, &
!$OMP            GradGradExner_XZ, GradGradExner_YZ, stab_XX, stab_YY, stab_ZZ, stab_XY, &
!$OMP            stab_XZ, stab_YZ, cubic_coeff_A, cubic_coeff_B, cubic_coeff_C,          &
!$OMP            comp_mat, work, eval_Re, eval_Im, schurvec) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      ! initialize fields
      IF (minute) THEN
        parcelfreq2(:,:,jb)   = 0._wp
        GradGradExner_XX(:,:) = 0._wp
        GradGradExner_YY(:,:) = 0._wp
        GradGradExner_ZZ(:,:) = 0._wp
        GradGradExner_XY(:,:) = 0._wp
        GradGradExner_XZ(:,:) = 0._wp
        GradGradExner_YZ(:,:) = 0._wp
      ENDIF

      CALL hessian_of_scalar_cc_inblk( scalar       = exner(:,:,:),          & !in
        &                              hessian_XX   = GradGradExner_XX(:,:), & !inout
        &                              hessian_YY   = GradGradExner_YY(:,:), & !inout
        &                              hessian_ZZ   = GradGradExner_ZZ(:,:), & !inout
        &                              hessian_XY   = GradGradExner_XY(:,:), & !inout
        &                              hessian_XZ   = GradGradExner_XZ(:,:), & !inout
        &                              hessian_YZ   = GradGradExner_YZ(:,:), & !inout
        &                              X            = X(:,:,:),              & !in
        &                              Y            = Y(:,:,:),              & !in
        &                              Z            = Z(:,:,:),              & !in
        &                              blk          = blk(:,:,jb),           & !in
        &                              idx          = idx(:,:,jb),           & !in
        &                              jb           = jb,                    & !in
        &                              i_startidx   = i_startidx,            & !in
        &                              i_endidx     = i_endidx,              & !in
        &                              i_startlev   = i_startlev,            & !in
        &                              i_endlev     = i_endlev,              & !in
        &                              nproma       = nproma,                & !in
        &                              nlev         = nlev,                  & !in
        &                              nstencil_max = nstencil_max,          & !in
        &                              opt_status   = istat                  ) !optout

      ierr_hessian(jb) = istat

      IF (switch_off_hor_percept) THEN
        
        ! E.g. for comparison with the square of the Brunt-Vaisala frequency:
        ! * parcel perceives only vertical variation of atmosphere
        ! * Coriolis effects are switched off
        ! * displacement is purely vertical D[r] = Dr[e_r], 
        !   with radial unit vector [e_r]
        ! These constraints reduce Eq. (14) to:
        !
        !    d**2Dr/dt**2 + ([e_r].{S}.[e_r])Dr = 0
        !
        ! => omega**2 = [e_r].{S}.[e_r]
        !
        DO jk =  i_startlev, i_endlev
          DO jc = i_startidx, i_endidx

            stab_XX = GradGradPhi_XX(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_XX(jc,jk)
            stab_YY = GradGradPhi_YY(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_YY(jc,jk)
            stab_ZZ = GradGradPhi_ZZ(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_ZZ(jc,jk)
            stab_XY = GradGradPhi_XY(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_XY(jc,jk)
            stab_XZ = GradGradPhi_XZ(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_XZ(jc,jk)
            stab_YZ = GradGradPhi_YZ(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_YZ(jc,jk)

            parcelfreq2(jc,jk,jb) = eR_X(jc,jb)**2 * stab_XX                      &
              &                   + eR_Y(jc,jb)**2 * stab_YY                      &
              &                   + eR_Z(jc,jb)**2 * stab_ZZ                      &
              &                   + 2._wp * ( eR_X(jc,jb) * eR_Y(jc,jb) * stab_XY &
              &                             + eR_X(jc,jb) * eR_Z(jc,jb) * stab_XZ &
              &                             + eR_Y(jc,jb) * eR_Z(jc,jb) * stab_YZ )
            
          ENDDO  !jc
        ENDDO  !jk
        
      ELSE ! The parcel perceives the atmospheric variation unrestrictedly

        DO jk =  i_startlev, i_endlev
          ! please see NEC user guide why we need to close 
          ! and reopen the jc-loop over and over again in the following:
          DO jc = i_startidx, i_endidx

            stab_XX = GradGradPhi_XX(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_XX(jc,jk)
            stab_YY = GradGradPhi_YY(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_YY(jc,jk)
            stab_ZZ = GradGradPhi_ZZ(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_ZZ(jc,jk)
            stab_XY = GradGradPhi_XY(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_XY(jc,jk)
            stab_XZ = GradGradPhi_XZ(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_XZ(jc,jk)
            stab_YZ = GradGradPhi_YZ(jc,jk,jb) + cpd * thetav(jc,jk,jb) * GradGradExner_YZ(jc,jk)

            ! coefficients of cubic equation (10):
            !
            ! * a = stab_x,x + stab_y,y + stab_z,z + 4 Omega**2
            !
            ! NOTE: we already add the negative sign required for Eq. (10)!
            cubic_coeff_A = -( stab_XX + stab_YY + stab_ZZ + earth_angular_velocity2_t_4 )
            !
            ! * b = stab_x,x stab_y,y - stab_x,y**2 
            !     + stab_y,y stab_z,z - stab_y,z**2 
            !     + stab_x,x stab_z,z - stab_x,z**2 
            !     + f_z**2 stab_z,z
            !
            cubic_coeff_B = stab_XX * stab_YY - stab_XY**2 &
              &           + stab_XX * stab_ZZ - stab_XZ**2 &
              &           + stab_YY * stab_ZZ - stab_YZ**2 &
              &           + earth_angular_velocity2_t_4 * stab_ZZ
            !
            ! * c = det(stab)
            !
            ! NOTE: we already add the negative sign required for Eq. (10)!
            cubic_coeff_C = -( stab_XX * ( stab_YY * stab_ZZ - stab_YZ**2        ) &
              &              - stab_XY * ( stab_XY * stab_ZZ - stab_YZ * stab_XZ ) &
              &              + stab_XZ * ( stab_XY * stab_YZ - stab_YY * stab_XZ ) )

            ! The following part is concerned with finding the roots of Eq. (10):
            !
            !     (omega**2)**3 -             a (omega**2)**2 +             b (omega**2) - c             = 0
            ! <=> (omega**2)**3 + cubic_coeff_A (omega**2)**2 + cubic_coeff_B (omega**2) + cubic_coeff_C = 0.
            !
            ! Subroutines or functions that solve a cubic equation do not exist in ICON as far as we know. 
            ! We make use of the following method:
            ! The problem of finding the roots of the cubic equation 
            ! is transformed into solving an eigenvalue problem whose characteristic polynomial 
            ! is equal to the cubic equation. For the latter a LAPACK routine exists.

            ! set up 3x3 companion matrix
            !
            !            ( 0 0 -cubic_coeff_C )
            ! comp_mat = ( 1 0 -cubic_coeff_B )
            !            ( 0 1 -cubic_coeff_A )
            !
            comp_mat(1,1,jc) = 0._wp
            comp_mat(2,1,jc) = 1._wp 
            comp_mat(3,1,jc) = 0._wp
            comp_mat(1,2,jc) = 0._wp
            comp_mat(2,2,jc) = 0._wp
            comp_mat(3,2,jc) = 1._wp 
            comp_mat(1,3,jc) = -cubic_coeff_C
            comp_mat(2,3,jc) = -cubic_coeff_B
            comp_mat(3,3,jc) = -cubic_coeff_A

          ENDDO  !jc

          DO jc = i_startidx, i_endidx

            ! Finding the roots of Eq. (10) is now equivalent 
            ! to solving the eigenvalue problem:
            !
            ! det(comp_mat - l * E) = 0, 
            !
            ! where E is the 3x3 identity matrix, and l = l_1,2,3 are the three eigenvalues 
            ! and at the same time the roots of Eq. (10).

            ! LAPACK routine to solve the eigenvalue problem:
            ! the companion matrix is a upper Hessenberg matrix, 
            ! and for the eigenvalues of such a kind of matrix a special routine does exist.
            ! Its documentation can be found on: http://www.netlib.org/lapack/explore-html/
            CALL DHSEQR( 'E',                  & != JOB   (in)    compute eigenvalues only
              &          'N',                  & != COMPZ (in)    no Schur vectors are computed
              &          3,                    & != N     (in)    comp_mat(N,N=3)
              &          1,                    & != ILO   (in)    should be 1 in our case
              &          3,                    & != IHI   (in)    should be N = 3 in our case
              &          comp_mat(1:3,1:3,jc), & != H     (inout) companion matrix comp_mat(N=3,N=3)
              &          3,                    & != LDH   (in)    comp_mat(N=LDH=3,N)
              &          eval_Re(1:3,jc),      & != WR    (out)   real part of eigenvalues eval_Re(3)
              &          eval_Im(1:3,jc),      & != WI    (out)   imaginary part of eigenvalues eval_Im(3)
              &          schurvec(1:1,1:3),    & != Z     (inout) Schur vectors schurvec(LDZ=1,3) (not referenced)
              &          1,                    & != LDZ   (in)    can be 1, because COMPZ = 'N'
              &          work(1:nwork),        & != WORK  (out)   work space: work(LWORK=nwork)
              &          nwork,                & != LWORK (in)
              &          istat                 ) != INFO  (out)   error indicator

            ! DHSEQR indicates errors by negative or positive non-zero values of istat, 
            ! we only store the information if any error occured to simplify matters
            ierr_dhseqr(jb) = ierr_dhseqr(jb) + ABS(istat)

          ENDDO  !jc

          DO jc = i_startidx, i_endidx

            ! Following Eq. (9) the air parcel oscillates with:
            !
            !     D[r](t) = <[r]> exp{-i omega t} = <[r]> exp{-i [Re(omega) + i Im(omega)] t}
            !             = <[r]> exp{-i Re(omega) t} exp{Im(omega) t}.
            !
            ! The most unstable/least stable mode 
            ! is the mode with the greatest imaginary part of omega.
            !
            ! The final result we are intereste in is:
            !
            !     parcelfreq2 = -[Im(omega)_max]**2.
            !
            ! Now, omega**2 = eval_Re + i eval_Im, 
            ! and  omega    = +-sqrt( omega**2 ).
            !
            ! The principal value for the root 
            ! of a complex number C = A + iB is:
            !
            !               ( |C| + A )                ( |C| - A )
            ! SQRT(C) = SQRT( ------- ) + i sgn(B) SQRT( ------- ),
            !               (    2    )                (    2    )
            !                                          -----------
            ! where sgn(B) = -1 for B < 0 and sgn(B) = 1 otherwise.
            ! The second root is the negative of the principal value.
            ! So [Im(omega)_max]**2 is equal to the maximum of the underlined term.

            parcelfreq2(jc,jk,jb) = -0.5_wp * MAXVAL(SQRT(eval_Re(1:3,jc)**2 + eval_Im(1:3,jc)**2) - eval_Re(1:3,jc))

          ENDDO  !jc
        ENDDO  !jk

      ENDIF  !IF (switch_off_hor_percept)

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL 
    
    ! check for error messages from computation of Hessian and LAPACK routine DHSEQR
    IF (ANY(ierr_hessian(:) /= SUCCESS)) CALL finish(TRIM(routine), 'Error in computation of Hessian')
    IF (ANY(ierr_dhseqr(:)  /= SUCCESS)) CALL finish(TRIM(routine), 'Error in LAPACK routine DHSEQR')

    !----------
    ! clean-up
    !----------

    IF (mode_hydrostatic) THEN
      DEALLOCATE( exner_hydro,  &
        &         thetav_hydro, &
        &         STAT=istat    )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of exner_hydro and thetav_hydro failed')
    ENDIF

    NULLIFY(blk, idx, eR_X, eR_Y, eR_Z, X, Y, Z, GradGradPhi_XX, GradGradPhi_YY, GradGradPhi_ZZ, &
      &     GradGradPhi_XY, GradGradPhi_XZ, GradGradPhi_YZ)

    IF (present_condensate_list) condensate_list => NULL()

    IF (lastcall) CALL destruct_sqrOfParcelFreq( p_patch     = p_patch,    & !inout
      &                                          opt_verbose = opt_verbose ) !optin

    IF (verbose) CALL message(TRIM(routine), &
      & 'Computation of square of parcel frequency finished on domain '//dom_str)

    IF(timer) CALL timer_stop(timer_opt_diag_atmo_parcelfreq2)

  END SUBROUTINE sqr_of_parcel_freq

  !=====================================================================================================

  !>
  !! Initialization subroutine for: sqr_of_parcel_freq
  !!
  SUBROUTINE setup_sqrOfParcelFreq( p_patch,      & !inout
    &                               p_int_state,  & !in
    &                               p_metrics,    & !in
    &                               stencil_type, & !in
    &                               opt_n_dom,    & !optin
    &                               opt_verbose,  & !optin
    &                               opt_minute    ) !optin

    ! in/out variables

    TYPE(t_patch),                INTENT(INOUT) :: p_patch
    TYPE(t_int_state),            INTENT(IN)    :: p_int_state
    TYPE(t_nh_metrics),           INTENT(IN)    :: p_metrics
    INTEGER,                      INTENT(IN)    :: stencil_type
    INTEGER,            OPTIONAL, INTENT(IN)    :: opt_n_dom
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_verbose
    LOGICAL,            OPTIONAL, INTENT(IN)    :: opt_minute

    ! local variables

    REAL(wp), ALLOCATABLE :: phi(:,:,:)

    INTEGER  :: ierr_hessian(p_patch%nblks_c)
    REAL(wp) :: radius, grav_radiusEarth2
    INTEGER  :: ndom, nlev, nblks_c, nstencil_max
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: jg, jb, jc, jk, jstencil
    INTEGER  :: istat
    LOGICAL  :: verbose, minute
    CHARACTER(LEN=2) :: dom_str

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':setup_sqrOfParcelFreq'

    !------------------------------------------------

    ! Please note that this subroutine may not be called 
    ! for all of the n_dom domains.

    ! domain index
    jg = p_patch%id

    ! transmit messages if required
    IF (PRESENT(opt_verbose)) THEN
      verbose = opt_verbose
    ELSE
      verbose = .FALSE.
    ENDIF

    dom_str = TRIM(int2string(jg))

    IF (verbose) CALL message(TRIM(routine), &
      & 'Setup of square of parcel frequency started on domain '//dom_str)

    ! number of domains
    IF (PRESENT(opt_n_dom)) THEN
      ndom = opt_n_dom
    ELSE
      ndom = -999
    ENDIF

    ! perform additional checks?
    IF (PRESENT(opt_minute)) THEN
      minute = opt_minute
    ELSE
      minute = .FALSE.
    ENDIF

    ! number of blocks
    nblks_c = p_patch%nblks_c

    ! number of vertical grid levels
    nlev = p_patch%nlev
    
    IF (.NOT. spt4airParcel%initializedSqrOfParcelFreq) THEN
      
      IF (.NOT. PRESENT(opt_n_dom)) CALL finish(TRIM(routine), &
        & 'opt_n_dom has to be present as safety indicator of first call')
      ALLOCATE( spt4sqrOfParcelFreq(ndom), &
        &       STAT=istat                 )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of spt4sqrOfParcelFreq failed')
      
      ! indicate setup
      spt4airParcel%initializedSqrOfParcelFreq = .TRUE.
      
    ENDIF
    
    IF (.NOT. spt4sqrOfParcelFreq(jg)%initialized) THEN
      
      ! which triangle neighborhood to use for the computation
      ! of Hessians of scalar fields?
      ! * 10 -> "edge stencil":   current triangle + nearest neighbors + their nearest neighbors
      ! * 13 -> "vertex stencil": current triangle + all surrounding triangles that share 
      !                           one edge or one vertex with the current triangle
      SELECT CASE(stencil_type)
      CASE(10)
        nstencil_max = 10
      CASE(13)
        nstencil_max = p_int_state%cell_environ%max_nmbr_nghbr_cells
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Invalid stencil type')
      END SELECT

      spt4sqrOfParcelFreq(jg)%nstencil_max = nstencil_max

      ALLOCATE( spt4sqrOfParcelFreq(jg)%trglNbhdBlk(nproma,nstencil_max,nblks_c), &
        &       spt4sqrOfParcelFreq(jg)%trglNbhdIdx(nproma,nstencil_max,nblks_c), &
        &       spt4sqrOfParcelFreq(jg)%X(nproma,nlev,nblks_c),                   &
        &       spt4sqrOfParcelFreq(jg)%Y(nproma,nlev,nblks_c),                   &
        &       spt4sqrOfParcelFreq(jg)%Z(nproma,nlev,nblks_c),                   &
        &       spt4sqrOfParcelFreq(jg)%GradGradPhi_XX(nproma,nlev,nblks_c),      &
        &       spt4sqrOfParcelFreq(jg)%GradGradPhi_YY(nproma,nlev,nblks_c),      &
        &       spt4sqrOfParcelFreq(jg)%GradGradPhi_ZZ(nproma,nlev,nblks_c),      &
        &       spt4sqrOfParcelFreq(jg)%GradGradPhi_XY(nproma,nlev,nblks_c),      &
        &       spt4sqrOfParcelFreq(jg)%GradGradPhi_XZ(nproma,nlev,nblks_c),      &
        &       spt4sqrOfParcelFreq(jg)%GradGradPhi_YZ(nproma,nlev,nblks_c),      &
        &       phi(nproma,nlev,nblks_c),                                         &
        &       STAT=istat                                                        )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation failed')

      ! some auxiliary quantities
      grav_radiusEarth2 = grav * grid_sphere_radius**2

      ! initialize error indicator for computation of Hessian 
      ierr_hessian(:) = SUCCESS

      ! initialize with zero
      IF (minute) THEN
!$OMP PARALLEL
        CALL init(spt4sqrOfParcelFreq(jg)%trglNbhdBlk(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%trglNbhdIdx(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%X(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%Y(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%Z(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%GradGradPhi_XX(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%GradGradPhi_YY(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%GradGradPhi_ZZ(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%GradGradPhi_XY(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%GradGradPhi_XZ(:,:,:))
        CALL init(spt4sqrOfParcelFreq(jg)%GradGradPhi_YZ(:,:,:))
        CALL init(phi(:,:,:))
!$OMP END PARALLEL
      ENDIF

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)
      
      ! loop over prognostic domain
      rl_start   = start_prog_cells
      rl_end     = end_prog_cells
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)
      
!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx, jstencil) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        IF (stencil_type == 10) THEN

          ! nstencil_max = 10 -> "edge stencil": 
          ! * current triangle + nearest neighbors + their nearest neighbors

          DO jstencil = 1, nstencil_max
            DO jc = i_startidx, i_endidx
              spt4sqrOfParcelFreq(jg)%trglNbhdBlk(jc,jstencil,jb) = p_int_state%rbf_c2grad_blk(jstencil,jc,jb)
              spt4sqrOfParcelFreq(jg)%trglNbhdIdx(jc,jstencil,jb) = p_int_state%rbf_c2grad_idx(jstencil,jc,jb)
            ENDDO  !jc
          ENDDO  !jstencil

        ELSEIF (stencil_type == 13) THEN

          ! 13 -> "vertex stencil": 
          ! * current triangle + all surrounding triangles that share 
          !   one edge or one vertex with the current triangle

          ! Please note: if the three vertices of the current triangle are hexagon points, 
          ! the number of triangles in the neighborhood is 13. If one of the vertices is 
          ! one of the 12 pentagon points on the globe, the number of triangles in the neighborhood is 12. 
          ! In the latter case p_int_state%cell_environ%blk/idx contain block and index 
          ! of the current cell for jstencil = nstencil_max = 13
          ! (see /src/shr_horizontal/mo_intp_rbf_coeffs: gen_index_list_radius).
          ! This "double counting" of the current triangle should do no harm 
          ! to the least-squares equation system which is set up and solved in hessian_of_scalar_cc_inblk below, 
          ! since the presence of two linearly dependent equations should not change its character.

          DO jstencil = 1, nstencil_max
            DO jc = i_startidx, i_endidx
              spt4sqrOfParcelFreq(jg)%trglNbhdBlk(jc,jstencil,jb) = p_int_state%cell_environ%blk(jc,jb,jstencil)
              spt4sqrOfParcelFreq(jg)%trglNbhdIdx(jc,jstencil,jb) = p_int_state%cell_environ%idx(jc,jb,jstencil)
            ENDDO  !jc
          ENDDO  !jstencil

        ENDIF  !IF (stencil_type == ...)

      ENDDO  !jb
!$OMP END DO
      
      ! loop over prognostic domain + first halo cell row 
      ! (here, the first halo cell row means triangles that share one edge 
      ! with triangles of the prognostic domain + triangles that share only one vertex 
      ! with triangles of the prognostic domain)
      ! + first boundary cell row (in case of nests or limited areas)
      rl_start   = start_bdry_lev1_cells
      rl_end     = end_halo_lev2_cells
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)
      
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, radius) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jk = 1, nlev
          
          DO jc = i_startidx, i_endidx

            ! Cartesian components of position:
            !
            ! (x, y, z) = radius * [ cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat) ]
            !
            ! Please note: the coordinate z we refer to in the following is NOT the common vertical coordinate!
            ! If we denote the latter with zeta, the relation is:
            !
            !     z = (radius_Earth + zeta) * sin(lat).
            !
            ! Important: the spherical geometry of the atmosphere imprints on the Cartesian coordinates 
            ! in case of the shallow-atmosphere mode, too!

            radius = grid_sphere_radius + p_metrics%z_mc(jc,jk,jb)

            spt4sqrOfParcelFreq(jg)%X(jc,jk,jb) = radius * p_patch%cells%cartesian_center(jc,jb)%x(1)
            spt4sqrOfParcelFreq(jg)%Y(jc,jk,jb) = radius * p_patch%cells%cartesian_center(jc,jb)%x(2)
            spt4sqrOfParcelFreq(jg)%Z(jc,jk,jb) = radius * p_patch%cells%cartesian_center(jc,jb)%x(3)

          ENDDO  !jc
        ENDDO  !jk

        IF (.NOT. upatmo_config(jg)%dyn%l_constgrav) THEN

          DO jk = 1, nlev          
            DO jc = i_startidx, i_endidx
 
              ! geopotential of the deep atmosphere:
              !
              ! phi = grav * (r_E / r)**2, 
              !
              ! where r_E denotes the radius of the Earth

              radius = grid_sphere_radius + p_metrics%z_mc(jc,jk,jb)

              phi(jc,jk,jb) = grav_radiusEarth2 / (grid_sphere_radius + p_metrics%z_mc(jc,jk,jb))**2

            ENDDO  !jc
          ENDDO  !jk

        ELSE

          DO jk = 1, nlev          
            DO jc = i_startidx, i_endidx
 
              ! geopotential of the shallow atmosphere: 
              ! please see NEC user guide why we cannot set 
              ! GradGradPhi_... = 0 directly
              phi(jc,jk,jb) = 0._wp

            ENDDO  !jc
          ENDDO  !jk

        ENDIF  !IF (.NOT. upatmo_config(jg)%dyn%l_constgrav)

      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      ! The six independent components of the Hessian of the geopotential ([G][G]phi)_Cart:
      !
      !                    ( GradGradPhi_x,x  GradGradPhi_x,y  GradGradPhi_x,z )
      ! ([G][G]phi)_Cart = ( GradGradPhi_x,y  GradGradPhi_y,y  GradGradPhi_y,z )
      !                    ( GradGradPhi_x,z  GradGradPhi_y,z  GradGradPhi_z,z )
      !
      ! Please note: in case of the shallow atmosphere 
      ! or in case of the deep-atmosphere with const. gravity [G][G]phi = 0 holds. 
      !
      ! The analytic expressions for the components of ([G][G]phi)_Cart read:
      !
      !                       grav r_E**2 ( r**2 - 3 x**2 )
      !     GradGradPhi_x,x = ----------------------------- 
      !                                  r**5            
      !
      ! and the like for GradGradPhi_y,y and GradGradPhi_z,z, 
      !
      !                         3 grav r_E**2 x y
      !     GradGradPhi_x,y = - -----------------
      !                              r**5            
      !
      ! and the like for GradGradPhi_x,y, GradGradPhi_x,z and  GradGradPhi_y,z.
      !
      ! However, we cannot compute the analytic expression for ([G][G]phi)_Cart, 
      ! please constult NEC user guide for reasons.
      
      ! loop over prognostic domain
      rl_start   = start_prog_cells
      rl_end     = end_prog_cells
      i_startblk = p_patch%cells%start_block(rl_start)
      i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, istat) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        CALL hessian_of_scalar_cc_inblk( scalar       = phi(:,:,:),                                     & !in
          &                              hessian_XX   = spt4sqrOfParcelFreq(jg)%GradGradPhi_XX(:,:,jb), & !inout
          &                              hessian_YY   = spt4sqrOfParcelFreq(jg)%GradGradPhi_YY(:,:,jb), & !inout
          &                              hessian_ZZ   = spt4sqrOfParcelFreq(jg)%GradGradPhi_ZZ(:,:,jb), & !inout
          &                              hessian_XY   = spt4sqrOfParcelFreq(jg)%GradGradPhi_XY(:,:,jb), & !inout
          &                              hessian_XZ   = spt4sqrOfParcelFreq(jg)%GradGradPhi_XZ(:,:,jb), & !inout
          &                              hessian_YZ   = spt4sqrOfParcelFreq(jg)%GradGradPhi_YZ(:,:,jb), & !inout
          &                              X            = spt4sqrOfParcelFreq(jg)%X(:,:,:),               & !in
          &                              Y            = spt4sqrOfParcelFreq(jg)%Y(:,:,:),               & !in
          &                              Z            = spt4sqrOfParcelFreq(jg)%Z(:,:,:),               & !in
          &                              blk          = spt4sqrOfParcelFreq(jg)%trglNbhdBlk(:,:,jb),    & !in
          &                              idx          = spt4sqrOfParcelFreq(jg)%trglNbhdIdx(:,:,jb),    & !in
          &                              jb           = jb,                                             & !in
          &                              i_startidx   = i_startidx,                                     & !in
          &                              i_endidx     = i_endidx,                                       & !in
          &                              i_startlev   = 1,                                              & !in
          &                              i_endlev     = nlev,                                           & !in
          &                              nproma       = nproma,                                         & !in
          &                              nlev         = nlev,                                           & !in
          &                              nstencil_max = nstencil_max,                                   & !in
          &                              opt_status   = istat                                           ) !optout

        ierr_hessian(jb) = istat

      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      ! check for error messages from computation of Hessian
      IF (ANY(ierr_hessian(:) /= SUCCESS)) CALL finish(TRIM(routine), 'Error in computation of Hessian')

      ! clean-up
      DEALLOCATE(phi, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of phi failed')

      ! indicate setup
      spt4sqrOfParcelFreq(jg)%initialized = .TRUE.
      
   ENDIF  !IF (.NOT. spt4sqrOfParcelFreq(jg)%initialized)
    
    IF (verbose) CALL message(TRIM(routine), &
      & 'Setup of square of parcel frequency finished on domain '//dom_str)

  END SUBROUTINE setup_sqrOfParcelFreq

  !=====================================================================================================

  !>
  !! Finalization subroutine for: sqr_of_parcel_freq
  !!
  SUBROUTINE destruct_sqrOfParcelFreq( p_patch,    & !inout
    &                                  opt_verbose ) !optin

    ! in/out variables

    TYPE(t_patch),           INTENT(INOUT) :: p_patch
    LOGICAL,       OPTIONAL, INTENT(IN)    :: opt_verbose

    ! local variables

    INTEGER  :: jg, istat
    LOGICAL  :: verbose
    CHARACTER(LEN=2) :: dom_str

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':destruct_sqrOfParcelFreq'
    
    !------------------------------------------------

    ! domain index
    jg = p_patch%id

    ! transmit messages if required
    IF (PRESENT(opt_verbose)) THEN
      verbose = opt_verbose
    ELSE
      verbose = .FALSE.
    ENDIF
    
    dom_str = TRIM(int2string(jg))
    
    IF (verbose) CALL message(TRIM(routine), &
      & 'Destruction of square of parcel frequency started on domain '//dom_str)
    
    IF (spt4airParcel%initializedSqrOfParcelFreq) THEN
      
      ! please note that the fields of spt4sqrOfParcelFreq(jg)
      ! may not have been allocated for all of the n_dom domains
      
      IF (spt4sqrOfParcelFreq(jg)%initialized) THEN

        DEALLOCATE( spt4sqrOfParcelFreq(jg)%trglNbhdBlk,    &
          &         spt4sqrOfParcelFreq(jg)%trglNbhdIdx,    &
          &         spt4sqrOfParcelFreq(jg)%X,              &
          &         spt4sqrOfParcelFreq(jg)%Y,              &
          &         spt4sqrOfParcelFreq(jg)%Z,              &
          &         spt4sqrOfParcelFreq(jg)%GradGradPhi_XX, &
          &         spt4sqrOfParcelFreq(jg)%GradGradPhi_YY, &
          &         spt4sqrOfParcelFreq(jg)%GradGradPhi_ZZ, &
          &         spt4sqrOfParcelFreq(jg)%GradGradPhi_XY, &
          &         spt4sqrOfParcelFreq(jg)%GradGradPhi_XZ, &
          &         spt4sqrOfParcelFreq(jg)%GradGradPhi_YZ, &
          &         STAT=istat                              )
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of spt4sqrOfParcelFreq%... failed')
        
        ! indicate destruction
        spt4sqrOfParcelFreq(jg)%initialized = .FALSE.
        
      ENDIF  !IF (spt4sqrOfParcelFreq(jg)%initialized)
      
      ! if this is the last call among all domains, 
      ! we deallocate spt4sqrOfParcelFreq itself
      IF (.NOT. ANY(spt4sqrOfParcelFreq(:)%initialized)) THEN
        
        DEALLOCATE(spt4sqrOfParcelFreq, STAT=istat)
        IF(istat /= SUCCESS) &
          & CALL finish(TRIM(routine), 'Deallocation of spt4sqrOfParcelFreq failed')
        
        ! indicate destruction
        spt4airParcel%initializedSqrOfParcelFreq = .FALSE.
        
      ENDIF
      
    ENDIF  !IF (spt4airParcel%initializedSqrOfParcelFreq)
    
    IF (verbose) CALL message(TRIM(routine), &
      & 'Destruction of square of parcel frequency finished on domain '//dom_str)

  END SUBROUTINE destruct_sqrOfParcelFreq

  !=====================================================================================================

  !>
  !! Compute Hessian of a scalar field in Cartesian coordinates
  !!
  !! Compute Hessian of a scalar field in Cartesian coordinates,
  !! for use within block loops.
  !!
  !! Method:
  !! The scalar field (denoted by s) is approximated by a Taylor expansion 
  !! up to second order in the vicinity of the cell center:
  !!
  !! s ~ s0 + D[r] * [G]s + (1/2) D[r] * ([G][G]s) * D[r],             (1)
  !!
  !! where s0 denotes the value at the cell center, D[r] is the distance vector 
  !! to the position where s is approximated, [G]s and [G][G]s are the gradient 
  !! and the Hessian of the scalar field at the cell center, and * denotes the scalar product. 
  !! s0, [G]s and [G][G]s are the unknowns, which are determined by a least-squares fit 
  !! of Eq. (1) to a prescribed cell neighborhood. 
  !! s0 is added to the list of unknowns just to make the fitting a bit more flexible.
  !! 
  SUBROUTINE hessian_of_scalar_cc_inblk( scalar,       & !in
    &                                    hessian_XX,   & !inout
    &                                    hessian_YY,   & !inout
    &                                    hessian_ZZ,   & !inout
    &                                    hessian_XY,   & !inout
    &                                    hessian_XZ,   & !inout
    &                                    hessian_YZ,   & !inout
    &                                    X,            & !in
    &                                    Y,            & !in
    &                                    Z,            & !in
    &                                    blk,          & !in
    &                                    idx,          & !in
    &                                    jb,           & !in
    &                                    i_startidx,   & !in
    &                                    i_endidx,     & !in
    &                                    i_startlev,   & !in
    &                                    i_endlev,     & !in
    &                                    nproma,       & !in
    &                                    nlev,         & !in
    &                                    nstencil_max, & !in
    &                                    opt_status    ) !optout

    ! in/out variables

    REAL(wp),           INTENT(IN)    :: scalar(:,:,:)   !< (nproma,nlev,nblks_c) scalar field
    ! Please see NEC user guide for reasons 
    ! why we need fields for each component separately:
    REAL(wp),           INTENT(INOUT) :: hessian_XX(:,:) !< (nproma,nlev) Hessian of scalar: xx-component
    REAL(wp),           INTENT(INOUT) :: hessian_YY(:,:) !< (nproma,nlev) Hessian of scalar: yy-component
    REAL(wp),           INTENT(INOUT) :: hessian_ZZ(:,:) !< (nproma,nlev) Hessian of scalar: zz-component
    REAL(wp),           INTENT(INOUT) :: hessian_XY(:,:) !< (nproma,nlev) Hessian of scalar: xy-component
    REAL(wp),           INTENT(INOUT) :: hessian_XZ(:,:) !< (nproma,nlev) Hessian of scalar: xz-component
    REAL(wp),           INTENT(INOUT) :: hessian_YZ(:,:) !< (nproma,nlev) Hessian of scalar: yz-component
    REAL(wp),           INTENT(IN)    :: X(:,:,:)        !< (nproma,nlev,nblks_c) x-coordinate of cell center
    REAL(wp),           INTENT(IN)    :: Y(:,:,:)        !< (nproma,nlev,nblks_c) y-coordinate of cell center
    REAL(wp),           INTENT(IN)    :: Z(:,:,:)        !< (nproma,nlev,nblks_c) z-coordinate of cell center
    INTEGER,            INTENT(IN)    :: blk(:,:)        !< (nproma,nstencil_max) triangle neighborhood stencil: blocks
    INTEGER,            INTENT(IN)    :: idx(:,:)        !< (nproma,nstencil_max) triangle neighborhood stencil: indices
    INTEGER,            INTENT(IN)    :: jb              !< block index
    INTEGER,            INTENT(IN)    :: i_startidx      !< start index for horizontal loop within block
    INTEGER,            INTENT(IN)    :: i_endidx        !< end index for horizontal loop within block
    INTEGER,            INTENT(IN)    :: i_startlev      !< start index for vertical loop
    INTEGER,            INTENT(IN)    :: i_endlev        !< end index for vertical loop
    INTEGER,            INTENT(IN)    :: nproma          !< number of grid columns in block
    INTEGER,            INTENT(IN)    :: nlev            !< number of full levels
    INTEGER,            INTENT(IN)    :: nstencil_max    !< max. number of triangles in stencil
    INTEGER,  OPTIONAL, INTENT(OUT)   :: opt_status      !< optional error code
    
    ! local variables
    
    INTEGER,  PARAMETER :: ncol  = 10
    INTEGER,  PARAMETER :: nwork = 700

    ! Why we have to allocate memory for a field with
    ! 39 x 10 x nproma x nelv elements?  => NEC user guide
    REAL(wp) :: mat(nstencil_max * 3, ncol, nproma, nlev)   !!! ~ 400 (nproma x nlev) fields !!! 
    REAL(wp) :: rhs(nstencil_max * 3, 1, nproma, nlev)      !!! ~ 40  (nproma x nlev) fields !!!
    REAL(wp) :: work(nwork)

    REAL(wp) :: dx, dy, dz
    INTEGER  :: startlev, endlev
    INTEGER  :: jc, jk, jstencil, jkk
    INTEGER  :: nrow_max, jrow
    INTEGER  :: info, error
    LOGICAL  :: present_status

    !------------------------------------------------

    !-------------
    ! preparation
    !-------------

    ! check presence of optional arguments
    present_status = PRESENT(opt_status)
 
    IF (present_status) opt_status = SUCCESS
    
    ! adjust start and end level if necessary
    IF (i_startlev == 1) THEN
      startlev = 2
    ELSE
      startlev = i_startlev
    ENDIF
    IF (i_endlev == nlev) THEN
      endlev = nlev - 1
    ELSE
      endlev = i_endlev
    ENDIF

    !------------------
    ! main computation
    !------------------ 

    ! Cartesian coordinates:
    ! The point of origin is in the center of Earth 
    ! and the positive z-axis crosses the North Pole.
    ! The relation to the geographical spherical coordinates is: 
    !
    !     x = r cos(lat) cos(lon),  y = r cos(lat) sin(lon),  z = r sin(lat).
    !
    ! Please note: the coordinate z we refer to in the following is NOT the common vertical coordinate!
    ! If we denote the latter with zeta, the relation is:
    !
    !     z =  r sin(lat) = (r_E + zeta) * sin(lat),
    !
    ! where r_E denotes the radius of the Earth.
    !
    ! The expansion of Eq. (1) in Cartesian coordinates reads:
    ! 
    ! s ~ s0 + Dx * ds/dx + Dy * ds/dy + Dz * ds/dz
    !   + (1/2) Dx**2 * d2s/dx2 + (1/2) Dy**2 * d2s/dy2 + (1/2) Dz**2 * d2s/dz2
    !   + Dx * Dy * d2s/dx/dy + Dx * Dz * d2s/dx/dz + Dy * Dz * d2s/dy/duz        (11)
    ! 
    ! The ten unknowns in Eq. (11) are ordered in the following way:
    !
    !  1: s0        (-> no output intended)
    !
    !  2: ds/dx     (1st component of gradient -> no output intended)
    !  3: ds/dy     (2nd component of gradient -> no output intended)
    !  4: ds/dz     (3rd component of gradient -> no output intended)
    !
    !  5: d2s/dx2   (1st component of Hessian)
    !  6: d2s/dy2   (2nd component of Hessian)
    !  7: d2s/dz2   (3rd component of Hessian)
    !  8: d2s/dx/dy (4th component of Hessian)
    !  9: d2s/dx/dz (5th component of Hessian)
    ! 10: d2s/dy/dz (6th component of Hessian)
    !
    ! Please note that s0 is included into the list of unknowns to make the fitting more flexible. 
    ! However, it should not be used. For instance, the least-squares fit does not guarantee s0 > 0 
    ! (the incorporation of inequality constraints would be too cumbersome).
    ! 
    ! Cell neighborhood stencil:
    ! * Horizontal:
    !   The triangle neighborhood should be composed of at least 10 triangles (9 neighbors).
    !   The simple 4-triangle neighborhood is insufficient. 
    !   This might be understood in the following way (without proof):
    !   If we assume a purely horizontal scalar field, 
    !   the Cartesian expansion of Eq. (1) would read:
    !
    !   s ~ s0 + Dx * ds/dx + Dy * ds/dy 
    !     + (1/2) Dx**2 * d2s/dx2 + (1/2) Dy**2 * d2s/dy2 + Dx * Dy * d2s/dx/dy.  (12) 
    !
    !   Eq. (12) has 6 unknowns (s0, ds/dx, ds/dy, d2s/dx2, d2s/dy2, d2s/dx/dy),
    !   so the equation system would be underdetermined, 
    !   if we would only use the 4-triangle neighborhood. 
    !   For the least-squares fit, however, we need an overdetermined system.
    !
    !   Please note: in case of nstencil_max = 13, the neighborhood is composed 
    !   of 13 triangles in the majority of cases. Only if one of the three vertices 
    !   of the central triangle is a pentagon point, the neighborhood only consists 
    !   of 12 neighbors. In this case the 13th entry of blk and idx contain 
    !   block and index of the central triangle. Double counting the equation 
    !   for the central cell in the equation system should, however, do no harm, 
    !   since the solver should actually "count" two linearly dependent equations as one.
    !   Even if the double counting should not happen to be fully neutral, 
    !   12 pentagon points in total on a global grid with thousands or millions of hexagon points
    !   do not justify a case differentiation for the problem at hand.
    !
    ! * Vertical:
    !   In the vertical direction, the neighborhood covers the current layer, 
    !   the layer below and the layer above.

    ! initialize error indicator for LAPACK routine DGELS
    error = SUCCESS

    ! max. number of equations in equation system
    nrow_max = 3 * nstencil_max

    ! initialize row index of equation system
    jrow = 0

    ! please see NEC user guide for reasons why the following challenging 
    ! loop nesting is necessary

    ! loop over triangles in neighborhood
    DO jstencil = 1, nstencil_max  ! 1 is current triangle itself, 
                                   ! 2 to nstencil_max are neighbors

      ! loop from layer above over current layer to layer below
      DO jkk = -1, 1

        ! update row index
        jrow = jrow + 1

#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = startlev, endlev  ! leave bottom and top layers off
#else
        DO jk = startlev, endlev
          DO jc = i_startidx, i_endidx
#endif
            ! distances to neighbors
            dx = X(idx(jc,jstencil),jk+jkk,blk(jc,jstencil)) - X(jc,jk,jb)
            dy = Y(idx(jc,jstencil),jk+jkk,blk(jc,jstencil)) - Y(jc,jk,jb)
            dz = Z(idx(jc,jstencil),jk+jkk,blk(jc,jstencil)) - Z(jc,jk,jb)
            
            !------------------------------------------------------------------------
            ! right-hand side of equation system: 
            ! scalar field in cell neighborhood
            ! (we subtract the scalar field of the central cell as a const. field)
            rhs(jrow,1,jc,jk) = scalar(idx(jc,jstencil),jk+jkk,blk(jc,jstencil)) - scalar(jc,jk,jb)
            !------------------------------------------------------------------------
            ! coefficients for the ten unknowns above:
            ! 1: s0 (or rather s0 - s(jc,jk,jb))
            ! => 1
            mat(jrow,1,jc,jk) = 1._wp
            !------------------------------------------------------------------------
            ! 2: ds/dx
            ! => Dx
            mat(jrow,2,jc,jk) = dx
            !------------------------------------------------------------------------
            ! 3: ds/dy
            ! => Dy
            mat(jrow,3,jc,jk) = dy
            !------------------------------------------------------------------------
            ! 4: ds/dz
            ! => Dz
            mat(jrow,4,jc,jk) = dz
            !------------------------------------------------------------------------
            ! 5: d2s/dx2
            ! => (1/2) Dx**2
            mat(jrow,5,jc,jk) = 0.5_wp * dx**2
            !------------------------------------------------------------------------
            ! 6: d2s/dy2
            ! => (1/2) Dy**2
            mat(jrow,6,jc,jk) = 0.5_wp * dy**2
            !------------------------------------------------------------------------
            ! 7: d2s/dz2
            ! => (1/2) Dz**2 
            mat(jrow,7,jc,jk) = 0.5_wp * dz**2
            !------------------------------------------------------------------------
            ! 8: d2s/dx/dy
            ! => Dx * Dy
            mat(jrow,8,jc,jk) = dx * dy
            !------------------------------------------------------------------------
            ! 9: d2s/dx/dz
            ! => Dx * Dz
            mat(jrow,9,jc,jk) = dx * dz
            !------------------------------------------------------------------------
            ! 10: d2s/dy/dz
            ! => Dx * Dz
            mat(jrow,10,jc,jk) = dy * dz
            !------------------------------------------------------------------------

          ENDDO  !jc/jk
        ENDDO  !jk/jc
        
      ENDDO  !jkk
    ENDDO  !jstencil

    DO jk = startlev, endlev
      DO jc = i_startidx, i_endidx

        ! Please note: solving for the Hessian is by far the most expensive part 
        ! of the computation of parcelfreq2. Although the LAPACK routine DGELS, 
        ! called in the following, is vectorized in some way, it seems that this way
        ! is not optimal for the problem at hand. 
        ! If we denote the least-squares equation system for the moment:
        !
        !     A * x = b,                                                          (21)
        !
        ! with matrix A having m x n elements, solution vector x having n elements 
        ! and right-hand side vector b having m elements, 
        ! we find for stencil 13: m = 13 * 3 = 39 and n = 10. 
        ! The costs for solving for x at each grid point are probably 
        ! of the order of the costs of solving one time step of the dynamical core 
        ! (at the minimum), all vectorization aside. 
        ! In the following we list some options to reduce the costs, 
        ! but none is easily implemented:
        !
        ! 1) Typically, algorithms that solve Eq. (21) directly, such as DGELS, 
        !    use some kind of QR-decomposition:
        !
        !      A * x = Q * R * x  | Q^T *                                         (22)
        ! =>   Q^T * Q * R * x = E * R * x = R * x = Q^T * b,                     (23)
        !
        !    with the orthogonal m x n matrix Q and the upper triangular n x n matrix R. 
        !    (.)^T denotes the transpose and E = Q^T * Q is the n x n identity matrix. 
        !    Eq. (23) can then be solved directly for x. 
        !    In the case at hand, A is time-independent and so are Q and R. 
        !    So we could save a large part of the costs, if we would do the QR-decomposition once 
        !    and store Q (or Q^T) and R for reuse. 
        !    However, if we estimate for the minimum momory consumption of Q and R 
        !    the memory consumption of A, we would have to store:
        !      39 x 10 x nproma x nlev x nblks_c
        !    elements, i.e. about 400 additional 3d fields! 
        !    This is out of the question, unfortunately.
        !
        !    BUT NOTE: currently a call of this subroutine allocates already 
        !              ~ 400 x nproma x nlev for mat!
        !
        ! 2) We could try to extract the algorithm of DGELS and tailor it 
        !    to our problem in order to gain a more optimal vectorization. 
        !    Unfortunately, the algorithm is quite complicated, DGELS makes much use 
        !    of nested subroutine calls. So the implementation of this option is not trivial.
        !
        ! 3) If we consider a quadratic invertible m x m matrix A = A0 + A', where A' is "small" 
        !    compared to A0, the inverse of A can be approximated by:
        !
        !      A^(-1) ~ A0^(-1) - A0^(-1) * A' * A0^(-1).                         (24)
        !
        !    Provided something comparable exists for m x n matrices, 
        !    we could proceed in the following way:
        !    * Do the QR-decomposition once for a grid-specific idealized nominal 
        !      pentagon and hexagon column and store Q and R, 
        !      i.e. determine A0^(-1) in terms of Eq. (24).
        !    * At each computation of parcelfreq2, determine the difference between 
        !      the real stencil of a grid cell and the idealized stencil, i.e. A'.
        !    * Assume that this difference is "small" and apply the analog to Eq. (24) 
        !      in order to obtain an approximate solution for x. 
        !    Above flat_height this yields probably a relatively accurate solution, 
        !    but below, where the grid layers are terrain-following, we have to expect 
        !    relatively large errors. 
        !
        ! 4) Probably the most preferable solution would be to dispense with the Eq. (21) altoghether 
        !    and compute the Hessian directly. This is relatively simple on a regular rectangular grid. 
        !    Unfortunately, on a triangular terrain-following grid this is extremely cumbersome 
        !    and of all options the most difficult to implement.

        ! Algorithms to solve the linear least-squares system: 
        ! mat * x = rhs for x directly (Gram-Schmidt, Householder, Givens etc.) 
        ! are extremely complicated. For this reason we resort to a LAPACK routine
        ! (the extraction and modification of its algorithm in order to inline it at this place 
        ! is beyond our means for the time being.)
        ! Its documentation can be found on:
        ! http://www.netlib.org/lapack/explore-html/ 
        ! In short the algorithm does a QR-decomposition of mat. 
        CALL DGELS( 'N',                          & != TRANS (in)    we solve mat * x = rhs, not mat**T * x = rhs
          &         nrow_max,                     & != M     (in)    number of rows of mat(M=nrow_max,N)
          &         ncol,                         & != N     (in)    number of columns of mat(M,N=ncol)
          &         1,                            & != NRHS  (in)    number of columns of rhs(M,NRHS=1)
          &         mat(1:nrow_max,1:ncol,jc,jk), & != A     (inout) on entry: the matrix mat(M=nrow_max,N=ncol)
                                                    !                on exit:  a representation of the QR-decomposition 
                                                    !                (which we will not use)
          &         nrow_max,                     & != LDA   (in)    here equal to M
          &         rhs(1:nrow_max,1:1,jc,jk),    & != B     (inout) on entry: the vector rhs(M=nrow_max,1)
                                                    !                on exit:  the solution x(N=ncol) = rhs(1:N=ncol,1)
          &         nrow_max,                     & != LDB   (in)    here equal to M
          &         work(1:nwork),                & != WORK  (out)   work space: work(LWORK=nwork)
          &         nwork,                        & != LWORK (in)
          &         info                          ) != INFO  (out)   error indicator

        ! DGELS indicates errors by negative or positive non-zero values of istat, 
        ! we only store the information if any error occured to simplify matters
        error = error + ABS(info)
          
        ! now, rhs should contain the solution:
        ! -----------------------
        !    i   |   rhs(i,1)
        ! -----------------------
        !    1   |   s0 (or s0 - s(jc,jk,jb))
        !    2   |   ds/dx
        !    3   |   ds/dy
        !    4   |   ds/dz
        !    5   |   d2s/dx2
        !    6   |   d2s/dy2
        !    7   |   d2s/dz2
        !    8   |   d2s/dx/dy
        !    9   |   d2s/dx/dz
        !   10   |   d2s/dy/dz
        ! -----------------------

      ENDDO  !jc

      DO jc = i_startidx, i_endidx
        ! the 6 independent components of the Hessian tensor of the scalar
        hessian_XX(jc,jk) = rhs(5, 1,jc,jk)
        hessian_YY(jc,jk) = rhs(6, 1,jc,jk)
        hessian_ZZ(jc,jk) = rhs(7, 1,jc,jk)
        hessian_XY(jc,jk) = rhs(8, 1,jc,jk)
        hessian_XZ(jc,jk) = rhs(9, 1,jc,jk)
        hessian_YZ(jc,jk) = rhs(10,1,jc,jk)
      ENDDO  !jc
    ENDDO  !jk

    ! simple constant extrapolation to model bottom and top
    IF (i_startlev == 1) THEN
      DO jc = i_startidx, i_endidx
        hessian_XX(jc,1) = hessian_XX(jc,2)
        hessian_YY(jc,1) = hessian_YY(jc,2)
        hessian_ZZ(jc,1) = hessian_ZZ(jc,2)
        hessian_XY(jc,1) = hessian_XY(jc,2)
        hessian_XZ(jc,1) = hessian_XZ(jc,2)
        hessian_YZ(jc,1) = hessian_YZ(jc,2)
      ENDDO  !jc
    ENDIF
    IF (i_endlev == nlev) THEN
      DO jc = i_startidx, i_endidx
        hessian_XX(jc,nlev) = hessian_XX(jc,nlev-1)
        hessian_YY(jc,nlev) = hessian_YY(jc,nlev-1)
        hessian_ZZ(jc,nlev) = hessian_ZZ(jc,nlev-1)
        hessian_XY(jc,nlev) = hessian_XY(jc,nlev-1)
        hessian_XZ(jc,nlev) = hessian_XZ(jc,nlev-1)
        hessian_YZ(jc,nlev) = hessian_YZ(jc,nlev-1)
      ENDDO  !jc
    ENDIF

    IF (error /= SUCCESS) THEN
      IF (present_status) opt_status = 101
    ENDIF

  end SUBROUTINE hessian_of_scalar_cc_inblk

END MODULE mo_diag_atmo_air_parcel

