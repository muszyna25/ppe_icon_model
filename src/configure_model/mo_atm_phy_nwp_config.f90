!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_atm_phy_nwp_config

  USE mo_kind,                ONLY: wp
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_parallel_config,     ONLY: nproma
  USE mo_io_units,            ONLY: filename_max
  USE mo_impl_constants,      ONLY: max_dom, MAX_CHAR_LENGTH, itconv, itccov,  &
    &                               itrad, itradheat, itsso, itgscp, itsatad,  &
    &                               itturb, itsfc, itgwd, itfastphy,           &
    &                               iphysproc, iphysproc_short, ismag, iedmf
  USE mo_math_constants,      ONLY: dbl_eps, pi_2, deg2rad
  USE mo_exception,           ONLY: message, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_radiation_config,    ONLY: irad_o3
  USE mo_les_config,          ONLY: configure_les
  USE mo_limarea_config,      ONLY: configure_latbc

  USE mo_name_list_output_config, ONLY: first_output_name_list, is_variable_in_output

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_phy_nwp_config, t_atm_phy_nwp_config, dt_phy
  PUBLIC :: configure_atm_phy_nwp
  PUBLIC :: lrtm_filename
  PUBLIC :: cldopt_filename
  PUBLIC :: ltuning_kessler, icpl_aero_conv, icpl_o3_tp, iprog_aero
  PUBLIC :: ltuning_ozone

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for atm dynamics
  !!--------------------------------------------------------------------------
  TYPE :: t_atm_phy_nwp_config

    ! namelist variables

    INTEGER ::  inwp_gscp        !> microphysics
    INTEGER ::  inwp_satad       !! saturation adjustment
    INTEGER ::  inwp_convection  !! convection
    LOGICAL ::  lshallowconv_only !! use shallow convection only
    LOGICAL ::  ldetrain_conv_prec !! detrain convective rain and snow
    INTEGER ::  inwp_radiation   !! radiation
    INTEGER ::  inwp_sso         !! sso
    INTEGER ::  inwp_gwd         !! non-orographic gravity wave drag
    INTEGER ::  inwp_cldcover    !! cloud cover
    INTEGER ::  inwp_turb        !! turbulence
    INTEGER ::  inwp_surface     !! surface including soil, ocean, ice,lake
    INTEGER  :: itype_z0         !! type of roughness length data
    REAL(wp) :: dt_conv          !> time step for convection
    REAL(wp) :: dt_ccov          !! time step for subscale cloud cover
    REAL(wp) :: dt_rad           !! "-"                     radiation
    REAL(wp) :: dt_sso           !! "-"  for subscale orographic gravity waves
    REAL(wp) :: dt_gwd           !! "-"  for subscale gravity waves
    REAL(wp) :: dt_fastphy       !! time step for fast physics processes
                                 !! microphysics, saturation adjustment, turbulence, 
                                 !! surface (in addition: update and radheat)
    ! hydci_pp                   
    REAL(wp) :: mu_rain          !! parameter in gamma distribution for rain
    REAL(wp) :: mu_snow          !! ...for snow
    REAL(wp) :: qi0, qc0

    INTEGER  :: icpl_aero_gscp     !! type of aerosol-microphysics coupling

    REAL(wp) :: ustart_raylfric    !! velocity at which extra Rayleigh friction starts
    REAL(wp) :: efdt_min_raylfric  !! e-folding time corresponding to maximum relaxation 
                                   !! coefficient
    LOGICAL  :: latm_above_top     !! use extra layer above model top for radiation 
                                   !! (reduced grid only)

    ! Derived variables

    LOGICAL :: lproc_on(iphysproc) !> contains information about status of 
                                   !! corresponding physical process
                                   !! ON: TRUE; OFF: FALSE
    LOGICAL :: lcalc_acc_avg       ! TRUE: calculate accumulated and averaged quantities

    LOGICAL :: lcalc_moist_integral_avg ! TRUE: calculate temporally averaged vertical integrals of moisture fields

    LOGICAL :: lcalc_extra_avg     ! TRUE: calculate aditional temporally averaged fields, which normally 
                                   !       are not computed in operational runs.
                                   !       lcalc_extra_avg is set to true automatically, if any of the 
                                   !       non-standard fields is specified in the output namelist.

    LOGICAL :: lcalc_dpsdt         ! TRUE: compute dpsdt for output even if a low message level (<= 10) is selected

    LOGICAL :: is_les_phy          !>TRUE is turbulence is 3D 
                                   !>FALSE otherwise

    INTEGER :: nclass_gscp         !> number of hydrometeor classes for 
                                   ! chosen grid scale microphysics

    ! Tuning variables

    REAL(wp), ALLOCATABLE :: fac_ozone(:)         ! vertical profile funtion for ozone tuning
    REAL(wp), ALLOCATABLE :: shapefunc_ozone(:,:) ! horizontal profile funtion for ozone tuning
    REAL(wp)              :: ozone_maxinc         ! maximum allowed change of ozone mixing ratio

  END TYPE t_atm_phy_nwp_config

  !>
  !!
  TYPE(t_atm_phy_nwp_config) :: atm_phy_nwp_config(max_dom) !< shape: (n_dom)

  !> NetCDF file containing longwave absorption coefficients and other data
  !> for RRTMG_LW k-distribution model ('rrtmg_lw.nc')
  CHARACTER(LEN=filename_max) :: lrtm_filename

  !> NetCDF file with RRTM Cloud Optical Properties for ECHAM6
  CHARACTER(LEN=filename_max) :: cldopt_filename

  INTEGER  :: icpl_aero_conv     !! type of coupling between aerosols and convection scheme
  INTEGER  :: iprog_aero         !! type of prognostic aerosol
  INTEGER  :: icpl_o3_tp         !! type of coupling between ozone and the tropopause

  REAL(wp) ::  &                       !> Field of calling-time interval (seconds) for
    &  dt_phy(max_dom,iphysproc_short) !! each domain and phys. process

  !!--------------------------------------------------------------------------
  !! Tuning parameters for physics
  !!--------------------------------------------------------------------------
  
  ! convection:
  ! GZ, 2013-09-13: tuning to reduce drizzle (may be overridden by icpl_aero_conv=1)
  LOGICAL,  PARAMETER :: ltuning_kessler  = .TRUE.

  ! ozone tuning if GEMS climatology is used:
  LOGICAL  :: ltuning_ozone 
  REAL(wp) :: tune_ozone_ztop
  REAL(wp) :: tune_ozone_zmid, tune_ozone_zmid2
  REAL(wp) :: tune_ozone_zbot
  REAL(wp) :: tune_ozone_fac
  REAL(wp) :: tune_ozone_lat
  REAL(wp) :: tune_ozone_maxinc
  INTEGER  :: ozone_shapemode

CONTAINS

  SUBROUTINE configure_atm_phy_nwp( n_dom, p_patch, dtime )
    !-------------------------------------------------------------------------
    !
    !>
    !! Setup NWP physics
    !!
    !! Read namelist for physics. Choose the physical package and subsequent
    !! parameters.
    !!
    !! @par Revision History
    !! Initial revision by Daniel Reinert, DWD (2010-10-06)
    !! revision for restructuring by Kristina Froehlich MPI-M (2011-07-12)
  
    TYPE(t_patch), TARGET,INTENT(IN) :: p_patch(:)
    INTEGER, INTENT(IN) :: n_dom
    REAL(wp),INTENT(IN) :: dtime
  
    INTEGER :: jg, jk, jk_shift, jb, jc
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &      routine = 'mo_atm_phy_nwp_config:configure_atm_phy_nwp'
    REAL(wp) :: z_mc_ref


    DO jg = 1,n_dom
      atm_phy_nwp_config(jg)%dt_fastphy = (dtime/2._wp**(p_patch(jg)%level &
          &                            -  p_patch(1)%level))            !seconds
    ENDDO

    dt_phy(:,:) = 0._wp

    DO jg = 1,n_dom

      ! initialize lproc_on
      atm_phy_nwp_config(jg)%lproc_on(1:iphysproc)    = .FALSE.


      ! Fill derived variable lproc_on
 
      IF ( atm_phy_nwp_config(jg)%inwp_satad > 0 )               &
        &  atm_phy_nwp_config(jg)%lproc_on(itsatad)   = .TRUE. 

      IF ( atm_phy_nwp_config(jg)%inwp_convection > 0 )          &
        &  atm_phy_nwp_config(jg)%lproc_on(itconv)    = .TRUE. 

      IF ( atm_phy_nwp_config(jg)%inwp_cldcover > 0 )            &
        &  atm_phy_nwp_config(jg)%lproc_on(itccov)    = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_radiation > 0 )           &
        &  atm_phy_nwp_config(jg)%lproc_on(itrad)     = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_sso > 0 )                 &
        &  atm_phy_nwp_config(jg)%lproc_on(itsso)     = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_gscp > 0 )                &
        &  atm_phy_nwp_config(jg)%lproc_on(itgscp)    = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_turb > 0 )                &
        &  atm_phy_nwp_config(jg)%lproc_on(itturb)    = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_radiation > 0 )           &
        &  atm_phy_nwp_config(jg)%lproc_on(itradheat) = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_surface > 0 )             &
        &  atm_phy_nwp_config(jg)%lproc_on(itsfc)     = .TRUE.

      IF ( atm_phy_nwp_config(jg)%inwp_gwd > 0 )                 &
        &  atm_phy_nwp_config(jg)%lproc_on(itgwd)     = .TRUE.


      ! Slow physics:
      ! currently for each domain the same time intervals are set
      !
      ! Fast physics:
      ! for each fast physics process the time interval is set equal to the 
      ! time interval for advection.


      ! Slow physics time steps

      dt_phy(jg,itconv) = atm_phy_nwp_config(jg)% dt_conv    ! sec

      dt_phy(jg,itsso)  = atm_phy_nwp_config(jg)% dt_sso     ! sec

      dt_phy(jg,itgwd)  = atm_phy_nwp_config(jg)% dt_gwd     ! sec

      dt_phy(jg,itrad)  = atm_phy_nwp_config(jg)% dt_rad     ! sec

      !> KF always call clouds after convection
      !! to ensure the proper output of Qx_tot
      IF ( atm_phy_nwp_config(jg)%lproc_on(itccov)   .AND.  &
        &  atm_phy_nwp_config(jg)%lproc_on(itconv) ) THEN
        atm_phy_nwp_config(jg)% dt_ccov = atm_phy_nwp_config(jg)% dt_conv
      ELSE
        atm_phy_nwp_config(jg)% dt_ccov = atm_phy_nwp_config(jg)% dt_fastphy ! sec
      ENDIF

      ! For EDMF DUALM cloud cover is called every turbulence time step
      IF ( atm_phy_nwp_config(jg)%inwp_turb == iedmf ) THEN
        atm_phy_nwp_config(jg)% dt_ccov = atm_phy_nwp_config(jg)% dt_fastphy ! sec
      ENDIF

      dt_phy(jg,itccov) = atm_phy_nwp_config(jg)% dt_ccov   ! sec

      ! Fast physics time step
      dt_phy(jg,itfastphy) = atm_phy_nwp_config(jg)%dt_fastphy ! sec



      IF (jg == 1) THEN
        ! issue a warning, if advective and convective timesteps are not synchronized
        !
        ! so far, only for Domain 1
        IF( MOD(atm_phy_nwp_config(jg)%dt_conv,dtime) > 10._wp*dbl_eps )  THEN
          WRITE(message_text,'(a,2F8.1)') &
            &'WARNING: convective timestep is not a multiple of advective timestep: ', &
            & dt_phy(jg,itconv), dtime
          CALL message(TRIM(routine), TRIM(message_text))
          WRITE(message_text,'(a,F8.1)') &
            &'implicit synchronization in time_ctrl_physics: dt_conv !=!', &
            & REAL((FLOOR(atm_phy_nwp_config(jg)%dt_conv/dtime) + 1),wp) &
            & * dtime
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF
      ENDIF  ! jg=1


    ENDDO  ! jg loop

    atm_phy_nwp_config(1:n_dom)%lcalc_dpsdt = &
      is_variable_in_output(first_output_name_list, var_name='ddt_pres_sfc')

    !Configure LES physics
    DO jg = 1, n_dom
    
      atm_phy_nwp_config(jg)%is_les_phy = .FALSE. 
    
      IF(atm_phy_nwp_config(jg)%inwp_turb==ismag)THEN
        CALL configure_les(jg,dtime)
        atm_phy_nwp_config(jg)%is_les_phy = .TRUE. 
      END IF 
    
      !convection should be turned off for LES
      IF(atm_phy_nwp_config(jg)%inwp_convection>0.AND.atm_phy_nwp_config(jg)%is_les_phy)THEN
        CALL message(TRIM(routine),'Turning off convection for LES!')
        atm_phy_nwp_config(jg)%inwp_convection = 0
      END IF

      !SSO should be turned off for LES
      IF(atm_phy_nwp_config(jg)%inwp_sso>0.AND.atm_phy_nwp_config(jg)%is_les_phy)THEN
        CALL message(TRIM(routine),'Turning off SSO scheme for LES!')
        atm_phy_nwp_config(jg)%inwp_sso = 0
      END IF

      !GWD should be turned off for LES
      IF(atm_phy_nwp_config(jg)%inwp_gwd>0.AND.atm_phy_nwp_config(jg)%is_les_phy)THEN
        CALL message(TRIM(routine),'Turning off GWD scheme for LES!')
        atm_phy_nwp_config(jg)%inwp_gwd = 0
      END IF
 
 
      ! Round up dt_rad to nearest advection timestep
      IF( MOD(dt_phy(jg,itrad),dtime)>10._wp*dbl_eps .AND. atm_phy_nwp_config(jg)%is_les_phy)THEN
        ! write warning only for global domain
        IF (jg==1) THEN
          WRITE(message_text,'(a,2F8.1)') &
            &'WARNING: radiation timestep is not a multiple of fastphy timestep: ', &
            & dt_phy(jg,itrad), dtime
          CALL message(TRIM(routine), TRIM(message_text))
          WRITE(message_text,'(a,F8.1)') &
            &'radiation time step is rounded up to next multiple: dt_rad !=!', &
            & REAL((FLOOR(dt_phy(jg,itrad)/dtime) + 1),wp) * dtime
          CALL message(TRIM(routine), TRIM(message_text))
        ENDIF
        dt_phy(jg,itrad) = REAL((FLOOR(dt_phy(jg,itrad)/dtime) + 1),wp) * dtime
      ENDIF

    END DO !ndom

    !Configure lateral boundary condition for limited area model
    IF(l_limited_area) THEN
      CALL configure_latbc()
    END IF

    ! Settings for ozone tuning, depending on option for ozone climatology
    SELECT CASE (irad_o3)
    CASE (7)  ! GEMS climatology
      CALL message(TRIM(routine), 'Use GEMS ozone climatology with tuning')
      ltuning_ozone     = .TRUE.
      tune_ozone_ztop   = 30000.0_wp
      tune_ozone_zmid2  = 15000.0_wp
      tune_ozone_zmid   = 15000.0_wp
      tune_ozone_zbot   = 10000.0_wp
      tune_ozone_fac    = 0.5_wp
      ozone_shapemode   = 1        ! tuning is applied at low latitudes only
      tune_ozone_lat    = 45._wp   ! tuning ends at 45 deg
      tune_ozone_maxinc = 2.e-6_wp ! maximum absolute change of O3 mixing ratio
                                   ! this value is about 12% of the climatological maximum in the tropics
    CASE (79,97) ! Blending between GEMS and MACC climatologies
      CALL message(TRIM(routine), 'Use blending between GEMS and MACC ozone climatologies with tuning')
      ltuning_ozone     = .TRUE.
      tune_ozone_ztop   = 29000.0_wp
      tune_ozone_zmid2  = 26000.0_wp
      tune_ozone_zmid   = 20000.0_wp
      tune_ozone_zbot   = 17000.0_wp
      tune_ozone_fac    = 0.15_wp
      ozone_shapemode   = 2
      tune_ozone_lat    = 30._wp
      tune_ozone_maxinc = 1.e-6_wp
    CASE DEFAULT
      ltuning_ozone     = .FALSE.
      tune_ozone_ztop   = 30000.0_wp
      tune_ozone_zmid2  = 15000.0_wp
      tune_ozone_zmid   = 15000.0_wp
      tune_ozone_zbot   = 10000.0_wp
      tune_ozone_fac    = 0._wp
      ozone_shapemode   = 1
      tune_ozone_lat    = -1._wp
      tune_ozone_maxinc = 0._wp
    END SELECT

    ! Ozone tuning function, applied to the ozone climatology as 
    ! o3clim_tuned = o3clim*(1.+fac_ozone*shapefunc_ozone)
    DO jg = 1, n_dom
      atm_phy_nwp_config(jg)%ozone_maxinc = tune_ozone_maxinc
      ALLOCATE(atm_phy_nwp_config(jg)%fac_ozone(p_patch(jg)%nlev), &
               atm_phy_nwp_config(jg)%shapefunc_ozone(nproma,p_patch(jg)%nblks_c) )
      ! Vertical profile function
      DO jk = 1,p_patch(jg)%nlev
        jk_shift = jk+p_patch(jg)%nshift_total
        z_mc_ref = 0.5_wp*(vct_a(jk_shift)+vct_a(jk_shift+1))
        IF ( z_mc_ref > tune_ozone_zbot .AND. z_mc_ref < tune_ozone_ztop ) THEN
          IF ( z_mc_ref < tune_ozone_zmid ) THEN
            atm_phy_nwp_config(jg)%fac_ozone(jk) = tune_ozone_fac *                        &
                 & sin((z_mc_ref-tune_ozone_zbot)/(tune_ozone_zmid-tune_ozone_zbot)*pi_2)**2
          ELSE IF ( z_mc_ref < tune_ozone_zmid2 ) THEN
            atm_phy_nwp_config(jg)%fac_ozone(jk) = tune_ozone_fac 
          ELSE
            atm_phy_nwp_config(jg)%fac_ozone(jk) = tune_ozone_fac *                        &
                 & cos((z_mc_ref-tune_ozone_zmid2)/(tune_ozone_ztop-tune_ozone_zmid2)*pi_2)**2
          END IF
        ELSE
          atm_phy_nwp_config(jg)%fac_ozone(jk) = 0.0_wp
        ENDIF
      ENDDO
      ! Horizontal profile function for fac_ozone
      DO jb = 1, p_patch(jg)%nblks_c
        DO jc = 1, nproma          
          IF (ozone_shapemode == 1 .AND. tune_ozone_lat > 0._wp) THEN
            IF (ABS(p_patch(jg)%cells%center(jc,jb)%lat) < tune_ozone_lat * deg2rad) THEN
              atm_phy_nwp_config(jg)%shapefunc_ozone(jc,jb) = &
                COS(p_patch(jg)%cells%center(jc,jb)%lat * 90._wp/tune_ozone_lat)**2
            ELSE
              atm_phy_nwp_config(jg)%shapefunc_ozone(jc,jb) = 0._wp
            END IF
          ELSE IF (ozone_shapemode == 2 .AND. tune_ozone_lat > 0._wp) THEN
            IF (ABS(p_patch(jg)%cells%center(jc,jb)%lat) < tune_ozone_lat * deg2rad) THEN
              atm_phy_nwp_config(jg)%shapefunc_ozone(jc,jb) = &
                1._wp - 0.75_wp*SQRT(COS(p_patch(jg)%cells%center(jc,jb)%lat * 90._wp/tune_ozone_lat))
            ELSE
              atm_phy_nwp_config(jg)%shapefunc_ozone(jc,jb) = 1._wp
            END IF
          ELSE
            atm_phy_nwp_config(jg)%shapefunc_ozone(jc,jb) = 1.0_wp
          END IF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE configure_atm_phy_nwp


END MODULE mo_atm_phy_nwp_config
