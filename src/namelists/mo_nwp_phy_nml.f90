!>
!!  Namelist for NWP physics
!!
!!  these Subroutines are called by control model and construct the
!!  physics composition
!!
!! @author <Kristina Froehlich, DWD>
!!
!!
!! @par Revision History
!! First implementation by Kristina Froehlich, DWD (2010-06-20>)
!!  namelist varaibles for calling frequency of fast physics schemes
!!   have been removed to ensure the high frequent calls
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_phy_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: max_dom, icosmo
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,      ONLY: use_restart_namelists, isRestart

  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,    &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config,                        &
    &                               config_lrtm_filename   => lrtm_filename,   &
    &                               config_cldopt_filename => cldopt_filename, &
    &                               config_icpl_aero_conv  => icpl_aero_conv,  &
    &                               config_iprog_aero      => iprog_aero,      &
    &                               config_icpl_o3_tp      => icpl_o3_tp

  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_cuparameters,        ONLY: icapdcycl

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_nwp_phy_namelist 

   !
   ! user defined calling intervals
   !
  REAL(wp) :: dt_conv(max_dom)   !> field element for convection
  REAL(wp) :: dt_rad(max_dom)    !! "-"                     radiation
  REAL(wp) :: dt_sso(max_dom)    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gwd(max_dom)    !! "-"  for subscale gravity waves

  ! switches defining physics packages
  INTEGER  :: inwp_convection(max_dom)    !! convection
  LOGICAL  :: lshallowconv_only(max_dom)  !! use shallow convection only
  LOGICAL  :: ldetrain_conv_prec(max_dom) !! detrain convective rain and snow
  INTEGER  :: inwp_cldcover(max_dom)      !! cloud cover
  INTEGER  :: inwp_radiation(max_dom)     !! radiation
  INTEGER  :: inwp_sso(max_dom)           !! sso
  INTEGER  :: inwp_gwd(max_dom)           !! non-orographic gravity wave drag
  INTEGER  :: inwp_gscp(max_dom)          !! microphysics
  INTEGER  :: inwp_satad(max_dom)         !! saturation adjustment
  INTEGER  :: inwp_turb(max_dom)          !! turbulence
  INTEGER  :: inwp_surface(max_dom)       !! surface including soil, ocean, ice,lake

  INTEGER  :: itype_z0           !! type of roughness length data
  INTEGER  :: icpl_aero_gscp     !! type of aerosol-microphysics coupling
  INTEGER  :: icpl_aero_conv     !! type of coupling between aerosols and convection scheme
  INTEGER  :: iprog_aero         !! type of prognostic aerosol
  INTEGER  :: icpl_o3_tp         !! type of ozone-tropopause coupling
  REAL(wp) :: qi0, qc0           !! variables for hydci_pp
  REAL(wp) :: ustart_raylfric    !! velocity at which extra Rayleigh friction starts
  REAL(wp) :: efdt_min_raylfric  !! e-folding time corresponding to maximum relaxation coefficient
  LOGICAL  :: latm_above_top(max_dom) !! use extra layer above model top for radiation (reduced grid only)
  ! parameter for cloud microphysics
  real(wp) :: mu_rain            !! shape parameter in gamma distribution for rain
  real(wp) :: mu_snow            !! ...for snow
  
  !> NetCDF file containing longwave absorption coefficients and other data
  !> for RRTMG_LW k-distribution model ('rrtmg_lw.nc')
  CHARACTER(LEN=filename_max) :: lrtm_filename

  !> NetCDF file with RRTM Cloud Optical Properties for ECHAM6
  CHARACTER(LEN=filename_max) :: cldopt_filename

  NAMELIST /nwp_phy_nml/ inwp_convection, inwp_cldcover,             &
    &                    inwp_radiation, inwp_sso, inwp_gwd,         &
    &                    inwp_gscp, inwp_satad,                      &
    &                    inwp_turb, inwp_surface,                    &
    &                    dt_conv, dt_rad, dt_sso, dt_gwd,            &
    &                    qi0, qc0, icpl_aero_gscp,                   &
    &                    ustart_raylfric, efdt_min_raylfric,         &
    &                    latm_above_top, itype_z0, mu_rain,          &
    &                    mu_snow, icapdcycl, icpl_aero_conv,         &
    &                    lrtm_filename, cldopt_filename, icpl_o3_tp, &
    &                    iprog_aero, lshallowconv_only, ldetrain_conv_prec

  LOGICAL :: l_nwp_phy_namelist_read = .false.
  
CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - performs sanity checks
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_nwp_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, jg
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = 'mo_nwp_phy_nml:read_nwp_phy_namelist'

    !-----------------------
    ! 0a. dummy settings; will be replaced with defaults after reading the namelist
    !     wherever nothing is specified explicitly   
    !-----------------------
    inwp_gscp(:)       = -1
    inwp_satad(:)      = -1
    inwp_convection(:) = -1
    inwp_radiation(:)  = -1
    inwp_sso(:)        = -1
    inwp_gwd(:)        = -1
    inwp_cldcover(:)   = -1
    inwp_turb(:)       = -1
    inwp_surface(:)    = -1

    lshallowconv_only(:)  = .FALSE.
    ldetrain_conv_prec(:) = .FALSE.

    dt_conv (:) = -999._wp
    dt_rad  (:) = -999._wp
    dt_sso  (:) = -999._wp
    dt_gwd  (:) = -999._wp

    lrtm_filename   = 'rrtmg_lw.nc'  
    cldopt_filename = 'ECHAM6_CldOptProps.nc'


    !------------------------------------------------------------------
    ! 0b. Real defaults for some other variables
    !------------------------------------------------------------------

    itype_z0 = 2  !  2 = land-cover related roughness lenght only (i.e. no orographic contrib)
    qi0      = 0.0_wp 
    qc0      = 0.0_wp 

    ! shape parameter for gamma distribution for rain and snow
    mu_rain = 0.0_wp
    mu_snow = 0.0_wp

    ustart_raylfric    = 160._wp
    efdt_min_raylfric  = 10800._wp

    latm_above_top(:)  = .FALSE.  ! no extra layer above model top for radiation computation

    ! CAPE correction to improve diurnal cycle of convection (moved from mo_cuparameters)
    icapdcycl = 0  ! 0= no CAPE diurnal cycle correction (IFS default prior to cy40r1, i.e. 2013-11-19)
                   ! 1=    CAPE - surface buoyancy flux (intermediate testing option)
                   ! 2=    CAPE - subcloud CAPE (IFS default starting with cy40r1)
                   ! 3=    Apply CAPE modification of (2) over land only, with additional restriction to the tropics

    ! coupling between aersols and cloud microphysics
    icpl_aero_gscp = 0  ! 0 = none
                        ! 1 = simple coupling with aerosol climatology disregarding the dependency of aerosol activation on vertical wind speed
                        ! 2 = more accurate coupling with aerosol climatology as a function of vertical wind speed

    ! coupling between aersols and convection scheme
    icpl_aero_conv = 0  ! 0 = none
                        ! 1 = specify thresholds (QC and cloud thickness) for precip initiation depending on aerosol climatology instead of land-sea mask

    ! type of prognostic aerosol
    iprog_aero = 0  ! 0 = pure climatology
                    ! 1 = very simple prognostic scheme based on advection of and relaxation towards climatology

    ! coupling between ozone and the tropopause
    icpl_o3_tp = 1      ! 0 = none
                        ! 1 = take climatological values from 100/350 hPa above/below the tropopause in the extratropics

    !------------------------------------------------------------------
    ! 1. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('nwp_phy_nml')
      READ(funit,NML=nwp_phy_nml)
      CALL close_tmpfile(funit)
      l_nwp_phy_namelist_read = .true.
    END IF

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_phy_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nwp_phy_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nwp_phy_nml)   ! write settings to temporary text file
      END IF
      l_nwp_phy_namelist_read = .true.
    END SELECT
    CALL close_nml

    !-----------------------
    ! 3. apply default settings where nothing is specified explicitly (except for restart)
    !-----------------------
    
    IF (.NOT. isRestart() .or. .NOT. l_nwp_phy_namelist_read) THEN

      ! 3a. Set default values for global domain where nothing at all has been specified

      ! Physics packages
      IF (inwp_gscp(1)       < 0) inwp_gscp(1)       = 1  !> 1 = hydci (COSMO-EU microphysics)
      IF (inwp_satad(1)      < 0) inwp_satad(1)      = 1  !> 1 = saturation adjustment on
      IF (inwp_convection(1) < 0) inwp_convection(1) = 1  !> 1 = Tiedtke/Bechthold convection
      IF (inwp_radiation(1)  < 0) inwp_radiation(1)  = 1  !> 1 = RRTM radiation
      IF (inwp_sso(1)        < 0) inwp_sso(1)        = 1  !> 1 = Lott and Miller scheme (COSMO)
      IF (inwp_gwd(1)        < 0) inwp_gwd(1)        = 1  !> 1 = Orr-Ern-Bechthold scheme (IFS)
      IF (inwp_cldcover(1)   < 0) inwp_cldcover(1)   = 1  !> 1 = diagnostic cloud cover (by Martin Koehler)
      IF (inwp_turb(1)       < 0) inwp_turb(1)       = icosmo  !> 1 = turbdiff (COSMO diffusion and transfer)
      IF (inwp_surface(1)    < 0) inwp_surface(1)    = 1  !> 1 = TERRA

      ! Time steps
      IF (dt_conv (1) < 0._wp) dt_conv (1) = 600._wp    !seconds
      IF (dt_sso  (1) < 0._wp) dt_sso  (1) = 1200._wp   !seconds
      IF (dt_gwd  (1) < 0._wp) dt_gwd  (1) = 1200._wp   !seconds
      IF (dt_rad  (1) < 0._wp) dt_rad  (1) = 1800._wp   !seconds

      
      ! 3b. Copy values of parent domain (in case of linear nesting) to nested domains where nothing has been specified

      DO jg = 2, max_dom

        ! Physics packages
        IF (inwp_gscp(jg)       < 0) inwp_gscp(jg)       = inwp_gscp(jg-1)
        IF (inwp_satad(jg)      < 0) inwp_satad(jg)      = inwp_satad(jg-1)
        IF (inwp_convection(jg) < 0) inwp_convection(jg) = inwp_convection(jg-1)
        IF (inwp_radiation(jg)  < 0) inwp_radiation(jg)  = inwp_radiation(jg-1)
        IF (inwp_sso(jg)        < 0) inwp_sso(jg)        = inwp_sso(jg-1)
        IF (inwp_gwd(jg)        < 0) inwp_gwd(jg)        = inwp_gwd(jg-1)
        IF (inwp_cldcover(jg)   < 0) inwp_cldcover(jg)   = inwp_cldcover(jg-1)
        IF (inwp_turb(jg)       < 0) inwp_turb(jg)       = inwp_turb(jg-1)
        IF (inwp_surface(jg)    < 0) inwp_surface(jg)    = inwp_surface(jg-1)

        ! Time steps
        IF (dt_conv (jg) < 0._wp) dt_conv (jg) = dt_conv (jg-1) 
        IF (dt_sso  (jg) < 0._wp) dt_sso  (jg) = dt_sso  (jg-1)
        IF (dt_gwd  (jg) < 0._wp) dt_gwd  (jg) = dt_gwd  (jg-1)
        IF (dt_rad  (jg) < 0._wp) dt_rad  (jg) = dt_rad  (jg-1)

      ENDDO

    ENDIF

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    
    ! check for valid parameters in namelists

    DO jg = 1, max_dom

      IF ( ALL((/0,1,2,3,4,5,6,9/) /= inwp_gscp(jg)) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting for inwp_gscp. Must be 0,1,2,3,4,5 or 9.')
      END IF
      
#ifndef __ICON_ART
    IF (inwp_gscp(jg) == 6) THEN
      CALL finish( TRIM(routine),'inwp_gscp == 6, but ICON was compiled without -D__ICON_ART')
    ENDIF
#endif
      
      IF (inwp_surface(jg) == 0 .AND. itype_z0 > 1) THEN
        CALL message(TRIM(routine), 'Warning: itype_z0 is reset to 1 because surface scheme is turned off')
        itype_z0 = 1
      ENDIF

      IF (icpl_aero_gscp > 0 .AND. inwp_gscp(jg) > 2) THEN
        CALL finish( TRIM(routine), 'Aerosol-microphysics coupling currently available only for inwp_gscp=1,2')
      ENDIF


      ! For backward compatibility, do not throw an error message, if inwp_turb=10,11 or 12 
      ! is chosen. reset inwp_turb to 1, instead
      IF ( ANY((/10,11,12/) == inwp_turb(jg)) ) THEN
        inwp_turb(jg) = 1
        WRITE(message_text,'(a,i2)') 'Reset inwp_turb to 1 for domain ', jg
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF

    ENDDO



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg=1,max_dom
      atm_phy_nwp_config(jg)%inwp_convection = inwp_convection(jg)
      atm_phy_nwp_config(jg)%inwp_cldcover   = inwp_cldcover(jg)
      atm_phy_nwp_config(jg)%inwp_radiation  = inwp_radiation(jg)
      atm_phy_nwp_config(jg)%inwp_sso        = inwp_sso(jg)
      atm_phy_nwp_config(jg)%inwp_gwd        = inwp_gwd(jg) 
      atm_phy_nwp_config(jg)%inwp_gscp       = inwp_gscp(jg)
      atm_phy_nwp_config(jg)%inwp_satad      = inwp_satad(jg)
      atm_phy_nwp_config(jg)%inwp_turb       = inwp_turb(jg)
      atm_phy_nwp_config(jg)%inwp_surface    = inwp_surface(jg)
      atm_phy_nwp_config(jg)%itype_z0        = itype_z0

      atm_phy_nwp_config(jg)%lshallowconv_only  = lshallowconv_only(jg)
      atm_phy_nwp_config(jg)%ldetrain_conv_prec = ldetrain_conv_prec(jg)

      atm_phy_nwp_config(jg)%dt_conv         = dt_conv (jg)
      atm_phy_nwp_config(jg)%dt_ccov         = dt_conv (jg)
      atm_phy_nwp_config(jg)%dt_rad          = dt_rad  (jg)
      atm_phy_nwp_config(jg)%dt_sso          = dt_sso  (jg)
      atm_phy_nwp_config(jg)%dt_gwd          = dt_gwd  (jg)
      atm_phy_nwp_config(jg)%qi0             = qi0 
      atm_phy_nwp_config(jg)%qc0             = qc0 
      atm_phy_nwp_config(jg)%ustart_raylfric = ustart_raylfric 
      atm_phy_nwp_config(jg)%efdt_min_raylfric = efdt_min_raylfric
      atm_phy_nwp_config(jg)%latm_above_top  = latm_above_top(jg)
      atm_phy_nwp_config(jg)%mu_rain         = mu_rain
      atm_phy_nwp_config(jg)%mu_snow         = mu_snow
      atm_phy_nwp_config(jg)%icpl_aero_gscp  = icpl_aero_gscp
    ENDDO

    config_lrtm_filename   = TRIM(lrtm_filename)
    config_cldopt_filename = TRIM(cldopt_filename)
    config_icpl_aero_conv  = icpl_aero_conv
    config_iprog_aero      = iprog_aero
    config_icpl_o3_tp      = icpl_o3_tp

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nwp_phy_nml)                    
      CALL store_and_close_namelist(funit, 'nwp_phy_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_phy_nml)

  END SUBROUTINE read_nwp_phy_namelist

END MODULE mo_nwp_phy_nml

