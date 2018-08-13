!>
!! Configuration of the parameterization for ...,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI (2017-04)
!!
!! Based on earlier codes of:
!!     ...
!!
!! References: 
!!     ...
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_rad_config

  USE mo_exception            ,ONLY: message, message_text, print_value, finish
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_run_config           ,ONLY: iqt, ico2, io3, ntracer, lart

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_rad_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_rad_config   !< allocate and initialize echam_rad_config
  PUBLIC ::    eval_echam_rad_config   !< evaluate echam_rad_config
  PUBLIC ::   print_echam_rad_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_rad'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the MPI physics package
  !!
  TYPE t_echam_rad_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! -- Spectral solar irradiance at 1 AU distance from the sun
     !
     INTEGER  :: isolrad        !< solar spectrum, constant or variable in time
     REAL(wp) :: fsolrad        !< scaling factor for solar irradiance
     !
     ! -- Orbit shape
     !
     LOGICAL  :: l_orbvsop87    !< .TRUE. for VSOP87 orbit, .FALSE. for Kepler orbit
     REAL(wp) :: cecc           !< Eccentricity  of  orbit
     REAL(wp) :: cobld          !< Obliquity of Earth axis [Deg]
     !
     ! -- Orbit motion
     !
     LOGICAL  :: lyr_perp       !< .FALSE.: transient, following vsop87
     !                          !  .TRUE. : orbit of year yr_perp of the vsop87 orbit is perpetuated
     INTEGER  :: yr_perp        !< year of vsop87 orbit to be perpetuated
     INTEGER  :: nmonth         !< nmonth=0    : Earth circles on orbit
     !                          !  nmonth=1-12 : Earth position on orbit fixed for month i
     LOGICAL  :: ldiur          !< .TRUE. : with diurnal cycle
     !                          !< .FALSE.: zonally averaged irradiation
     LOGICAL  :: l_sph_symm_irr !< .TRUE. for globally averaged irradiation (RCE)
     !                          !< .FALSE. for lat (lon) dependent irradiation
     !
     !                          ! PROVISIONAL - ONLY BEST METHOD WILL BE KEPT ("0" or "3")
     INTEGER  :: icosmu0        !< selects method for the definition of cosmu0_rt in the extended
     !                          !  sunlit area, as needed if solar fluxes are adjusted to the
     !                          !  current time between radiation time steps.
     !                          !  0: no adjustment, the original cosmu0 is used for the rad. transfer
     !                          !     Has small effects on land temperture (less smooth intraday time series)
     !                          !  1: MAX(0.1,cosmu0), as used in ECHAM6 and icon-aes-1.0 and -1.1.
     !                          !  2: (cosmu0+dcosmu0)/(1+dcosmu0), dcosmu0 = SIN(dmu0), dmu0=pi*dt_rad/86400s
     !                          !     DO NOT USE! Strong effects on MA temp. and wind and land surface temp.
     !                          !  3: 0.5*SIN(dmu0)*(1+(pi/2-mu0)/dmu0), dmu0=pi*dt_rad/86400s
     !                          !     Has small effects on the MA temp. and wind and the land surface temp.
     !                          !  4: sin(mu0s)*(pi/2+dmu0-mu0), , dmu0=pi*dt_rad/86400s, mu0s = tangent point
     !                          !     Has moderate effects on the MA temp. and wind and the land surface temp.
     !
     ! --- Radiative agents
     !     irad_x=0 : radiation uses tracer x = 0
     !     irad_x=1 : radiation uses tracer x from a tracer variable
     !     irad_x>1 : radiation uses tracer x following external specifications of various kinds:
     !                - globally constant  or spatially varying
     !                - constant in time, constant annual cycle, or transient
     !
     INTEGER  :: irad_h2o       !< water vapor, clouds and ice for radiation
     INTEGER  :: irad_co2       !< CO2
     INTEGER  :: irad_ch4       !< CH4
     INTEGER  :: irad_n2o       !< N2O
     INTEGER  :: irad_o3        !< O3
     INTEGER  :: irad_o2        !< O2
     INTEGER  :: irad_cfc11     !< CFC 11
     INTEGER  :: irad_cfc12     !< CFC 12
     INTEGER  :: irad_aero      !< aerosols
     !
     ! --- Greenhouse gases scenario (read from file)
     !     ighg = 0 : select default gas volume mixing ratios - 1990 values (CMIP5)
     !     ighg = 1 : transient CMIP5 scenario from file
     !
     INTEGER  :: ighg
     !
     ! --- Volume mixing ratios - 1990 values (CMIP5)
     !
     REAL(wp) :: vmr_co2
     REAL(wp) :: vmr_n2o
     REAL(wp) :: vmr_o2
     REAL(wp) :: vmr_ch4
     REAL(wp) :: vmr_cfc11
     REAL(wp) :: vmr_cfc12
     !
     ! --- Scaling factor for mixing ratios
     !
     REAL(wp) :: frad_h2o
     REAL(wp) :: frad_co2
     REAL(wp) :: frad_n2o
     REAL(wp) :: frad_o3
     REAL(wp) :: frad_o2
     REAL(wp) :: frad_ch4
     REAL(wp) :: frad_cfc
     !
  END TYPE t_echam_rad_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_rad_config), TARGET :: echam_rad_config(max_dom)

CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_rad_config
    !
    ! ECHAM radiation configuration
    ! -----------------------------
    !
    echam_rad_config(:)% isolrad        = 0
    echam_rad_config(:)% fsolrad        = 1.0_wp
    !
    echam_rad_config(:)% l_orbvsop87    = .TRUE.
    echam_rad_config(:)% cecc           =  0.016715_wp
    echam_rad_config(:)% cobld          =  23.44100_wp
    !
    echam_rad_config(:)% lyr_perp       = .FALSE.
    echam_rad_config(:)% yr_perp        = -99999
    echam_rad_config(:)% nmonth         =  0   
    echam_rad_config(:)% ldiur          = .TRUE.
    echam_rad_config(:)% l_sph_symm_irr = .FALSE.
    !
    echam_rad_config(:)% icosmu0        = 3
    !
    echam_rad_config(:)% irad_h2o       = 1
    echam_rad_config(:)% irad_co2       = 2
    echam_rad_config(:)% irad_ch4       = 3
    echam_rad_config(:)% irad_n2o       = 3
    echam_rad_config(:)% irad_o3        = 0
    echam_rad_config(:)% irad_o2        = 2
    echam_rad_config(:)% irad_cfc11     = 2
    echam_rad_config(:)% irad_cfc12     = 2
    echam_rad_config(:)% irad_aero      = 2
    !
    echam_rad_config(:)% ighg           = 0
    !
    ! Default volume mixing ratios: 1990 values (CMIP5)
    echam_rad_config(:)% vmr_co2        =  348.0e-06_wp
    echam_rad_config(:)% vmr_ch4        = 1650.0e-09_wp
    echam_rad_config(:)% vmr_n2o        =  306.0e-09_wp
    echam_rad_config(:)% vmr_o2         =    0.20946_wp
    echam_rad_config(:)% vmr_cfc11      =  214.5e-12_wp
    echam_rad_config(:)% vmr_cfc12      =  371.1e-12_wp
    !
    echam_rad_config(:)% frad_h2o       = 1.0_wp
    echam_rad_config(:)% frad_co2       = 1.0_wp
    echam_rad_config(:)% frad_ch4       = 1.0_wp
    echam_rad_config(:)% frad_n2o       = 1.0_wp
    echam_rad_config(:)% frad_o3        = 1.0_wp
    echam_rad_config(:)% frad_o2        = 1.0_wp
    echam_rad_config(:)% frad_cfc       = 1.0_wp
    !
  END SUBROUTINE init_echam_rad_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_rad_config
    !
    CHARACTER(LEN=*), PARAMETER :: routine = 'eval_echam_rad_config'

    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg

    ! Shortcuts to components of echam_rad_config
    !
    LOGICAL , POINTER :: l_orbvsop87, l_sph_symm_irr, ldiur, lyr_perp
    INTEGER , POINTER :: isolrad, nmonth, yr_perp
    REAL(wp), POINTER :: fsolrad, cecc, cobld
    INTEGER , POINTER :: irad_h2o, irad_co2, irad_ch4, irad_n2o, irad_o3, irad_o2, irad_cfc11, irad_cfc12, irad_aero
    REAL(wp), POINTER ::            vmr_co2,  vmr_ch4,  vmr_n2o,           vmr_o2,  vmr_cfc11,  vmr_cfc12
    REAL(wp), POINTER :: frad_h2o, frad_co2, frad_ch4, frad_n2o, frad_o3, frad_o2, frad_cfc

    CALL message    ('','')
    CALL message    ('','------------------------------------------------------------------------')
    CALL message    ('','')
    CALL message    ('','Effective input to the radiation')
    CALL message    ('','================================')
    CALL message    ('','')

    DO jg = 1,n_dom
       !
       isolrad    => echam_rad_config(jg)% isolrad
       fsolrad    => echam_rad_config(jg)% fsolrad
       !
       l_orbvsop87=> echam_rad_config(jg)% l_orbvsop87
       cecc       => echam_rad_config(jg)% cecc
       cobld      => echam_rad_config(jg)% cobld
       !
       lyr_perp   => echam_rad_config(jg)% lyr_perp
       yr_perp    => echam_rad_config(jg)% yr_perp
       nmonth     => echam_rad_config(jg)% nmonth
       !
       ldiur      => echam_rad_config(jg)% ldiur
       l_sph_symm_irr => echam_rad_config(jg)% l_sph_symm_irr
       !
       irad_h2o   => echam_rad_config(jg)% irad_h2o
       irad_co2   => echam_rad_config(jg)% irad_co2
       irad_ch4   => echam_rad_config(jg)% irad_ch4
       irad_n2o   => echam_rad_config(jg)% irad_n2o
       irad_o3    => echam_rad_config(jg)% irad_o3
       irad_o2    => echam_rad_config(jg)% irad_o2
       irad_cfc11 => echam_rad_config(jg)% irad_cfc11
       irad_cfc12 => echam_rad_config(jg)% irad_cfc12
       irad_aero  => echam_rad_config(jg)% irad_aero
       !
       vmr_co2    => echam_rad_config(jg)% vmr_co2
       vmr_ch4    => echam_rad_config(jg)% vmr_ch4
       vmr_n2o    => echam_rad_config(jg)% vmr_n2o
       vmr_o2     => echam_rad_config(jg)% vmr_o2
       vmr_cfc11  => echam_rad_config(jg)% vmr_cfc11
       vmr_cfc12  => echam_rad_config(jg)% vmr_cfc12
       !
       frad_h2o   => echam_rad_config(jg)% frad_h2o
       frad_co2   => echam_rad_config(jg)% frad_co2
       frad_ch4   => echam_rad_config(jg)% frad_ch4
       frad_n2o   => echam_rad_config(jg)% frad_n2o
       frad_o3    => echam_rad_config(jg)% frad_o3
       frad_o2    => echam_rad_config(jg)% frad_o2
       frad_cfc   => echam_rad_config(jg)% frad_cfc
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message   ('','For domain '//cg)
       CALL message   ('','------------')
       CALL message   ('','')
       CALL message   ('','Solar spectrum')
       CALL message   ('','--------------')
       CALL message   ('','')
       !
       SELECT CASE (isolrad)
       CASE (0) 
          CALL message('','SRTM default solar spectrum')
       CASE (1) 
          CALL message('','Time dependent solar spectrum from file')
       CASE (2) 
          CALL message('','Average 1844-1856 of transient CMIP5 solar')
       CASE (3) 
          CALL message('','Average 1979-1988 of transient CMIP5 solar spectrum')
       CASE (4)
          CALL message('','Solar flux for RCE simulations with diurnal cycle')
       CASE (5)
          CALL message('','Solar flux for RCE simulations without diurnal cycle')
       CASE (6)
          CALL message('','Average 1850-1873 of transient CMIP6 solar')
       CASE default 
          WRITE (message_text, '(a,i0,a)') &
               'ERROR: isolrad = ', isolrad, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       IF (fsolrad >= 0.0_wp) THEN
          CALL print_value('Factor applied to solar spectrum',fsolrad)
       ELSE
          CALL finish(routine,'ERROR: Negative fsolrad is not allowed')
       END IF
       !
       CALL message   ('','')
       CALL message   ('','Earth orbit')
       CALL message   ('','-----------')
       !
       IF (l_orbvsop87) THEN
          CALL message('','The VSOP87 orbit is used')
          IF (lyr_perp) THEN
             CALL print_value('The VSOP87 orbit is perpetuated for the year',yr_perp)
          END IF
       ELSE
          CALL message('','The Kepler orbit is used')
          CALL print_value('Eccentricity of Earth orbit   =',cecc)
          CALL print_value('Obliquity of Earth axis [Deg] =',cobld)
       END IF
       !
       SELECT CASE (nmonth)
       CASE(0)
          CALL message('','Earth cycles on orbit --> annual cycle in solar irradiation at TOA')
       CASE(1:12)
          CALL print_value('Earth rests at orbit position of the mid-time of month', nmonth)
       CASE default
          WRITE (message_text, '(a,i0,a)') &
               'ERROR: nmonth=', nmonth, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       CALL message   ('','')
       CALL message   ('','Solar irradiation at the top of the atmosphere')
       CALL message   ('','----------------------------------------------')
       !
       IF (.NOT. ldiur) THEN
          CALL message('','Zonally averaged solar irradiation at TOA --> diurnal cycle off')
       ENDIF
       !
       IF (l_sph_symm_irr) THEN
          CALL message('','Globally uniform irradiation')
       ENDIF
       !
       CALL message   ('','')
       CALL message   ('','Sources of volume/mass mixing ratios used in radiation')
       CALL message   ('','------------------------------------------------------')
       !
       ! --- Check  H2O
       !
       SELECT CASE (irad_h2o)
       CASE(0)
          CALL message('','No H2O(gas,liq,ice) in radiation')
       CASE(1)
          CALL message('','H2O   (gas,liq,ice) mass mixing ratios from tracer fields')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_h2o   =', irad_h2o, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check  CO2
       ! 
       SELECT CASE (irad_co2)
       CASE(0)
          CALL message('','No CO2 in radiation')
       CASE(1)
          IF ( iqt <= ico2 .AND. ico2 <= ntracer .AND. .NOT. lart) THEN
             CALL message('','CO2   mass mixing ratio from tracer field')
          ELSE
             CALL finish(routine,'ERROR: irad_co2 = 1 (CO2 tracer in radiation) is not '// &
                  &              'a valid choice because no CO2 tracer is available')
          END IF
       CASE(2)
          WRITE (message_text, '(a,e16.8)') &
               'CO2   volume mixing ratio =', vmr_co2
          CALL message('',message_text)
       CASE(4)
          CALL message('','CO2   volume mixing ratio from ghg scenario file')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_co2   = ', irad_co2, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check CH4
       ! 
       SELECT CASE (irad_ch4)
       CASE(0)
          CALL message('','No CH4 in radiation')
       CASE(2)
          WRITE (message_text, '(a,e16.8)') &
                          'CH4   volume mixing ratio =', vmr_ch4
          CALL message('',message_text)
       CASE(3)
          WRITE (message_text, '(a,e16.8)') &
                          'CH4   tanh-profile with surface volume mixing ratio =', vmr_ch4
          CALL message('',message_text)
       CASE(4)
          CALL message('','CH4   tanh-profile with surface volume mixing ratio from ghg scenario file')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_ch4   =', irad_ch4, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check N2O
       ! 
       SELECT CASE (irad_n2o)
       CASE(0)
          CALL message('','No N2O in radiation')
       CASE(2)
          WRITE (message_text, '(a,e16.8)') &
                          'N2O   volume mixing ratio =', vmr_n2o
          CALL message('',message_text)
       CASE(3)
          WRITE (message_text, '(a,e16.8)') &
                          'N2O   tanh-profile with surface volume mixing ratio from echam_rad_config =', vmr_n2o
          CALL message('',message_text)
       CASE(4)
          CALL message('','N2O   tanh-profile with surface volume mixing ratio from ghg scenario file')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_n2o   =',irad_n2o,' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check CFCs
       ! 
       SELECT CASE (irad_cfc11)
       CASE(0)
          CALL message('','No CFC11 in radiation')
       CASE(2)
          WRITE (message_text, '(a,e16.8)') &
                          'CFC11 volume mixing ratio =', vmr_cfc11
          CALL message('',message_text)
       CASE(4)
          CALL message('','CFC11 volume mixing ratio from ghg scenario file')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_cfc11 =', irad_cfc11, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT

       SELECT CASE (irad_cfc12)
       CASE(0)
          CALL message('','No CFC12 in radiation')
       CASE(2)
          WRITE (message_text, '(a,e16.8)') &
                          'CFC12 volume mixing ratio =', vmr_cfc12
          CALL message('',message_text)
       CASE(4)
          CALL message('','CFC12 volume mixing ratio from ghg scenario file')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_cfc12 =', irad_cfc12, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check O3
       ! 
       SELECT CASE (irad_o3)
       CASE(0)
          CALL message('','No O3 in radiation')
       CASE(1)
          IF ( iqt <= io3  .AND. io3  <= ntracer .AND. .NOT. lart) THEN
             CALL message('','O3    mass mixing ratio from tracer field')
          ELSE
             CALL finish(routine,'ERROR: irad_o3  = 1 (O3  tracer in radiation) is not '// &
                  &              'a valid choice because no O3  tracer is available '   // &
                  &              ' or lart = .TRUE.')
          END IF
       CASE(2)
          CALL message('','O3    clim. annual cycle 3-dim. volume mixing ratio from file')
       CASE(4)
          CALL message('','O3    constant-in-time 3-dim. volume mixing ratio from file')
       CASE(8)
          CALL message('','O3    transient 3-dim. volume mixing ratio from file')
       CASE(10)
          CALL message('','O3    from ART')
       CASE default
          WRITE (message_text, '(a,i0,a)') 'ERROR: irad_o3    =', irad_o3, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check O2
       ! 
       SELECT CASE (irad_o2)
       CASE(0)
          CALL message('','No O2  in radiation')
       CASE(2)
          WRITE (message_text, '(a,e16.8)') &
               'O2    volume mixing ratio =', vmr_o2
          CALL message('',message_text)
       CASE default
          WRITE (message_text, '(a,i0,a)') &
               'ERROR: irad_o2    =', irad_o2, ' is not supported'
          CALL finish(routine,message_text)
       END SELECT
       !
       ! --- Check aerosol
       ! 
       SELECT CASE (irad_aero)
       CASE(0)
          CALL message('','No aerosol in radiation')
       END SELECT
       !
       !
       ! --- Check scaling factors
       !
       CALL message('','')
       CALL message('','Multiplication factors applied in radiation to vol./mass mixing ratio sources')
       CALL message('','-----------------------------------------------------------------------------')
       !
       IF (frad_h2o >= 0.0_wp) THEN
          CALL print_value('H2O(gas,liq,ice): frad_h2o =',frad_h2o)
       ELSE
          CALL finish(routine,'ERROR: Negative frad_h2o is not allowed')
       END IF
       !
       IF (frad_co2 >= 0.0_wp) THEN
          CALL print_value('CO2             : frad_co2 =',frad_co2)
       ELSE
          CALL finish(routine,'ERROR: Negative frad_co2 is not allowed')
       END IF
       !
       IF (frad_ch4 >= 0.0_wp) THEN
          CALL print_value('CH4             : frad_ch4 =',frad_ch4)
       ELSE
          CALL finish(routine,'ERROR: Negative frad_ch4 is not allowed')
       END IF
       !
       IF (frad_n2o >= 0.0_wp) THEN
          CALL print_value('N2O             : frad_n2o =',frad_n2o)
       ELSE
          CALL finish(routine,'ERROR: Negative frad_n2o is not allowed')
       END IF
       !
       IF (frad_o3  >= 0.0_wp) THEN
          CALL print_value('O3              : frad_o3  =',frad_o3 )
       ELSE
          CALL finish(routine,'ERROR: Negative frad_o3  is not allowed')
       END IF
       !
       IF (frad_o2  >= 0.0_wp) THEN
          CALL print_value('O2              : frad_o2  =',frad_o2 )
       ELSE
          CALL finish(routine,'ERROR: Negative frad_o2  is not allowed')
       END IF
       !
       IF (frad_cfc >= 0.0_wp) THEN
          CALL print_value('CFC11 and CFC12 : frad_cfc =',frad_cfc)
       ELSE
          CALL finish(routine,'ERROR: Negative frad_cfc is not allowed')
       END IF
       !
       CALL message('','')
       !
    END DO ! jg loop
    !
    CALL message('','')
    !
  END SUBROUTINE eval_echam_rad_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_rad_config

    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg

    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM radiation configuration')
    CALL message    ('','=============================')
    CALL message    ('','')

    IF (n_dom > 1) THEN
       CALL message    ('','-----------------------------------------------------------------------------------')
       CALL message    ('','!! WARNING: The current radiation uses the the setup of domain 1 for all domains !!')
       CALL message    ('','-----------------------------------------------------------------------------------')
    END IF
    
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% isolrad       ',echam_rad_config(jg)% isolrad       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% fsolrad       ',echam_rad_config(jg)% fsolrad       )
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% l_orbvsop87   ',echam_rad_config(jg)% l_orbvsop87   )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% cecc          ',echam_rad_config(jg)% cecc          )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% cobld         ',echam_rad_config(jg)% cobld         )
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% lyr_perp      ',echam_rad_config(jg)% lyr_perp      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% yr_perp       ',echam_rad_config(jg)% yr_perp       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% nmonth        ',echam_rad_config(jg)% nmonth        )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% ldiur         ',echam_rad_config(jg)% ldiur         )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% l_sph_symm_irr',echam_rad_config(jg)% l_sph_symm_irr)
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_h2o      ',echam_rad_config(jg)% irad_h2o      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_co2      ',echam_rad_config(jg)% irad_co2      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_ch4      ',echam_rad_config(jg)% irad_ch4      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_n2o      ',echam_rad_config(jg)% irad_n2o      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_o3       ',echam_rad_config(jg)% irad_o3       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_o2       ',echam_rad_config(jg)% irad_o2       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_cfc11    ',echam_rad_config(jg)% irad_cfc11    )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_cfc12    ',echam_rad_config(jg)% irad_cfc12    )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% irad_aero     ',echam_rad_config(jg)% irad_aero     )
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% ighg          ',echam_rad_config(jg)% ighg          )
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% vmr_co2       ',echam_rad_config(jg)% vmr_co2       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% vmr_ch4       ',echam_rad_config(jg)% vmr_ch4       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% vmr_n2o       ',echam_rad_config(jg)% vmr_n2o       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% vmr_o2        ',echam_rad_config(jg)% vmr_o2        )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% vmr_cfc11     ',echam_rad_config(jg)% vmr_cfc11     )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% vmr_cfc12     ',echam_rad_config(jg)% vmr_cfc12     )
       CALL message    ('','')
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_h2o      ',echam_rad_config(jg)% frad_h2o      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_co2      ',echam_rad_config(jg)% frad_co2      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_ch4      ',echam_rad_config(jg)% frad_ch4      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_n2o      ',echam_rad_config(jg)% frad_n2o      )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_o3       ',echam_rad_config(jg)% frad_o3       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_o2       ',echam_rad_config(jg)% frad_o2       )
       CALL print_value('    echam_rad_config('//TRIM(cg)//')% frad_cfc      ',echam_rad_config(jg)% frad_cfc      )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_rad_config

  !----

END MODULE mo_echam_rad_config
