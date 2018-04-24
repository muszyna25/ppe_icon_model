!>
!! Configuration of the parameterization for cloud microphysics,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI-M (2017-12)
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
MODULE mo_echam_cld_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_physical_constants   ,ONLY: tmelt

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_cld_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_cld_config   !< allocate and initialize echam_cld_config
!!$  PUBLIC ::    eval_echam_cld_config   !< evaluate echam_cld_config
  PUBLIC ::   print_echam_cld_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_cld'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM microphysics
  !!
  TYPE t_echam_cld_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! vertical loop limit
     INTEGER  :: jks      !          start index of jk loops in "cover" and "cloud"
     !
     ! general thresholds
     REAL(wp) :: ccwmin   ! [kg/kg]  cloud water and ice minimum mass mixing ratio for cover>0
     REAL(wp) :: cqtmin   ! [kg/kg]  cloud water/ice minimum for microphysical processes
     !
     ! freezing/deposition/sublimation
     REAL(wp) :: cthomi   ! [K]      maximum temperature for homogeneous freezing
     REAL(wp) :: csecfrl  ! [kg/kg]  minimum in-cloud water mass mixing ratio in mixed phase clouds
     !
     ! warm clouds
     REAL(wp) :: ccraut   !          coefficient of autoconversion of cloud droplets to rain
     REAL(wp) :: ccracl   !          coefficient of accretion  of cloud droplets by falling rain
     REAL(wp) :: cauloc   !          coefficient of local rainwater production by autoconversion
     REAL(wp) :: clmin    !          minimum for cauloc*dz/5000
     REAL(wp) :: clmax    !          maximum for cauloc*dz/5000
     !
     ! cold clouds
     REAL(wp) :: cvtfall  !          coefficient of sedimentation velocity of cloud ice
     REAL(wp) :: ceffmin  ! [1e-6m]  min effective radius for ice cloud
     REAL(wp) :: ceffmax  ! [1e-6m]  max effective radius for ice cloud
     REAL(wp) :: crhoi    ! [kg/kg]  density of cloud ice
     REAL(wp) :: crhosno  ! [kg/kg]  bulk density of snow
     REAL(wp) :: cn0s     !          intercept parameter for snow size distribution
     REAL(wp) :: ccsaut   !          coefficient of autoconversion of cloud ice to snow
     REAL(wp) :: ccsacl   !          coefficient of accretion of cloud droplets by falling snow
     !
     ! cloud droplet number concentration
     REAL(wp) :: cn1lnd   ! [1e6/m3] over land, p <= 100 hPa
     REAL(wp) :: cn2lnd   ! [1e6/m3] over land, p >= 800 hPa
     REAL(wp) :: cn1sea   ! [1e6/m3] over sea , p <= 100 hPa
     REAL(wp) :: cn2sea   ! [1e6/m3] over sea , p >= 800 hPa
     !
     ! cloud optics
     REAL(wp) :: cinhomi  !          ice    cloud inhomogeneity factor
     !
     !                    !          liquid cloud inhomogeneity factor:
     REAL(wp) :: cinhoml1 !          - ktype = 0 = stratiform clouds
     REAL(wp) :: cinhoml2 !          - ktype = 4 = shallow conv. (cf. clwprat)
     REAL(wp) :: cinhoml3 !          - ktype = 1 = deep convection             and
     !                    !            ktype = 2 = shallow conv. (cf. clwprat) and
     !                    !            ktype = 3 = mid-level conv.
     REAL(wp) :: clwprat  !          critical ratio of cloud liq.+ice paths below and above
     !                    !          the top of shallow convection
     !                    !          for ratio > clwprat -> change ktype from 2 to 4
     !
     ! cloud cover
     REAL(wp) :: crs      !          critical relative humidity at surface
     REAL(wp) :: crt      !          critical relative humidity aloft
     INTEGER  :: nex      !          transition parameter for critical relative humidity profile
     INTEGER  :: jbmin    !          index of highest level for search of top level of inversion layer over sea (ca. 2 km)
     INTEGER  :: jbmax    !          index of bottom level of inversion layer over sea
     REAL(wp) :: cinv     !          fraction of dry adiabatic lapse rate for search of top level of inversion layer over sea
     REAL(wp) :: csatsc   !          minimum effective saturation for cloud cover below an invesion layer over sea
     !
     ! tropopause diagnostics
!!$     REAL(wp) :: cptop    ! [Pa]     pressure of highest level for tropopause calculation
!!$     REAL(wp) :: cpbot    ! [Pa]     pressure of lowest  level for tropopause calculation
     INTEGER  :: ncctop   !          index of highest level for tropopause calculation
     INTEGER  :: nccbot   !          index of lowest  level for tropopause calculation
     !
  END TYPE t_echam_cld_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_cld_config), TARGET :: echam_cld_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_cld_config
    !
    ! ECHAM cloud microphyiscs configuration
    ! --------------------------------------
    !
    ! general thresholds
    echam_cld_config(:)% jks      = 15             ! L47: jks=15 is at ca. 30km
    echam_cld_config(:)% ccwmin   = 1.e-7_wp
    echam_cld_config(:)% cqtmin   = 1.e-12_wp
    !
    ! freezing/deposition/sublimation
    echam_cld_config(:)% cthomi   = tmelt-35.0_wp
    echam_cld_config(:)% csecfrl  = 5.0e-6_wp
    !
    ! warm clouds
    echam_cld_config(:)% ccraut   = 15.0_wp
    echam_cld_config(:)% ccracl   =  6.0_wp
    echam_cld_config(:)% cauloc   = 10.0_wp
    echam_cld_config(:)% clmin    = 0.0_wp
    echam_cld_config(:)% clmax    = 0.5_wp
    !
    ! cold clouds
    echam_cld_config(:)% cvtfall  = 2.5_wp
    echam_cld_config(:)% ceffmin  = 10.0_wp
    echam_cld_config(:)% ceffmax  = 150.0_wp
    echam_cld_config(:)% crhoi    = 500.0_wp
    echam_cld_config(:)% crhosno  = 100.0_wp
    echam_cld_config(:)% cn0s     = 3.e6_wp
    echam_cld_config(:)% ccsaut   = 95.0_wp
    echam_cld_config(:)% ccsacl   = 0.10_wp
    !
    ! cloud droplet number concentration
    echam_cld_config(:)% cn1lnd   =  20._wp
    echam_cld_config(:)% cn2lnd   = 180._wp
    echam_cld_config(:)% cn1sea   =  20._wp
    echam_cld_config(:)% cn2sea   =  80._wp
    !
    ! cloud optics
    echam_cld_config(:)% cinhomi  = 0.80_wp
    echam_cld_config(:)% cinhoml1 = 0.80_wp
    echam_cld_config(:)% cinhoml2 = 0.40_wp
    echam_cld_config(:)% cinhoml3 = 0.80_wp
    echam_cld_config(:)% clwprat  = 4.0_wp
    !
    ! cloud cover
    echam_cld_config(:)% crs      = 0.968_wp
    echam_cld_config(:)% crt      = 0.8_wp
    echam_cld_config(:)% nex      = 2
    echam_cld_config(:)% jbmin    = 40
    echam_cld_config(:)% jbmax    = 45
    echam_cld_config(:)% cinv     = 0.25_wp
    echam_cld_config(:)% csatsc   = 0.7_wp
    !
    ! tropopause diagnostics
!!$    echam_cld_config(:)% cptop  = 1000.0_wp
!!$    echam_cld_config(:)% cpbot  = 50000.0_wp
    echam_cld_config(:)% ncctop   = 13
    echam_cld_config(:)% nccbot   = 35
    !
  END SUBROUTINE init_echam_cld_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_echam_cld_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_echam_cld_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_cld_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM cloud microphyiscs configuration')
    CALL message    ('','======================================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% jks      ',echam_cld_config(jg)% jks     )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccwmin   ',echam_cld_config(jg)% ccwmin  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cqtmin   ',echam_cld_config(jg)% cqtmin  )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cthomi   ',echam_cld_config(jg)% cthomi  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% csecfrl  ',echam_cld_config(jg)% csecfrl )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccraut   ',echam_cld_config(jg)% ccraut  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccracl   ',echam_cld_config(jg)% ccracl  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cauloc   ',echam_cld_config(jg)% cauloc  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% clmin    ',echam_cld_config(jg)% clmin   )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% clmax    ',echam_cld_config(jg)% clmax   )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cvtfall  ',echam_cld_config(jg)% cvtfall )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ceffmin  ',echam_cld_config(jg)% ceffmin )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ceffmax  ',echam_cld_config(jg)% ceffmax )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% crhoi    ',echam_cld_config(jg)% crhoi   )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% crhosno  ',echam_cld_config(jg)% crhosno )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cn0s     ',echam_cld_config(jg)% cn0s  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccsaut   ',echam_cld_config(jg)% ccsaut  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccsacl   ',echam_cld_config(jg)% ccsacl  )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cn1lnd   ',echam_cld_config(jg)% cn1lnd  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cn2lnd   ',echam_cld_config(jg)% cn2lnd  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cn1sea   ',echam_cld_config(jg)% cn1sea  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cn2sea   ',echam_cld_config(jg)% cn2sea  )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cinhomi  ',echam_cld_config(jg)% cinhomi )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cinhoml1 ',echam_cld_config(jg)% cinhoml1)
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cinhoml2 ',echam_cld_config(jg)% cinhoml2)
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cinhoml3 ',echam_cld_config(jg)% cinhoml3)
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% clwprat  ',echam_cld_config(jg)% clwprat )
       CALL message    ('','')
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% crs      ',echam_cld_config(jg)% crs     )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% crt      ',echam_cld_config(jg)% crt     )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% nex      ',echam_cld_config(jg)% nex     )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% jbmin    ',echam_cld_config(jg)% jbmin   )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% jbmax    ',echam_cld_config(jg)% jbmax   )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cinv     ',echam_cld_config(jg)% cinv    )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% csatsc   ',echam_cld_config(jg)% csatsc  )
       CALL message    ('','')
!!$       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cptop    ',echam_cld_config(jg)% cptop   )
!!$       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cpbot    ',echam_cld_config(jg)% cpbot   )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ncctop   ',echam_cld_config(jg)% ncctop  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% nccbot   ',echam_cld_config(jg)% nccbot  )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_cld_config

  !----

END MODULE mo_echam_cld_config
