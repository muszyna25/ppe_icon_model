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

  USE mo_exception            ,ONLY: finish, message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_vertical_coord_table ,ONLY: vct_a
  USE mo_physical_constants   ,ONLY: tmelt
  USE mo_echam_phy_config     ,ONLY: echam_phy_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_cld_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_cld_config   !< allocate and initialize echam_cld_config
  PUBLIC ::    eval_echam_cld_config   !< evaluate echam_cld_config
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
     REAL(wp) :: zmaxcld  ! [m]      maximum height for cloud microphysics calculation
     INTEGER  :: jkscld   !          vertical start index for cloud microphysics calculation
     !                               diagnosed in eval_echam_cld_config
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
     ! cloud liquid/ice path ratio
     REAL(wp) :: clwprat  !          critical ratio of cloud liq.+ice paths below and above
     !                    !          the top of shallow convection
     !                    !          for ratio > clwprat -> change ktype from 2 to 4
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
    ! vertical range
    echam_cld_config(:)% zmaxcld  = echam_phy_config(:)% zmaxcloudy
    !
    ! general thresholds
    echam_cld_config(:)% ccwmin   = 1.e-7_wp
    echam_cld_config(:)% cqtmin   = 1.e-12_wp
    !
    ! freezing/deposition/sublimation
    echam_cld_config(:)% cthomi   = tmelt-35.0_wp
    echam_cld_config(:)% csecfrl  = 1.5e-5_wp
    !
    ! warm clouds
    echam_cld_config(:)% ccraut   =  2.0_wp
    echam_cld_config(:)% ccracl   =  6.0_wp
    echam_cld_config(:)% cauloc   =  1.0_wp
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
    echam_cld_config(:)% ccsaut   =  2.0_wp
    echam_cld_config(:)% ccsacl   = 0.10_wp
    !
    ! cloud liquid/ice path ratio
    echam_cld_config(:)% clwprat  = 4.0_wp
    !
  END SUBROUTINE init_echam_cld_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_cld_config
    !
    INTEGER           :: jg, jk, klev
    CHARACTER(LEN=2)  :: cg
    !
    klev = SIZE(vct_a)-1
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       ! diagnose jkscld
       echam_cld_config(jg)% jkscld = 1
       !
       DO jk = 1,klev
          !
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > echam_cld_config(jg)% zmaxcld) THEN
             echam_cld_config(jg)% jkscld = echam_cld_config(jg)% jkscld + 1
          ELSE
             EXIT
          END IF
          !
       END DO
       !
    END DO
    !
  END SUBROUTINE eval_echam_cld_config

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
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% zmaxcld  ',echam_cld_config(jg)% zmaxcld )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% jkscld   ',echam_cld_config(jg)% jkscld  )
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
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% cn0s     ',echam_cld_config(jg)% cn0s    )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccsaut   ',echam_cld_config(jg)% ccsaut  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% ccsacl   ',echam_cld_config(jg)% ccsacl  )
       CALL print_value('    echam_cld_config('//TRIM(cg)//')% clwprat  ',echam_cld_config(jg)% clwprat )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_cld_config

  !----

END MODULE mo_echam_cld_config
