!>
!! Configuration of the parameterization for cloud cover,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI-M (2019-06)
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
MODULE mo_echam_cov_config

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
  PUBLIC ::         echam_cov_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_cov_config   !< allocate and initialize echam_cov_config
  PUBLIC ::    eval_echam_cov_config   !< evaluate echam_cov_config
  PUBLIC ::   print_echam_cov_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_cov'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM microphysics
  !!
  TYPE t_echam_cov_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! vertical loop limit
     REAL(wp) :: zmaxcov  !          maximum height (m) for cloud cover calculation
     INTEGER  :: jkscov   !          vertical start index for cloud cover calculation
     !                               diagnosed in eval_echam_cov_config
     ! cloud cover
     INTEGER  :: icov     !          cloud cover scheme
     !                               1: fractional cloud cover dependent on relative humidity
     !                               2:     0/1    cloud cover dependent on relative humidity
     ! icov=1 and icov=2:
     REAL(wp) :: csat     !          relative humidity for 100% cloud cover
     !
     ! icov=1 only:
     REAL(wp) :: crs      !          critical relative humidity at surface
     REAL(wp) :: crt      !          critical relative humidity aloft
     INTEGER  :: nex      !          transition parameter for critical relative humidity profile
     REAL(wp) :: zmaxinv  !          maximum height (m) above sea level for search of inversion layer
     INTEGER  :: jksinv   !          vertical start index for search of inversion layer
     !                               diagnosed in eval_echam_cov_config
     REAL(wp) :: zmininv  !          minimum height (m) above sea level for search of inversion layer
     INTEGER  :: jkeinv   !          vertical end index for search of inversion layer
     !                               diagnosed in eval_echam_cov_config
     REAL(wp) :: cinv     !          fraction of dry adiabatic lapse rate for search of top level of inversion layer over sea
     REAL(wp) :: csatsc   !          lower limit of scaling factor for saturation mixing ratio in layer below inversion
     !                               (csatsc=1 defaults to the standard scheme without accounting of inversion layers)
     !
  END TYPE t_echam_cov_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_cov_config), TARGET :: echam_cov_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_cov_config
    !
    ! ECHAM cloud cover configuration
    ! -------------------------------
    !
    echam_cov_config(:)% zmaxcov  = echam_phy_config(:)% zmaxcloudy
    echam_cov_config(:)% icov     = 1
    echam_cov_config(:)% csat     = 1.0_wp
    echam_cov_config(:)% crs      = 0.968_wp
    echam_cov_config(:)% crt      = 0.8_wp
    echam_cov_config(:)% nex      = 2
    echam_cov_config(:)% zmaxinv  = 2000.0_wp
    echam_cov_config(:)% zmininv  =  200.0_wp
    echam_cov_config(:)% cinv     = 0.25_wp
    echam_cov_config(:)% csatsc   = 0.7_wp
    !
  END SUBROUTINE init_echam_cov_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_cov_config
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
       ! diagnose jkscov
       echam_cov_config(jg)% jkscov = 1
       DO jk = 1,klev
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > echam_cov_config(jg)% zmaxcov) THEN
             echam_cov_config(jg)% jkscov = echam_cov_config(jg)% jkscov + 1
          ELSE
             EXIT
          END IF
       END DO
       !
       SELECT CASE (echam_cov_config(jg)% icov)
       CASE (1)
          !
          ! diagnose jksinv
          echam_cov_config(jg)% jksinv = 1
          DO jk = 1,klev
             IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > echam_cov_config(jg)% zmaxinv) THEN
                echam_cov_config(jg)% jksinv = echam_cov_config(jg)% jksinv + 1
             ELSE
                EXIT
             END IF
          END DO
          !
          ! diagnose jkeinv
          echam_cov_config(jg)% jkeinv = klev
          DO jk = klev,1,-1
             IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp < echam_cov_config(jg)% zmininv) THEN
                echam_cov_config(jg)% jkeinv = echam_cov_config(jg)% jkeinv - 1
             ELSE
                EXIT
             END IF
          END DO
          !
          ! check that crs and crt are smaller than csat 
          IF (echam_cov_config(jg)% crs >= echam_cov_config(jg)% csat) THEN
             CALL finish('eval_echam_cov_config', &
                  &      'echam_cov_config('//TRIM(cg)//')% crs >= echam_phy_config('//TRIM(cg)//')% csat is not allowed')
          END IF
          !
          IF (echam_cov_config(jg)% crt >= echam_cov_config(jg)% csat) THEN
             CALL finish('eval_echam_cov_config', &
                  &      'echam_cov_config('//TRIM(cg)//')% crt >= echam_phy_config('//TRIM(cg)//')% csat is not allowed')
          END IF
          !
       CASE (2)
          !
       CASE DEFAULT
          !
          CALL finish('eval_echam_cov_config', &
                  &   'echam_cov_config('//TRIM(cg)//')% icov /= 1 or 2 is not allowed')
       END SELECT
       !
    END DO
    !
  END SUBROUTINE eval_echam_cov_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_cov_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM cloud cover configuration')
    CALL message    ('','===============================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_cov_config('//TRIM(cg)//')% zmaxcov  ',echam_cov_config(jg)% zmaxcov )
       CALL print_value('    echam_cov_config('//TRIM(cg)//')% jkscov   ',echam_cov_config(jg)% jkscov  )
       CALL print_value('    echam_cov_config('//TRIM(cg)//')% icov     ',echam_cov_config(jg)% icov    )
       SELECT CASE (echam_cov_config(jg)% icov)
       CASE(1)
          CALL message    ('','---      --> use the fractional cloud cover scheme')
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% csat     ',echam_cov_config(jg)% csat    )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% crs      ',echam_cov_config(jg)% crs     )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% crt      ',echam_cov_config(jg)% crt     )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% nex      ',echam_cov_config(jg)% nex     )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% zmaxinv  ',echam_cov_config(jg)% zmaxinv )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% jksinv   ',echam_cov_config(jg)% jksinv  )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% zmininv  ',echam_cov_config(jg)% zmininv )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% jkeinv   ',echam_cov_config(jg)% jkeinv  )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% cinv     ',echam_cov_config(jg)% cinv    )
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% csatsc   ',echam_cov_config(jg)% csatsc  )
       CASE(2)
          CALL message    ('','---      --> use the 0/1 cloud cover scheme')
          CALL print_value('    echam_cov_config('//TRIM(cg)//')% csat     ',echam_cov_config(jg)% csat    )
       END SELECT
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_cov_config

  !----

END MODULE mo_echam_cov_config
