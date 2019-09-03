!>
!! Configuration of the parameterization for NWP graupel microphysics,
!! that is used in the ECHAM physics package.
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!! First version by Monika Esch, MPI-M (2018-06)
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
MODULE mo_echam_mig_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_mig_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_mig_config   !< allocate and initialize echam_mig_config
!!$  PUBLIC ::    eval_echam_mig_config   !< evaluate echam_mig_config
  PUBLIC ::   print_echam_mig_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_mig'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM microphysics
  !!
  TYPE t_echam_mig_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! copies of the according nwp settings
     !
     ! thresholds
     REAL(wp) :: qi0_nwp        ! cloud ice threshold for autoconversion
     REAL(wp) :: qc0_nwp        ! cloud water threshold for autoconversion
     !
     ! grid scale microphysics
     REAL(wp) :: tune_zceff_min
     REAL(wp) :: tune_v0snow    ! previous ICON value was 20
     REAL(wp) :: tune_zvz0i     ! original value of Heymsfield+Donner 1990: 3.29
     REAL(wp) :: mu_rain        ! COSMO_EU default
     !
     !
  END TYPE t_echam_mig_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_mig_config), TARGET :: echam_mig_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_mig_config
    !
    ! Graupel microphyiscs configuration
    ! --------------------------------------
    !
    ! general thresholds
    echam_mig_config(:)% qi0_nwp         = 0.0_wp
    echam_mig_config(:)% qc0_nwp         = 0.0_wp
    !
    ! grid scale microphysics
    echam_mig_config(:)% tune_zceff_min  = 0.075_wp
    echam_mig_config(:)% tune_v0snow     = 25.0_wp      ! previous ICON value was 20
    echam_mig_config(:)% tune_zvz0i      = 1.25_wp      ! original value of Heymsfield+Donner 1990: 3.29
    echam_mig_config(:)% mu_rain         = 0.0_wp       ! COSMO_EU default
    !
  END SUBROUTINE init_echam_mig_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_echam_mig_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_echam_mig_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_mig_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','Graupel microphyiscs configuration')
    CALL message    ('','==================================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_mig_config('//TRIM(cg)//')% qi0_nwp        ',echam_mig_config(jg)% qi0_nwp  )
       CALL print_value('    echam_mig_config('//TRIM(cg)//')% qc0_nwp        ',echam_mig_config(jg)% qc0_nwp  )
       CALL message    ('','')
       CALL print_value('    echam_mig_config('//TRIM(cg)//')% tune_zceff_min ',echam_mig_config(jg)% tune_zceff_min  )
       CALL print_value('    echam_mig_config('//TRIM(cg)//')% tune_v0snow    ',echam_mig_config(jg)% tune_v0snow )
       CALL print_value('    echam_mig_config('//TRIM(cg)//')% tune_zvz0i     ',echam_mig_config(jg)% tune_zvz0i  )
       CALL print_value('    echam_mig_config('//TRIM(cg)//')% mu_rain        ',echam_mig_config(jg)% mu_rain  )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_mig_config

  !----

END MODULE mo_echam_mig_config
