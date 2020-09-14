!>
!! Configuration of the upper and lower levels for WMO_tropopause calculation
!! that is used in the ECHAM physics package.
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!! First version by Monika Esch, MPI-M (2019-11)
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
MODULE mo_echam_wmo_config

  USE mo_exception            ,ONLY: finish, message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_vertical_coord_table ,ONLY: vct_a

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_wmo_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_wmo_config   !< allocate and initialize echam_wmo_config
  PUBLIC ::    eval_echam_wmo_config   !< evaluate echam_wmo_config
  PUBLIC ::   print_echam_wmo_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_wmo'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM microphysics
  !!
  TYPE t_echam_wmo_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! vertical loop limits
     REAL(wp) :: zmaxwmo  !          maximum height (m) for tropopause calculation
     REAL(wp) :: zminwmo  !          minimum height (m) for tropopause calculation
     INTEGER  :: jkswmo   !          vertical start index for tropopause calculation
     INTEGER  :: jkewmo   !          vertical end   index for tropopause calculation
     !                               diagnosed in eval_echam_wmo_config
  END TYPE t_echam_wmo_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_wmo_config), TARGET :: echam_wmo_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_wmo_config
    !
    ! WMO_tropopause configuration
    ! ----------------------------
    !
    echam_wmo_config(:)% zmaxwmo  = 38000.
    echam_wmo_config(:)% zminwmo  =  5000.
    !
  END SUBROUTINE init_echam_wmo_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_wmo_config
    !
    INTEGER           :: jg, jk, klev
    !
    klev = SIZE(vct_a)-1
    !
    DO jg = 1,n_dom
       !
       ! diagnose jkswmo = start index
       !
       echam_wmo_config(jg)% jkswmo = 1
       DO jk = 1,klev
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > echam_wmo_config(jg)% zmaxwmo) THEN
             echam_wmo_config(jg)% jkswmo = jk + 1
          ELSE
             EXIT
          END IF
       END DO
       !
       ! diagnose jkewmo = end index
       !
       echam_wmo_config(jg)% jkewmo = echam_wmo_config(jg)% jkswmo
       DO jk = echam_wmo_config(jg)% jkswmo,klev
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > echam_wmo_config(jg)% zminwmo) THEN
             echam_wmo_config(jg)% jkewmo = jk + 1
          ELSE
             EXIT
          END IF
       END DO
       !
    END DO
    !
  END SUBROUTINE eval_echam_wmo_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_wmo_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','WMO tropopause configuration')
    CALL message    ('','============================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_wmo_config('//TRIM(cg)//')% zmaxwmo ',echam_wmo_config(jg)% zmaxwmo )
       CALL print_value('    echam_wmo_config('//TRIM(cg)//')% zminwmo ',echam_wmo_config(jg)% zminwmo )
       CALL print_value('    echam_wmo_config('//TRIM(cg)//')% jkswmo  ',echam_wmo_config(jg)% jkswmo  )
       CALL print_value('    echam_wmo_config('//TRIM(cg)//')% jkewmo  ',echam_wmo_config(jg)% jkewmo    )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_wmo_config

  !----

END MODULE mo_echam_wmo_config
