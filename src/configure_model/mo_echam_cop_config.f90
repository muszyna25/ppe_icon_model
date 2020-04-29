!>
!! Configuration of the parameterization for cloud optical properties,
!! that is used in the ECHAM physics package.
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!! First version by Monika Esch, MPI-M (2019-10)
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
MODULE mo_echam_cop_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_echam_phy_config     ,ONLY: echam_phy_tc, dt_zero
  USE mtime,                   ONLY: OPERATOR(>)

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_cop_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_cop_config   !< allocate and initialize echam_cop_config
!!$  PUBLIC ::    eval_echam_cop_config   !< evaluate echam_cop_config
  PUBLIC ::   print_echam_cop_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_cop'
  
  !>
  !! Configuration type containing parameters for the configuration of the cloud optical properties
  !!
  TYPE t_echam_cop_config
     !
     ! configuration parameters
     ! ------------------------
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
  END TYPE t_echam_cop_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_cop_config), TARGET :: echam_cop_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_cop_config

    !
    ! echam cloud microphysics active
    !
    ! cloud optical properties configuration
    ! --------------------------------------
    !
    ! cloud droplet number concentration
    echam_cop_config(:)% cn1lnd   =  20._wp
    echam_cop_config(:)% cn2lnd   = 180._wp
    echam_cop_config(:)% cn1sea   =  20._wp
    echam_cop_config(:)% cn2sea   =  80._wp
    !
    ! cloud optics
    echam_cop_config(:)% cinhomi  = 0.80_wp
    echam_cop_config(:)% cinhoml1 = 0.80_wp
    echam_cop_config(:)% cinhoml2 = 0.40_wp
    echam_cop_config(:)% cinhoml3 = 0.80_wp
    !
    !
  END SUBROUTINE init_echam_cop_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_echam_cop_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_echam_cop_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_cop_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','cloud optical properties configuration')
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
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cn1lnd   ',echam_cop_config(jg)% cn1lnd  )
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cn2lnd   ',echam_cop_config(jg)% cn2lnd  )
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cn1sea   ',echam_cop_config(jg)% cn1sea  )
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cn2sea   ',echam_cop_config(jg)% cn2sea  )
       CALL message    ('','')
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cinhomi  ',echam_cop_config(jg)% cinhomi )
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cinhoml1 ',echam_cop_config(jg)% cinhoml1)
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cinhoml2 ',echam_cop_config(jg)% cinhoml2)
       CALL print_value('    echam_cop_config('//TRIM(cg)//')% cinhoml3 ',echam_cop_config(jg)% cinhoml3)
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_echam_cop_config

  !----

END MODULE mo_echam_cop_config
