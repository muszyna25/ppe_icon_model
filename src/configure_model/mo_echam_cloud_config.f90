!>
!! Configuration of the ECHAM cumulus cloud scheme. 
!! Includes switches and tuning parameters, as well as 
!! other control variables.
!!≤
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI (2011-07)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_cloud_config

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: print_value, message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_echam_cloud_config, echam_cloud_config
  PUBLIC :: configure_echam_cloud

  !>
  !! Derived type containing main swithes for configuring 
  !! the cumulus cloud scheme of ECHAM
  !!
  TYPE t_echam_cloud_config

     ! Namelist variables

     REAL(wp) :: cthomi
     REAL(wp) :: cn0s
     REAL(wp) :: crhoi
     REAL(wp) :: crhosno
     REAL(wp) :: ccsaut
     REAL(wp) :: clmax
     REAL(wp) :: clmin
     REAL(wp) :: ceffmax    ! max eff.radius for ice cloud

     REAL(wp) :: ccsacl
     REAL(wp) :: ccracl
     REAL(wp) :: ccraut
     REAL(wp) :: ceffmin    ! min eff.radius for ice cloud
     REAL(wp) :: ccwmin     ! cloud water limit for cover>0
     REAL(wp) :: cinv       ! fraction of dry adiabatic lapse rate
     REAL(wp) :: cauloc
     REAL(wp) :: cqtmin     ! total water minimum

     REAL(wp) :: cn1lnd
     REAL(wp) :: cn2lnd
     REAL(wp) :: cn1sea
     REAL(wp) :: cn2sea

     REAL(wp) :: cinhomi
     REAL(wp) :: cinhoml1
     REAL(wp) :: cinhoml2
     REAL(wp) :: cinhoml3

     REAL(wp) :: csecfrl
     REAL(wp) :: crs     ! Critical relative humidity at surface
     REAL(wp) :: crt     ! Critical relative humidity aloft
     REAL(wp) :: cvtfall
     REAL(wp) :: clwprat
     REAL(wp) :: csatsc
     INTEGER  :: nex     ! Transition parameter for critical relative humidity profile
     INTEGER  :: nadd

     REAL(wp) :: cptop      ! min. pressure level for cond.
     REAL(wp) :: cpbot      ! max. pressure level for tropopause calc.

     INTEGER  :: ncctop     ! max. level for condensation
     INTEGER  :: nccbot     ! lowest level for tropopause calculation
     INTEGER  :: jbmin      ! highest inversion level
     INTEGER  :: jbmax      ! lowest inversion level

  END TYPE t_echam_cloud_config

  !>
  !! The configuration state (variable).
  !! So far we have not yet tried to use different configurations for different
  !! domains (grid levels) in experiments with local grid refinement (nesting),
  !! thus the variable is declared as a scalar. Later it might be changed into
  !! an array of shape (/n_dom/) or (/MAX_DOM/).
  !!
  TYPE(t_echam_cloud_config), TARGET :: echam_cloud_config

CONTAINS
  !---------------------------------------------------------------------------
  !>
  !!
  !! Assign value to derived variables in echam_cloud_config.
  !!
  !! @Revision history
  !! Adapted from ECHAM6 by Hui Wan (MPI-M, 2010-2011)
  !!
  SUBROUTINE configure_echam_cloud( )

!!$    CHARACTER(LEN=*),PARAMETER :: &
!!$             routine = 'mo_echam_cloud_config:configure_echam_cloud'

    !------------------------------------------------------------------------
    ! Print the configuration on stdio
    !
    CALL message('','')
    CALL message('','------- configuration of the ECHAM cloud scheme --------')
    CALL message('','')
    CALL print_value(' cthomi  ', echam_cloud_config% cthomi  )
    CALL print_value(' cn0s    ', echam_cloud_config% cn0s  )
    CALL print_value(' crhoi   ', echam_cloud_config% crhoi   )
    CALL print_value(' crhosno ', echam_cloud_config% crhosno )
    CALL print_value(' ccsaut  ', echam_cloud_config% ccsaut  )
    CALL print_value(' clmax   ', echam_cloud_config% clmax   )
    CALL print_value(' clmin   ', echam_cloud_config% clmin   )
    CALL print_value(' ceffmax ', echam_cloud_config% ceffmax )
    CALL message('','')
    CALL print_value(' ccsacl  ', echam_cloud_config% ccsacl  )
    CALL print_value(' ccracl  ', echam_cloud_config% ccracl  )
    CALL print_value(' ccraut  ', echam_cloud_config% ccraut  )
    CALL print_value(' ceffmin ', echam_cloud_config% ceffmin )
    CALL print_value(' ccwmin  ', echam_cloud_config% ccwmin  )
    CALL print_value(' cinv    ', echam_cloud_config% cinv    )
    CALL print_value(' cauloc  ', echam_cloud_config% cauloc  )
    CALL print_value(' cqtmin  ', echam_cloud_config% cqtmin  )
    CALL message('','')
    CALL print_value(' cn1lnd  ', echam_cloud_config% cn1lnd  )
    CALL print_value(' cn2lnd  ', echam_cloud_config% cn2lnd  )
    CALL print_value(' cn1sea  ', echam_cloud_config% cn1sea  )
    CALL print_value(' cn2sea  ', echam_cloud_config% cn2sea  )
    CALL message('','')
    CALL print_value(' cinhomi ', echam_cloud_config% cinhomi )
    CALL print_value(' cinhoml1', echam_cloud_config% cinhoml1)
    CALL print_value(' cinhoml2', echam_cloud_config% cinhoml2)
    CALL print_value(' cinhoml3', echam_cloud_config% cinhoml3)
    CALL message('','')
    CALL print_value(' csecfrl ', echam_cloud_config% csecfrl )
    CALL print_value(' crs     ', echam_cloud_config% crs     )
    CALL print_value(' crt     ', echam_cloud_config% crt     )
    CALL print_value(' cvtfall ', echam_cloud_config% cvtfall )
    CALL print_value(' clwprat ', echam_cloud_config% clwprat )
    CALL print_value(' csatsc  ', echam_cloud_config% csatsc  )
    CALL print_value(' nex     ', echam_cloud_config% nex     )
    CALL print_value(' nadd    ', echam_cloud_config% nadd    )
    CALL message('','')
    CALL print_value(' cptop   ', echam_cloud_config% cptop   )
    CALL print_value(' cpbot   ', echam_cloud_config% cpbot   )
    CALL message('','')
    CALL print_value(' ncctop  ', echam_cloud_config% ncctop  )
    CALL print_value(' nccbot  ', echam_cloud_config% nccbot  )
    CALL print_value(' jbmin   ', echam_cloud_config% jbmin   )
    CALL print_value(' jbmax   ', echam_cloud_config% jbmax   )
    CALL message('','')

  END SUBROUTINE configure_echam_cloud
  !-------------

END MODULE mo_echam_cloud_config
