!>
!! Contains the tunable parameters for the cloud microphysics.
!! Subroutines for intializing these values are also included.
!!
!! @par Revision History
!!  Originally the module "mo_cloud" from ECHAM, by A. Tompkins (2000-07)
!!  Transfered to and modified for ICON by Hui Wan, MPI (2010-07)
!!  Tompkins removed again and adapted to ECHAM6.3 by M. Esch (2015-05)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_cloud_params

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: tmelt
  USE mo_exception,          ONLY: print_value

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_cloud_params'

  PUBLIC :: sucloud

  PUBLIC :: ceffmax
  PUBLIC :: cthomi,cn0s,crhoi,crhosno
  PUBLIC :: ccsaut
  PUBLIC :: clmin,clmax
  PUBLIC :: clwprat
  PUBLIC :: crs,crt,nex,nadd
  PUBLIC :: cqtmin
  PUBLIC :: cptop, cpbot, ncctop, nccbot
  PUBLIC :: jbmin, jbmax
  PUBLIC :: lonacc

  PUBLIC :: csatsc, cinv, ceffmin, csecfrl, cvtfall, ccwmin
  PUBLIC :: ccraut, ccsacl, ccracl, cauloc

  !----------------------------------------
  ! default values for cloud microphysics
  !----------------------------------------

  REAL(wp),    PARAMETER :: cthomi  = tmelt-35.0_wp
  REAL(wp),    PARAMETER :: cn0s    = 3.e6_wp
  REAL(wp),    PARAMETER :: crhoi   = 500.0_wp
  REAL(wp),    PARAMETER :: crhosno = 100.0_wp
  REAL(wp),    PARAMETER :: ccsaut  = 95.0_wp
  REAL(wp),    PARAMETER :: clmax   = 0.5_wp
  REAL(wp),    PARAMETER :: clmin   = 0.0_wp
  REAL(wp),    PARAMETER :: ceffmax = 150.0_wp   ! max eff.radius for ice cloud
  LOGICAL,     PARAMETER :: lonacc  = .TRUE.

  REAL(wp),    PARAMETER :: ccsacl  = 0.10_wp
  REAL(wp),    PARAMETER :: ccracl  = 12.0_wp
  REAL(wp),    PARAMETER :: ccraut  = 20.0_wp
  REAL(wp),    PARAMETER :: ceffmin = 10.0_wp    ! min eff.radius for ice cloud
  REAL(wp),    PARAMETER :: ccwmin  = 1.e-7_wp   ! cloud water limit for cover>0
  REAL(wp),    PARAMETER :: cinv    = 0.25_wp    ! fraction of dry adiabatic lapse rate
  REAL(wp),    PARAMETER :: cauloc  = 10.0_wp
  REAL(wp),    PARAMETER :: cqtmin = 1.e-12_wp   ! total water minimum
  !-----------------------------------------------------------------------------------
  ! Define parameters depending on resolution set in subroutine sucloud of this module
  !----------------------------------------------------------------------------------- 
  REAL(wp) ::    csecfrl
  REAL(wp) ::    crs            ! Critical relative humidity at surface
  REAL(wp) ::    crt            ! Critical relative humidity aloft
  REAL(wp) ::    cvtfall
  REAL(wp) ::    clwprat
  REAL(wp) ::    csatsc
  INTEGER  ::    nex            ! Transition parameter for critical relative humidity profile
  INTEGER  ::    nadd

  !---------------------------------------
  ! default values for cloud cover scheme
  !---------------------------------------

  REAL(wp),    PARAMETER :: cptop  = 1000.0_wp   ! min. pressure level for cond.
  REAL(wp),    PARAMETER :: cpbot  = 50000.0_wp  ! max. pressure level for tropopause calc.

  !-------------------------------------------------------------
  ! parameters initialized in subroutine sucloud of this module
  !-------------------------------------------------------------

  INTEGER            :: ncctop           ! max. level for condensation
  INTEGER            :: nccbot           ! lowest level for tropopause calculation
  INTEGER            :: jbmin            ! highest inversion level
  INTEGER            :: jbmax            ! lowest inversion level

CONTAINS
  !>
  !!
  !! @author E. Roeckner, MPI, October 2001
  !!
  !! Defines highest level *ncctop* where condensation is allowed.
  !!
  !! This routine is called from *iniphy* in ECHAM
  !! This routine is called from *mo_echam_phy_init* in ICON
  !!
  SUBROUTINE sucloud ( )

    ! local variables
    !
    ! Calculate values for ECHAM using vct
    ! Preset values for ICON (first attempt): Use values
    !  of ECHAM L47
    !
    jbmin=40
    jbmax=45
    ncctop=13
    nccbot=35
    !
    CALL print_value ('highest inversion level           : jbmin  = ', jbmin )
    CALL print_value ('lowest  inversion level           : jbmax  = ', jbmax )
    CALL print_value ('highest level for condensation    : ncctop = ', ncctop)
    CALL print_value ('lowest  level for tropopause calc.: nccbot = ', nccbot)
    !
    ! -- set resolution dependent parameters
    !
    crs     = 0.975_wp
    crt     = 0.75_wp
    cvtfall = 2.5_wp
    csecfrl = 5.e-6_wp
    clwprat = 4.0_wp
    csatsc  = 0.7_wp
    nex     = 2
    nadd    = 0

END SUBROUTINE sucloud

END MODULE mo_echam_cloud_params
