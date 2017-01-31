!>
!! @brief control variables for bgc modules
!!
!! Definition of variables
!!
MODULE mo_control_bgc

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

!  LOGICAL     :: lspinbgc = .FALSE.       !<  Switch for BGC spinup (alkalinity correction)
!  REAL(wp)    :: fspinbgc = 1._wp         !<  Factor for manual adjustment of alkalinity



  REAL(wp) :: deltacalc = 0._wp  !<  Correction for calcite pool [kMol/s]
  REAL(wp) :: deltaorg  = 0._wp  !<  Correction for organic carbon pool [kMol/s]
  REAL(wp) :: deltasil  = 0._wp  !<  Correction for silicate pool [kMol/s]


  ! Control variables

  REAL(wp) :: dtbgc              !<  time step length [sec].
  REAL(wp) :: inv_dtbgc          !<  inverse time step length [sec^-1].
  REAL(wp) :: dtb                !<  time step length [days].
  REAL(wp) :: inv_dtb            !<  inverse time step length [days].
  INTEGER  :: ndtdaybgc          !<  time steps per day.

  INTEGER  :: ldtbgc             !<  time step number from bgc restart file
  INTEGER  :: ldtrunbgc          !<  actual time steps of run.

  INTEGER  :: icyclibgc          !<  switch for cyclicity.
  INTEGER  :: ndtrunbgc          !<  total no. of time steps of run.

  INTEGER  :: bgcstartyear       !<  year of ocean restart file
  INTEGER  :: bgcstartmonth      !<  month of ocean restart file
  INTEGER  :: bgcstartday        !<  day of ocean restart file

  INTEGER  :: bgc_zlevs           !<  time step number from bgc restart file
  INTEGER  :: bgc_nproma          !<  actual time steps of run.

  REAL(wp) :: rmasks = 0.0_wp     !<  value at wet cells in sediment.
  REAL(wp) :: rmasko = 0.0_wp     !<  value at wet cells in ocean.    


  INTEGER:: bgc_gin, bgc_arctic, bgc_lab, bgc_natl, bgc_atl, bgc_tatl, bgc_tropac,  bgc_land, bgc_ind, bgc_soce, bgc_npac, bgc_carb
  
END MODULE mo_control_bgc
