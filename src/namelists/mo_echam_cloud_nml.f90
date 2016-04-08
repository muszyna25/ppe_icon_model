!>
!! Namelist for configuring cloud parameterization in
!! the ECHAM physics package
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_cloud_nml

  USE mo_echam_cloud_config,  ONLY: echam_cloud_config
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_physical_constants,  ONLY: tmelt

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_echam_cloud_namelist

  !--------------------------------------------------------------
  ! Namelist variables 
  !
  REAL(wp) :: cthomi
  REAL(wp) :: cn0s
  REAL(wp) :: crhoi
  REAL(wp) :: crhosno
  REAL(wp) :: ccsaut
  REAL(wp) :: clmax
  REAL(wp) :: clmin
  REAL(wp) :: ceffmax    ! max eff.radius for ice cloud
  LOGICAL  :: lonacc

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
  !
  NAMELIST /echam_cloud_nml/                                &
       & cthomi  , cn0s    , crhoi   , crhosno , ccsaut   , &
       & clmax   , clmin   , ceffmax , lonacc  ,            &
       & ccsacl  , ccracl  , ccraut  , ceffmin , ccwmin   , &
       & cinv    , cauloc  , cqtmin  ,                      &
       & cn1lnd  , cn2lnd  , cn1sea  , cn2sea  ,            &
       & cinhomi , cinhoml1, cinhoml2, cinhoml3,            &
       & csecfrl , crs     , crt     , cvtfall , clwprat  , &
       & csatsc  , nex     , nadd    , cptop   , cpbot    , &
       & ncctop  , nccbot  , jbmin   , jbmax

CONTAINS
  !>
  !! Read the cloud namelist
  !!
  SUBROUTINE read_echam_cloud_namelist( filename )

    CHARACTER(LEN=*),INTENT(IN) :: filename
!!$    CHARACTER(LEN=*),PARAMETER :: &
!!$    routine = 'mo_echam_cloud_nml:read_echam_cloud_namelist'

    INTEGER  :: ist, funit, iunit

    !------------------------------------------------------------
    ! Set default values
    !
    cthomi  = tmelt-35.0_wp
    cn0s    = 3.e6_wp
    crhoi   = 500.0_wp
    crhosno = 100.0_wp
    ccsaut  = 95.0_wp
    clmax   = 0.5_wp
    clmin   = 0.0_wp
    ceffmax = 150.0_wp   ! max eff.radius for ice cloud
    lonacc  = .TRUE.

    ccsacl  = 0.10_wp
    ccracl  = 12.0_wp
    ccraut  = 20.0_wp
    ceffmin = 10.0_wp    ! min eff.radius for ice cloud
    ccwmin  = 1.e-7_wp   ! cloud water limit for cover>0
    cinv    = 0.25_wp    ! fraction of dry adiabatic lapse rate
    cauloc  = 10.0_wp
    cqtmin  = 1.e-12_wp  ! total water minimum

    cn1lnd  =  20._wp
    cn2lnd  = 120._wp
    cn1sea  =  20._wp
    cn2sea  =  40._wp

    cinhomi = 0.80_wp
    cinhoml1= 0.60_wp
    cinhoml2= 0.40_wp
    cinhoml3= 0.60_wp

    csecfrl = 1.0e-5_wp
    crs     = 0.999_wp   ! Critical relative humidity at surface
    crt     = 0.90_wp    ! Critical relative humidity aloft
    cvtfall = 3.0_wp
    clwprat = 4.0_wp
    csatsc  = 0.7_wp
    nex     = 2          ! Transition parameter for critical relative humidity profile
    nadd    = 0

    cptop  = 1000.0_wp   ! min. pressure level for cond.
    cpbot  = 50000.0_wp  ! max. pressure level for tropopause calc.

    ncctop = 13          ! max. level for condensation
    nccbot = 35          ! lowest level for tropopause calculation
    jbmin  = 40          ! highest inversion level
    jbmax  = 45          ! lowest inversion level

    !-------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('echam_cloud_nml')
      READ(funit,NML=echam_cloud_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !
    CALL open_nml(TRIM(filename))
    CALL position_nml('echam_cloud_nml',STATUS=ist)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, echam_cloud_nml)      ! write defaults to temporary text file
    END IF
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, echam_cloud_nml)       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, echam_cloud_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------
    ! Sanity check
    !

    !-----------------------------------------------------
    ! Store the namelist for restart
    !
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=echam_cloud_nml)
      CALL store_and_close_namelist(funit, 'echam_cloud_nml')
    ENDIF

    !-----------------------------------------------------
    ! Fill configuration state
    !
    echam_cloud_config% cthomi  = cthomi
    echam_cloud_config% cn0s    = cn0s
    echam_cloud_config% crhoi   = crhoi
    echam_cloud_config% crhosno = crhosno
    echam_cloud_config% ccsaut  = ccsaut
    echam_cloud_config% clmax   = clmax
    echam_cloud_config% clmin   = clmin
    echam_cloud_config% ceffmax = ceffmax
    echam_cloud_config% lonacc  = lonacc
    !
    echam_cloud_config% ccsacl  = ccsacl
    echam_cloud_config% ccracl  = ccracl
    echam_cloud_config% ccraut  = ccraut
    echam_cloud_config% ceffmin = ceffmin
    echam_cloud_config% ccwmin  = ccwmin
    echam_cloud_config% cinv    = cinv
    echam_cloud_config% cauloc  = cauloc
    echam_cloud_config% cqtmin  = cqtmin
    !
    echam_cloud_config% cn1lnd  = cn1lnd
    echam_cloud_config% cn2lnd  = cn2lnd
    echam_cloud_config% cn1sea  = cn1sea
    echam_cloud_config% cn2sea  = cn2sea
    !
    echam_cloud_config% cinhomi = cinhomi
    echam_cloud_config% cinhoml1= cinhoml1
    echam_cloud_config% cinhoml2= cinhoml2
    echam_cloud_config% cinhoml3= cinhoml3
    !
    echam_cloud_config% csecfrl = csecfrl
    echam_cloud_config% crs     = crs
    echam_cloud_config% crt     = crt
    echam_cloud_config% cvtfall = cvtfall
    echam_cloud_config% clwprat = clwprat
    echam_cloud_config% csatsc  = csatsc
    echam_cloud_config% nex     = nex
    echam_cloud_config% nadd    = nadd
    !
    echam_cloud_config% cptop   = cptop
    echam_cloud_config% cpbot   = cpbot
    !
    echam_cloud_config% ncctop  = ncctop
    echam_cloud_config% nccbot  = nccbot
    echam_cloud_config% jbmin   = jbmin
    echam_cloud_config% jbmax   = jbmax

  END SUBROUTINE read_echam_cloud_namelist
  !-------------

END MODULE mo_echam_cloud_nml
