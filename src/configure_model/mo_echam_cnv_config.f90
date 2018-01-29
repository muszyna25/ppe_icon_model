!>
!! Configuration of the parameterization for moist convection,
!! that is used in the ECHAM physics package.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!! First version by Marco Giorgetta, MPI-M (2017-11)
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
MODULE mo_echam_cnv_config

  USE mo_exception            ,ONLY: message, print_value, finish
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom, success
  USE mo_grid_config          ,ONLY: n_dom
  USE mo_run_config           ,ONLY: nlev
  USE mo_vertical_coord_table ,ONLY: vct_a
  USE mo_physical_constants   ,ONLY: grav, p0sl_bg, p0ref

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         echam_cnv_config   !< user specified configuration parameters
  PUBLIC ::    init_echam_cnv_config   !< allocate and initialize echam_cnv_config
  PUBLIC ::    eval_echam_cnv_config   !< evaluate echam_cnv_config
  PUBLIC ::   print_echam_cnv_config   !< print out
  PUBLIC ::   alloc_echam_cnv_config   !< allocate
  PUBLIC :: dealloc_echam_cnv_config   !< deallocate

  ! variables
  PUBLIC ::                  cevapcu   !< evaporation coefficient (n_dom,nlev)

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'echam_cnv'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the ECHAM convection
  !!
  TYPE t_echam_cnv_config
     !
     ! configuration parameters
     ! ------------------------
     !
     LOGICAL  :: lmfpen    !< true when penetrative convection is switched on
     LOGICAL  :: lmfmid    !< true when midlevel    convection is switched on
     LOGICAL  :: lmfdd     !< true when cumulus downdraft      is switched on
     LOGICAL  :: lmfdudv   !< true when cumulus friction       is switched on
     !
     REAL(wp) :: entrmid   !< average entrainment rate for midlevel convection
     REAL(wp) :: entrscv   !< average entrainment rate for shallow convection
     REAL(wp) :: entrpen   !< average entrainment rate for penetrative convection
     REAL(wp) :: entrdd    !< average entrainment rate for cumulus downdrafts
     !
     REAL(wp) :: cprcon    !< coefficient for determining conversion from cloud water to rain
     REAL(wp) :: cmfctop   !< fractional convective mass flux across the top of cloud
     REAL(wp) :: cmfdeps   !< fractional convective mass flux for downdrafts at lfs
     !
     REAL(wp) :: cminbuoy  !< minimum excess buoyancy
     REAL(wp) :: cmaxbuoy  !< maximum excess buoyancy
     REAL(wp) :: cbfac     !< factor for std dev of virtual pot temp
     REAL(wp) :: centrmax  !< maximum entrainment/detrainment rate
     !
     REAL(wp) :: dlev_land !< minimum pressure thickness of clouds for precipitation over land
     REAL(wp) :: dlev_ocean!< minimum pressure thickness of clouds for precipitation over ocean
     !
     REAL(wp) :: cmftau    !< characteristic adjustment time scale (s)
     !
     REAL(wp) :: cmfcmin   !< minimum massflux value (for safety)
     REAL(wp) :: cmfcmax   !< maximum massflux value allowed for updrafts etc
     !
     INTEGER  :: nmctop    !< max. level for cloud base of mid level conv.
     !
  END TYPE t_echam_cnv_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_echam_cnv_config), TARGET :: echam_cnv_config(max_dom)

  !>
  !! Evaporation coefficient for kuo0, for multiple domains/grids.
  !!
  REAL(wp),    ALLOCATABLE, TARGET :: cevapcu(:,:)  !< with dimensions (nlev,n_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_echam_cnv_config
    !
    ! ECHAM convection configuration
    ! ------------------------------
    !
    echam_cnv_config(:)% lmfmid     = .TRUE.
    echam_cnv_config(:)% lmfpen     = .TRUE.
    echam_cnv_config(:)% lmfdd      = .TRUE.
    echam_cnv_config(:)% lmfdudv    = .TRUE.
    !
    echam_cnv_config(:)% entrscv    = 3.0e-3_wp
    echam_cnv_config(:)% entrmid    = 2.0e-4_wp
    echam_cnv_config(:)% entrpen    = 2.0e-4_wp
    echam_cnv_config(:)% entrdd     = 4.0e-4_wp
    !
    echam_cnv_config(:)% cprcon     = 2.5e-4_wp
    echam_cnv_config(:)% cmfctop    = 0.2_wp
    echam_cnv_config(:)% cmfdeps    = 0.3_wp
    !
    echam_cnv_config(:)% cminbuoy   = 0.2_wp
    echam_cnv_config(:)% cmaxbuoy   = 1.0_wp
    echam_cnv_config(:)% cbfac      = 1.0_wp
    echam_cnv_config(:)% centrmax   = 3.0e-4_wp
    !
    echam_cnv_config(:)% dlev_land  = 0.0_wp
    echam_cnv_config(:)% dlev_ocean = 0.0_wp
    !
    echam_cnv_config(:)% cmftau     = 3600.0_wp
    !
    echam_cnv_config(:)% cmfcmin    = 1.0e-10_wp
    echam_cnv_config(:)% cmfcmax    = 1.0_wp
    !
    echam_cnv_config(:)% nmctop     = nlev/3
    !
  END SUBROUTINE init_echam_cnv_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_echam_cnv_config
    !
    CHARACTER(LEN=*),PARAMETER :: routine = 'eval_echam_cnv_config'
    !
    REAL(wp) :: zp(nlev), zph(nlev+1), zeta(nlev), ztmp
    INTEGER  :: jk,nk
    !
    !------------------------------------------------------------------------
    ! Determine highest level *nmctop* for cloud base of midlevel convection
    ! search highest lower half level below 300 hPa9000 m height
    !
    DO jk = 1, nlev
      nk = jk
      IF(vct_a(jk+1) < 9000.0_wp) EXIT
    END DO
    echam_cnv_config(:)%nmctop = nk
    !
    !------------------------------------------------------------------------
    ! Set evaporation coefficient for kuo0
    !
    ! half level pressure assuming 7500 m scale height
    DO jk=1,nlev+1
       zph(jk) = p0sl_bg*EXP(-vct_a(jk)/7500._wp) 
    END DO
    !
    ! full level pressure and eta
    DO jk = 1, nlev
      zp(jk)   = (zph(jk)+zph(jk+1))*0.5_wp
      zeta(jk) = zp(jk)/p0ref
    END DO
    !
    DO jk = 1,nlev
       ztmp          = 1.E3_wp/(38.3_wp*0.293_wp)*SQRT(zeta(jk))
       cevapcu(jk,:) = 1.93E-6_wp*261._wp*SQRT(ztmp)*0.5_wp/grav
    END DO
    !
  END SUBROUTINE eval_echam_cnv_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_echam_cnv_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','ECHAM convection configuration')
    CALL message    ('','==============================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% lmfmid    ', echam_cnv_config(jg)% lmfmid    )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% lmfpen    ', echam_cnv_config(jg)% lmfpen    )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% lmfdd     ', echam_cnv_config(jg)% lmfdd     )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% lmfdudv   ', echam_cnv_config(jg)% lmfdudv   )
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% entrscv   ', echam_cnv_config(jg)% entrscv   )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% entrmid   ', echam_cnv_config(jg)% entrmid   )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% entrpen   ', echam_cnv_config(jg)% entrpen   )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% entrdd    ', echam_cnv_config(jg)% entrdd    )
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cprcon    ', echam_cnv_config(jg)% cprcon    )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cmfctop   ', echam_cnv_config(jg)% cmfctop   )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cmfdeps   ', echam_cnv_config(jg)% cmfdeps   )
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cminbuoy  ', echam_cnv_config(jg)% cminbuoy  )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cmaxbuoy  ', echam_cnv_config(jg)% cmaxbuoy  )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cbfac     ', echam_cnv_config(jg)% cbfac     )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% centrmax  ', echam_cnv_config(jg)% centrmax  )
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% dlev_land ', echam_cnv_config(jg)% dlev_land )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% dlev_ocean', echam_cnv_config(jg)% dlev_ocean)
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cmftau    ', echam_cnv_config(jg)% cmftau    )
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cmfcmin   ', echam_cnv_config(jg)% cmfcmin   )
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% cmfcmax   ', echam_cnv_config(jg)% cmfcmax   )
       CALL message('','')
       CALL print_value('    echam_cnv_config('//TRIM(cg)//')% nmctop    ', echam_cnv_config(jg)% nmctop    )
       !
    END DO
    !
  END SUBROUTINE print_echam_cnv_config

  !----

  !>
  !! Allocate memory
  !!
  SUBROUTINE alloc_echam_cnv_config
    !
    CHARACTER(LEN=*),PARAMETER :: routine = 'alloc_echam_cnv_config'
    !
    INTEGER :: istat
    !
    ALLOCATE( cevapcu(nlev,n_dom),STAT=istat )
    IF (istat/=SUCCESS) CALL finish(TRIM(routine),'allocation of cevapcu failed')
    !
  END SUBROUTINE alloc_echam_cnv_config

  !----

  !>
  !! Deallocate memory
  !!
  SUBROUTINE dealloc_echam_cnv_config
    !
    CHARACTER(LEN=*),PARAMETER :: routine = 'dealloc_echam_cnv_config'
    !
    INTEGER :: istat
    !
    DEALLOCATE( cevapcu,STAT=istat )
    IF (istat/=SUCCESS) CALL finish(TRIM(routine),'deallocation of cevapcu failed')
    !
  END SUBROUTINE dealloc_echam_cnv_config

  !----

END MODULE mo_echam_cnv_config
