!>
!! Configuration of the ECHAM cumulus convection scheme. 
!! Includes switches and tuning parameters, as well as 
!! other control variables.
!!
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
MODULE mo_echam_conv_config

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, print_value, message, message_text
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_physical_constants, ONLY: grav

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_echam_conv_config, echam_conv_config
  PUBLIC :: configure_echam_convection, cleanup_echam_convection

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  !>
  !! Derived type containing main swithes for configuring 
  !! the cumulus convection scheme of ECHAM
  !!
  TYPE t_echam_conv_config

    ! Namelist variables

    INTEGER :: iconv     !< 1,2,3 for different convection schemes
    INTEGER :: ncvmicro  !< 0 or 1. Scheme for convective microphysics
  
    LOGICAL :: lmfpen    !< true when penetrative convection is switched on
!    LOGICAL :: lmfmid    !< true when midlevel    convection is switched on
!    LOGICAL :: lmfdd     !< true when cumulus downdraft      is switched on
!    LOGICAL :: lmfdudv   !< true when cumulus friction       is switched on

!    REAL(wp) :: dlev     !< "zdlev" in subroutine "cuasc". 
                         !< Critical thickness (unit: Pa) necessary for 
                         !< the onset of convective precipitation
  
!    REAL(wp) :: cmftau   !< characteristic adjustment time scale
                         !< (replaces "ztau" in "cumastr"
!    REAL(wp) :: cmfctop  !< fractional convective mass flux across the 
                         !< top of cloud
!    REAL(wp) :: cprcon   !< coefficient for determining conversion
                         !< from cloud water to rain
  
!    REAL(wp) :: cminbuoy !< minimum excess buoyancy
!    REAL(wp) :: entrpen  !< entrainment rate for penetrative convection
 
    ! Currently unused namelist variables 
    !INTEGER :: nauto        !< autoconversion scheme. 1 or 2.
    !LOGICAL :: lconvmassfix !< aerosol mass fixer in convection
    !LOGICAL :: lmfscv       !< true when shallow convection is switched on

    ! Derived variables

    INTEGER :: nmctop    !< max. level for cloud base of mid level conv.

    REAL(wp),ALLOCATABLE :: cevapcu(:)  !< evaporation coefficient for kuo0
                                        !< In ECHAM6 it is declared in 
                                        !< mo_physc2, allocated in subroutine 
                                        !< alloc_mods, and initialized in 
                                        !< subroutine iniphy.
  
  END TYPE t_echam_conv_config

  !>
  !! The configuration state (variable).
  !! So far we have not yet tried to use different configurations for different
  !! domains (grid levels) in experiments with local grid refinement (nesting),
  !! thus the variable is declared as a scalar. Later it might be changed into
  !! an array of shape (/n_dom/) or (/MAX_DOM/).
  !!
  TYPE(t_echam_conv_config) :: echam_conv_config

CONTAINS
  !---------------------------------------------------------------------------
  !>
  !!
  !! Assign value to derived variables in echam_conv_config.
  !!
  !! @Revision history
  !! Adapted from ECHAM6 by Hui Wan (MPI-M, 2010-2011)
  !!
  SUBROUTINE configure_echam_convection( nlev, vct_a, vct_b, ceta )

    INTEGER, INTENT(IN) :: nlev
    REAL(WP),INTENT(IN) :: vct_a(nlev+1)
    REAL(WP),INTENT(IN) :: vct_b(nlev+1)
    REAL(WP),INTENT(IN) :: ceta (nlev)

    REAL(wp) :: zp(nlev), zph(nlev+1), ztmp
    INTEGER  :: jk, istat

    CHARACTER(LEN=*),PARAMETER :: &
             routine = 'mo_echam_conv_config:config_echam_convection'

    !------------------------------------------------------------------------
    ! Print the configuration on stdio
    !------------------------------------------------------------------------
    CALL message('','')
    CALL message('','------- configuration of the ECHAM convection scheme --------')

    SELECT CASE (echam_conv_config% iconv)
    CASE(1); CALL message('','--- iconv = 1 -> Convection: Nordeng (default)')
    CASE(2); CALL message('','--- iconv = 2 -> Convection: Tiedtke')
    CASE(3); CALL message('','--- iconv = 3 -> Convection: Hybrid')
    CASE default
      WRITE(message_text,'(a,i0,a)') 'iconv = ',echam_conv_config% iconv, &
                                     ' is not supported'
      CALL finish(TRIM(routine),message_text)
    END SELECT

    SELECT CASE(echam_conv_config% ncvmicro)
    CASE (0); CALL message('','--- ncvmicro = 0')
    CASE DEFAULT
      CALL finish(TRIM(routine),'ncvmicro > 0 not yet supported in ICON')
    END SELECT

    CALL print_value(' lmfpen   ', echam_conv_config% lmfpen)

!    CALL print_value(' cmftau   ', echam_conv_config% cmftau)
!    CALL print_value(' cmfctop  ', echam_conv_config% cmfctop)
!    CALL print_value(' cprcon   ', echam_conv_config% cprcon)
!    CALL print_value(' cminbuoy ', echam_conv_config% cminbuoy)
!    CALL print_value(' entrpen  ', echam_conv_config% entrpen)
!    CALL print_value(' dlev     ', echam_conv_config% dlev)

    CALL message('','---------------------------')
    CALL message('','')

    !------------------------------------------------------------------------
    ! Determine highest level *nmctop* for cloud base of midlevel convection
    ! assuming nmctop=9 (300 hPa) for the standard 19 level model
    !------------------------------------------------------------------------
    ! Compute half level pressure values, assuming 101320 Pa surface pressure

    DO jk=1,nlev+1
      zph(jk) = vct_a(jk) + vct_b(jk)*101320.0_wp
    END DO

    ! Compute full level pressure

    DO jk = 1, nlev
      zp(jk) = (zph(jk)+zph(jk+1))*0.5_wp
    END DO

    ! Search for 300 hPa level

    DO jk = 1, nlev
      echam_conv_config%nmctop = jk
      IF(zp(jk).GE.30000.0_wp) EXIT
    END DO

    CALL message('','')
    CALL print_value('lowest model level for cloud base '//&
                     'of mid level convection: nmctop = ', &
                     echam_conv_config%nmctop)
    CALL message('','')

    !--------------------------------------
    ! Set evaporation coefficient for kuo0
    !--------------------------------------

    ALLOCATE( echam_conv_config%cevapcu(nlev),STAT=istat )
    IF (istat/=SUCCESS) CALL finish(TRIM(routine),'allocation of cevapcu failed')

    DO jk = 1,nlev
       ztmp = 1.E3_wp/(38.3_wp*0.293_wp)*SQRT(ceta(jk))
       echam_conv_config%cevapcu(jk) = 1.93E-6_wp*261._wp*SQRT(ztmp)*0.5_wp/grav
    END DO

  END SUBROUTINE configure_echam_convection
  !------------
  !>
  !!
  !! Deallocate memory
  !!
  SUBROUTINE cleanup_echam_convection

    INTEGER :: istat
    CHARACTER(LEN=*),PARAMETER :: &
             routine = 'mo_echam_conv_config:cleanup_echam_convection'

    DEALLOCATE( echam_conv_config%cevapcu,STAT=istat )
    IF (istat/=SUCCESS) CALL finish(TRIM(routine),'deallocation of cevapcu failed')

  END SUBROUTINE cleanup_echam_convection
  !-------------

END MODULE mo_echam_conv_config
