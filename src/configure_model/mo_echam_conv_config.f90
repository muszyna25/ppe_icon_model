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
  USE mo_exception,          ONLY: finish, print_value, message
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_physical_constants, ONLY: grav, p0sl_bg, p0ref

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_echam_conv_config, echam_conv_config
  PUBLIC :: configure_echam_convection, cleanup_echam_convection

  !>
  !! Derived type containing main swithes for configuring 
  !! the cumulus convection scheme of ECHAM
  !!
  TYPE t_echam_conv_config

    ! Namelist variables

    LOGICAL  :: lmfpen    !< true when penetrative convection is switched on
    LOGICAL  :: lmfmid    !< true when midlevel    convection is switched on
    LOGICAL  :: lmfdd     !< true when cumulus downdraft      is switched on
    LOGICAL  :: lmfdudv   !< true when cumulus friction       is switched on

    REAL(wp) :: entrmid   !< average entrainment rate for midlevel convection
    REAL(wp) :: entrscv   !< average entrainment rate for shallow convection
    REAL(wp) :: entrpen   !< average entrainment rate for penetrative convection
    REAL(wp) :: entrdd    !< average entrainment rate for cumulus downdrafts

    REAL(wp) :: cprcon    !< coefficient for determining conversion from cloud water to rain
    REAL(wp) :: cmfctop   !< fractional convective mass flux across the top of cloud
    REAL(wp) :: cmfdeps   !< fractional convective mass flux for downdrafts at lfs

    REAL(wp) :: cminbuoy  !< minimum excess buoyancy
    REAL(wp) :: cmaxbuoy  !< maximum excess buoyancy
    REAL(wp) :: cbfac     !< factor for std dev of virtual pot temp
    REAL(wp) :: centrmax  !< maximum entrainment/detrainment rate

    REAL(wp) :: cmftau    !< characteristic adjustment time scale (s)


    ! Variables set by routine configure_echam_convection

    REAL(wp) :: cmfcmin   !< minimum massflux value (for safety)
    REAL(wp) :: cmfcmax   !< maximum massflux value allowed for updrafts etc

    INTEGER  :: nmctop    !< max. level for cloud base of mid level conv.

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
  TYPE(t_echam_conv_config), TARGET :: echam_conv_config

CONTAINS
  !---------------------------------------------------------------------------
  !>
  !!
  !! Assign value to derived variables in echam_conv_config.
  !!
  !! @Revision history
  !! Adapted from ECHAM6 by Hui Wan (MPI-M, 2010-2011)
  !!
  SUBROUTINE configure_echam_convection( nlev, vct_a, vct_b )

    INTEGER, INTENT(IN) :: nlev
    REAL(WP),INTENT(IN) :: vct_a(nlev+1)
    REAL(WP),INTENT(IN) :: vct_b(nlev+1)

    REAL(wp) :: zp(nlev), zph(nlev+1), ztmp, zeta(nlev)
    INTEGER  :: jk, istat

    CHARACTER(LEN=*),PARAMETER :: &
             routine = 'mo_echam_conv_config:config_echam_convection'

    !------------------------------------------------------------------------
    ! Set lower and upper fractional limits of the convective mass flux
    !
    echam_conv_config%cmfcmin = 1.e-10_wp
    echam_conv_config%cmfcmax = 1.0_wp

    !------------------------------------------------------------------------
    ! Determine highest level *nmctop* for cloud base of midlevel convection
    ! - highest level below 300 hPa
    !
    ! Use the vertical coordinate tables that describe the hybrid sigma 
    ! coordinate as follows:
    ! - Pressure sigma coord.: ph(jk) = a(jk) + p_sfc*b(jk) [Pa]
    ! - Height   sigma coord.: zh(jk) = a(jk) + z_sfc*b(jk) [m]
    !                                 = a(jk) for a flat surface at z_sfc = 0 m
    !
    ! The vertical arrays are ordered from the top of the model (tom) to the surface (sfc):
    ! - ph(1) = p_tom, ph(nlev+1) = p_sfc
    ! - zh(1) = z_tom, zh(nlev+1) = z_sfc
    !
    ! The a(:) coefficients at the top of the model, where the sigma portion b(:) is zero,
    ! can be used to distinguish the pressure and height sigma grid:
    ! - pressure sigma grid : a(1) < a(2)
    ! - height   sigma grid : a(1) > a(2)
    !
    IF (vct_a(1) < vct_a(2)) THEN
      !
      ! pressure sigma grid
      !
      ! half level pressure
      DO jk=1,nlev+1
        zph(jk) = vct_a(jk) + vct_b(jk)*p0sl_bg
      END DO
      !
    ELSE
      !
      ! height sigma grid
      !
      ! half level pressure assuming 7500 m scale height
      DO jk=1,nlev+1
        zph(jk) = p0sl_bg*EXP(-vct_a(jk)/7500._wp) 
      END DO
      !
    END IF
    !
    ! full level pressure and eta
    DO jk = 1, nlev
      zp(jk)   = (zph(jk)+zph(jk+1))*0.5_wp
      zeta(jk) = zp(jk)/p0ref
    END DO
    !
    ! search for 300 hPa level
    DO jk = 1, nlev
      echam_conv_config%nmctop = jk
      IF(zp(jk) > 30000.0_wp) EXIT
    END DO

    !------------------------------------------------------------------------
    ! Set evaporation coefficient for kuo0
    !
    ALLOCATE( echam_conv_config%cevapcu(nlev),STAT=istat )
    IF (istat/=SUCCESS) CALL finish(TRIM(routine),'allocation of cevapcu failed')
    !
    DO jk = 1,nlev
       ztmp = 1.E3_wp/(38.3_wp*0.293_wp)*SQRT(zeta(jk))
       echam_conv_config%cevapcu(jk) = 1.93E-6_wp*261._wp*SQRT(ztmp)*0.5_wp/grav
    END DO

    !------------------------------------------------------------------------
    ! Print the configuration on stdio
    !
    CALL message('','')
    CALL message('','------- configuration of the ECHAM convection scheme --------')
    CALL message('','')
    CALL print_value(' lmfmid   ', echam_conv_config% lmfmid  )
    CALL print_value(' lmfpen   ', echam_conv_config% lmfpen  )
    CALL print_value(' lmfdd    ', echam_conv_config% lmfdd   )
    CALL print_value(' lmfdudv  ', echam_conv_config% lmfdudv )
    CALL message('','')
    CALL print_value(' entrscv  ', echam_conv_config% entrscv )
    CALL print_value(' entrmid  ', echam_conv_config% entrmid )
    CALL print_value(' entrpen  ', echam_conv_config% entrpen )
    CALL print_value(' entrdd   ', echam_conv_config% entrdd  )
    CALL message('','')
    CALL print_value(' cprcon   ', echam_conv_config% cprcon  )
    CALL print_value(' cmfctop  ', echam_conv_config% cmfctop )
    CALL print_value(' cmfdeps  ', echam_conv_config% cmfdeps )
    CALL message('','')
    CALL print_value(' cminbuoy ', echam_conv_config% cminbuoy)
    CALL print_value(' cmaxbuoy ', echam_conv_config% cmaxbuoy)
    CALL print_value(' cbfac    ', echam_conv_config% cbfac   )
    CALL print_value(' centrmax ', echam_conv_config% centrmax)
    CALL message('','')
    CALL print_value(' cmftau   ', echam_conv_config% cmftau  )
    CALL message('','')
    CALL message('','---------------------------')
    CALL message('','')
    CALL print_value(' nmctop   ', echam_conv_config% nmctop )
    CALL print_value(' cmfcmin  ', echam_conv_config% cmfcmin )
    CALL print_value(' cmfcmax  ', echam_conv_config% cmfcmax )
    CALL message('','')

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
