!>
!! @brief Parameters for the boundary layer parameterization. 
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-08-31)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_echam_vdiff_params

  USE mo_kind,              ONLY: wp
 !USE mo_impl_constants,    ONLY: SUCCESS
  USE mo_exception,         ONLY: print_value, message
  USE mo_namelist,          ONLY: position_nml, POSITIONED
  USE mo_io_units,          ONLY: nnml
  USE mo_physical_constants,ONLY: grav

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: clam, ckap, cchar, cb, cc, cfreec, cgam     !< parameters
  PUBLIC :: eps_shear, eps_corio, tke_min, z0m_min      !< parameters
  PUBLIC :: z0m_ice, z0m_oce
  PUBLIC :: chneu, shn, smn, da1, custf, cwstf          !< parameters
  PUBLIC :: cons2, cons25, cons5                        !< parameters
  PUBLIC :: cvdifts, tpfac1, tpfac2, tpfac3, tpfac4     !< parameters
  PUBLIC :: itop, itopp1, ibl, iblm1, iblmin, iblmax    !< parameters
  PUBLIC :: echam_vdiff_ctl                             !< namelist
  PUBLIC :: setup_vdiff, init_vdiff_params              !< subroutine

  !--------------------
  ! Namelist variables
  !--------------------

  LOGICAL,PUBLIC :: lsfc_mom_flux   !< switch on/off surface momentum flux
  LOGICAL,PUBLIC :: lsfc_heat_flux  !< switch on/off surface heat flux
                                    !< (sensible AND latent)

  NAMELIST/echam_vdiff_ctl/ lsfc_mom_flux, lsfc_heat_flux

  !-------------------
  ! Module parameters
  !-------------------

  REAL(wp),PARAMETER :: clam    = 150._wp      !< asymptotic mixing length for momentum
  REAL(wp),PARAMETER :: ckap    = 0.4_wp       !< karman constant.
  REAL(wp),PARAMETER :: cchar   = 0.018_wp     !< charnock constant.
  REAL(wp),PARAMETER :: cb      = 5._wp        !< stability parameter near neutrality.
  REAL(wp),PARAMETER :: cc      = 5._wp        !< stability parameter for unstable cases.
! REAL(wp),PARAMETER :: cd      = 5._wp        !< stability parameter for stable cases.
  REAL(wp),PARAMETER :: cfreec  = 0.001_wp     !< free convection parameter
  REAL(wp),PARAMETER :: cgam    = 1.25_wp      !< free convection parameter

  REAL(wp),PARAMETER :: eps_shear = 1.e-5_wp   !< zepshr in sbr. vdiff of ECHAM6
  REAL(wp),PARAMETER :: eps_corio = 5.e-5_wp   !< zepcor in sbr. vdiff of ECHAM6
! REAL(wp),PARAMETER :: tke_min   = 1.e-10_wp  !< ztkemin in sbr. vdiff of ECHAM6
  REAL(wp)           :: tke_min   = 1.e-10_wp  !< ztkemin in sbr. vdiff of ECHAM6
! REAL(wp),PARAMETER :: tke_min   = 1.e-4_wp   !< ztkemin in sbr. vdiff of ECHAM6
  REAL(wp),PARAMETER :: z0m_min   = 1.5e-5_wp  !< zepzzo in sbr. vdiff of ECHAM5
  REAL(wp),PARAMETER :: z0m_ice   = 1e-3_wp    !< cz0ice in sbr. vdiff of ECHAM5
  REAL(wp),PARAMETER :: z0m_oce   = 5e-4_wp    !< see mo_surface_ocean.f90 of ECHAM6

  REAL(wp),PARAMETER :: chneu  = 0.3_wp

  REAL(wp),PARAMETER :: cons2  = 0.5_wp*ckap/grav
  REAL(wp),PARAMETER :: cons25 = cons2/clam
  REAL(wp),PARAMETER :: cons5  = 3._wp*cb*cc*grav**2

  ! Parameters related to time step weighting in *rhs* of *vdiff* and *scv*

  REAL(wp),PARAMETER :: cvdifts = 1.5_wp
  REAL(wp),PARAMETER :: tpfac1  = cvdifts
  REAL(wp),PARAMETER :: tpfac2  = 1._wp / tpfac1
  REAL(wp),PARAMETER :: tpfac3  = 1._wp - tpfac2
  REAL(wp),PARAMETER :: tpfac4  = 1._wp + tpfac3

  !-------------------
  ! Module variables
  !-------------------

  REAL(wp) :: shn
  REAL(wp) :: smn
  REAL(wp) :: da1

  INTEGER :: itop
  INTEGER :: itopp1
  INTEGER :: ibl
  INTEGER :: iblm1
  INTEGER :: iblmin
  INTEGER :: iblmax

  ! For the surface scheme of ECHAM5

  REAL(wp) :: custf, cwstf

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  SUBROUTINE setup_vdiff

    INTEGER :: ist

    ! Set default values

    lsfc_mom_flux  = .TRUE.
    lsfc_heat_flux = .TRUE.

    ! Read namelist (every CPU does this)

    CALL position_nml('echam_vdiff_ctl',status=ist)
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, echam_vdiff_ctl)
    END SELECT

    ! Check validity; send values to stdout

    CALL message('','')
    CALL message('','------- namelist echam_vdiff_ctl --------')

    CALL print_value(' lsfc_mom_flux  ',lsfc_mom_flux)
    CALL print_value(' lsfc_heat_flux ',lsfc_heat_flux)

    CALL message('','---------------------------')
    CALL message('','')

  END SUBROUTINE setup_vdiff
  !-------------
  !>
  !!
  SUBROUTINE init_vdiff_params( klev,klevp1,kvclev,vct )

    INTEGER,INTENT(IN)  :: klev, klevp1, kvclev
    REAL(wp),INTENT(IN) :: vct(2*kvclev)

    INTEGER  :: jk
    REAL(wp) :: zph(klevp1), zp(klev), zh(klev)

    !------------------

    shn    = 2.22_wp*0.22_wp*SQRT(2._wp)
    smn    = shn*1.24_wp*2.37_wp/3.69_wp
    da1    = 1._wp/smn**3
    custf  = 1._wp/smn**2
    cwstf  = 0.2_wp

    itop   = 1
    itopp1 = itop + 1

    ibl    = klev - 1    ! CO2 flux is distributed into lowest layers klev to ibl
    iblm1  = ibl - 1

    ! CO2 mixing in PBL.
    ! Compute lowest (klev-2, approx. 500m in L19/31) and highest
    !   (highest below 2km) level ktop.
    ! CO2 is ideally mixed between bottom and itop.

    iblmin = klev - 2

    ! Search for the highest model level below 2000 m:
    ! First compute the half level pressure values,
    ! assuming 101320 Pa surface pressure

    DO jk=1,klevp1
      zph(jk)=vct(jk)+vct(jk+kvclev)*101320.0_wp
    END DO

    ! Then the full level pressure

    DO jk = 1,klev
      zp(jk)=(zph(jk)+zph(jk+1))*0.5_wp
    END DO

    ! Compute the elevation with respect to the Earth's surface

    DO jk = 1,klev
      zh(jk)=(zph(klevp1)-zp(jk))/(grav*1.25_wp)
    END DO

    ! Search for highest level below 2000m

    DO jk = 1,klev
      iblmax=jk
      IF(zh(jk).LT.2000.0_wp) EXIT
    END DO

  END SUBROUTINE init_vdiff_params
  !-------------

END MODULE mo_echam_vdiff_params

