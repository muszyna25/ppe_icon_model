!>
!! This module prepares aerosol climatologies in a format that can be used by ecRad
!!
!! @author Daniel Rieger, Deutscher Wetterdienst, Offenbach
!!
!! @par Revision History
!! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (YYYY-MM-DD)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
#if defined __xlC__
@PROCESS SPILL(1058)
#endif
MODULE mo_nwp_ecrad_prep_aerosol

  USE mo_kind,                   ONLY: wp
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: t_ecrad_aerosol_type, t_ecrad_conf
#endif
  USE mo_impl_constants,         ONLY: max_char_length
  
  USE mo_aerosol_util,           ONLY: zaea_rrtm,zaes_rrtm,zaeg_rrtm

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_prep_aerosol'

#ifdef __ECRAD
  PUBLIC :: nwp_ecrad_prep_aerosol

INTERFACE nwp_ecrad_prep_aerosol
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_constant
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_tegen_tanre
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_art
END INTERFACE nwp_ecrad_prep_aerosol
  
  
CONTAINS

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_constant
  !! Prepare aerosol from constant values. If these optional values are not passed,
  !! the corresponding field in ecRad is set to zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-15)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_constant ( ecrad_conf, ecrad_aerosol,                &
    &                                          od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw)
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf                        !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    REAL(wp), INTENT(in), OPTIONAL :: &
      &  od_lw, ssa_lw, g_lw,   & !< Optical depth, single scattering albedo, assymetry factor long wave
      &  od_sw, ssa_sw, g_sw      !< Optical depth, single scattering albedo, assymetry factor short wave

    IF (ecrad_conf%do_lw) THEN
      ecrad_aerosol%od_lw(:,:,:)  = 0._wp
      ecrad_aerosol%ssa_lw(:,:,:) = 0._wp
      ecrad_aerosol%g_lw(:,:,:)   = 0._wp
      IF ( PRESENT(od_lw) )  ecrad_aerosol%od_lw(:,:,:)  = od_lw
      IF ( PRESENT(ssa_lw) ) ecrad_aerosol%ssa_lw(:,:,:) = ssa_lw
      IF ( PRESENT(g_lw) )   ecrad_aerosol%g_lw(:,:,:)   = g_lw
    ENDIF

    IF (ecrad_conf%do_sw) THEN
      ecrad_aerosol%od_sw(:,:,:)  = 0._wp
      ecrad_aerosol%ssa_sw(:,:,:) = 0._wp
      ecrad_aerosol%g_sw(:,:,:)   = 0._wp
      IF ( PRESENT(od_sw) )  ecrad_aerosol%od_sw(:,:,:)  = od_sw
      IF ( PRESENT(ssa_sw) ) ecrad_aerosol%ssa_sw(:,:,:) = ssa_sw
      IF ( PRESENT(g_sw) )   ecrad_aerosol%g_sw(:,:,:)   = g_sw
    ENDIF

  END SUBROUTINE nwp_ecrad_prep_aerosol_constant
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_tegen_tanre
  !! Prepare aerosol from Tanre and Tegen climatology for ecRad. All the necessary 
  !! information on the vertical and spatial distribution of the aerosol is given by the
  !! preprocessed fields zaeq1 - zaeq5. Code taken and adapted from rrtm, module mo_radiation
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2018-12-07)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_tegen_tanre ( slev, nlev, i_startidx, i_endidx,       &
    &                                             zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,      &
    &                                             ecrad_conf, ecrad_aerosol )
    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx     !< Start and end index of horizontal loop
    REAL(wp), INTENT(in)     :: &
      &  zaeq1(:,:),            & !< Tegen optical thicknesses         1: continental
      &  zaeq2(:,:),            & !<   relative to 550 nm, including   2: maritime
      &  zaeq3(:,:),            & !<   a vertical profile              3: desert
      &  zaeq4(:,:),            & !<   for 5 different                 4: urban
      &  zaeq5(:,:)               !<   aerosol species.                5: stratospheric background
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf                        !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
! Local variables
    REAL(wp)                 :: &
      &  tau_abs, tau_sca         !< Absorption and scattering optical depth
    CHARACTER(len=max_char_length), PARAMETER :: &
      &  routine = modname//'::nwp_ecrad_prep_aerosol_tegen' 
    INTEGER                  :: &
      &  jc, jk, jband,         & !< Loop indices
      &  jband_rrtm               !< Band index in rrtm container (for shortwave: shifted by nbands_lw)

! LONGWAVE
    IF (ecrad_conf%do_lw) THEN
      DO jband = 1, ecrad_conf%n_bands_lw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ! LW optical thickness
            ecrad_aerosol%od_lw (jband,jk,jc) =  zaeq1(jc,jk) * zaea_rrtm(jband,1) &
              &                                + zaeq2(jc,jk) * zaea_rrtm(jband,2) &
              &                                + zaeq3(jc,jk) * zaea_rrtm(jband,3) &
              &                                + zaeq4(jc,jk) * zaea_rrtm(jband,4) &
              &                                + zaeq5(jc,jk) * zaea_rrtm(jband,5)
            ! No scattering at aerosol in longwave
            ecrad_aerosol%ssa_lw(jband,jk,jc) = 0._wp
            ecrad_aerosol%g_lw  (jband,jk,jc) = 0._wp
          ENDDO ! jc
        ENDDO ! jk
      ENDDO ! jband
    ENDIF

! SHORTWAVE
    IF (ecrad_conf%do_sw) THEN
      DO jband = 1, ecrad_conf%n_bands_sw
        jband_rrtm = ecrad_conf%n_bands_lw + jband
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ! SW absorption optical depth
            tau_abs =  zaeq1(jc,jk) * zaea_rrtm(jband_rrtm,1) &
              &      + zaeq2(jc,jk) * zaea_rrtm(jband_rrtm,2) &
              &      + zaeq3(jc,jk) * zaea_rrtm(jband_rrtm,3) &
              &      + zaeq4(jc,jk) * zaea_rrtm(jband_rrtm,4) &
              &      + zaeq5(jc,jk) * zaea_rrtm(jband_rrtm,5)
            ! SW scattering optical depth
            tau_sca =  zaeq1(jc,jk) * zaes_rrtm(jband_rrtm,1) &
              &      + zaeq2(jc,jk) * zaes_rrtm(jband_rrtm,2) &
              &      + zaeq3(jc,jk) * zaes_rrtm(jband_rrtm,3) &
              &      + zaeq4(jc,jk) * zaes_rrtm(jband_rrtm,4) &
              &      + zaeq5(jc,jk) * zaes_rrtm(jband_rrtm,5)
    
            ! Total optical depth for band jband
            ecrad_aerosol%od_sw (jband,jk,jc) = tau_abs + tau_sca
    
            ! Bulk SW single scattering albedo for band jband
            ecrad_aerosol%ssa_sw(jband,jk,jc) = tau_sca / ( tau_abs + tau_sca )
    
            ! Bulk SW asymmetry factor
            ecrad_aerosol%g_sw  (jband,jk,jc) =                                        &
              & (   zaeq1(jc,jk) * zaes_rrtm(jband_rrtm,1) * zaeg_rrtm(jband_rrtm,1)   &
              &   + zaeq2(jc,jk) * zaes_rrtm(jband_rrtm,2) * zaeg_rrtm(jband_rrtm,2)   &
              &   + zaeq3(jc,jk) * zaes_rrtm(jband_rrtm,3) * zaeg_rrtm(jband_rrtm,3)   &
              &   + zaeq4(jc,jk) * zaes_rrtm(jband_rrtm,4) * zaeg_rrtm(jband_rrtm,4)   &
              &   + zaeq5(jc,jk) * zaes_rrtm(jband_rrtm,5) * zaeg_rrtm(jband_rrtm,5) ) / tau_sca
          ENDDO ! jc
        ENDDO ! jk
      ENDDO ! jband
    ENDIF

  END SUBROUTINE nwp_ecrad_prep_aerosol_tegen_tanre
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_art
  !! Plugin aerosol from ART
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2019-05-15)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_art ( slev, nlev, i_startidx, i_endidx, jb, jg,   &
    &                                     zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,          &
    &                                     ecrad_conf, ecrad_aerosol )
    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx,  & !< Start and end index of horizontal loop
      &  jb, jg                   !< Block and domain index
    REAL(wp), INTENT(in)     :: &
      &  zaeq1(:,:),            & !< Tegen optical thicknesses         1: continental
      &  zaeq2(:,:),            & !<   relative to 550 nm, including   2: maritime
      &  zaeq3(:,:),            & !<   a vertical profile              3: desert
      &  zaeq4(:,:),            & !<   for 5 different                 4: urban
      &  zaeq5(:,:)               !<   aerosol species.                5: stratospheric background
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf                        !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)

    ! CALL art interface: This has to be modified somehow as we are here inside the jb loop in
    !                     contrast to the call inside mo_radiation
    !                     Longwave scattering might also be a point to consider here
    !CALL art_rad_aero_interface(zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,   & !< Tegen aerosol
    !  &                         zaea_rrtm,zaes_rrtm,zaeg_rrtm,   & !< Tegen coefficients for wavelength bands
    !  &                         jg, jb, slev, nlev,              & !< Indices domain, block, level
    !  &                         i_startidx, i_endidx,            & !< Indices nproma loop
    !  &                         ecrad_conf%n_bands_lw,           & !< Number of SW bands
    !  &                         ecrad_conf%n_bands_sw,           & !< Number of LW bands
    !  &                         ecrad_aerosol%od_lw,             & !< OUT: Optical depth LW
    !  &                         ecrad_aerosol%od_sw,             & !< OUT: Optical depth SW
    !  &                         ecrad_aerosol%ssa_sw,            & !< OUT: Single scattering albedo SW
    !  &                         ecrad_aerosol%g_sw)                !< OUT: Assymetry parameter SW
  END SUBROUTINE nwp_ecrad_prep_aerosol_art
#endif
END MODULE mo_nwp_ecrad_prep_aerosol
