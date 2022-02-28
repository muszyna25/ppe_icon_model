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
  USE mo_exception,              ONLY: finish
#ifdef __ECRAD
  USE mo_ecrad,                  ONLY: t_ecrad_aerosol_type, t_ecrad_conf, t_opt_ptrs
#endif

  USE mo_aerosol_util,           ONLY: tegen_scal_factors

  USE mo_art_radiation_interface, ONLY: art_rad_aero_interface

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_prep_aerosol'

#ifdef __ECRAD
  PUBLIC :: nwp_ecrad_prep_aerosol

INTERFACE nwp_ecrad_prep_aerosol
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_constant
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_tegen
  MODULE PROCEDURE nwp_ecrad_prep_aerosol_td
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
    &                                          od_lw, ssa_lw, g_lw, od_sw, ssa_sw, g_sw, use_acc)
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf                        !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    REAL(wp), INTENT(in), OPTIONAL :: &
      &  od_lw, ssa_lw, g_lw,   & !< Optical depth, single scattering albedo, assymetry factor long wave
      &  od_sw, ssa_sw, g_sw      !< Optical depth, single scattering albedo, assymetry factor short wave
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
    LOGICAL                  :: lacc

    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if

#ifdef _OPENACC
    IF (lacc) CALL finish('nwp_ecrad_prep_aerosol_constant',' not ported to gpu')
#endif

    IF (ecrad_conf%do_lw) THEN
      ecrad_aerosol%od_lw(:,:,:)  = 0._wp
      ecrad_aerosol%ssa_lw(:,:,:) = 0._wp
      ecrad_aerosol%g_lw(:,:,:)   = 0._wp      
      IF ( PRESENT(od_lw) ) THEN 
        ecrad_aerosol%od_lw(:,:,:)  = od_lw
      ENDIF
      IF ( PRESENT(ssa_lw) ) THEN 
        ecrad_aerosol%ssa_lw(:,:,:) = ssa_lw
      ENDIF
      IF ( PRESENT(g_lw) ) THEN  
        ecrad_aerosol%g_lw(:,:,:)   = g_lw
      ENDIF
    ENDIF

    IF (ecrad_conf%do_sw) THEN
      ecrad_aerosol%od_sw(:,:,:)  = 0._wp
      ecrad_aerosol%ssa_sw(:,:,:) = 0._wp
      ecrad_aerosol%g_sw(:,:,:)   = 0._wp
      IF ( PRESENT(od_sw) ) THEN 
        ecrad_aerosol%od_sw(:,:,:)  = od_sw
      ENDIF
      IF ( PRESENT(ssa_sw) ) THEN 
        ecrad_aerosol%ssa_sw(:,:,:) = ssa_sw
      ENDIF
      IF ( PRESENT(g_sw) ) THEN  
        ecrad_aerosol%g_sw(:,:,:)   = g_sw
      ENDIF
    ENDIF

  END SUBROUTINE nwp_ecrad_prep_aerosol_constant
  !---------------------------------------------------------------------------------------



  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_tegen
  !! Prepare aerosol from Tegen climatology for ecRad. All the necessary
  !! information on the vertical and spatial distribution of the aerosol is given by the
  !! preprocessed fields zaeq1 - zaeq5. Code taken and adapted from rrtm, module mo_radiation
  !!
  !! @par Revision History
  !! Initial release by Daniel Rieger, Deutscher Wetterdienst, Offenbach (2021-10-20)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_tegen ( slev, nlev, i_startidx, i_endidx,       &
    &                                       zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,      &
    &                                       ecrad_conf, ecrad_aerosol, use_acc )
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
      &  ecrad_conf               !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol            !< ecRad aerosol information (input)
    LOGICAL, INTENT(in), OPTIONAL :: use_acc
! Local variables
    REAL(wp)                 :: &
      &  tau_abs, tau_sca         !< Absorption and scattering optical depth
    REAL(wp), POINTER        :: &
      &  scal_abs(:,:),         & !< Scaling factor absorption
      &  scal_sct(:,:),         & !< Scaling factor scattering
      &  scal_asy(:,:)            !< Scaling factor asymmetry
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//'::nwp_ecrad_prep_aerosol_tegen' 
    INTEGER                  :: &
      &  jc, jk, jband,         & !< Loop indices
      &  jband_shift              !< Band index in container (for shortwave: shifted by n_bands_lw)
    LOGICAL                  :: lacc

    if (present(use_acc)) then
      lacc = use_acc
    else
      lacc = .false.
    end if

    scal_abs => tegen_scal_factors%absorption
    scal_sct => tegen_scal_factors%scattering
    scal_asy => tegen_scal_factors%asymmetry

    !$ACC DATA PRESENT(ecrad_conf, ecrad_aerosol, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, scal_abs, scal_sct, scal_asy) IF(lacc)

! LONGWAVE
    IF (ecrad_conf%do_lw) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3) 
      DO jk = slev, nlev
!NEC$ nointerchange
!NEC$ nounroll
        DO jband = 1, ecrad_conf%n_bands_lw
          DO jc = i_startidx, i_endidx
            ! LW optical thickness
            ecrad_aerosol%od_lw (jband,jk,jc) =  zaeq1(jc,jk) * scal_abs(jband,1) &
              &                                + zaeq2(jc,jk) * scal_abs(jband,2) &
              &                                + zaeq3(jc,jk) * scal_abs(jband,3) &
              &                                + zaeq4(jc,jk) * scal_abs(jband,4) &
              &                                + zaeq5(jc,jk) * scal_abs(jband,5)
            ! No scattering at aerosol in longwave
            ecrad_aerosol%ssa_lw(jband,jk,jc) = 0._wp
            ecrad_aerosol%g_lw  (jband,jk,jc) = 0._wp
          ENDDO ! jc
        ENDDO ! jband
      ENDDO ! jk
      !$ACC END PARALLEL
    ENDIF

! SHORTWAVE
    IF (ecrad_conf%do_sw) THEN
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) IF(lacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(3) PRIVATE(jband_shift, tau_abs, tau_sca)
      DO jk = slev, nlev
!NEC$ nointerchange
        DO jband = 1, ecrad_conf%n_bands_sw
          DO jc = i_startidx, i_endidx

            jband_shift = ecrad_conf%n_bands_lw + jband

            ! SW absorption optical depth
            tau_abs =  zaeq1(jc,jk) * scal_abs(jband_shift,1) &
              &      + zaeq2(jc,jk) * scal_abs(jband_shift,2) &
              &      + zaeq3(jc,jk) * scal_abs(jband_shift,3) &
              &      + zaeq4(jc,jk) * scal_abs(jband_shift,4) &
              &      + zaeq5(jc,jk) * scal_abs(jband_shift,5)
            ! SW scattering optical depth
            tau_sca =  zaeq1(jc,jk) * scal_sct(jband_shift,1) &
              &      + zaeq2(jc,jk) * scal_sct(jband_shift,2) &
              &      + zaeq3(jc,jk) * scal_sct(jband_shift,3) &
              &      + zaeq4(jc,jk) * scal_sct(jband_shift,4) &
              &      + zaeq5(jc,jk) * scal_sct(jband_shift,5)
    
            ! Total optical depth for band jband
            ecrad_aerosol%od_sw (jband,jk,jc) = tau_abs + tau_sca
            
            ! Bulk SW single scattering albedo for band jband
            ecrad_aerosol%ssa_sw(jband,jk,jc) = tau_sca / ( tau_abs + tau_sca )
    
            ! Bulk SW asymmetry factor
            ecrad_aerosol%g_sw  (jband,jk,jc) =                                        &
              & (   zaeq1(jc,jk) * scal_sct(jband_shift,1) * scal_asy(jband_shift,1)   &
              &   + zaeq2(jc,jk) * scal_sct(jband_shift,2) * scal_asy(jband_shift,2)   &
              &   + zaeq3(jc,jk) * scal_sct(jband_shift,3) * scal_asy(jband_shift,3)   &
              &   + zaeq4(jc,jk) * scal_sct(jband_shift,4) * scal_asy(jband_shift,4)   &
              &   + zaeq5(jc,jk) * scal_sct(jband_shift,5) * scal_asy(jband_shift,5) ) / tau_sca
          ENDDO ! jc
        ENDDO ! jband
      ENDDO ! jk
      !$ACC END PARALLEL
    ENDIF

    !$ACC END DATA

  END SUBROUTINE nwp_ecrad_prep_aerosol_tegen

  !---------------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_td
  !! Time-dependent aerosol
  !! @par Revision History
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_td (slev, nlev, i_startidx, i_endidx, &
    &                                   opt_ptrs_lw, opt_ptrs_sw,         &
    &                                   ecrad_conf, ecrad_aerosol)

    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx     !< Start and end index of horizontal loop

    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf               !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol            !< ecRad aerosol information (input)

    INTEGER                  :: &
      &  jc, jk, jband            !< Loop indices

    TYPE(t_opt_ptrs), DIMENSION(ecrad_conf%n_bands_lw), INTENT(in):: opt_ptrs_lw
    TYPE(t_opt_ptrs), DIMENSION(ecrad_conf%n_bands_sw), INTENT(in):: opt_ptrs_sw

   ! LONGWAVE
    IF (ecrad_conf%do_lw) THEN
      DO jband = 1, ecrad_conf%n_bands_lw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ! LW optical thickness
            ecrad_aerosol%od_lw  (jband,jk,jc) = opt_ptrs_lw(jband)%ptr_od(jc,jk)
            ! No scattering at aerosol in longwave
            ecrad_aerosol%ssa_lw (jband,jk,jc) = 0._wp
            ecrad_aerosol%g_lw   (jband,jk,jc) = 0._wp
          ENDDO ! jc
        ENDDO   ! jk
      ENDDO     ! jband
    ENDIF

   !SHORTWAVE
    IF (ecrad_conf%do_sw) THEN
      DO jband = 1, ecrad_conf%n_bands_sw
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ecrad_aerosol%od_sw  (jband,jk,jc) = opt_ptrs_sw(jband)%ptr_od  (jc,jk)
            ecrad_aerosol%ssa_sw (jband,jk,jc) = opt_ptrs_sw(jband)%ptr_ssa (jc,jk)
            ecrad_aerosol%g_sw   (jband,jk,jc) = opt_ptrs_sw(jband)%ptr_g   (jc,jk)
          ENDDO ! jc
        ENDDO   ! jk
      ENDDO     ! jband
    ENDIF
  END SUBROUTINE nwp_ecrad_prep_aerosol_td
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  !! SUBROUTINE nwp_ecrad_prep_aerosol_art
  !! Preparation of aerosol optical properties from ART for ecRad. The Subroutine reshapes the ART
  !! output variables and reverses them vertically. 
  !!
  !! @par Revision History
  !! Initial release by Anika Rohde, KIT, Karlsruhe (2021-09-15)
  !!
  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_ecrad_prep_aerosol_art ( slev, nlev, i_startidx, i_endidx, jb, jg,   &
    &                                     nproma, zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,  &
    &                                     ecrad_conf, ecrad_aerosol )
    INTEGER, INTENT(in)      :: &
      &  slev, nlev,            & !< Start and end index of vertical loop
      &  i_startidx, i_endidx,  & !< Start and end index of horizontal loop
      &  jb, jg, nproma           !< Block and domain index
    REAL(wp), INTENT(in)     :: &
      &  zaeq1(:,:),            & !< Tegen optical thicknesses       1: continental
      &  zaeq2(:,:),            & !< relative to 550 nm, including   2: maritime
      &  zaeq3(:,:),            & !< a vertical profile              3: desert
      &  zaeq4(:,:),            & !< for 5 different                 4: urban
      &  zaeq5(:,:)               !< aerosol species.                5: stratospheric background
    TYPE(t_ecrad_conf),        INTENT(in)    :: &
      &  ecrad_conf                        !< ecRad configuration object
    TYPE(t_ecrad_aerosol_type),INTENT(inout) :: &
      &  ecrad_aerosol                     !< ecRad aerosol information (input)
    ! Local variables
    REAL(wp)                                                :: &
      &  tau_lw_art(nproma, nlev, ecrad_conf%n_bands_lw),     &  !< AOD LW (vertically reversed)      
      &  tau_sw_art(nproma, nlev, ecrad_conf%n_bands_sw),     &  !< AOD SW (vertically reversed)      
      &  ssa_sw_art(nproma, nlev, ecrad_conf%n_bands_sw),     &  !< SSA SW (vertically reversed)       
      &  g_sw_art(nproma,   nlev, ecrad_conf%n_bands_sw)         !< Assymetry parameter SW (vertically reversed)                  
    INTEGER                  :: &
      &  jc, jk, jband, jkb          !< Loop indices
  
    CALL art_rad_aero_interface(zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,   & !< Tegen aerosol
      &                         tegen_scal_factors%absorption,   & !
      &                         tegen_scal_factors%scattering,   & !< Tegen coefficients
      &                         tegen_scal_factors%asymmetry,    & !
      &                         jg, jb, slev, nlev,              & !< Indices domain, block, level
      &                         i_startidx, i_endidx,            & !< Indices nproma loop
      &                         ecrad_conf%n_bands_lw,           & !< Number of SW bands
      &                         ecrad_conf%n_bands_sw,           & !< Number of LW bands
      &                         tau_lw_art,                      & !< OUT: Optical depth LW
      &                         tau_sw_art,                      & !< OUT: Optical depth SW
      &                         ssa_sw_art,                      & !< OUT: SSA SW
      &                         g_sw_art)                          !< OUT: Assymetry parameter SW

    ! ART gives vertically reversed optical properties
    ! Furthermore, the shape of ecrad variables (bands, nlev, nproma) 
    ! and art variables (nproma, nlev, bands) differ

    DO jc = i_startidx, i_endidx 
      DO jk = slev, nlev
        jkb = nlev+1-jk
        DO jband = 1, ecrad_conf%n_bands_lw
          ecrad_aerosol%od_lw(jband, jk, jc)  = tau_lw_art(jc, jkb, jband)
        ENDDO
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx
      DO jk = slev, nlev
        jkb = nlev+1-jk
        DO jband = 1, ecrad_conf%n_bands_sw
          ecrad_aerosol%od_sw(jband, jk, jc)  = tau_sw_art(jc, jkb, jband)
          ecrad_aerosol%ssa_sw(jband, jk, jc) = ssa_sw_art(jc, jkb, jband)
          ecrad_aerosol%g_sw(jband, jk, jc)   = g_sw_art(jc, jkb, jband)
        ENDDO
      ENDDO
    ENDDO 

    ecrad_aerosol%ssa_lw(:,:,:) = 0._wp
    ecrad_aerosol%g_lw  (:,:,:) = 0._wp
   
  END SUBROUTINE nwp_ecrad_prep_aerosol_art
#endif
END MODULE mo_nwp_ecrad_prep_aerosol
