!>
!! This module contains routines for the vertical interpolation of
!! surface/soil fields provided by external analyses to the ICON grid
!!
!! Soil moisture is read from IFS2ICON as soil moisture index and then
!! converted back to soil moisture mass [m] using TERRA soil types.
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! Guenther Zaengl, DWD (2011-07-29):                    first version
!! Martin Koehler and Juergen Helmert, DWD (2011-08-12): soil moisture index conversion
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
MODULE mo_nwp_sfc_interp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma 
  USE mo_prepicon_config,     ONLY: nlevsoil_in, nlev_in
  USE mo_prepicon_types,      ONLY: t_prepicon_state
  USE mo_lnd_nwp_config,      ONLY: nlev_soil
  USE mo_impl_constants,      ONLY: zml_soil
  USE mo_physical_constants,  ONLY: grav
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_exception,           ONLY: message, message_text, finish

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: process_sfcfields

CONTAINS


  !-------------
  !>
  !! SUBROUTINE process_sfcfields
  !! Routine to convert surface fields interpolated horizontally by IFS2ICON
  !! to the ICON prognostic variables. Important ingredients are 
  !! - height adjustment of temperatures (partly done)
  !! - vertical interpolation of soil temperature and moisture 
  !! - conversion of soil moisture information
  !! - height adjustment of snow cover information (not yet done)
  !!
  !! Other open items - just to document them somewhere
  !! - IFS2ICON needs to take into account land-sea-mask information for horizontal
  !!   interpolation of surface fields
  !! - And, the most complicated problem, lakes/islands not present at all in the
  !!   source data, is not yet addressed at all here!

   SUBROUTINE process_sfcfields(p_patch, prepicon)


    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_prepicon_state), INTENT(INOUT)    :: prepicon

    ! LOCAL VARIABLES

    INTEGER  :: jg, jb, jk, jc, jk1, idx0(nlev_soil-1)
    INTEGER  :: nlen, nlev
    ! Soil layer depths in IFS
    REAL(wp) :: zsoil_ifs(4)=(/ 0.07_wp,0.21_wp,0.72_wp,1.89_wp/)
    REAL(wp) :: dzsoil_icon(8)=(/ 0.01_wp,0.02_wp,0.06_wp,0.18_wp,&
      & 0.54_wp,1.62_wp,4.86_wp,14.58_wp/)
    ! Standard atmosphere vertical temperature gradient
    REAL(wp) :: dtdz_clim = -6.5e-3_wp

    REAL(wp) :: tcorr1(nproma),tcorr2(nproma),wfac,wfac_vintp(nlev_soil-1),wfac_snow,snowdep

!-------------------------------------------------------------------------

    nlev = p_patch%nlev
    jg   = p_patch%id


    ! Vertical interpolation indices and weights
    DO jk = 1, nlev_soil-1
      IF (zml_soil(jk) < zsoil_ifs(1)) THEN
        idx0(jk)       = 0
        wfac_vintp(jk) = 1._wp - zml_soil(jk)/zsoil_ifs(1)
      ELSE IF (zml_soil(jk) > zsoil_ifs(nlevsoil_in)) THEN
        idx0(jk)       = nlevsoil_in
        wfac_vintp(jk) = 1._wp - (zml_soil(jk)-zsoil_ifs(nlevsoil_in))/&
                                 (zml_soil(8) -zsoil_ifs(nlevsoil_in))
      ELSE
        DO jk1 = 1, nlevsoil_in-1
          IF (zml_soil(jk) > zsoil_ifs(jk1) .AND. zml_soil(jk) <= zsoil_ifs(jk1+1)) THEN
            idx0(jk)       = jk1
            wfac_vintp(jk) = (zsoil_ifs(jk1+1)-zml_soil(jk)) / &
                             (zsoil_ifs(jk1+1)-zsoil_ifs(jk1))
          ENDIF
        ENDDO
      ENDIF
    ENDDO

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,wfac,tcorr1,tcorr2,snowdep,wfac_snow)
    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      ! 2D fields that do not require height adjustment
      ! these fields are simply copied
      ! Note: land-sea-mask is currently unused
      DO jc = 1, nlen
        prepicon%sfc%skinres(jc,jb) = prepicon%sfc_in%skinres(jc,jb)
        prepicon%sfc%ls_mask(jc,jb) = prepicon%sfc_in%ls_mask(jc,jb)
        IF (prepicon%sfc_in%seaice(jc,jb) >= 0._wp) THEN
          prepicon%sfc%seaice(jc,jb)  = prepicon%sfc_in%seaice(jc,jb) 
        ELSE
          prepicon%sfc%seaice(jc,jb)  = -999.9_wp 
        ENDIF
        IF (prepicon%sfc_in%sst(jc,jb) > 10._wp) THEN
          prepicon%sfc%sst(jc,jb)  = prepicon%sfc_in%sst(jc,jb) 
        ELSE
          prepicon%sfc%sst(jc,jb)  = -999.9_wp 
        ENDIF
      ENDDO

      ! 2D fields that require height adjustment
      ! unfortunately, skin temperature is the only variable for which it is
      ! intuitively clear what to do ...
      DO jc = 1, nlen
        ! Adjust skin temperature with the difference between the atmospheric
        ! temperatures at the lowest model level
        prepicon%sfc%tskin(jc,jb)      = prepicon%sfc_in%tskin(jc,jb) +       &
          (prepicon%atm%temp(jc,nlev,jb) - prepicon%atm_in%temp(jc,nlev_in,jb))

        ! Height adjustment for snow variables is not yet implemented
        prepicon%sfc%tsnow(jc,jb)    = prepicon%sfc_in%tsnow(jc,jb) 
        prepicon%sfc%snowweq(jc,jb)  = prepicon%sfc_in%snowweq(jc,jb) 
        prepicon%sfc%snowdens(jc,jb) = prepicon%sfc_in%snowdens(jc,jb) 
        prepicon%sfc%snowalb(jc,jb)  = prepicon%sfc_in%snowalb(jc,jb) 
      ENDDO

      ! Height adjustment of soil temperatures
      DO jc = 1, nlen
        ! correction for ground level
        tcorr1(jc) = prepicon%atm%temp(jc,nlev,jb) - prepicon%atm_in%temp(jc,nlev_in,jb)
        ! climatological correction for deep soil levels
        tcorr2(jc) = dtdz_clim*(prepicon%topography_c(jc,jb)-prepicon%sfc_in%phi(jc,jb)/grav)
      ENDDO

      DO jk = 1, nlevsoil_in
        wfac = REAL(jk,wp)/REAL(nlevsoil_in,wp) ! weight for climatological correction
        DO jc = 1, nlen
          prepicon%sfc_in%tsoil(jc,jb,jk) = prepicon%sfc_in%tsoil(jc,jb,jk) + &
            wfac*tcorr2(jc) + (1._wp-wfac)*tcorr1(jc)
        ENDDO
      ENDDO

      ! Fill extra levels of incoming data to simplify vertical interpolation
      DO jc = 1, nlen
        ! factor depending on snow depth to weight initialization of soil top temperature
        ! between skin temperature and soil level 1 temperature
        ! For a snow depth of more than 25 cm, it is assumed that the soil top temperature
        ! is the same as the soil level 1 temperature
        snowdep = 1000._wp*prepicon%sfc_in%snowweq(jc,jb) / & ! snow depth in m
                  MAX(25._wp,prepicon%sfc_in%snowdens(jc,jb))
        wfac_snow = SQRT(MIN(1._wp,4._wp*snowdep))
        prepicon%sfc_in%tsoil(jc,jb,0) = (1._wp-wfac_snow)*prepicon%sfc%tskin(jc,jb) + &
                                          wfac_snow*prepicon%sfc_in%tsoil(jc,jb,1) ! already height-adjusted

        prepicon%sfc%tsoil(jc,jb,0)    = prepicon%sfc_in%tsoil(jc,jb,0) ! copy soil-top temperature

        prepicon%sfc_in%wsoil(jc,jb,0) = prepicon%sfc_in%wsoil(jc,jb,1) ! no-gradient condition for moisture

        ! outgoing tsoil(nlev_soil) has been initialized with the external parameter field t_cl before
        prepicon%sfc_in%tsoil(jc,jb,nlevsoil_in+1) = prepicon%sfc%tsoil(jc,jb,nlev_soil)
        ! assume no-gradient condition for soil moisture
        prepicon%sfc_in%wsoil(jc,jb,nlevsoil_in+1) = prepicon%sfc_in%wsoil(jc,jb,nlevsoil_in)
      ENDDO

      ! Vertical interpolation of multi-layer soil fields from IFS levels to TERRA levels
      DO jk = 1, nlev_soil-1
        DO jc = 1, nlen
          prepicon%sfc%tsoil(jc,jb,jk) = wfac_vintp(jk) *prepicon%sfc_in%tsoil(jc,jb,idx0(jk))+ &
                                  (1._wp-wfac_vintp(jk))*prepicon%sfc_in%tsoil(jc,jb,idx0(jk)+1)
          prepicon%sfc%wsoil(jc,jb,jk) = wfac_vintp(jk) *prepicon%sfc_in%wsoil(jc,jb,idx0(jk))+ &
                                  (1._wp-wfac_vintp(jk))*prepicon%sfc_in%wsoil(jc,jb,idx0(jk)+1)
        ENDDO
      ENDDO


      ! Conversion of IFS soil moisture index (vertically interpolated) into TERRA soil moisture [m]
      !   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
      !   safety: min=air dryness point, max=pore volume
      DO jk = 1, nlev_soil-1
        DO jc = 1, nlen
          IF(ext_data(jg)%atm%soiltyp(jc,jb) == 3) prepicon%sfc%wsoil(jc,jb,jk) = &
            & dzsoil_icon(jk) * &
            & MIN(0.364_wp,&
            & MAX((prepicon%sfc%wsoil(jc,jb,jk)*(0.196_wp - 0.042_wp) + 0.042_wp),0.012_wp))

          IF(ext_data(jg)%atm%soiltyp(jc,jb) == 4) prepicon%sfc%wsoil(jc,jb,jk) = &
            & dzsoil_icon(jk) * &
            & MIN(0.445_wp,&
            & MAX((prepicon%sfc%wsoil(jc,jb,jk)*(0.26_wp  - 0.1_wp  ) + 0.1_wp)  ,0.03_wp ))

          IF(ext_data(jg)%atm%soiltyp(jc,jb) == 5) prepicon%sfc%wsoil(jc,jb,jk) = &
            & dzsoil_icon(jk) * &
            & MIN(0.455_wp,&
            & MAX((prepicon%sfc%wsoil(jc,jb,jk)*(0.34_wp  - 0.11_wp ) + 0.11_wp) ,0.035_wp))

          IF(ext_data(jg)%atm%soiltyp(jc,jb) == 6) prepicon%sfc%wsoil(jc,jb,jk) = &
            & dzsoil_icon(jk) * &
            & MIN(0.475_wp,&
            & MAX((prepicon%sfc%wsoil(jc,jb,jk)*(0.37_wp  - 0.185_wp) + 0.185_wp),0.06_wp ))

          IF(ext_data(jg)%atm%soiltyp(jc,jb) == 7) prepicon%sfc%wsoil(jc,jb,jk) = &
            & dzsoil_icon(jk) * &
            & MIN(0.507_wp,&
            & MAX((prepicon%sfc%wsoil(jc,jb,jk)*(0.463_wp - 0.257_wp) + 0.257_wp),0.065_wp))

          IF(ext_data(jg)%atm%soiltyp(jc,jb) == 8) prepicon%sfc%wsoil(jc,jb,jk) = &
            & dzsoil_icon(jk) * &
            & MIN(0.863_wp,&
            & MAX((prepicon%sfc%wsoil(jc,jb,jk)*(0.763_wp - 0.265_wp) + 0.265_wp),0.098_wp))

          ! We need to catch problematic coast cases: IFS-ocean, ICON-land - &
          ! for moisture and temperature 
          prepicon%sfc%wsoil(jc,jb,jk) =  &
            &               MIN(dzsoil_icon(jk), MAX(0.0_wp, prepicon%sfc%wsoil(jc,jb,jk)))
        ENDDO
      ENDDO

      ! assume no-gradient condition for soil moisture reservoir layer
      DO jc = 1, nlen
        prepicon%sfc%wsoil(jc,jb,nlev_soil) = prepicon%sfc%wsoil(jc,jb,nlev_soil-1)*          &
                                                dzsoil_icon(nlev_soil)/dzsoil_icon(nlev_soil-1)
      ENDDO

    ENDDO
!$OMP END DO 
!$OMP END PARALLEL


  END SUBROUTINE process_sfcfields

END MODULE mo_nwp_sfc_interp
