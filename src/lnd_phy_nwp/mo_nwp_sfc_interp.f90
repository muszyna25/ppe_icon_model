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
  USE mo_initicon_config,     ONLY: nlevsoil_in, nlev_in
  USE mo_nh_initicon_types,   ONLY: t_initicon_state
  USE mo_lnd_nwp_config,      ONLY: nlev_soil
  USE mo_impl_constants,      ONLY: zml_soil, dzsoil_icon => dzsoil
  USE mo_physical_constants,  ONLY: grav, dtdz_standardatm
  USE mo_phyparam_soil,       ONLY: cporv
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_exception,           ONLY: message, message_text, finish

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: process_sfcfields
  PUBLIC :: smi_to_sm_mass

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
  !!   interpolation of surface fields (already available)
  !! - And, the most complicated problem, lakes/islands not present at all in the
  !!   source data, is not yet addressed at all here!

   SUBROUTINE process_sfcfields(p_patch, initicon)


    TYPE(t_patch),          INTENT(IN)       :: p_patch
    TYPE(t_initicon_state), INTENT(INOUT)    :: initicon

    ! LOCAL VARIABLES
    CHARACTER(LEN=*), PARAMETER       :: routine = 'process_sfcfields'

    INTEGER  :: jg, jb, jk, jc, jk1, idx0(nlev_soil-1)
    INTEGER  :: nlen, nlev

    REAL(wp) :: tcorr1(nproma),tcorr2(nproma),wfac,wfac_vintp(nlev_soil-1),wfac_snow,snowdep

    REAL(wp) :: zwsoil(nproma)        ! local soil moisture field
    INTEGER  :: i_count, ic           ! counter
    LOGICAL  :: lerr                  ! error flag

    ! Soil layer depths in IFS
    REAL(wp), PARAMETER :: zsoil_ifs(4)=(/ 0.07_wp,0.21_wp,0.72_wp,1.89_wp/)

!-------------------------------------------------------------------------

    IF (nlev_in == 0) THEN
      CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF

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
!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,wfac,tcorr1,tcorr2,snowdep,wfac_snow,ic,i_count,zwsoil,lerr)
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
        initicon%sfc%skinres(jc,jb) = initicon%sfc_in%skinres(jc,jb)
        initicon%sfc%ls_mask(jc,jb) = initicon%sfc_in%ls_mask(jc,jb)
        IF (initicon%sfc_in%seaice(jc,jb) >= 0._wp) THEN
          initicon%sfc%seaice(jc,jb)  = initicon%sfc_in%seaice(jc,jb) 
        ELSE
          initicon%sfc%seaice(jc,jb)  = -999.9_wp 
        ENDIF
        ! Remark (GZ): the second condition is a workaround until a proper treatment of missing values
        !              becomes available in prep_icon
        IF (initicon%sfc_in%sst(jc,jb) > 10._wp .AND. initicon%sfc_in%sst(jc,jb) < 305._wp) THEN
          initicon%sfc%sst(jc,jb)  = initicon%sfc_in%sst(jc,jb) 
        ELSE
          initicon%sfc%sst(jc,jb)  = -999.9_wp 
        ENDIF
      ENDDO

      ! 2D fields that require height adjustment
      ! unfortunately, skin temperature is the only variable for which it is
      ! intuitively clear what to do ...
      DO jc = 1, nlen
        ! Adjust skin temperature with the difference between the atmospheric
        ! temperatures at the lowest model level
        initicon%sfc%tskin(jc,jb)      = initicon%sfc_in%tskin(jc,jb) +       &
          (initicon%atm%temp(jc,nlev,jb) - initicon%atm_in%temp(jc,nlev_in,jb))

        ! Height adjustment for snow variables is not yet implemented
        initicon%sfc%tsnow(jc,jb)    = initicon%sfc_in%tsnow(jc,jb) 
        initicon%sfc%snowweq(jc,jb)  = initicon%sfc_in%snowweq(jc,jb) 
        initicon%sfc%snowdens(jc,jb) = initicon%sfc_in%snowdens(jc,jb) 
        initicon%sfc%snowalb(jc,jb)  = initicon%sfc_in%snowalb(jc,jb) 
      ENDDO

      ! Height adjustment of soil temperatures
      DO jc = 1, nlen
        ! correction for ground level
        tcorr1(jc) = initicon%atm%temp(jc,nlev,jb) - initicon%atm_in%temp(jc,nlev_in,jb)
        ! climatological correction for deep soil levels
        tcorr2(jc) = dtdz_standardatm &
          &        * (initicon%topography_c(jc,jb)-initicon%sfc_in%phi(jc,jb)/grav)
      ENDDO

      DO jk = 1, nlevsoil_in
        wfac = REAL(jk,wp)/REAL(nlevsoil_in,wp) ! weight for climatological correction
        DO jc = 1, nlen
          initicon%sfc_in%tsoil(jc,jb,jk) = initicon%sfc_in%tsoil(jc,jb,jk) + &
            wfac*tcorr2(jc) + (1._wp-wfac)*tcorr1(jc)
        ENDDO
      ENDDO

      ! Fill extra levels of incoming data to simplify vertical interpolation
      DO jc = 1, nlen
        ! factor depending on snow depth to weight initialization of soil top temperature
        ! between skin temperature and soil level 1 temperature
        ! For a snow depth of more than 25 cm, it is assumed that the soil top temperature
        ! is the same as the soil level 1 temperature
        snowdep = 1000._wp*initicon%sfc_in%snowweq(jc,jb) / & ! snow depth in m
                  MAX(25._wp,initicon%sfc_in%snowdens(jc,jb))
        wfac_snow = SQRT(MIN(1._wp,4._wp*snowdep))
        initicon%sfc_in%tsoil(jc,jb,0) = (1._wp-wfac_snow)*initicon%sfc%tskin(jc,jb) + &
                                          wfac_snow*initicon%sfc_in%tsoil(jc,jb,1) ! already height-adjusted

        initicon%sfc%tsoil(jc,jb,0)    = initicon%sfc_in%tsoil(jc,jb,0) ! copy soil-top temperature

        initicon%sfc_in%wsoil(jc,jb,0) = initicon%sfc_in%wsoil(jc,jb,1) ! no-gradient condition for moisture

        ! outgoing tsoil(nlev_soil) has been initialized with the external parameter field t_cl before
        initicon%sfc_in%tsoil(jc,jb,nlevsoil_in+1) = initicon%sfc%tsoil(jc,jb,nlev_soil)
        ! assume no-gradient condition for soil moisture
        initicon%sfc_in%wsoil(jc,jb,nlevsoil_in+1) = initicon%sfc_in%wsoil(jc,jb,nlevsoil_in)
      ENDDO

      ! Vertical interpolation of multi-layer soil fields from IFS levels to TERRA levels
      DO jk = 1, nlev_soil-1
        DO jc = 1, nlen
          initicon%sfc%tsoil(jc,jb,jk) = wfac_vintp(jk) *initicon%sfc_in%tsoil(jc,jb,idx0(jk))+ &
                                  (1._wp-wfac_vintp(jk))*initicon%sfc_in%tsoil(jc,jb,idx0(jk)+1)
          initicon%sfc%wsoil(jc,jb,jk) = wfac_vintp(jk) *initicon%sfc_in%wsoil(jc,jb,idx0(jk))+ &
                                  (1._wp-wfac_vintp(jk))*initicon%sfc_in%wsoil(jc,jb,idx0(jk)+1)
        ENDDO
      ENDDO


      ! (re)-initialize error flag
      lerr = .FALSE.

      ! Conversion of IFS soil moisture index (vertically interpolated) into TERRA soil moisture [m]
      !   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
      !   safety: min=air dryness point, max=pore volume
      DO jk = 1, nlev_soil-1

        ! loop over target (ICON) land points only
        i_count = ext_data(jg)%atm%lp_count(jb)

        zwsoil(:) = 0._wp

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data(jg)%atm%idx_lst_lp(ic,jb)

          SELECT CASE(ext_data(jg)%atm%soiltyp(jc,jb))
            CASE (1,2)  !ice,rock
            ! set wsoil to 0 for ice and rock
            zwsoil(jc) = 0._wp

            CASE(3)  !sand
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.364_wp, &
            & MAX((initicon%sfc%wsoil(jc,jb,jk)*(0.196_wp - 0.042_wp) + 0.042_wp),0.012_wp))

            CASE(4)  !sandyloam
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.445_wp, &
            & MAX((initicon%sfc%wsoil(jc,jb,jk)*(0.26_wp  - 0.1_wp  ) + 0.1_wp)  ,0.03_wp ))

            CASE(5)  !loam
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.455_wp, &
            & MAX((initicon%sfc%wsoil(jc,jb,jk)*(0.34_wp  - 0.11_wp ) + 0.11_wp) ,0.035_wp))

            CASE(6)  !clayloam
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.475_wp, &
            & MAX((initicon%sfc%wsoil(jc,jb,jk)*(0.37_wp  - 0.185_wp) + 0.185_wp),0.06_wp ))

            CASE(7)  !clay
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.507_wp, &
            & MAX((initicon%sfc%wsoil(jc,jb,jk)*(0.463_wp - 0.257_wp) + 0.257_wp),0.065_wp))

            CASE(8)  !peat
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.863_wp, &
            & MAX((initicon%sfc%wsoil(jc,jb,jk)*(0.763_wp - 0.265_wp) + 0.265_wp),0.098_wp))

            CASE(9,10)!sea water, sea ice
            ! ERROR landpoint has soiltype sea water or sea ice
            lerr = .TRUE.
          END SELECT

          ! We need to catch problematic coast cases: ICON-land but IFS-ocean &
          ! for moisture (and in principle for temperature as well)
          ! 
          ! The following criterion is probably too hard. ls_mask(jc,jb) < 0.5_wp does 
          ! not necessarily mean that the stencil used for the interpolation onto 
          ! this target point did not contain any IFS-land point. So we may throw away 
          ! some valid points.
          ! However, this criterion is save in the sense that we will not miss any 
          ! problematic point. For a less severe but nevertheless save criterion, 
          ! information about the type of interpolation (CONSERVATIVE, RBF, NN) would be 
          ! necessary. This information is not readily available within ICON.  
          IF ( (initicon%sfc%ls_mask(jc,jb)  < 0.5_wp  ) .AND.   & ! ICON-land but IFS-ocean
            &  (initicon%sfc%wsoil(jc,jb,jk) < -999._wp) )  THEN   ! check for missing value
            ! set dummy value (50% of pore volume)
            zwsoil(jc) = 0.5_wp * cporv(ext_data(jg)%atm%soiltyp(jc,jb)) * dzsoil_icon(jk)
          ENDIF

        ENDDO  ! ic
        ! overwrite wsoil
        initicon%sfc%wsoil(1:nlen,jb,jk) = zwsoil(1:nlen)

      ENDDO  ! jk

      ! assume no-gradient condition for soil moisture reservoir layer
      DO jc = 1, nlen
        initicon%sfc%wsoil(jc,jb,nlev_soil) = initicon%sfc%wsoil(jc,jb,nlev_soil-1)*          &
                                                dzsoil_icon(nlev_soil)/dzsoil_icon(nlev_soil-1)
      ENDDO

      IF (lerr) THEN
        CALL finish(routine, "Landpoint has invalid soiltype (sea water or sea ice)")
      ENDIF

    ENDDO  ! jb
!$OMP END DO 
!$OMP END PARALLEL


  END SUBROUTINE process_sfcfields
  !-------------
  !>
  !! SUBROUTINE smi_to_sm_mass
  !!
  !! Conversion of soil moisture index into TERRA soil moisture [m]
  !!   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
  !!   safety: min=air dryness point, max=pore volume
  !!
  !! Required input: initicon state
  !! Output is written on fields of land state
  !!
  !! @par Revision History
  !! Initial version by P Ripodas, DWD(2013-05)
  !! - extracted from process_sfcfields
  !! Modification by Daniel Reinert, DWD (2013-10-17)
  !! - updated soil moisture initialization according to process_sfcfields
  !
  SUBROUTINE smi_to_sm_mass(p_patch, wsoil)


    TYPE(t_patch), INTENT(IN)    :: p_patch
    REAL(wp),      INTENT(INOUT) :: wsoil(:,:,:)!soil moisture index in
                                                !soil moisture mass out [mm H2O]
                                                ! (nproma,nlev_soil,nblks)
    ! LOCAL VARIABLES
    !
    REAL(wp) :: zwsoil(nproma)        ! local soil moisture field
    INTEGER  :: i_count, ic           ! counter
    INTEGER  :: jg, jb, jk, jc, nlen
    LOGICAL  :: lerr                  ! error flag

    CHARACTER(LEN=*), PARAMETER       :: routine = 'mo_nwp_sfc_interp:smi_to_sm_mass'
!-------------

    jg = p_patch%id

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen,lerr,ic,i_count,zwsoil)
    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF


      ! (re)-initialize error flag
      lerr = .FALSE.

      ! Conversion of IFS soil moisture index (vertically interpolated) into TERRA soil moisture [m]
      !   soil moisture index = (soil moisture - wilting point) / (field capacity - wilting point)
      !   safety: min=air dryness point, max=pore volume
      DO jk = 1, nlev_soil-1

        ! loop over target (ICON) land points only
        i_count = ext_data(jg)%atm%lp_count(jb)

        zwsoil(:) = 0._wp

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, i_count
          jc = ext_data(jg)%atm%idx_lst_lp(ic,jb)

          SELECT CASE(ext_data(jg)%atm%soiltyp(jc,jb))
            CASE (1,2)  !ice,rock
            ! set wsoil to 0 for ice and rock
            zwsoil(jc) = 0._wp

            CASE(3)  !sand
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.364_wp, &
            & MAX((wsoil(jc,jb,jk)*(0.196_wp - 0.042_wp) + 0.042_wp),0.012_wp))

            CASE(4)  !sandyloam
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.445_wp, &
            & MAX((wsoil(jc,jb,jk)*(0.26_wp  - 0.1_wp  ) + 0.1_wp)  ,0.03_wp ))

            CASE(5)  !loam
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.455_wp, &
            & MAX((wsoil(jc,jb,jk)*(0.34_wp  - 0.11_wp ) + 0.11_wp) ,0.035_wp))

            CASE(6)  !clayloam
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.475_wp, &
            & MAX((wsoil(jc,jb,jk)*(0.37_wp  - 0.185_wp) + 0.185_wp),0.06_wp ))

            CASE(7)  !clay
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.507_wp, &
            & MAX((wsoil(jc,jb,jk)*(0.463_wp - 0.257_wp) + 0.257_wp),0.065_wp))

            CASE(8)  !peat
            zwsoil(jc) = dzsoil_icon(jk) * MIN(0.863_wp, &
            & MAX((wsoil(jc,jb,jk)*(0.763_wp - 0.265_wp) + 0.265_wp),0.098_wp))

            CASE(9,10)!sea water, sea ice
            ! ERROR landpoint has soiltype sea water or sea ice
            lerr = .TRUE.
          END SELECT

          ! We need to catch problematic coast cases: ICON-land but IFS/GME-ocean &
          ! for moisture (and in principle for temperature as well)
          ! 
          ! can we do better than this??
          IF ( wsoil(jc,jb,jk) <= -999._wp )  THEN   ! check for missing value
            ! set dummy value (50% of pore volume)
            zwsoil(jc) = 0.5_wp * cporv(ext_data(jg)%atm%soiltyp(jc,jb)) * dzsoil_icon(jk)
          ENDIF

        ENDDO  ! ic
        ! overwrite wsoil
        wsoil(1:nlen,jb,jk) = zwsoil(1:nlen)
      ENDDO  ! jk


      ! assume no-gradient condition for soil moisture reservoir layer
      DO jc = 1, nlen
        wsoil(jc,nlev_soil,jb) = wsoil(jc,nlev_soil-1,jb)*          &
                                        dzsoil_icon(nlev_soil)/dzsoil_icon(nlev_soil-1)
      ENDDO

      IF (lerr) THEN
        CALL finish(routine, "Landpoint has invalid soiltype (sea water or sea ice)")
      ENDIF

    ENDDO  ! jb
!$OMP END DO 
!$OMP END PARALLEL

  END SUBROUTINE smi_to_sm_mass


END MODULE mo_nwp_sfc_interp
