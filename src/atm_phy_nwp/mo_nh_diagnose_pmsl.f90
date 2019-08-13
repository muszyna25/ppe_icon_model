!! Contains utilities for diagnose mean sea level pressure.
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

MODULE mo_nh_diagnose_pmsl

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, finish
  USE mo_impl_constants,      ONLY: PRES_MSL_METHOD_IFS
  USE mo_physical_constants,  ONLY: rd, grav, dtdz_standardatm
  USE mo_parallel_config,     ONLY: nproma
  USE mo_initicon_config,     ONLY: zpbl1, zpbl2

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_diagnose_pmsl'

  ! Artificial limits on temperature profile used for extrapolation below the ground
  REAL(wp), PARAMETER :: t_low  = 255.0_wp
  REAL(wp), PARAMETER :: t_high = 290.5_wp

  ! (note: for the 'diagnose_pmsl...' routines no deep-atmosphere copy exists, 
  ! because the error from grav=const is assumed to be negligibly small for them)
  PUBLIC :: diagnose_pmsl
  PUBLIC :: diagnose_pmsl_gme
  PUBLIC :: diagnose_pmsl_ifs


CONTAINS

  !------------------------------------------------------------------
  !> SUBROUTINE diagnose_pmsl_gme
  !  Compute mean sea level pressure from the surface pressure.
  !
  ! @par Revision History
  ! Initial version by F. Prill, DWD (2012-07-16)
  !
  ! Based on GME's "pp_extrap.f90" which is a modified GM routine by
  ! J.Haseler (ECMWF).
  !
  ! Cf. also
  ! Trenberth, K.E., Berry, J.C., Buja, L.E. 1993:
  ! "Vertical Interpolation and Truncation of Model-coordinate Data."
  !
  ! We assume an atmosphere in hydrostatic equilibrium and having
  ! a constant lapse rate.
  !
  ! Method
  !
  ! For temperature, geopotential and mean sea-level pressure, the
  ! surface temperature, "ztstar", and the msl temperature, "ztmsl"
  ! are calculated, assuming a lapse rate of 6.5 degrees *C/km.  If
  ! both "ztstar" and "ztmsl" are greater than 290.5, "ztmsl" is
  ! corrected according to the formula:
  !
  !    ztmsl = 290.5-0.005*(ztstar-290.5)**2
  !
  ! Comment on this method (from Trenberth et al.):
  ! "The above method seems to work reasonably well as long as the
  ! surface height is less than about 2000m elevation. Problems remain
  ! especially in the Himalayan-Tibetan Plateau region. At very high
  ! elevations, the critical dependence on the local surface
  ! temperature that is used to set the temperature of the whole
  ! artificial below-ground column becomes amplified and leads to
  ! considerable noise in the mean sea level pressure."
  !
  SUBROUTINE diagnose_pmsl_gme(pres3d_in, pres_sfc_in, temp3d_in, z3d_in, pres_out, &
    &                          nblks, npromz, nlevs_in)

    ! Input fields
    REAL(wp), INTENT(IN)  :: pres3d_in   (:,:,:) ! pressure field of input data
    REAL(wp), INTENT(IN)  :: pres_sfc_in (:,:)   ! surface pressure field (input data)
    REAL(wp), INTENT(IN)  :: temp3d_in   (:,:,:) ! temperature of input data
    REAL(wp), INTENT(IN)  :: z3d_in      (:,:,:) ! 3D height coordinate field of input data

    ! Output
    REAL(wp), INTENT(OUT) :: pres_out   (:,:)    ! mean sea level pressure field

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! local variables
!!$    CHARACTER(*), PARAMETER :: routine = TRIM("mo_nh_vert_interp:diagnose_pmsl_gme")
    REAL(wp), PARAMETER :: zthird = 1._wp/3._wp, &
      &                    zlapse = 6.5E-3_wp
    INTEGER             :: nlev, nlevp1, nlen, jb, jc
    REAL(wp)            :: zrg, ztstar, ztmsl, zalph, zprt, zprtal, &
      &                    geop_sfc, temp_in, pres_in, pres_sfc

    !  INITIALISATION
    !  --------------

    nlev   = nlevs_in
    nlevp1 = nlev + 1
    zrg    = 1.0_wp/grav

    ! TEMPERATURE-RELATED FIELDS
    ! --------------------------

    ! initialize output field
    pres_out(:,nblks) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,geop_sfc,temp_in,pres_in,pres_sfc,ztstar,ztmsl,zalph,zprt,zprtal) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
      ENDIF

      DO jc = 1, nlen

        geop_sfc = z3d_in(jc,nlevp1,jb)*grav    ! surface geopotential
        temp_in  = temp3d_in(jc,nlev,jb)        ! temperature at lowest (full) model level
        pres_in  = pres3d_in(jc,nlev,jb)        ! pressure at lowest (full) level
        pres_sfc = pres_sfc_in(jc,jb)           ! surface pressure

        ! Surface pressure is given on half levels, temperature is given
        ! on full levels. Therefore extrapolate temperature with
        ! zlapse=6.5E-3 to model surface.
        !
        ! This formula follows from a Taylor series expansion of
        !   temp/ztstar  =  (pres_sfc/pres_in)**(R*lapse/grav)
        ! for small increments of p:
        ztstar = ( 1._wp + zlapse*rd*zrg*(pres_sfc/pres_in-1._wp) )*temp_in

        ! slight modifications for exceptionally warm or cold
        ! conditions:
        IF ( ztstar < t_low) THEN
          ztstar = 0.5_wp*(t_low+ztstar)
        END IF

        ztmsl = ztstar + zlapse*zrg*geop_sfc

        IF ( ztmsl > t_high) THEN
          IF ( ztstar > t_high) THEN
            ztstar = 0.5_wp*(t_high+ztstar)
            ztmsl  = ztstar
          ELSE
            ztmsl  = t_high
          END IF
        END IF

        IF ( ABS(ztmsl-ztstar) < 1.E-6_wp) THEN
          zalph = 0._wp
        ELSE IF ( ABS(zrg*geop_sfc) > 1.E-4_wp) THEN
          zalph= rd*(ztmsl-ztstar)/geop_sfc
        ELSE
          zalph= rd*zlapse*zrg
        END IF

        ! MEAN SEA LEVEL PRESSURE
        ! -----------------------

        ! The formula
        !   pres_msl = pres_sfc ( 1 + (zalph*geop_sfc)/(rd*ztstar) )**(1/zalph)
        ! is approximated (using the log-power series) as follows:
        zprt   = geop_sfc/(rd*ztstar)
        zprtal = zprt*zalph
        pres_out(jc,jb) = pres_sfc*EXP( zprt * ( 1.0_wp-zprtal*(0.5_wp-zthird*zprtal) ) )

      END DO ! jc
    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE diagnose_pmsl_gme



  !------------------------------------------------------------------
  !> SUBROUTINE diagnose_pmsl_ifs
  !  Compute mean sea level pressure from the surface pressure.
  !
  ! @par Revision History
  ! Initial version by Guenther Zaengl, DWD (2013-04-11)
  !
  ! Based on IFS's "ppmer.f90"
  !
  ! The method is basically the same as for the GME routine above, with the
  ! exception that the surface temperature used for extrapolating the pressure
  ! is based on the temperature at 150 m AGL rather than that at the lowest
  ! model level. Moreover, the artificial corrections applied to extremely cold
  ! temperatures are slightly different (remark: they are not self-consistent
  ! for geopotential and sea-level pressure, and are again different to the
  ! temperatures used for temperature extrapolation)
  !
  !
  SUBROUTINE diagnose_pmsl_ifs(pres_sfc_in, temp3d_in, z3d_in, pmsl_out,      &
    &                          nblks, npromz, nlevs_in, wfac_extrap, kextrap, &
    &                          zextrap, method )

    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_sfc_in (:,:)   ! surface pressure field (input data)
    REAL(wp), INTENT(IN)  :: temp3d_in   (:,:,:) ! temperature of input data
    REAL(wp), INTENT(IN)  :: z3d_in      (:,:,:) ! 3D height coordinate field of input data

    ! Output
    REAL(wp), INTENT(OUT) :: pmsl_out   (:,:)    ! mean sea level pressure field

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Coefficients
    INTEGER , INTENT(IN) :: kextrap(:,:)     ! index of model level immediately above (by default) 150 m AGL
    REAL(wp), INTENT(IN) :: wfac_extrap(:,:) ! corresponding interpolation coefficient
    REAL(wp), INTENT(IN) :: zextrap(:,:)     ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! Method (PRES_MSL_METHOD_IFS or PRES_MSL_METHOD_IFS_CORR)
    INTEGER , INTENT(IN) :: method           ! valid input: PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":diagnose_pmsl_ifs"

    INTEGER             :: nlev, nlevp1, nlen, jb, jc

    REAL(wp)            :: temp_extrap, vtgrad, dtdz_thresh
    REAL(wp), DIMENSION(nproma) :: tsfc, tmsl, tsfc_mod, tmsl_mod

!-------------------------------------------------------------------------

    nlev   = nlevs_in
    nlevp1 = nlev + 1

    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    dtdz_thresh = 1.e-4_wp ! 0.1 K/km


    ! sanity check
    !
    IF (method < 3) THEN
      WRITE(message_text,'(a)') 'invalid value for itype_pres_msl! must be >= 3'
      CALL finish(routine, message_text)
    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jc,temp_extrap,vtgrad,tsfc,tmsl,tsfc_mod,tmsl_mod) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pmsl_out(nlen+1:nproma,jb)  = 0.0_wp
      ENDIF

      DO jc = 1, nlen

        ! Auxiliary temperature at 150 m AGL from which the sfc temp is extrapolated
        temp_extrap = wfac_extrap(jc,jb) *temp3d_in(jc,kextrap(jc,jb),jb  ) + &
               (1._wp-wfac_extrap(jc,jb))*temp3d_in(jc,kextrap(jc,jb)+1,jb)

        ! extrapolated surface temperature (tstar in IFS nomenclature)
        tsfc(jc) = temp_extrap - dtdz_standardatm*zextrap(jc,jb)

        ! extrapolated sea-level temperature (t0 in IFS nomenclature)
        tmsl(jc) = tsfc(jc) - dtdz_standardatm*z3d_in(jc,nlevp1,jb)

        ! artificial modifications in the presence of extreme temperatures
        IF (tsfc(jc) < t_low) THEN
          tsfc_mod(jc) = 0.5_wp*(tsfc(jc)+t_low)
        ELSE IF (tsfc(jc) < t_high) THEN
          tsfc_mod(jc) = tsfc(jc)
        ELSE
          tsfc_mod(jc) = 0.5_wp*(tsfc(jc)+t_high)
        ENDIF

        tmsl_mod(jc) = tsfc_mod(jc) - dtdz_standardatm*z3d_in(jc,nlevp1,jb)

        IF (tsfc_mod(jc) < t_high .AND. tmsl_mod(jc) > t_high) THEN
          tmsl_mod(jc) = t_high
        ELSE IF (tsfc_mod(jc) >= t_high .AND. tmsl_mod(jc) > tsfc_mod(jc)) THEN
          tmsl_mod(jc) = tsfc_mod(jc)
        ELSE
          IF (method == PRES_MSL_METHOD_IFS) THEN
            tmsl_mod(jc) = tmsl(jc) ! this is missing in the geopotential computation!!!
          ENDIF
          ! tmsl_mod unmodified for method == PRES_MSL_METHOD_IFS_CORR
        ENDIF

        vtgrad = (tsfc_mod(jc) - tmsl_mod(jc))/SIGN(MAX(1.e-4_wp,ABS(z3d_in(jc,nlevp1,jb))),z3d_in(jc,nlevp1,jb))

        IF (ABS(vtgrad) > dtdz_thresh) THEN
          pmsl_out(jc,jb) = pres_sfc_in(jc,jb)*EXP(-grav/(rd*vtgrad)* &
                            LOG(tmsl_mod(jc)/tsfc_mod(jc)) )
        ELSE
          pmsl_out(jc,jb) = pres_sfc_in(jc,jb)*EXP(grav*z3d_in(jc,nlevp1,jb)/ &
                           (rd*0.5_wp*(tsfc_mod(jc)+tmsl_mod(jc))) )
        ENDIF


      END DO ! jc
    END DO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE diagnose_pmsl_ifs


  !-------------
  !>
  !! SUBROUTINE diagnose_pmsl
  !! Diagnoses sea-level pressure
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate field
  !! Output: sea-level pressure
  !!
  !! Method: stepwise analytical integration of the hydrostatic equation
  !!  step 1: from surface level of input data 500 m (zpbl1) upward
  !!  step 2: from there to zpbl1 meters above sea level
  !!  step 3: from there down to sea level
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2012-07-17)
  !!
  !!
  !!
  SUBROUTINE diagnose_pmsl(pres_in, tempv_in, z3d_in, pmsl_out, nblks, npromz, nlevs_in, &
                           wfacpbl1, kpbl1, wfacpbl2, kpbl2, z_target, extrapol_dist)


    ! Input fields
    REAL(wp), INTENT(IN)  :: pres_in  (:,:,:) ! pressure field on main levels
    REAL(wp), INTENT(IN)  :: tempv_in (:,:,:) ! virtual temperature on main levels
    REAL(wp), INTENT(IN)  :: z3d_in   (:,:,:) ! 3D height coordinate field of main levels

    ! Output
    REAL(wp), INTENT(OUT) :: pmsl_out (:,:) ! sea-level pressure

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Coefficients
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    REAL(wp), INTENT(IN) :: z_target       ! target height of pressure diagnosis (usually 0 m = sea level)
    REAL(wp), INTENT(IN) :: extrapol_dist  ! maximum extrapolation distance up to which the local vertical
                                           ! temperature gradient is used

    ! LOCAL VARIABLES

    INTEGER  :: jb, jc
    INTEGER  :: nlen
    REAL(wp), DIMENSION(nproma) :: tempv1, tempv2, vtgrad_pbl_in, vtgrad_up, zdiff_inout, &
                                   sfc_inv, tempv_out_pbl1, vtgrad_up_out, tempv_out,     &
                                   vtgrad_pbl_out, p_pbl1, p_pbl1_out
    REAL(wp) :: dtvdz_thresh

!-------------------------------------------------------------------------

    ! Threshold for switching between analytical formulas for constant temperature and
    ! constant vertical gradient of temperature, respectively
    dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jc,nlen,tempv1,tempv2,vtgrad_pbl_in,vtgrad_up,zdiff_inout,sfc_inv, &
!$OMP            tempv_out_pbl1,vtgrad_up_out,tempv_out,vtgrad_pbl_out,p_pbl1,         &
!$OMP            p_pbl1_out) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pmsl_out(nlen+1:nproma,jb)  = 0.0_wp
      ENDIF

      DO jc = 1, nlen
        ! Virtual temperature at height zpbl1
        tempv1(jc) = wfacpbl1(jc,jb) *tempv_in(jc,kpbl1(jc,jb),jb  ) + &
              (1._wp-wfacpbl1(jc,jb))*tempv_in(jc,kpbl1(jc,jb)+1,jb)

        ! Virtual temperature at height zpbl2
        tempv2(jc) = wfacpbl2(jc,jb) *tempv_in(jc,kpbl2(jc,jb),jb  ) + &
              (1._wp-wfacpbl2(jc,jb))*tempv_in(jc,kpbl2(jc,jb)+1,jb)

        ! Vertical gradient between surface level zpbl1
        vtgrad_pbl_in(jc) = (tempv1(jc) - tempv_in(jc,nlevs_in,jb))/zpbl1

        ! Vertical gradient between zpbl1 and zpbl2
        vtgrad_up(jc) = (tempv2(jc) - tempv1(jc))/(zpbl2 - zpbl1)

        ! Set reasonably restrictive limits
        vtgrad_up(jc) = MAX(vtgrad_up(jc),-8.0e-3_wp)
        vtgrad_up(jc) = MIN(vtgrad_up(jc),-5.0e-3_wp)

        ! height distance between lowest input level
        ! (negative if model grid point lies above z_target)
        zdiff_inout(jc) = z_target - z3d_in(jc,nlevs_in,jb)

        ! "surface inversion", defined by the difference between the extrapolated
        ! extrapolated temperature from above and the original input temperature
        sfc_inv(jc) = tempv1(jc)-vtgrad_up(jc)*zpbl1 - tempv_in(jc,nlevs_in,jb)

        ! Reduction of the surface inversion depending on the extrapolation
        ! distance. The surface inversion is fully restored for extrapolation distances
        ! up to zpbl1 and disregarded for distances larger than 3*zpbl1
        IF (ABS(zdiff_inout(jc)) > 3._wp*zpbl1) THEN
          sfc_inv(jc) = 0._wp
        ELSE IF (ABS(zdiff_inout(jc)) > zpbl1) THEN
          sfc_inv(jc) = sfc_inv(jc)*(1._wp - (ABS(zdiff_inout(jc))-zpbl1)/(2._wp*zpbl1))
        ENDIF

        ! Linear extrapolation to zpbl1 above target height using the upper temperature gradient
        ! extrapol_dist is the maximum extrapolation distance up to which the local
        ! vertical gradient is considered; for larger distance, the standard atmosphere
        ! vertical gradient is taken
        IF (zdiff_inout(jc) > extrapol_dist) THEN
          tempv_out_pbl1(jc) = tempv1(jc) + zdiff_inout(jc)*vtgrad_up(jc)
          vtgrad_up_out(jc)  = vtgrad_up(jc)
        ELSE
          tempv_out_pbl1(jc) = tempv1(jc) + extrapol_dist*vtgrad_up(jc) +   &
                            (zdiff_inout(jc)-extrapol_dist)*dtdz_standardatm

          ! Artificial limitation analogous to GME method
          IF (tempv_out_pbl1(jc) < 255.65_wp) THEN
            tempv_out_pbl1(jc) = 0.5_wp*(tempv_out_pbl1(jc)+255.65_wp)
          ELSE IF (tempv_out_pbl1(jc) > 290.65_wp) THEN
            tempv_out_pbl1(jc) = 0.5_wp*(tempv_out_pbl1(jc)+290.65_wp)
            tempv_out_pbl1(jc) = MIN(298.15_wp,tempv_out_pbl1(jc))
            tempv_out_pbl1(jc) = MAX(248.15_wp,tempv_out_pbl1(jc))
          ENDIF

          ! Averaged vertical temperature gradient
          vtgrad_up_out(jc) = (tempv_out_pbl1(jc) - tempv1(jc))/zdiff_inout(jc)
        ENDIF

        ! Final extrapolation of temperature to target height, including restored surface inversion
        tempv_out(jc) = tempv_out_pbl1(jc) - vtgrad_up_out(jc)*zpbl1 - sfc_inv(jc)

        ! Boundary-layer vertical gradient of extrapolation profile
        vtgrad_pbl_out(jc) = (tempv_out_pbl1(jc)-tempv_out(jc))/zpbl1


        ! Three-step vertical integration of pressure
        ! step 1: from surface level of input data to zpbl1 (to get rid of surface inversion)

        IF (ABS(vtgrad_pbl_in(jc)) > dtvdz_thresh) THEN
          p_pbl1(jc) = pres_in(jc,nlevs_in,jb)*EXP(-grav/(rd*vtgrad_pbl_in(jc))* &
                       LOG(tempv1(jc)/tempv_in(jc,nlevs_in,jb)) )
        ELSE
          p_pbl1(jc) = pres_in(jc,nlevs_in,jb)*EXP(-grav*zpbl1 / &
                      (rd*0.5_wp*(tempv1(jc)+tempv_in(jc,nlevs_in,jb))) )
        ENDIF

        ! step 2: from zpbl1 of input data to zpbl1 above target level
        IF (ABS(vtgrad_up_out(jc)) > dtvdz_thresh) THEN
          p_pbl1_out(jc) = p_pbl1(jc)*EXP(-grav/(rd*vtgrad_up_out(jc))* &
                           LOG(tempv_out_pbl1(jc)/tempv1(jc)) )
        ELSE
          p_pbl1_out(jc) = p_pbl1(jc)*EXP(-grav*zdiff_inout(jc)/ &
                          (rd*0.5_wp*(tempv1(jc)+tempv_out_pbl1(jc))) )
        ENDIF

        ! step 3: from zpbl1 above target level to target level
        IF (ABS(vtgrad_pbl_out(jc)) > dtvdz_thresh) THEN
          pmsl_out(jc,jb) = p_pbl1_out(jc)*EXP(-grav/(rd*vtgrad_pbl_out(jc))* &
                            LOG(tempv_out(jc)/tempv_out_pbl1(jc)) )
        ELSE
          pmsl_out(jc,jb) = p_pbl1_out(jc)*EXP(grav*zpbl1/ &
                           (rd*0.5_wp*(tempv_out(jc)+tempv_out_pbl1(jc))) )
        ENDIF

      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE diagnose_pmsl


END MODULE  mo_nh_diagnose_pmsl

