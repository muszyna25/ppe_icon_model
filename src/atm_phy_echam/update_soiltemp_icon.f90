SUBROUTINE update_soiltemp_icon(nidx,nsoil                   &
                         , delta_time                        &
                         , time_steps_soil                   &
                         , pts                               &
!!$ TR                         , psn                               &
                         , psodif, prgcgn                    &
                         , pgrndc, pgrndd                    &
                         , ptsoil                            &
                         , pgrndcapc, pgrndhflx              &
!!$ TR                         , lmask, ldglac)
                         )
!
!   AUTHOR:  FREDERIC HOURDIN     30/01/92
!
!            ADAPTED TO THE LMD-GCM BY JAN POLCHER  26/02/92
!            ADAPTED TO THE ECHAM-GCM BY JAN-PETER SCHULZ, MPI  03/02/96
!
!            J.-P. SCHULZ   MPI - OCTOBER 1997 :
!               ROUTINE USED FOR IMPLEMENTATION OF AN IMPLICIT
!               COUPLING BETWEEN LAND SURFACE AND ATMOSPHERE IN THE
!               ECHAM4 GCM.
!            U.SCHLESE DKRZ - NOVEMBER 1999  MODIFIED FOR ECHAM5
!            U.Schlese DKRZ - February 2000  new soil temperatures
!            L Kornblueh, MPI, January 2003, removed MERGE
!            
!            Adapted to JSBACH by Thomas Raddatz, Mai 2004
!            Adapted to ICON  by Thomas Raddatz, Sep 2011
!   OBJECTIVE:  COMPUTATION OF:
!               THE GROUND TEMPERATURE EVOLUTION
!               THE GROUND SPECIFIC HEAT "CAPCAL"
!               THE SURFACE DIFFUSIVE FLUX FROM GROUND "F0"
!
!
!   METHOD:  IMPLICIT TIME INTEGRATION
!
!   CONSECUTIVES GROUND TEMPERATURES ARE RELATED BY:
!           T(K+1) = C(K) + D(K)*T(K)  (1)
!   THE COEFFICIENTS C (=pgrndc) AND D (=pgrndd) ARE COMPUTED AT THE
!   T-DT TIME-STEP.
!   ROUTINE STRUCTURE:
!   1)NEW TEMPERATURES ARE COMPUTED  USING (1)
!   2)C AND D COEFFICIENTS ARE COMPUTED FROM THE NEW TEMPERATURE
!     PROFILE FOR THE T+DT TIME-STEP
!   3)THE COEFFICIENTS A AND B ARE COMPUTED WHERE THE DIFFUSIVE
!     FLUXES AT THE T+DT TIME-STEP IS GIVEN BY
!            FDIFF = A + B TS(T+DT)
!     OR     FDIFF = F0 + CAPCAL (TS(T+DT)-TS(T))/DT
!            WITH F0 = A + B (TS(T))
!                 CAPCAL = B*DT
!
!     ------------------------------------------------------------------
!
!   DECLARATIONS:
!
!!$ TR USE mo_jsbach_constants   , ONLY: RhoH2O
!!$ TR USE mo_time_control       , ONLY: lstart
USE mo_kind               , ONLY: wp

IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  ARGUMENTS
!
  INTEGER, Intent(in)  ::  nidx                 !! length of the vector
  INTEGER, Intent(in)  ::  nsoil                !! number of soil layers (fixed to 5)
  REAL(wp), Intent(in)     ::  delta_time       !! time step
  REAL(wp), Intent(in)     ::  time_steps_soil(nidx)   !! number of time steps since initialisation of the soil
  REAL(wp), Intent(in)     ::  pts(nidx)            !! surface temperature at top of soil [K]
!!$ TR  REAL(wp), Intent(in)     ::  psn(nidx)            !! equivalent snow depth [m water]
  REAL(wp), Intent(in)     ::  psodif(nidx)         !! soil temperature diffusivity [m^2/s]
  REAL(wp), Intent(in)     ::  prgcgn(nidx)         !! soil heat capacity [J/m^3K]
  REAL(wp), Intent(inout)  ::  pgrndc(nidx,nsoil)   !!
  REAL(wp), Intent(inout)  ::  pgrndd(nidx,nsoil)   !!
  REAL(wp), Intent(inout)  ::  ptsoil(nidx,nsoil)   !! soil temperature [K]
  REAL(wp), Intent(out)    ::  pgrndcapc(nidx)      !!
  REAL(wp), Intent(out)    ::  pgrndhflx(nidx)      !! ground heat flux
!!$ TR  LOGICAL, Intent(in)  ::  lmask(nidx)
!!$ TR  LOGICAL, Intent(in)  ::  ldglac(nidx)         !! glacier mask
!
!     ------------------------------------------------------------------
!
!  local Variables
!
  INTEGER :: jk
  REAL(wp) :: zso_cond(nidx), zso_capa(nidx)
  REAL(wp) :: z1(nidx)
  REAL(wp) :: zd1(nsoil)
  REAL(wp) :: zdz1(nidx,nsoil),   zdz2(nidx,nsoil)
  REAL(wp) :: zkappa(nidx,nsoil), zcapa(nidx,nsoil)
  REAL(wp) :: zsnow_h(nidx), zx1(nidx), zx2(nidx)
  REAL(wp) :: cdel(nsoil),cmid(nsoil)
  REAL(wp) :: zrici, zdifiz, zsn_cond, zsn_dens, zsn_capa
!
!     ------------------------------------------------------------------
!
!*    1.  SPECIFYING THE DEPTHS OF THE TEMPERATURE LEVELS.
!
!*    1.1 SOME CONSTANTS.
!
  zrici = 2.09e+06_wp                                 !! volumetric heat capacity of ice [j/m**3/k]
  zdifiz = 12.e-07_wp                                 !! temperature diffusivity of ice  [m**2/s]
  zsn_cond = 0.31_wp                                  !! snow thermal conductivity [j/s/m/k]
  zsn_dens = 330.0_wp                                 !! snow density              [kg/m**3]
  zsn_capa = 634500.0_wp                              !! snow  heat capacity   [j/m**3/k]
  cdel = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)         !! thicknesses of soil layers [m]
  cmid = (/0.0325_wp,0.192_wp,0.7755_wp,2.683_wp,6.984_wp/)       !! depth of mids of soil layers [m]
!
!*    1.2 COMPUTING SOME USEFUL CONSTANTS.
!
  DO jk = 1,nsoil-1
     zd1(jk) = 1._wp / (cmid(jk+1) - cmid(jk))
  END DO
!
!*    1.3 COMPUTE OF THE SOIL THERMAL CONDUCTIVITY [J/S/M/K] FROM
!*        THE SOIL TEMPERATURE DIFFUSIVITY [M**2/S].
!
!!$ TR  WHERE (lmask(:) .AND. ldglac(:))
!!$ TR    zso_capa(:) = zrici
!!$ TR    zso_cond(:) = zso_capa(:) * zdifiz
!!$ TR  ELSEWHERE (lmask(:))
  zso_capa(:) = prgcgn(:)
  zso_cond(:) = zso_capa(:) * psodif(:)
!
!*    1.4 PRE-SET THERMAL CONDUCTIVITY AT ALL LEVELS.
!
!!$ TR  DO jk = 1,nsoil
!!$ TR     WHERE (lmask)
!!$ TR        zkappa(:,jk) = zso_cond(:)
!!$ TR        zcapa(:,jk)  = zso_capa(:)
!!$ TR     END WHERE
!!$ TR  END DO

  DO jk = 1,nsoil
     zkappa(:,jk) = zso_cond(:)
     zcapa(:,jk)  = zso_capa(:)
  END DO

!
!   --------------------------------------------------------------
!   COMPUTATION OF THE GROUND TEMPERATURES USING THE CGRD AND DGRD
!   COEFFICIENTS COMPUTED AT THE PREVIOUS TIME-STEP
!   --------------------------------------------------------------
!
!!$ TR  IF(.NOT.lstart) THEN
!   Upper layer
!
!!$ TR    ptsoil(:,1) = pts(:)
    WHERE (time_steps_soil(:) > 0.5_wp) ptsoil(:,1) = pts(:)
!
!   Deeper layers
!
    DO jk = 1,nsoil-1
!!$ TR      WHERE(lmask) ptsoil(:,jk+1) = pgrndc(:,jk) + pgrndd(:,jk) * ptsoil(:,jk)
      WHERE (time_steps_soil(:) > 0.5_wp) ptsoil(:,jk+1) = pgrndc(:,jk) + &
        pgrndd(:,jk) * ptsoil(:,jk)
    END DO
!!$ TR  END IF
!
!   ---------------------------------------------------------------
!   COMPUTATION OF THE CGRD AND DGRD COEFFICIENTS FOR THE NEXT STEP
!   ---------------------------------------------------------------
!
!!$ TR  WHERE (lmask) 
!!$ TR    zsnow_h(:) = psn(:) * RhoH2O / zsn_dens
     zsnow_h(:) = 0._wp
     
!
!*       Special treatment for first layer
!
     WHERE ( zsnow_h(:) > cmid(2) )
        zcapa(:,1) = zsn_capa
        zkappa(:,1) = zsn_cond
     ELSEWHERE( zsnow_h(:) > 0.0_wp .AND. zsnow_h(:) <= cmid(2) )
        zx1 = zsnow_h(:) / cmid(2)
        zx2 = ( cmid(2) - zsnow_h(:)) / cmid(2)
        zcapa(:,1) = zx1 * zsn_capa + zx2 * zso_capa(:)
        zkappa(:,1) = 1._wp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
     ELSEWHERE
        zcapa(:,1) = zso_capa(:)
        zkappa(:,1) = zso_cond(:)
     ENDWHERE
!!$ TR  END WHERE
!
  DO jk = 2, nsoil - 2
!!$ TR    WHERE (lmask)
       WHERE ( zsnow_h(:) > cmid(jk+1) )
          zcapa(:,jk) = zsn_capa
          zkappa(:,jk) = zsn_cond
       ELSEWHERE ( zsnow_h(:) > cmid(jk) .AND. zsnow_h(:) <= cmid(jk+1) )
          zx1 = (zsnow_h(:) - cmid(jk)) * zd1(jk)
          zx2 = ( cmid(jk+1) - zsnow_h(:)) * zd1(jk)
          zcapa(:,jk) = zx1 * zsn_capa + zx2 * zso_capa(:)
          zkappa(:,jk) = 1._wp / ( zx1 / zsn_cond + zx2 / zso_cond(:) )
       ELSEWHERE
          zcapa(:,jk) = zso_capa(:)
          zkappa(:,jk) = zso_cond(:)
       ENDWHERE
!!$ TR    END WHERE
  END DO
!
  DO jk=1,nsoil
!!$ TR    WHERE (lmask) zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
     zdz2(:,jk) = zcapa(:,jk) * cdel(jk) / delta_time
  END DO
!
  DO jk=1,nsoil-1
!!$ TR    WHERE (lmask) zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
     zdz1(:,jk) = zd1(jk) * zkappa(:,jk)
  END DO
!
!!$ TR  WHERE (lmask)
     z1(:) = zdz2(:,nsoil) + zdz1(:,nsoil-1)
     pgrndc(:,nsoil-1) = zdz2(:,nsoil) * ptsoil(:,nsoil) / z1(:)
     pgrndd(:,nsoil-1) = zdz1(:,nsoil-1) / z1(:)
!!$ TR  END WHERE
!
  DO jk=nsoil-1,2,-1
!!$ TR     WHERE (lmask)
        z1(:) = 1._wp / (zdz2(:,jk) + zdz1(:,jk-1) + zdz1(:,jk) * (1._wp - pgrndd(:,jk)))
        pgrndc(:,jk-1) = (ptsoil(:,jk) * zdz2(:,jk) + zdz1(:,jk) * pgrndc(:,jk)) * z1(:)
        pgrndd(:,jk-1) = zdz1(:,jk-1) * z1(:)
!!$ TR     END WHERE
  END DO
!
!   ---------------------------------------------------------
!   COMPUTATION OF THE SURFACE DIFFUSIVE FLUX FROM GROUND AND
!   CALORIFIC CAPACITY OF THE GROUND:
!   ---------------------------------------------------------
!
!!$ TR  WHERE (lmask)
     pgrndhflx(:) = zdz1(:,1) * (pgrndc(:,1) + (pgrndd(:,1) - 1._wp) * ptsoil(:,1))
     pgrndcapc(:) = (zdz2(:,1) * delta_time + delta_time * (1._wp - pgrndd(:,1)) * zdz1(:,1))
!!$ TR  ELSEWHERE
!!$ TR     pgrndhflx = 0._wp
!!$ TR     pgrndcapc = 0._wp
!!$ TR  END WHERE

END SUBROUTINE update_soiltemp_icon
