SUBROUTINE update_surf_down ( kbdim, delta_time, time_steps_soil  &
     , pqm1, ptsl, pspeed10, ptm1  &
     , pwl, pcvw, pwlmx &
     , psn, pcvs, psnc,       pgld             &
     , pgrndcapc &
!!$     , paphm1, paphms1, ptte             &     ! ATTENTION!
     , pevapl,     pevapot, pevwsd                    &
     , prain, psnow                            &
     , lmask, lpglac       &
     , palac &
     , psnacl, psnmel, progl &
     , pmlres, tte_corr &
!++mgs : add instantaneous snow melt for drydep calculation
     , psmelt   &
!--mgs

     ) 
  !
  !     ------------------------------------------------------------------
  ! ptsl          Sfc temp [K] at timestep t+dt (unfiltered)
  ! pspeed10      Wind speed at 10m height [m/s] ... from 'vdiff'
  ! ptm1          Air temp [K] at timestep t-dt (filtered)
  ! pwl           Water content [m] in skin reservoir
  !               (vegetation and bare soil)
  ! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
  ! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
  !               of vegetation index and leaf area index)
  ! psn           Snow depth [m water equivalent] at the ground
  ! pcvs          Fractional snow cover (function of psn in *physc*)
  ! psnc          Snow depth [m water equivalent] at the canopy
  ! pgld          Glacier depth (including snow) [m water equivalent]
  ! pgrndcapc     Heat capacity of the uppermost soil layer [j/m**2/K]
  ! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
  ! pevapot       Potential evaporation/sublimation [kg/m**2/s]
  ! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
  ! prain         Totalrainfall [kg/m**2/s]
  ! psnow         Total snowfall [kg/m**2/s]
  ! lpglac        Logical glacier mask
  ! palac         Precipitation minus sublimation at glacier points [m water equivalent]
  ! psnacl        Snow budget at non-glacier points [kg/m**2] (accumul.)
  ! psnmel        Snow melt [kg/m**2] (accumulated for diagnostics)
  ! progl         Glacier runoff (rain+snow/ice melt) [kg/m**2] (accumul.)
  !
  ! The following local variables represent the respective fluxes
  ! integrated over one timestep (delta_time) and divided by the density of
  ! water (rhoh2o). Units are m water equivalent.
  !
  ! zraind        Total rain
  ! zsnowd        Total snowfall
  ! zevttd        Total evaporation
  ! zevsnd        Sublimation from snow
  ! zevwld        Evaporation from the skin reservoir
  ! zrogl         Runoff at glacier points (rain and melt, but no calving)
  ! zsnmel        Snow/ice melt at land and glacier points
  ! zsncmelt      Snow melt in the canopy
  ! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
  ! pmlres        Residual melt water available for infiltration into the
  !               non-frozen soil after being intercepted by the
  !               skin reservoir
  ! tte_corr      ???
  ! psmelt        instantaneous snow melt rate (for dry deposition)
  !
  !       Rest folgt spaeter ....
  !
  !     *SURF* - UPDATES LAND VALUES OF TEMPERATURE, MOISTURE AND SNOW.
  !              CALCULATE FLUXES OF TOTAL RAIN, TOTAL SNOW AND EVAPO-
  !              RATION FROM THE THREE RESERVOIRS (SN, WS, WL)
  !              CONVERT FLUXES (KG/M**2*S) TO CHANGES OF WATER LEVELS (M)
  !              DURING TIMESTEP DELTA_TIME.
  !
  !     J.F.GELEYN     E.C.M.W.F.     08/06/82.
  !     MODIFIED BY
  !     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
  !     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
  !     J.-P. SCHULZ   MPI - 1997 : IMPLEMENTATION OF IMPLICIT
  !                                 COUPLING BETWEEN LAND SURFACE
  !                                 AND ATMOSPHERE.
  !     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
  !     MODIFIED BY M. ESCH        MPI - APR  1999
  !     MODIFIED BY E. ROECKNER    MPI - JAN  2001
  !     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
  !     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
  !                                                for snow changed
  !     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
  !
  !     MODIFICATION
  !
  !     PURPOSE
  !
  !     INTERFACE.
  !
  !          *SURF* IS CALLED FROM *PHYSC*.
  !
  !     METHOD.
  !
  !     EXTERNALS.
  !
  !          NONE.
  !
  !     REFERENCE.
  !
  !          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
  !     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.
  !
!!$ TR  USE mo_param_switches ! only lsurf is used
!!$ TR  USE mo_physc2
!!$ TR  USE mo_constants
!!$ TR  USE mo_vegetation ! only cvinter is used
!!$ TR  USE mo_radiation
  USE mo_physical_constants, ONLY: alf, cpd, vtmpc2, rhoh2o, tmelt, g=>grav
  USE mo_kind, ONLY: dp=>wp
  !
  IMPLICIT NONE

  INTEGER,  INTENT(in) :: kbdim
  REAL(dp), INTENT(in) :: delta_time
  REAL(dp), INTENT(in) :: time_steps_soil(kbdim)
  LOGICAL,  INTENT(in) :: lmask(kbdim)
  LOGICAL,  INTENT(in) :: lpglac(kbdim)
  REAL(dp), INTENT(in) :: pqm1(kbdim)
  REAL(dp), INTENT(in) :: pspeed10(kbdim)
  REAL(dp), INTENT(in) :: ptm1(kbdim)
  REAL(dp), INTENT(in) :: pcvw(kbdim)
  REAL(dp), INTENT(in) :: pwlmx(kbdim)
  REAL(dp), INTENT(in) :: pcvs(kbdim)
  REAL(dp), INTENT(in) :: pgrndcapc(kbdim)
  REAL(dp), INTENT(in) :: pevapl(kbdim)
  REAL(dp), INTENT(in) :: pevapot(kbdim)
  REAL(dp), INTENT(in) :: prain(kbdim)
  REAL(dp), INTENT(in) :: psnow(kbdim)

  REAL(dp), INTENT(inout) :: ptsl(kbdim)
  REAL(dp), INTENT(inout) :: pwl(kbdim)
  REAL(dp), INTENT(inout) :: psn(kbdim)
  REAL(dp), INTENT(inout) :: psnc(kbdim)
  REAL(dp), INTENT(inout) :: pgld(kbdim)
  REAL(dp), INTENT(inout) :: psnacl(kbdim)
  REAL(dp), INTENT(inout) :: psnmel(kbdim)
  REAL(dp), INTENT(inout) :: progl(kbdim)

  REAL(dp), INTENT(out)   :: pevwsd(kbdim)
  REAL(dp), INTENT(out)   :: palac(kbdim)
  REAL(dp), INTENT(inout) :: tte_corr(kbdim)
  REAL(dp), INTENT(out)   :: pmlres(kbdim)
  REAL(dp), INTENT(out)   :: psmelt(kbdim)   !++mgs
  !
  !  local arrays
  !
  INTEGER :: jl
  REAL(dp) ::                                                             &
       zraind (kbdim),       zsnowd(kbdim),         zevttd(kbdim)        &
       , zevsnd(kbdim),        zevwld(kbdim)        &
       , zrogl(kbdim)         &
       , zsnmel(kbdim),        zsn(kbdim)                                 &
       ,        zsncmelt(kbdim)                             
!!$       , zdp(kbdim),          
  REAL(dp) :: zlfdcp(kbdim)
  !
  !  local scalars
  !
  REAL(dp) ::                                                             &
       zsmelt, zsnmlt,      &
       zmprcp,    &
       zc2, zc3, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac
  !
  !      Parameters
  !

  REAL(dp), PARAMETER :: cvinter = 0.25_dp  !!$ TR for JSBACH testing
  REAL(dp), PARAMETER :: cwlmax = 2.E-4_dp  !!$ TR for JSBACH testing

  zdtime = delta_time
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zsncfac=rhoh2o*g/zdtime
  !
  !     ------------------------------------------------------------------
  !
  !*    1.     Convert water fluxes to [m water equivalent * timestep]
  !
  DO 110 jl=1,kbdim
     palac(jl)=0._dp
     zrogl(jl)=0._dp
     pmlres(jl)=0._dp
     zsn(jl)=0._dp
     zsnmel(jl)=0._dp
     zsncmelt(jl)=0._dp
     zraind(jl)=prain(jl)             *zdtime/rhoh2o
     zsnowd(jl)=psnow(jl)             *zdtime/rhoh2o
     IF(lmask(jl)) THEN
        zrcp=1._dp/(cpd*(1._dp+vtmpc2*MAX(0.0_dp,pqm1(jl))))
        zlfdcp(jl)=alf*zrcp
        zevttd(jl)=pevapl(jl)                        *zdtime/rhoh2o
        zevsnd(jl)=pcvs(jl)*pevapot(jl)              *zdtime/rhoh2o
        zevwld(jl)=(1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
        pevwsd(jl)=zevttd(jl)-zevsnd(jl)
     ELSE
        pevwsd(jl) = 0._dp ! Sain initialization for diagnostics
     END IF
110 END DO
  !
!!$ TR  IF(lsurf) THEN
     !
     !     ------------------------------------------------------------------
     !
     !*    2.     Budgets of snow (canopy and ground) and glaciers
     !
     !*    2.1    Snow changes in the canopy (interception of snowfall,
     !            sublimation, melting, unloading due to wind)
     !
     DO 210 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           zsn(jl)=zsnowd(jl)+zevsnd(jl)
           zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
           zmprcp=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
           zsncp=psnc(jl)+zmprcp
           zsnowd(jl)=zsnowd(jl)-zmprcp
           psnc(jl)=MIN(MAX(0._dp,zsncp+zevsnd(jl)),zsncmax)
           zevsnd(jl)=zevsnd(jl)-(psnc(jl)-zsncp)
           zexpt=MAX(0._dp,ptm1(jl)+3._dp-tmelt)*zdtime/zc2
           zexpw=pspeed10(jl)*zdtime/zc3
           zsncmelt(jl)=psnc(jl)*(1._dp-EXP(-zexpt))
           psnc(jl)=psnc(jl)-zsncmelt(jl)
           zsncwind=psnc(jl)*(1._dp-EXP(-zexpw))
           psnc(jl)=psnc(jl)-zsncwind
           zsnowd(jl)=zsnowd(jl)+zsncwind
           tte_corr(jl)=zsncmelt(jl)*zsncfac*zlfdcp(jl)
           !
           !   pwl(jl)=pwl(jl)+zsncmelt(jl) see section 2.5
           !
        ELSE
           psnc(jl)=0._dp
        END IF
210  END DO
     !
     !*    2.2    Snowfall and sublimation on land (excluding glaciers)
     !
     DO 220 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
           IF (psn(jl).LT.0._dp) THEN
              pevwsd(jl)=pevwsd(jl)+psn(jl)
              psn(jl)=0._dp
           END IF
        ELSE
           psn(jl)=0._dp
        END IF
220  END DO
     !
     !*    2.3    Snowfall and sublimation on glaciers and diagnostics
     !
     DO 230 jl=1,kbdim
        IF (lpglac(jl)) THEN
           pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
           palac(jl)=zraind(jl)+zsnowd(jl)+zevttd(jl)
           zrogl(jl)=zraind(jl)
        END IF
230  END DO
     !
     !*    2.4    Snow and glacier melt
     !
!!$ TR     IF (.NOT. lstart) THEN
        DO 240 jl=1,kbdim
           IF (lmask(jl) .AND. time_steps_soil(jl) > 0.5_dp) THEN
              IF (ptsl(jl).GT.tmelt) THEN
                 IF (lpglac(jl)) THEN
                    zsnmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
                    pgld(jl)=pgld(jl)-zsnmel(jl)
                    zrogl(jl)=zrogl(jl)+zsnmel(jl)
                    ptsl(jl)=tmelt
                 ELSE IF (psn(jl).GT.0._dp) THEN
                    zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
                    zsnmel(jl)=MIN(psn(jl),zsmelt)
                    ptsl(jl)=ptsl(jl)-zsnmel(jl)*alf*rhoh2o/pgrndcapc(jl)
                    psn(jl)=MAX(psn(jl)-zsnmel(jl),0._dp)
                 END IF
              END IF
           END IF
240     END DO
!!$ TR     END IF
     !
     !*    2.5    Snow budget and meltwater (glacier-free land only)
     !
     DO 250 jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
           pwl(jl)=pwl(jl)+zsncmelt(jl)                  ! Add melt water on canopy to skin reservoir
           zsnmlt=zsnmel(jl)+MAX(0._dp,pwl(jl)-pwlmx(jl))   ! Excess water on canopy drips to ground and
           pwl(jl)=MIN(pwlmx(jl),pwl(jl))                ! skin reservoir on canopy is set to maximum
           pmlres(jl)=zsnmlt
           zsnmel(jl)=zsnmel(jl)+zsncmelt(jl)
           zsn(jl)=zsn(jl)-zsnmel(jl)
        END IF
250  END DO
     !
     ! Accumulate water fluxes
     ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
     !       accumulation - multiply by zdtime
      DO 601 jl=1,kbdim
         IF (lmask(jl)) THEN
            IF (.NOT. lpglac(jl)) psnacl(jl) = psnacl(jl)  +zsn(jl) * rhoh2o
            psnmel(jl) = psnmel(jl)  +zsnmel(jl) * rhoh2o
            IF (lpglac(jl)) progl(jl)  = progl(jl)   +zrogl(jl)  *rhoh2o
         END IF
 601  END DO
!++mgs
!--- Included for dry deposition ---------------------------------------
      psmelt(1:kbdim) = zsnmel(1:kbdim) * rhoh2o
!--mgs

!!$ TR  END IF
  !
  RETURN

END SUBROUTINE update_surf_down
