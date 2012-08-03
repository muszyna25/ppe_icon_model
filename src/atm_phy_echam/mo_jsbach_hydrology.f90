MODULE mo_jsbach_hydrology

  USE mo_kind, ONLY: dp=>wp
  USE mo_physical_constants, ONLY: alf, cpd, vtmpc2, rhoh2o, tmelt, g=>grav

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_surf_down, update_surf_up

CONTAINS
  !>
  !!

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

SUBROUTINE update_surf_up ( kbdim, itile, delta_time, nsoil &
     , pwl, pcvw, pwlmx &
     , pws, pwsmx &
     , ptsoil, pcvs          &
     , porostd &
     , pevapl,     pevapot, pevwsd                     &
     , prain                           &
     , lmask, lpglac &
     , prunoff, pdrain &
     , pros_hd,    pdrain_hd                   &
     , pmlres &
     , zevwld &
!     Soil moisture layers, evap. fluxes over land
     , zdel, wsi, zwsfc           &
     , dzr, dzs, fksat , fmpot , bclapp, vpor     &
     , etrans, redevap                    &
     , vfc, vpwp, spsi                           &
     , jllog                                     &
     )
!
!     ------------------------------------------------------------------
! pwl           Water content [m] in skin reservoir
!               (vegetation and bare soil)
! pwlmx         Skin reservoir [m] (calculated in *vdiff* as a function
!               of vegetation index and leaf area index)
! pcvw          Skin reservoir fraction (= pwl/pwlmx, see *vdiff*)
! pws           Soil water content [m]
! pwsmx         Water holding capacity [m] of the soil
! ptsoil        Temperature [K] in upper soil layer
! pcvs          Fractional snow cover (function of psn in *physc*)
! porostd       Subgrid standard diviation [m] used in runoff scheme
! pevapl        Total evapotranspiration, including sublimation [kg/m**2/s]
! etrans        Transpiration [kg/m**2/s]
! pevapot       Potential evaporation/sublimation [kg/m**2/s]
! pevwsd        Evapotranspiration without sublimation and evaporation from interception reservoir [kg/m**2]
!               BUT in this routine it is [m], such as zevwld.
! prain         Totalrainfall [kg/m**2/s]
! prunoff       Total runoff [kg/m**2] at non-glacier points (accumul.)
! pdrain        Drainage at non-glacier points [kg/m**2] (accumul.)
! pros_hd       Runoff for HD-Model (does NOT include drainage) [m]
! pdrain_hd     Drainage for HD-Model [m]
! lpglac        Logical glacier mask
!
! The following local variables represent the respective fluxes
! integrated over one timestep (delta_time) and divided by the density of
! water (rhoh2o). Units are m water equivalent.
!
! zraind        Total rain [m] 
! zevwld        Evaporation from the skin reservoir [m]
! zros          Total runoff (including drainage) at non-glacier points
! zdrain        Drainage at non-glacier points
! zsn           Snow budget at non-glacier points (snowfall-subl-melt)
! pmlres        Residual melt water available for infiltration into the
!               non-frozen soil after being intercepted by the
!               skin reservoir
!  
!  ******** 5 layer scheme
!     ***  kdeep = Number of soil layers = 5
!     ***   zdel = Thickness of the 5 soil layers [m]
!     ***  dzrsi = Rooted depth per soil layer (until rooting depth DZR) [m]
!     ***   dzsi = Depth that can be saturated with water per soil layer
!     ***          (until bedrock DZS) [m]
!     ***    dzr = Rooting depth in [m]
!     ***    dzs = Soil depth derived from textures in [m]
!     ***  fksat = Saturated hydraulic conductivity: unit [m/s]
!     ***  fmpot = Matrix potential [m]
!     *** bclapp = Exponent B in Clapp and Hornberger 
!     ***   vpor = Volumetric soil porosity [m/m]   following Cosby et al.
!     ***    vfc = Volumetric soil field capacity [m/m]   
!     ***   vpwp = Volumetric soil wilting point [m/m]   
!     ***   spsi = Soil Pore size distribution index    
!     ***  zwsat = MAXIMALE SPEICHERGEHALT DER SOIL LAYER I (POROSITAET) [M]
!     ***          WICHTIG FUER AUFNAHME DER INFILTRATION IN DEN BODENSCHICHTEN.
!     ***          IST DEFINIERT BIS ZUR TIEFE DES BODENS dzs, WO DER BEDROCK BEGINNT.
!     ***  zwsfc = Field capacity of soil layer I (analog zu ZWSAT) [m]
!     ***  zwpwp = Soil wilting point of soil layer I (analog to ZWSAT) [m]
!     ** aebsoil = Bare soil evaporation in [m], accumulated over TWODT/2
!     ** aetrans = Transpiration in [m] accumulated over one time step 
!     ***          (calculated from etrans which is in [kg/m**2/s])
!     *** redevap  = Diagnostic desired reduction of evapotranspiration
!     ***            due to limited storage in soil layers [kg/m**2/s]
!     *** wsi(I) = Soil moisture content in I. soil layer [m]
!     ***          Is updated with the fluxes calculated in SOILHYD
!     *** wsim1m(i) = Soil Moisture of previous time step in i. Soil layer [m] 
!     *** ZDRAIN = Drainage: Volumen in time step: unit [m]
!
!     *** jllog          relative Index of grid box for log output between 1 and kbdim
!
!       The rest is to be added later on ....
!
!     *SURF* - Updates land values of temperature, moisture and snow.
!              Calculate fluxes of total rain, total snow and evapo-
!              ration from the three reservoirs (SN, WS, WL)
!              convert fluxes (kg/m**2*s) to changes of water levels (m)
!              during timestep delta_time.
!
!     J.F.GELEYN     E.C.M.W.F.     08/06/82.
!     MODIFIED BY
!     C.A.BLONDIN    E.C.M.W.F.    18/12/86.
!     MODIFIED BY L.DUMENIL      MET.INST.HH     20/05/88
!                 J.-P. SCHULZ   MPI - 1997 : Implementation of implicit coupling between 
!                                             land surface and atmosphere.
!     MODIFIED BY E. ROECKNER    MPI - SEPT 1998
!     MODIFIED BY M. ESCH        MPI - APR  1999
!     MODIFIED BY E. ROECKNER    MPI - JAN  2001
!     MODIFIED BY I. Kirchner    MPI - MARCH 2001 date/time control
!     MODIFIED BY E. ROECKNER    MPI - SEP  2002 interception reservoir 
!                                                for snow changed
!     MODIFIED BY L. KORNBLUEH   MPI - JAN  2003 removed MERGE
!     MODIFIED BY ??? R. SCHNUR?        JSBACH Version
!
!     MODIFIED BY S. HAGEMANN    MPI - February  2009 
!     *** Implementation of 5 layer soil hydrology scheme
!
!     MODIFICATION
!
!     PURPOSE
!
!     INTERFACE.
!
!          *SURF* is called from *PHYSC*.
!
!     METHOD.
!
!     EXTERNALS.
!
!          NONE.
!
!     REFERENCE.
!
!          See soil processes' part of the model's documentation for
!     details about the mathematics of this routine.
!
!!$ TR  USE mo_control,        ONLY: ngl
!!$ TR  USE mo_param_switches, ONLY: lsurf
!!$ TR  USE mo_vegetation ! only cvinter is used
!!$ TR  USE mo_radiation
!!$ TR#ifdef STANDALONE
!!$ TR  USE mo_jsbach_comm_to_echam5mods, ONLY: nlat
!!$ TR#endif

  INTEGER,  INTENT(in)    :: kbdim
  INTEGER,  INTENT(in)    :: itile
  REAL(dp), INTENT(in)    :: delta_time
  INTEGER,  INTENT(in)    :: nsoil ! also used as switch for 5 layer calculation
  LOGICAL,  INTENT(in)    :: lmask(kbdim)
  LOGICAL,  INTENT(in)    :: lpglac(kbdim)
  REAL(dp), INTENT(in)    :: pcvw(kbdim)
  REAL(dp), INTENT(in)    :: pwlmx(kbdim)
  REAL(dp), INTENT(in)    :: pwsmx(kbdim)
  REAL(dp), INTENT(in)    :: ptsoil(kbdim) 
  REAL(dp), INTENT(in)    :: pcvs(kbdim)
  REAL(dp), INTENT(in)    :: porostd(kbdim)
  REAL(dp), INTENT(in)    :: pevapl(kbdim)      !! not used
  REAL(dp), INTENT(in)    :: pevapot(kbdim)
  REAL(dp), INTENT(in)    :: prain(kbdim)
  REAL(dp), INTENT(in)    :: pmlres(kbdim)

  REAL(dp), INTENT(inout) :: pwl(kbdim)
  REAL(dp), INTENT(inout) :: pws(kbdim)
  REAL(dp), INTENT(inout) :: pevwsd(kbdim)
  REAL(dp), INTENT(inout) :: prunoff(kbdim)
  REAL(dp), INTENT(inout) :: pdrain(kbdim)

  REAL(dp), INTENT(inout) :: pros_hd(kbdim)
  REAL(dp), INTENT(inout) :: pdrain_hd(kbdim)
  REAL(dp), INTENT(out)   :: zevwld(kbdim) ! From optional part
!
! Soil moisture layers, evap. fluxes over land.
  REAL(dp), INTENT(in)   , optional :: zdel(nsoil)
  REAL(dp), INTENT(inout), optional :: wsi(kbdim,nsoil)
  REAL(dp), INTENT(inout), optional :: zwsfc(kbdim,nsoil)
  REAL(dp), INTENT(in)   , optional :: dzr(kbdim)
  REAL(dp), INTENT(in)   , optional :: dzs(kbdim)
  REAL(dp), INTENT(in)   , optional :: fksat(kbdim)
  REAL(dp), INTENT(in)   , optional :: fmpot(kbdim)
  REAL(dp), INTENT(in)   , optional :: bclapp(kbdim)
  REAL(dp), INTENT(in)   , optional :: vpor(kbdim)
  REAL(dp), INTENT(in)   , optional :: etrans(kbdim)
  REAL(dp), INTENT(out)  , optional :: redevap(kbdim) ! reduction of evapotranspiration from SOILCHANGE
  REAL(dp), INTENT(in)   , optional :: vfc(kbdim)
  REAL(dp), INTENT(in)   , optional :: vpwp(kbdim)
  REAL(dp), INTENT(in)   , optional :: spsi(kbdim)
  INTEGER,  INTENT(inout), optional :: jllog      ! gridbox number for logoutput in SOILHYD
!
!  local variables
!
  INTEGER :: jl,jk
  REAL(dp) ::                                                             &
       zraind (kbdim)      &
     , zros(kbdim),          zdrain(kbdim)

  REAL(dp) ::                                                             &
       zorvari, zorvars, zdrmin, zdrmax, zdrexp, zsmelt,              &
       zmprcp, zwlp, zwptr, zwdtr, zwslim, zconw2, zconw3, zconw4,    &
       zroeff, zbws, zb1, zbm, zconw1, zlyeps, zvol, zprfl,           &
       zlysic, zwsup, zsncp, zexpt, zexpw, zsncmax, zsncwind
  REAL(dp) :: zdtime, zrcp, zsncfac

  ! Variables necessary for evaporation fluxes and 5 soil layers
  REAL(dp), PARAMETER :: cvinter = 0.25  !!$ TR for JSBACH testing
  REAL(dp), PARAMETER :: zeps  = 1.e-10_dp
  REAL(dp), PARAMETER :: zfak1 = 1._dp           ! =1., Dummy for SOILCHANGE (not used)
  INTEGER,  PARAMETER :: kb1   = 1               ! = 1, Digit for 1st field value for transfer to subroutine
  INTEGER,  PARAMETER :: ilog  = 0               ! Switch for logoutput in SOILHYD
!!!!  PARAMETER (jllog=918, ilog=3)        ! No logoutput for ilog=0

  REAL(dp) ::                                                         &
       dzrsi(kbdim, nsoil), dzsi(kbdim, nsoil),                       &
       zwsat(kbdim,nsoil),                        &
       zwpwp(kbdim,nsoil),                                            &
       wsim1m(kbdim, nsoil), aebsoil(kbdim), aetrans(kbdim)

  REAL(dp) :: zinfil(kbdim)    ! Infiltration for transfer to SOILCHANGE - was previously a local scalar 
  INTEGER  :: isch             ! Switch to control SOILCHANGE

!
!      Parameters
!
  zdtime  = delta_time
  zorvari = 100._dp
!!$ TR#ifndef STANDALONE
!!$ TR for JSBACH testing (ngl is number of latitutdes)
!!$ TR  zorvars = 1000._dp*64._dp/ngl
  zorvars = 1000._dp*64._dp/180._dp
!!$ TR#else
!!$  zorvars = 1000._dp*64._dp/nlat
!!$ TR#endif
  zdrmin  = 0.001_dp/(3600._dp*1000._dp)
  zdrmax  = 0.1_dp/(3600._dp*1000._dp)
  zdrexp  = 1.5_dp
  zsncfac = rhoh2o*g/zdtime

!     ------------------------------------------------------------------
!
!*    1.     Convert water fluxes to [m water equivalent * timestep]

  zros   = 0._dp
  zdrain = 0._dp
  zraind = prain * zdtime / rhoh2o
  zevwld = 0._dp
  DO jl=1,kbdim
    IF (lmask(jl)) THEN
!!$      zros(jl) = 0._dp
!!$      zdrain(jl) = 0._dp
!!$      zraind(jl) = prain(jl) * zdtime/rhoh2o
      zevwld(jl) = (1._dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o
    END IF
  ENDDO

! Initialization of working arrays for 5 soil layer scheme
  IF (nsoil == 5) THEN
    CALL SoilDef(kb1, kbdim, kbdim, nsoil,               &
                 dzr, dzs, vpor, vfc, vpwp, lmask, zdel, &
                 dzrsi, dzsi, zwsat, zwsfc, zwpwp)

!   *** Unit conversion for Transpiration [kg/m**2/s] --> [m]
    aetrans(1:kbdim) = etrans(1:kbdim)  * zdtime / rhoh2o
    redevap(1:kbdim) = 0._dp

  ENDIF   

!!$ TR  IF (lsurf) THEN

!     ------------------------------------------------------------------
!*  4.     Water budget
!
!*  4.1    Skin reservoir (vegetation and bare soil)
!
    DO jl=1,kbdim
      IF (lmask(jl) .AND. .NOT. lpglac(jl)) THEN
!
!*      4.1.1  Interception of rain
!
        zmprcp     = MIN(zraind(jl)*cvinter,pwlmx(jl)-pwl(jl))
        zwlp       = pwl(jl)    + zmprcp
        zraind(jl) = zraind(jl) - zmprcp
!
!*      4.1.2  Evaporation or dew collection
!
        pwl(jl)    = MIN(MAX(0._dp,zwlp+zevwld(jl)),pwlmx(jl))
        pevwsd(jl) = pevwsd(jl)-(pwl(jl)-zwlp)
!HAG
!       Get sure that ESKIN comprises only amounts that really change the skin reservoir
        zevwld(jl) = pwl(jl)-zwlp ! This is supposedly the correct way
!        if (nsoil==5) zevwld(jl) = pwl(jl)-zwlp ! Use this for compatibility with older jsbach revisions
      ELSE
        pwl(jl) = 0._dp
      END IF
    ENDDO

!   Calculation of Bare Soil Evaporation
!   Unit of pevwsd is [m]
    IF (nsoil == 5) aebsoil(1:kbdim) = pevwsd(1:kbdim) - aetrans(1:kbdim)
!
!*  4.2    Soil reservoir
!
    DO jl=1,kbdim

      zinfil(jl) = 0._dp
      IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
        zwslim = 0.05_dp*pwsmx(jl)

        zroeff = MAX(0._dp, porostd(jl)-zorvari)   &
                     /(porostd(jl)+zorvars)
        zbws   = MAX(MIN(zroeff,0.5_dp),0.01_dp)
        zb1    = 1._dp+zbws
        zbm    = 1._dp/zb1
        zconw1 = pwsmx(jl)*zb1
        zlyeps = 0._dp
        zvol   = 0._dp

!*    4.2.1  Surface runoff, infiltration and evaporation from soil

        ! Positive Enoskin_nosnow is added to throughfall if the old bucket is used
        IF (nsoil == 1) THEN

          IF (pevwsd(jl) >= 0.0_dp) THEN
            zprfl = pmlres(jl)+zraind(jl)+pevwsd(jl)
          ELSE
            pws(jl) = pws(jl)    + pevwsd(jl)
            zprfl   = pmlres(jl) + zraind(jl)
          END IF

        ELSE

!         ******* Positive EBSoil and Etrans are added to throughfall separately
!         ******* Negative fluxes change soil moisture in SOILCHANGE.
          zprfl = pmlres(jl)+zraind(jl)
          IF (aebsoil(jl) >= 0.0_dp) zprfl = zprfl + aebsoil(jl)
          IF (aetrans(jl) >= 0.0_dp) zprfl = zprfl + aetrans(jl)

        ENDIF

        IF (ptsoil(jl) < tmelt) THEN
          zros(jl) = zprfl
        ELSE
          IF (zprfl > 0._dp .AND. pws(jl) > zwslim) THEN
            IF (pws(jl) > pwsmx(jl)) THEN
              zlyeps = pws(jl)-pwsmx(jl)
            ELSE
              zlyeps = 0._dp
            END IF
            zlysic     = MIN(1._dp,(pws(jl)-zlyeps)/pwsmx(jl))
            zvol       = (1._dp-zlysic)**zbm-zprfl/zconw1
            zros(jl)   = zprfl-(pwsmx(jl)-pws(jl))
            IF (zvol > 0._dp) zros(jl) = zros(jl)+pwsmx(jl)*zvol**zb1
            zros(jl)   = MAX(zros(jl),0._dp)
            zinfil(jl) = zprfl-zros(jl)
          ELSE
            zros(jl)   = 0._dp
            zinfil(jl) = zprfl
          END IF
          ! In 5 layer scheme, infiltration changes soil moisture in SOILCHANGE.
          IF (nsoil == 1) pws(jl) = pws(jl)+zinfil(jl)
        END IF
!
!*    4.2.2  Drainage and total runoff
!
      ENDIF
    ENDDO
!
!   ****** Drainage if the old bucket is used
    IF (nsoil == 1) THEN

      DO jl=1,kbdim

        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN

!         *** Term shifting moved from the beginning of the original soil hydrology loop
          zwdtr  = 0.90_dp*pwsmx(jl)
          zwslim = 0.05_dp*pwsmx(jl)
          zconw2 = pwsmx(jl)-zwdtr
          zconw3 = zdrmax-zdrmin

          IF (pws(jl) <= zwslim) THEN
            zdrain(jl) = 0._dp
          ELSE
            IF (ptsoil(jl) > tmelt) THEN
              zdrain(jl) = zdrmin*MAX(pws(jl)/pwsmx(jl),0.7_dp)
              IF (pws(jl) > zwdtr) THEN
                 zdrain(jl) = zdrain(jl)+zconw3*((pws(jl)-zwdtr)/zconw2)**zdrexp
              END IF
              zdrain(jl) = zdrain(jl)*zdtime
              zdrain(jl) = MIN(zdrain(jl),pws(jl)-zwslim)
              pws(jl)    = pws(jl)-zdrain(jl)
            ELSE
              zdrain(jl) = 0._dp
            END IF
          END IF
          zwsup    = MAX(pws(jl)-pwsmx(jl),0._dp)
          pws(jl)  = pws(jl)-zwsup
          zros(jl) = zros(jl)+zdrain(jl)+zwsup
        ELSE
           pws(jl) = 0._dp
        END IF
      ENDDO

    ELSE ! nsoil == 1
      !  Change of soil moisture in 5 layers
      !  change belongs to one time step = zdtime

      ! Soil moisture of previous time step: wsim1m
      wsim1m(:,:) = wsi(:,:)

      ! Zero initialization for drainage field before SOILCHANGE
      zdrain(1:kbdim) = 0._dp

      ! Soilchange expects aetrans and aebsoil in [m]
      isch = 1     ! Changes by evaporation fluxes
                   ! Soilchange only uses negative evaporation fluxes as positive fluxes,
                   ! e.g. by dew formation, are part of the throughfall zprfl.
      CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
            isch, zdtime, zfak1, lmask,             &
            wsi , wsim1m,                           &
            dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
            zinfil, zdrain, aetrans, aebsoil,       &
            redevap)
      ! Unit of zinfil is [m] but soilchange expects m/s
      ! For numerical reasons, half of the water is infiltrated 
      ! before SOILHYD, half if it after
      zinfil(1:kbdim) = zinfil(1:kbdim) / zdtime / 2.

      isch = 2     ! Changes by infiltration

      CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
            isch, zdtime, zfak1, lmask,             &
            wsi , wsim1m,                           &
            dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
            zinfil, zdrain, aetrans, aebsoil,       &
            redevap)

      ! Drainage calculation [m] with 5 layer soil scheme
      CALL SoilHyd(kb1, kbdim, kbdim, nsoil,      &
           zdtime, lmask, ilog, jllog,            &
           wsi, dzsi, zwsat, zwsfc,               &
           fksat, vpor, bclapp, fmpot,            &     
           zdrain, spsi, zwpwp, zdel              &
           )

      ! Second half of the water is infiltrated after SOILHYD
      isch = 2     ! Changes by infiltration
      CALL SoilChange(kb1, kbdim, kbdim, nsoil,     &
            isch, zdtime, zfak1, lmask,             &
            wsi , wsim1m,                           &
            dzr, dzrsi,    dzsi, zwsat, zwsfc,      &
            zinfil, zdrain, aetrans, aebsoil,       &
            redevap)

      ! Does overflow occurs due to infiltration in the layers? --> added to drainage
      DO jk=1, nsoil
        DO jl=kb1, kbdim
          IF (dzsi(jl,1) > 0) THEN
            IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN
              zdrain(jl) = zdrain(jl) + wsi(jl,jk) - zwsfc(jl,jk)
              wsi(jl,jk) = zwsfc(jl,jk)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      ! WS must be re-derived from wsi after SOILHYD
      pws(1:kbdim) = 0.0_dp
      DO jk=1, nsoil
        DO jl=kb1, kbdim
          IF (lmask(jl)) THEN
            IF (dzrsi(jl,jk) >= zdel(jk)) THEN
              pws(jl) = pws(jl) + wsi(jl,jk)
            ELSE IF (dzrsi(jl,jk) > 0._dp) THEN
              pws(jl) = pws(jl) + wsi(jl,jk) * dzrsi(jl,jk)/dzsi(jl,jk)
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO jl=1,kbdim
        IF (lmask(jl) .AND. .NOT.lpglac(jl)) THEN
          zros(jl) = zros(jl)+zdrain(jl)
        ENDIF
      ENDDO

    ENDIF ! end of nsoil if clause.
!
!*  4.2.3  Runoff and drainage for the HD-Model [m]
!
    DO jl=1,kbdim
      IF (lmask(jl)) THEN
        pros_hd(jl)   = zros(jl)-zdrain(jl)
        pdrain_hd(jl) = zdrain(jl)
      END IF
    END DO
!
!     ------------------------------------------------------------------
!
!*    6.     Accumulate fluxes for diagnostics
!
!     6.1    Water fluxes
!
    ! Note: conversion from m water equivalent to kg/m^2s - multiply by rhoh2o / zdtime
    !       accumulation - multiply by zdtime
    DO jl=1,kbdim
      IF (lmask(jl)) THEN
        prunoff(jl) = prunoff(jl) +zros(jl)   *rhoh2o
        pdrain(jl)  = pdrain(jl)  +zdrain(jl) *rhoh2o
      END IF
    END DO
!
!     ------------------------------------------------------------------
!!$ TR  END IF

END SUBROUTINE update_surf_up

SUBROUTINE soilhyd(kidia, kfdia, klon, kdeep, zdiagt, loland, ilog, jllog,  &
    wsi, dzsi, zwsat, zwsfc, fksat, vpor, bclapp, fmpot,  &
    zdrain, spsi, zwpwp, zdel )
!     **************************************************************************

!     ******* SOIL HYDROLOGY - Routine to be build into REMO
!     *******   to calculate soil hydrology for 5 soil layers
!     ******* SOIL HYDROLOGY = percolation and vertical diffusion of water
!     ******* Type of diffusiont:  Richtmeyr/Morton

!     *** PROGRAMMING AND DEVELOPMENT: STEFAN HAGEMANN
!     ****
!     *** VS. 1.0 - Test version for soil processes without run-off and others
!              percolation dependent of IPERC type
!              REMO calculations in [m], but output is in [mm]/output interval.
!     ***
!     *** VS. 1.1 - Version for REMO
!     ***           Units must be made consistent.

!     *** VS. 1.2 - Limitation to IPERC=2, IDIFF=2
!     ***           --> alternatives removed

!     *** VS. 1.3 - Last loop splitted and loop order changed to get longest loop inside

!     *** VS. 1.4 - September 2007 - HAG
!     ***           Loop from Numerical Recipes changed to get longest loop inside.
!     ***               Therefore were ZTRI(KDEEP) and ZDC changed into fields
!     ***               ZTRI(KLON,KDEEP) and ZDC(KLON,KDEEP).
!     ***           Inconsistency for percolation overflow of layers eliminated
!     ***               now overflow exists if WSI > ZWSFC (instead of ZWSAT as used before)
!     ***           Test for 1E-10. was changed due to the calculation accuracy of the NEC
!     ***               at WSMX-1.E-10 < ZEPS (1E-15)
!     ***
!     ***
!     *** Fluxes are always calculated to be valid two time steps ahead.
!     *** Therefore the fluxes must be divided by 2 for accumulation.
!     *** The calculations themselves uses volumes, not fluxes.
!     *** The flux is calculated so that later the unit of the flux is the accumulated volume
!     *** per time step.
!     ***
!     *** VS. 1.5 - December 2007 - HAG
!     ***           Use of ZWSFC instead of ZWSAT for percolation

!     ***           Juli 2008 - S. HAGEMANN (R. PODZUN)
!     ***           Removal of all occurances of ZVFC in soilhyd.f and surf.f
!     ***           UNSAFE! Use may cause problems. on nec no obvious errors were caused.

!     *** VS. 1.6 - March 2009 - HAG
!     ***           Introduction of maximum percolation rate and solving the percolation
!     ***           equation with midpoint method instead of euler method.

!     *** VS. 1.7 - March-May 2009 - HAG
!     ***           Implementation of van genuchetn equation for percolation.
!     ***           preliminary using iperc to allow switching to previuos
!     ***               Clapp & Hornberger method.
!     ***           Allow ECHAM drainage for lowest layer above bedrock (IDRECH).
!     ***           Allow reduction of hydraulic conductivity with depth for percolation.

!     *** VS. 1.8 - October 2009 - HAG
!     ***           Percolation overflow is not allowed anymore from layers below the surface
!     ***               layer. This has led to a too strong drying
!     ***               of soil layers in wet areas, especially from the second lowest layer
!     ***               above the bedrock. Now it assumed that water piles upwards from
!     ***               SATURATEDsoil5layer.f90 layers, very similar to a rising water table. 
!     ***               Whether this may actually be used to describe a water table, has to be investigated.
!     ***               But this can be also done offline.

!     *** VS. 1.9 - June 2010 - HAG
!     ***            Initialization of zdrain with zero must happen in calling routine before
!     ***                soilchange is called.
!
!     *** vs. 1.10 - Jay 2011 - HAG
!     ***            Code converted using TO_F90 by Alan Miller
!     ***            plus REAL --> REAL(dp) and fine tuning afterwards
!
!     SOILHYD
!     ***  FKSAT = Saturated hydraulic conductivity: unit [m/s] -> Must be checked in ECMGBDT.F
!     ***  FMPOT = Matrix potential [m]
!     *** BCLAPP = Exponent B in Clapp and Hornberger
!     ***   VPOR = Volumetric soil porosity [m/m]
!     ***   SPSI = Soil pore size distribution index used in Van Genuchten

!     *** ZDIAGT = Time step: unit [s]   ##### org. was [day]

!     *** ZDRAIN = Drainage: Volumen during time step: unit [m]
!     ***  ZWSAT = Maximum storage in soil layer I (porosity) [m]
!     ***          Important for the take up of the infiltration in the soil layers
!     ***          Defined down to the depth of the soil - that is where the bedrock begins.
!     ***  ZWSFC = Field capacity for soil layer I (analog to ZWSAT) [m]
!     ***  ZWPWP = Texture depending wilting point of the soil layer I (analog to ZWSAT) [m]
!     ***
!     ***  ZPERC = Downward flux from soil layer I: First in m/day, then in [m].
!     ***          First because of gravitational drainage.
!     ***  ZDIFF = Diffusion of water between soil layer i and layer i+1 [m^2/day]
!     ***   ZMVG = Exponent M in Van Genuchten scheme = SPSI / (SPSI+1)
!     ***   ZDEL = Thickness of the 5 soil layers [m]
!     *** ZKREDU = Reduction factor of hydraulic conductivity with depth for percolation
!     ***  IREDU = Reduction of conductivity with depth yes/no (1/0)
!     ***
!     ***  NSOIL - Number of soil layers
!     *** WSI(I) = Soil moisture content in I. soil layer [m]
!     ***          Updated through the fluxes calculated in SOILHYD
!     ***  ILOG =  1 --> Log output to StdOut with keword SOILFLUX for grid cell JLLOG
!     ***                5 * DIFF. FLUXES, 5 PERC. FLUXES, ZDRAIN
!     ***          2 --> Log output to StdOut with keyword SOILHYD1 for grid cell JLLOG

!     ***  IPERC = Type of percolation
!     ***          1  ECHAM4 drainage (removed, not available)
!     ***          2  Clapp & Hornberger
!     ***          3  Van Genuchten methode 2, after these
!     ***
!     ***  IDRECH = ECHAM drainage
!     ***          0   No, drainage is percolation from lowest layer above bedrock.
!     ***          1   Yes, for lowest layer above bedrock (replaces percolation)
!     ***          2   For all layers above bedrock (in addition to percolation, except
!     ***              for lowest layer above bedrock where it replaces percolation.

INTEGER,  INTENT(IN)                     :: kidia
INTEGER,  INTENT(IN)                     :: kfdia
INTEGER,  INTENT(IN)                     :: klon
INTEGER,  INTENT(IN)                     :: kdeep
REAL(dp), INTENT(IN)                     :: zdiagt
LOGICAL,  INTENT(IN)                     :: loland(klon)
INTEGER,  INTENT(IN)                     :: ilog
INTEGER,  INTENT(INOUT)                  :: jllog
REAL(dp), INTENT(INOUT)                  :: wsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: dzsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: zwsat(klon,kdeep)
REAL(dp), INTENT(IN)                     :: zwsfc(klon,kdeep)
REAL(dp), INTENT(IN)                     :: fksat(klon)
REAL(dp), INTENT(IN)                     :: vpor(klon)
REAL(dp), INTENT(IN)                     :: bclapp(klon)
REAL(dp), INTENT(IN)                     :: fmpot(klon)
REAL(dp), INTENT(INOUT)                  :: zdrain(klon)
REAL(dp), INTENT(IN)                     :: spsi(klon)
REAL(dp), INTENT(IN)                     :: zwpwp(klon,kdeep)
REAL(dp), INTENT(IN)                     :: zdel(kdeep)

REAL(dp), PARAMETER :: zeps   = 1.e-15_dp
INTEGER,  PARAMETER :: iperc  = 3
INTEGER,  PARAMETER :: idrech = 2
INTEGER,  PARAMETER :: iredu  = 0
REAL(dp), PARAMETER :: credu  = 1.0_dp

REAL(dp) :: zdum, zdum2, zdum3
REAL(dp) :: zperc(klon,kdeep), zdum4, zctsat
REAL(dp) :: zmvg(klon), zwrel
REAL(dp) :: zdiflog(10), wslog(kdeep)
!     *** ECHAM-Drainage-Parameter
REAL(dp) :: zdrmin, zdrmax, zdrexp, zwdmin, zwdcri

!     *** Diffusion
REAL(dp) :: zda(klon,kdeep), zdb(klon,kdeep), zdiff(klon,kdeep-1)
REAL(dp) :: ztri(klon,kdeep), zdc(klon,kdeep)

REAL(dp) :: zpmax, zp2         ! Maximum percolation rate
REAL(dp) :: zkredu(kdeep)      ! KSAT reduction factor with depth

!     *** Dummy variables, running indices, etc.
INTEGER  :: i, jl

  zwdmin = 0.05_dp
  zwdcri = 0.9_dp
  zdrmin = 0.001_dp / (3600._dp*1000._dp)
  zdrmax = 0.1_dp / (3600._dp*1000._dp)
  zdrexp = 1.5_dp
  IF (iredu == 1) THEN
    zkredu(1) = (1._dp + EXP(-credu*zdel(1)) )/2._dp
    DO i=2, kdeep
      zkredu(i) = (EXP(-credu*SUM(zdel(1:i-1))) +  &
          EXP(-credu*SUM(zdel(1:i))) )/2._dp
    END DO
  ELSE
    zkredu(1:kdeep) = 1._dp
  END IF

!     ***** Initialization of local fields
  zperc(kidia:kfdia,1:kdeep) = 0._dp
  zmvg(kidia:kfdia) = spsi(kidia:kfdia) / (spsi(kidia:kfdia)+1._dp)

!     *** Percolation = Gravitational drainage - calculated from bottom upward.
!     *** Maximum percoplation due to numeric instability at large time steps
!     *** and small spatial distances.
  IF (iperc == 2) THEN          ! Clapp & Hornberger-formulation
    DO i=1, kdeep
      DO jl=kidia, kfdia
!cc        IF (DZSI(JL,I).GT.0) THEN
        IF (zwsfc(jl,i)-1.e-10_dp > zeps) THEN
          
!         *** MAXIMUM PERCOLATION RATE ZPMAX
          zdum = ( zwsfc(jl,i) -  &
              fksat(jl)*zkredu(i) * zdiagt - zwpwp(jl,i) ) /  &
              ( zwsfc(jl,i) - zwpwp(jl,i))
          IF (zdum > zeps) THEN
            zp2 = fksat(jl)*zkredu(i) * zdum ** (2._dp*bclapp(jl)+3._dp)
          ELSE
            zp2=0._dp
          END IF
          zpmax  = MAX( (fksat(jl)*zkredu(i) + zp2)/2._dp , 0._dp)
          zctsat = MIN( MAX((wsi(jl,i)-zwsfc(jl,i))/zdiagt ,0._dp)  &
              / (fksat(jl)*zkredu(i)) /zdiagt , 1._dp)
          zpmax  = zctsat * fksat(jl)*zkredu(i) + (1._dp-zctsat) * zpmax
          
          zwrel  = MAX(wsi(jl,i)-zwpwp(jl,i), 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          zwrel  = MIN(zwrel, 1._dp)
          
          IF (zwrel-1.e-10_dp > zeps) THEN
            zp2   = fksat(jl)*zkredu(i) * zdiagt / 2._dp * zwrel ** (2._dp*bclapp(jl)+3._dp)
            zdum2 = MAX( MAX(wsi(jl,i)-zp2,0._dp)  -zwpwp(jl,i), 0._dp)/  &
                (zwsfc(jl,i)-zwpwp(jl,i))
            zdum2 = MIN(zdum2, 1._dp)
            zperc(jl,i) = fksat(jl)*zkredu(i) * zdum2 ** (2._dp*bclapp(jl)+3._dp)
            zperc(jl,i) = MIN(zperc(jl,i) , zpmax )
          ELSE
            zperc(jl,i) = 0._dp
          END IF
        END IF
      END DO
    END DO
  ELSE IF (iperc == 3) THEN          ! Van Genuchten-formulation
    DO i = 1, kdeep
      DO jl=kidia, kfdia
        IF (zwsfc(jl,i)-1.e-10_dp > zeps) THEN
          
!         *** MAXIMUM PERCOLATION RATE ZPMAX
          zdum = MAX(zwsfc(jl,i)-zwpwp(jl,i) -  &
              fksat(jl)*zkredu(i)*zdiagt, 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          IF (zdum > zeps) THEN
            zp2 = fksat(jl)*zkredu(i) * SQRT(zdum) * ( 1._dp - ( 1._dp -  &
                zdum**(1._dp/zmvg(jl)) )**zmvg(jl) )**2._dp
          ELSE
            zp2=0._dp
          END IF
          zpmax = MAX( (fksat(jl)*zkredu(i) + zp2)/2._dp , 0._dp)
          zctsat = MIN( MAX((wsi(jl,i)-zwsfc(jl,i))/zdiagt ,0._dp)  &
              / (fksat(jl)*zkredu(i)) /zdiagt , 1._dp)
          zpmax = zctsat * fksat(jl)*zkredu(i) + (1._dp-zctsat) * zpmax
          
!         *** Midpoint-Method
          
          zwrel = MAX(wsi(jl,i)-zwpwp(jl,i), 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          zwrel = MIN(zwrel, 1._dp)
          
          zdum2 = fksat(jl)*zkredu(i) * SQRT(zwrel) * ( 1._dp - ( 1._dp -  &
              zwrel**(1._dp/zmvg(jl)) )**zmvg(jl) )**2._dp
          zdum3 = wsi(jl,i) - MIN(zdum2,zpmax) * zdiagt / 2._dp
          
          zdum4 = MAX(zdum3-zwpwp(jl,i), 0._dp) / (zwsfc(jl,i)-zwpwp(jl,i))
          zdum4 = MIN(zdum4, 1._dp)
          
          zperc(jl,i) = fksat(jl)*zkredu(i) * SQRT(zdum4) * ( 1._dp -  &
              ( 1._dp - zdum4**(1._dp/zmvg(jl)) )**zmvg(jl) )**2._dp
          
          zperc(jl,i) = MIN(zperc(jl,i) , zpmax)
        ELSE
          zperc(jl,i) = 0._dp
        END IF
      END DO
    END DO
  END IF

! *** Using ECHAM4 drainage calculation? for the lowest layer above bedrock
! *** this replaces the percolation value.

  IF (idrech == 1 .OR. idrech == 2) THEN
    DO i = 1, kdeep-1
      DO jl=kidia, kfdia
        IF (dzsi(jl,1) > 0) THEN
          IF (dzsi(jl,i) > 0) THEN
            IF (dzsi(jl,i+1) <= 0) THEN           ! *** Lowest layer above bedrock
              
!cc              ZDUM = ZWSFC(JL,I) * ZWDMIN
              zdum = zwpwp(jl,i)
              zdum2 = zwsfc(jl,i) * zwdcri
              IF (wsi(jl,i) <= zdum) THEN
                zperc(jl,i) = 0._dp
              ELSE
                zdum3 = zdrmin * (wsi(jl,i)-zdum) / (zwsfc(jl,i)-zdum)
                IF (wsi(jl,i) > zdum2) THEN
                  zdum4 = (zdrmax-zdrmin)*( (wsi(jl,i)-zdum2)  &
                      / (zwsat(jl,i) - zdum2) )**zdrexp
                  zdum3 = zdum3 + MIN(zdum4, (wsi(jl,i)-zdum2)/zdiagt)
                END IF
                zperc(jl,i) = MIN(zdum3, (wsi(jl,i)-zdum)/zdiagt)
              END IF
            ELSE IF (idrech == 2) THEN            ! *** Other layers for IDRECH=2
              
!             *** Drainage is calculated with percolation removed from WSI
              zp2 = wsi(jl,i) - zperc(jl,i) * zdiagt
              zdum = zwpwp(jl,i)
              zdum2 = zwsfc(jl,i) * zwdcri
              IF (zp2 > zdum) THEN         ! Drainage occurs only above wilting point
                zdum3 = zdrmin * (zp2-zdum) / (zwsfc(jl,i)-zdum)
                IF (zp2 > zdum2) THEN
                  zdum4 = (zdrmax-zdrmin)*( (zp2-zdum2)  &
                      / (zwsat(jl,i) - zdum2) )**zdrexp
                  zdum3 = zdum3 + MIN(zdum4, (zp2-zdum2)/zdiagt)
                END IF
                
!               *** Except for the lowest layer above bedrock, drainage is removed from
!               *** soil moisture layer wsi before diffusion occurs
!               *** here, drainage is used as accumulated value over ZDIAGT
                zdum4 = MIN(zdum3 * zdiagt, (zp2-zdum))
                wsi(jl,i) = wsi(jl,i) - zdum4
                zdrain(jl) = zdrain(jl) + zdum4
              END IF
            END IF                                 ! Lowest layer / IDRECH=2
          END IF
        END IF
      END DO
    END DO
    
!   *** Lowest layer (=5)
    DO jl=kidia, kfdia
      IF (dzsi(jl,kdeep) > 0) THEN
!cc          ZDUM = ZWSFC(JL,KDEEP) * ZWDMIN
        zdum = zwpwp(jl,kdeep)
        zdum2 = zwsfc(jl,kdeep) * zwdcri
        IF (wsi(jl,kdeep) <= zdum) THEN
          zperc(jl,kdeep) = 0._dp
        ELSE
          zdum3 = zdrmin * (wsi(jl,kdeep)-zdum) / (zwsfc(jl,kdeep)-zdum)
          IF (wsi(jl,kdeep) > zdum2) THEN
            zdum4 = (zdrmax-zdrmin)*( (wsi(jl,kdeep)-zdum2)  &
                / (zwsat(jl,kdeep) - zdum2) )**zdrexp
            zdum3= zdum3 + MIN(zdum4, (wsi(jl,kdeep)-zdum2)/zdiagt)
          END IF
          zperc(jl,kdeep) = MIN(zdum3, (wsi(jl,kdeep)-zdum)/zdiagt)
        END IF
      END IF
    END DO
  END IF

! *** Original unit of percolation fluxes (=unit of FKSAT) was [m/day]
! *** Now it's accumulatied [m] over ZDIAGE (unit of ZDIAGT=[s]
! *** --> Take time step into account : * ZDIAGT
  zperc(kidia:kfdia,1:kdeep) = zperc(kidia:kfdia,1:kdeep) * zdiagt

! *** Diffusion of water *******************************************

! *** Numerical recipes/Richtmyer & -Morton diffusion  *************

  DO i = 1, kdeep-1
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        
!       ***  Soil moisture diffusivity [m^2/day]
!       ***  method c: diffusivity of both levels is averaged
!       ***            to calculate diffusivity between layers
        IF (wsi(jl,i+1) > 0) THEN
          zdum = wsi(jl,i+1)/dzsi(jl,i+1)
          
!         *** Calculating the diffusivity of layer i+1
          zdum4 = zdum / vpor(jl)
          IF (zdum4-1.e-10_dp > zeps) THEN
            zdum2 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                zdum * zdum4 ** (bclapp(jl)+3.)
          ELSE
            zdum2 = 0._dp
          END IF
        ELSE
          zdum2=0._dp
        END IF
        IF (wsi(jl,i) > 0._dp) THEN
          zdum = wsi(jl,i)/dzsi(jl,i)
          
!         *** Calculating the diffusivity of layer i
          zdum4 = zdum / vpor(jl)
          IF (zdum4-1.e-10_dp > zeps) THEN
            zdum3 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                zdum *  zdum4 ** (bclapp(jl)+3._dp)
          ELSE
            zdum3 = 0._dp
          END IF
        ELSE
          zdum3=0._dp
        END IF
        
!       *** Calculating the diffusivity at the bottom of layer i
        zdiff(jl,i) = (zdum2+zdum3)*0.5_dp
      END IF
    END DO
  END DO

! *** Calculation of diffusion coefficients
! *** Deepest layer KDEEP = 5
  DO jl=kidia, kfdia
    zda(jl, kdeep) = 0._dp
    IF (dzsi(jl,kdeep) > 0._dp) THEN
      
!     *** [ZDIAGT] = s
      zdb(jl, kdeep) = zdiff(jl,kdeep-1) *zdiagt /dzsi(jl,kdeep)  &
          / (dzsi(jl,kdeep)+dzsi(jl,kdeep-1)) * 2._dp
    ELSE
      zdb(jl, kdeep) = 0._dp
    END IF
  END DO

! *** Layers i = 4,3,2
  DO i=kdeep-1, 2,-1
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        IF (dzsi(jl,i+1) > 0._dp) THEN
          zda(jl,i) = zdiff(jl,i) * zdiagt / dzsi(jl,i) /  &
              (dzsi(jl,i)+dzsi(jl,i+1)) * 2._dp
        ELSE
          zda(jl,i) = 0._dp
        END IF
        zdb(jl,i) = zdiff(jl,i-1) * zdiagt / dzsi(jl,i) /  &
            (dzsi(jl,i)+dzsi(jl,i-1)) * 2._dp
      ELSE
        zda(jl, i) = 0._dp
        zdb(jl, i) = 0._dp
      END IF
    END DO
  END DO

! *** Layer 1
  DO jl=kidia, kfdia
    zdb(jl, 1) = 0._dp
    IF (dzsi(jl,1) > 0) THEN
      IF (dzsi(jl,2) > 0) THEN
        zda(jl,1) = zdiff(jl,1) * zdiagt / dzsi(jl,1) /  &
            (dzsi(jl,1)+dzsi(jl,2)) * 2._dp
      ELSE
        zda(jl,1) = 0._dp
      END IF
    ELSE
      zda(jl, 1) = 0._dp
    END IF
  END DO

! *** Temporary storage for log output
  IF (ilog == 1 .OR. ilog == 3) wslog(1:kdeep) = wsi(jllog, 1:kdeep)

! *** Routine TRIDIAG from Numerical Recipes, p. 43
! *** -ZDA = CI, -ZDB=AI, ZDA+ZDB+1=BI, WSI(T)=RI, WSI(T+1)=UI
  DO jl=kidia, kfdia
    IF (dzsi(jl,1) > 0) THEN
      zdc(jl,1) = zda(jl,1)+zdb(jl,1)+1._dp
      wsi(jl,1) = wsi(jl,1) / zdc(jl,1)
    END IF
  END DO
! *** Decomposition and forward substitution
  DO i=2, kdeep
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        
        IF (dzsi(jl,i) > 0) THEN
          ztri(jl,i) = -zda(jl,i-1) / zdc(jl,i-1)
          zdc(jl,i) = zda(jl,i)+zdb(jl,i)+1._dp + zdb(jl,i)*ztri(jl,i)
          wsi(jl,i) = (wsi(jl,i)/dzsi(jl,i) + zdb(jl,i)  &
              * wsi(jl,i-1)/dzsi(jl,i-1) ) / zdc(jl,i) * dzsi(jl,i)
        END IF
      END IF
    END DO
  END DO

  DO i=kdeep-1,1,-1       ! Backsubstitution
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        
        IF (dzsi(jl,i+1) > 0) THEN
          wsi(jl,i) = wsi(jl,i) - ztri(jl,i+1)*wsi(jl,i+1)*  &
              (dzsi(jl,i)/dzsi(jl,i+1))
        END IF
      END IF
    END DO
  END DO

! *** Log output of fluxes ?
  IF (ilog == 1 .OR. ilog == 3) THEN
    DO i=1, kdeep
      zdiflog(i) = wslog(i) - wsi(jllog, i)
    END DO
  END IF

! *** If layer overflow occurs after diffusion, it is added to drainage.
  DO i = 1, kdeep
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        IF (wsi(jl,i) > zwsfc(jl,i)) THEN
          zdrain(jl) = zdrain(jl) + wsi(jl,i) - zwsfc(jl,i)
          wsi(jl,i) = zwsfc(jl,i)
        END IF
      END IF
    END DO
  END DO

  IF (ilog == 2 .OR. ilog == 3) THEN
    WRITE(*,'(A10,1X, 11(E16.9E2, 1X) )')  &
        "SOILHYD1: ", fksat(jllog)*1000._dp, wsi(jllog,1:kdeep)*1000._dp,  &
        zwsfc(jllog,1:kdeep)*1000._dp
  END IF

! ******* NEW SOIL LAYER SOIL MOISTURE CONTENTS *********************

! *** Percolation per time step must be smaller than the usable soil 
! *** layer contents.
  DO i = 1, kdeep
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0 .AND. zperc(jl,i) > 0) THEN
        zdum = zwdmin*zwsfc(jl,i)
        IF (wsi(jl,i)-zperc(jl,i) < zdum) THEN
          IF (wsi(jl,i) > zdum) THEN
            zperc(jl,i) = wsi(jl,i) - zdum
            wsi(jl,i) = zdum
          ELSE
            zperc(jl,i) = 0._dp
          END IF
        ELSE
          wsi(jl,i) = wsi(jl,i) - zperc(jl,i)
        END IF
      END IF
    END DO
  END DO

! *** Percolation added to downward layer or drainage
! *** if lowest layer or bedrock is reached.
  DO i = 1, kdeep-1
    DO jl=kidia, kfdia
      IF (dzsi(jl,1) > 0) THEN
        IF (dzsi(jl,i+1) > 0) THEN
          wsi(jl,i+1) = wsi(jl,i+1) + zperc(jl,i)
        ELSE
          zdrain(jl) = zdrain(jl) + zperc(jl,i)
          zperc(jl,i) = 0._dp
        END IF
      END IF
    END DO
  END DO

  DO jl=kidia, kfdia
    IF (dzsi(jl,1) > 0) THEN
      zdrain(jl) = zdrain(jl) + zperc(jl,kdeep)
    END IF
  END DO

! *** Percolation overflow is not allowed from layers below the surface layer.
! *** Instead it is assumed that water piles upwards from saturated layers.
  DO i = kdeep,2, -1
    DO jl=kidia, kfdia
      IF (dzsi(jl,i) > 0) THEN
        IF (wsi(jl,i) > zwsfc(jl,i)) THEN
          zperc(jl,i-1) = MAX( zperc(jl,i-1) - wsi(jl,i) + zwsfc(jl,i), 0._dp)
          wsi(jl,i-1) = wsi(jl,i-1) + wsi(jl,i) - zwsfc(jl,i)
          wsi(jl,i) = zwsfc(jl,i)
        END IF
      END IF
    END DO
  END DO

! *** Overflow from the surface layer enters drainage. this should not happen,
! *** so this is just for security.
  DO jl=kidia, kfdia
    IF (dzsi(jl,1) > 0) THEN
      IF (wsi(jl,1) > zwsfc(jl,1)) THEN
        zdrain(jl) = zdrain(jl) + wsi(jl,1) - zwsfc(jl,1)
        wsi(jl,1) = zwsfc(jl,1)
      END IF
    END IF
  END DO

! *** Log output of the fluxes ?
  IF (ilog == 1 .OR. ilog == 3) THEN
    DO i=1, kdeep
      zdiflog(kdeep+i) = zperc(jllog,i)
    END DO
    WRITE(*, '(A10,1X,11(E15.6E3,1X) )') &
       "SOILFLUX",zdiflog(1:10)*1000._dp, zdrain(jllog)*1000._dp
  END IF

! *** Log output of WSI!
  IF (ilog == 2 .OR. ilog == 3) THEN
    WRITE(*,'(A10,1X,10(E16.9E2,1X) )' ) &
       "SOILHYD2: ", wsi(jllog,1:kdeep)*1000._dp, dzsi(jllog,1:kdeep)*1000._dp
  END IF

END SUBROUTINE soilhyd

SUBROUTINE soildef (kidia, kfdia, klon, kdeep,  &
    dzr, dzs, vpor, vfc, vpwp, loland, zdel, dzrsi, dzsi, zwsat, zwsfc, zwpwp)
!     ***********************************************************************

!     *** Routine defining various fields for the 5-layer soil scheme.

!     ******** Version 1.0 - Juli 2006
!              Programming and development: Stefan Hagemann
!     ******** Version 1.1 - Juni 2007
!              Optimization : Ralf Podzun
!     ******** Version 1.2 - September 2007
!              A loop optimized : Stefan Hagemann
!     ******** Version 1.3 - Oktober 2007
!              Minor technical correction without effect: Stefan Hagemann
!     ******** Version 1.4 - January 2008
!              Minor technical correction without effect: Stefan Hagemann
!     ******** Version 1.5 - February 2008 - Stefan Hagemann
!              Minor technical correction that may effect only 1-2 gridboxes
!              where DZR=0. (near lakes)
!     ******** Version 1.6 - March 2009 - Stefan Hagemann
!              Volumetric soil field capacity vfc and wilting point vpwp are given
!              via input list and do not need to be computed in soildef.
!              new output: ZWPWP, WSMX eliminated.

!     ***   ZDEL = Thickness of the 5 soil layers [m]
!     ***  DZRSI = Rooted depth per soil layer (until rooting depth DZR) [m]
!     ***   DZSI = Depth that can be saturated with water per soil layer
!     ***          (until bedrock DZS) [m]
!     ***    DZR = Rooting depth in [m]
!     ***    DZS = Soil depth derived from textures in [m]
!     ***   VPOR = Volumetric porosity [m water / m soil]
!     ***    VFC = Volumetric field capacity [m water / m soil]
!     ***   VPWP = Volumetric permanent wilting point [m water / m soil]
!     ***  KDEEP = Number of soil layers
!     ***  ZWSAT = Maximum storage capacity of soil layer I (porosity) [m]
!     ***          Important for the uptake of infiltration in the soil layers.
!     ***          Defined to the depth of the soil, that is, where the bedrock begins.
!     ***  ZWSFC = Field capacity of soil I (analog to ZWSAT) [m]
!     ***  ZWPWP = Texture depending wilting point of soil layer I (analog to ZWSAT) [m]
!
!
INTEGER,  INTENT(IN)                     :: kidia
INTEGER,  INTENT(IN)                     :: kfdia
INTEGER,  INTENT(IN)                     :: klon
INTEGER,  INTENT(IN)                     :: kdeep

REAL(dp), INTENT(IN)                     :: dzr(klon)
REAL(dp), INTENT(IN)                     :: dzs(klon)
REAL(dp), INTENT(IN)                     :: vpor(klon)
REAL(dp), INTENT(IN)                     :: vfc(klon)
REAL(dp), INTENT(IN)                     :: vpwp(klon)
LOGICAL,  INTENT(IN)                     :: loland(klon)
REAL(dp), INTENT(IN)                     :: zdel(kdeep)

REAL(dp), INTENT(OUT)                    :: dzrsi(klon, kdeep)
REAL(dp), INTENT(OUT)                    :: dzsi(klon, kdeep)
REAL(dp), INTENT(OUT)                    :: zwsat(klon,kdeep)
REAL(dp), INTENT(OUT)                    :: zwsfc(klon,kdeep)
REAL(dp), INTENT(OUT)                    :: zwpwp(klon,kdeep)


REAL(dp), PARAMETER :: zeps = 1.e-15_dp

LOGICAL  :: loflag1(klon), loflag2(klon)

INTEGER  :: i, jl
REAL(dp) :: zdels(kdeep), zdelsm1(kdeep)

! ************ Field initialization -
  DO  jl=1,klon
    loflag1(jl) = .true.
    loflag2(jl) = .true.
  END DO

  dzsi (:,:) = 0._dp
  dzrsi(:,:) = 0._dp
!  CALL setra(dzsi,klon*kdeep,0._dp)
!  CALL setra(dzrsi,klon*kdeep,0._dp)

!CDIR NOVECTOR
  DO  i=2,kdeep
    zdels(i)   = sum(zdel(1:i))
    zdelsm1(i) = sum(zdel(1:i-1))
  END DO

  zdels(1)   = zdel(1)
  zdelsm1(1) = 0._dp

  DO  i=1,kdeep
    DO  jl=kidia, kfdia
      IF (loland(jl).AND.loflag1(jl)) THEN
        IF (dzr(jl) > zdels(i)) THEN
          dzrsi(jl,i) = zdel(i)
        ELSE
          dzrsi(jl,i) = dzr(jl) - zdelsm1(i)
          loflag1(jl) = .false.
        END IF
      END IF
    END DO
  END DO

! *** Saturatable depth
  DO  i=1,kdeep
    DO  jl=kidia, kfdia
      IF (loland(jl) .AND. loflag2(jl)) THEN
        IF (dzs(jl) > zdels(i)) THEN
          dzsi(jl,i) = zdel(i)
        ELSE
          dzsi(jl,i)  = dzs(jl) - zdelsm1(i)
          loflag2(jl) = .false.
        END IF
      END IF
    END DO
  END DO

  DO  i=1,kdeep
    DO  jl=1,klon
      zwsat(jl,i) = vpor(jl) * dzsi(jl,i)
      zwsfc(jl,i) = vfc(jl)  * dzsi(jl,i)
      zwpwp(jl,i) = vpwp(jl) * dzsi(jl,i)
    END DO
  END DO

END SUBROUTINE soildef


SUBROUTINE soilchange (kidia, kfdia, klon, kdeep,  isch ,  &
    dt, zdiags, loland, &
!         *** - Soil moisture layers,  
    wsi   , wsim1m, dzr, dzrsi  , dzsi, zwsat, zwsfc, &
!         ***   Vertical water fluxes  -- evap. fluxes over land
    finfil, drain, aetrans, aebsoil,                  &
!         ***   Reduction of evaporation if necessary 
    redevap)

!     Intermediate routine for changing 5 soil moisture layers due to
!     evaporation fluxes, infiltration and drainage that are calculated
!     with the old bucket scheme. The soil moisture changes are calculated
!     for one time step, not two as previously done for ws using the
!     3 time step scheme

!     *** VS. 1.0 -- Juli 2006 -- Stefan Hagemann

!     *** VS 1.1 August 2006 -- Stefan Hagemann
!     *** Security check moved to ISCH =1 loop

!     *** VS 1.2 -- January 2007 -- Stefan Hagemann
!     *** IWORK renamed to ISCH .

!     *** VS 1.3 -- September 2007 -- Stefan Hagemann
!     *** Small error in transpiration loop eliminated that occured if transpiration
!     ***   was taken from a layer where WSI < wilting point. The error lead to the
!     ***   the result that the WSI of the corresponding layer was set to the wilting point.
!     ***   The error only causes problems in dry region gridboxes with small WSMX that are
!     ***   close to the bare soil evaporation water depth of 10 cm, and probably
!     ***   only if small time steps are used, e.g. at 0.25 degree.
!     ***   resolution.
!     *** Introduction of zeps as small increment to be used instead of zero in some
!     ***   comparison-if statements. small flux: 0.01 mm/day = 1.16e-10 m/s
!     ***   ZEPS = 1e-15 << against a small flux.
!     *** If WSI < WWILT for one layer, transpiration amount will be reduced
!     ***   instead of re-distribution to other layers.
!     ***   Thus, AETRANS is changed in SOILCHANGE, but also EVAP and EVAPL
!     ***   must be changed.
!     *** Bare soil evaporation will only be taken from the most upper layer.
!     ***   overshooting amounts will be disregarded and AEBSOIL will be reduced.
!     ***   thus, AEBSOIL is changed in SOILCHANGE, but also EVAP and EVAPL
!     ***   must be changed.
!     ***  --> New array submitted to routine: REDEVAP
!     ***
!     *** VS 1.4 -- March 2008 -- Stefan Hagemann
!     *** Reduction of evaporation is taken out again as there was no feedback
!     *** (technical) to the atmosphere. An implementation of this feedback would
!     *** require more complex changes to VDIFF and comprehensive testing. In
!     *** my opinion a reduction of transpiration in VDIFF would be desirable.
!     *** Field REDEVAP is kept for diagnostic purposes.
!     ***
!     *** VS 1.5 - October 2009 Stefan Hagemann
!     ***   Small error in soil moisture change for REDEVAP not zero was corrected.
!     ***   This has led to errors in dry areas.
!     ***   Infiltration fills the layers from above until field capacity. Initially
!     ***   it was filled up until saturation but this caused (likely) unreal(dp)istic
!     ***   high layer overflow values after diffusion(and thus drainage) in dry areas
!     ***   as congo.
!     ***
!     *** VS 1.6 - February 2010 Stefan Hagemann
!     *** Adding infiltration security check if soil for some reason
!     *** cannot take it all

!     ******** Input fields
!     ***AEBSOIL = Bare soil evaporation in [m], accumulated over TWODT/2
!     ***AETRANS = Transpiration in [m], accumulated over TWODT/2
!     ***  DRAIN = Drainage in [m], accumulated over TWODT
!     ***          --> Conversion to TWODT/2 needed
!     *** FINFIL = Infiltration in [m/s] --> multiplication with DT required.

!     ***    DZR = Rooting depth in [m]
!     ***  DZRSI = Rooted depth per soil layer (until rooting depth DZR)
!     ***   DZSI = Depth that can be saturated with water per soil layer
!     ***          (until bedrock DZS)
!     ***  ZWSAT = Maximum storage capacity of soil layer I (porosity) [m]
!     ***          Important for the uptake of infiltration in the soil layers.
!     ***          Defined to the depth of the soil, that is, where the bedrock begins.
!     ***  ZWSFC = Field capacity of soil I (analog to ZWSAT) [m]

!     ***    DT = time step, usually [s]
!     *** ISCH  = Switch which values will change the soil moisture
!     ***         0 = No change: set WSI = WSIM1M
!     ***         1 = change by evaporation fluxes
!     ***         2 = change by infiltration
!     ***         3 = change by infiltration and drainage
!     ---------------------

INTEGER,  INTENT(IN)                     :: kidia
INTEGER,  INTENT(IN)                     :: kfdia
INTEGER,  INTENT(IN)                     :: klon
INTEGER,  INTENT(IN)                     :: kdeep
INTEGER,  INTENT(IN)                     :: isch
REAL(dp), INTENT(IN)                     :: dt
REAL(dp), INTENT(IN)                     :: zdiags
LOGICAL,  INTENT(IN)                     :: loland(klon) ! -- Land and ice mask      --
REAL(dp), INTENT(INOUT)                  :: wsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: wsim1m(klon, kdeep)
REAL(dp), INTENT(IN)                     :: dzr(klon)
REAL(dp), INTENT(IN)                     :: dzrsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: dzsi(klon, kdeep)
REAL(dp), INTENT(IN)                     :: zwsat(klon,kdeep)
REAL(dp), INTENT(IN)                     :: zwsfc(klon,kdeep)
REAL(dp), INTENT(IN)                     :: finfil(klon)
REAL(dp), INTENT(INOUT)                  :: drain(klon)
REAL(dp), INTENT(IN)                     :: aetrans(klon)
REAL(dp), INTENT(IN)                     :: aebsoil(klon)
REAL(dp), INTENT(INOUT)                  :: redevap(klon)

REAL(dp), PARAMETER :: zeps = 1.e-15_dp

INTEGER :: jl, jk

REAL(dp) :: zwdmin, zwilt
REAL(dp) :: zrebsoil(klon), zretrans(klon)

REAL(dp) :: zdum, zdum2, zrel
REAL(dp) :: zover(klon)

! *** Constants 
  zwdmin = 0.05_dp
  zwilt  = 0.35_dp ! Read WAVA later on

! ******* No change = set WSI to WSIM1M
  IF (isch == 0) THEN
    wsi(kidia:kfdia,1:kdeep) = wsim1m(kidia:kfdia,1:kdeep)
    RETURN
  END IF

! ******* Change by evaporation fluxes
  IF (isch == 1) THEN
    zrebsoil(kidia:kfdia) = 0._dp
    zretrans(kidia:kfdia) = 0._dp
    zover(kidia:kfdia) = aebsoil(kidia:kfdia)

!   *** Bare soil evaporation loop
    DO  jl=kidia, kfdia
      IF (loland(jl)) THEN
!        *** TAU is packed in SURF in precipitation
!        ***     --> Is contained in the infiltration!
        IF (zover(jl) < -zeps) THEN
          wsi(jl,1) = wsi(jl,1) + zover(jl)
          IF (wsi(jl,1) < 0._dp) THEN
            zover(jl) = wsi(jl,1)
            wsi(jl,1) = 0._dp
          ELSE
            zover(jl)=0._dp
          END IF
        END IF
      END IF
    END DO
    
!       *** Reduction of AEBSOIL necessary? --> diagnostic
    DO jl=kidia, kfdia
      IF (loland(jl)) THEN
        IF (zover(jl) < -zeps) THEN
          zrebsoil(jl) = zover(jl)
!!!            AEBSOIL(JL) = MIN(AEBSOIL(JL) - ZOVER(JL), 0.)
        END IF
      END IF
    END DO

!   *** transpiration loop
    zover(kidia:kfdia) = 0._dp
    DO  jk=1, kdeep
      DO  jl=kidia, kfdia
        IF (loland(jl)) THEN
          IF (aetrans(jl) < -zeps) THEN
            
!           *** Distribute transpiration over rooting depth
!           *** for this, calculate rooted water depth 
            IF (dzr(jl) <= zeps) THEN
              
!             *** Special case: sealed grid cell
              IF (jk == 1) zretrans(jl) = zretrans(jl) + aetrans(jl)
              
            ELSE IF (dzrsi(jl,jk) > 0) THEN
!             *** Note that transpiration should only access water until the
!             *** rooting depth --> pay regard where the relative term of the
!             *** rooted depth within the layer has to be applied
              zrel = MIN(dzrsi(jl,jk)/dzsi(jl,jk), 1._dp)
              
              zdum2 = zwsfc(jl,jk)*zwilt * zrel
              IF (wsi(jl,jk)*zrel >= zdum2) THEN
                zdum = wsi(jl,jk) * zrel + aetrans(jl) * dzrsi(jl,jk)/dzr(jl)
                IF (zdum < zdum2) THEN
                  zover(jl) = zdum - zdum2
                  wsi(jl,jk) = zdum2 + wsi(jl,jk) * MAX(1._dp-zrel, 0._dp)
                ELSE
                  wsi(jl,jk) = zdum + wsi(jl,jk) * MAX(1._dp-zrel, 0._dp)
                  zover(jl)=0._dp
                END IF
              ELSE
                zover(jl) = aetrans(jl) * dzrsi(jl,jk)/dzr(jl)
              END IF
              zretrans(jl) = zretrans(jl) + zover(jl)
            END IF
          END IF
        END IF
      END DO
    END DO

!   *** Adding ZRETRANS and ZREBSOIL to REDEVAP
    DO jl=kidia, kfdia
      IF (loland(jl)) THEN
!CC          IF (zretrans(JL).LT.-zeps) THEN
!CC            aetrans(JL) = MIN(aetrans(JL) - redevap(JL), 0.)
!CC          ENDIF
        redevap(jl) = redevap(jl) + zrebsoil(jl) + zretrans(jl)
      END IF
    END DO
    
!   *** If redevap not zero, soil moisture must be changed
!   *** to close the land surface water balance.
!   *** Reduction will be relatively distributed to the wetness
!   *** of each layer.
    DO  jl=kidia, kfdia
      IF (loland(jl)) THEN
        IF (redevap(jl) < -zeps) THEN
          zdum2 = SUM(wsi(jl,:))
          IF (zdum2 > 0._dp) THEN
            DO jk=1, kdeep
              zdum = wsi(jl,jk) + redevap(jl) * wsi(jl,jk)/zdum2
              wsi(jl,jk) = MAX(zdum, 0._dp)
            END DO
          END IF
        END IF
      END IF
    END DO
  END IF

! ******* Change by infiltration
  IF (isch == 2 .OR. isch == 3) THEN
    zover(kidia:kfdia) = finfil(kidia:kfdia)*dt
    
    DO  jk = 1, kdeep
      DO  jl=kidia, kfdia
        IF (loland(jl)) THEN
          
!         *** infiltration-loop for 5 layer routine
          wsi(jl,jk) = wsi(jl,jk) + zover(jl)
          IF (wsi(jl,jk) > zwsfc(jl,jk)) THEN
            zover(jl)  = wsi(jl,jk) - zwsfc(jl,jk)
            wsi(jl,jk) = zwsfc(jl,jk)
          ELSE
            zover(jl) = 0._dp
          END IF
        END IF
      END DO
    END DO
  
!   *** If infiltration cannot fully enter the soil
    DO jl=kidia, kfdia
      IF (loland(jl)) THEN
        IF (zover(jl) > 0) drain(jl) = drain(jl)+zover(jl)
      END IF
    END DO
  END IF

! ******* Change by drainage (unit m in SURF), used only for testing
  IF (isch == 3) THEN
    DO  jl=kidia, kfdia
      IF (loland(jl).AND.drain(jl) > 0._dp) THEN
        zover(jl) = -drain(jl) * zdiags
        DO jk=kdeep,1,-1
!         *** Start saturated(?) depth
          zdum = zwsfc(jl,jk)*zwdmin
          IF (zwsat(jl,jk) > 0 .AND. wsi(jl,jk)-zdum > 0) THEN
            wsi(jl,jk) =  wsi(jl,jk) + zover(jl)
            IF (wsi(jl,jk)-zdum < 0) THEN
              zover(jl) = wsi(jl,jk)-zdum
              wsi(jl,jk) = zdum
            ELSE
              zover(jl) =0._dp
              EXIT      ! Exit soil layer loop
            END IF
          END IF
        END DO
      END IF
    END DO
  END IF

! ****** Avoid negative values - shouldn't happen
  wsi(kidia:kfdia,1:kdeep) = MAX( wsi(kidia:kfdia,1:kdeep), 0._dp)

END SUBROUTINE soilchange

SUBROUTINE setra ( rar, idim, rwert )

!**** SETRA    -   Up: Set value for a REAL(dp)-field
!**
!**   Call     :   CALL SETRA ( RAR, IDIM, RWERT )
!**   Purpose  :   Fill REAL(dp)-field with given value
!**   Input parameters: 
!**                RAR  -  REAL(dp)-field to fill
!**                IDIM -  Number of elements to fill
!**                RWERT-  REAL(dp)-value to fill RAR with
!**
!**   Author   : Volker Wohlgemuth,F3
!
!
INTEGER, INTENT(IN)                      :: idim
REAL(dp), INTENT(IN)                     :: rwert
REAL(dp), INTENT(OUT)                    :: rar(idim)

INTEGER ii

  DO  ii = 1,idim
    rar(ii) = rwert
  END DO

END SUBROUTINE setra

END MODULE mo_jsbach_hydrology
