!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model, according to 
!! Semtner's Zero-Layer Model (1976)
!!
!! @author Achim Randelhoff
!! 
!! @par Revision History
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
MODULE mo_sea_ice_zerolayer

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, ki, ks, Tf,  &
    &                               mu,mus,ci, alf, I_0, alv, clw,    &
    &                               cpd, zemiss_def, rd, stbo, tmelt   
  USE mo_ocean_nml,           ONLY: no_tracer
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, hci_layer, use_constant_tfreez
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state 
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean
  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=12)           :: str_module    = 'SeaIceZeroLy'  ! Output of module for 1 line debug

  PUBLIC :: set_ice_temp_zerolayer
  PUBLIC :: set_ice_temp_zero_nogradients
  PUBLIC :: ice_growth_zerolayer

CONTAINS

  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! set_ice_temp_zerolayer:: calculate new ice + snow temperatures based on:
  !!    Semtner, Albert J., 1976: A Model for the Thermodynamic Growth of Sea
  !!    Ice in Numerical Investigations of Climate. J. Phys. Oceanogr., 6,
  !!    379--389.  doi:
  !!    http://dx.doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2
  !!   (Appendix)
  !!
  !! This function changes:
  !! ice % Tsurf    the new surface temperature   for each ice category     [deg C]
  !! ice % Qbot     Heat flux available for freezing/melting at ice bottom  [W/m^2]
  !! ice % Qtop     Heat flux available for melting at ice surface          [W/m^2]
  !!
  !!           all "dtime" in this function are atmospheric time step
  !! @par Revision History
  !! Initial release by Achim Randelhoff

  SUBROUTINE set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, i_therm_model, pdtime, &
            &   Tsurf,          & 
            &   hi,             & 
            &   hs,             & 
            &   Qtop,           & 
            &   Qbot,           & 
            &   SWnet,          & 
            &   nonsolar,       & 
            &   dnonsolardT,    &
            &   Tfw,            &
            &   doy)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice, i_therm_model
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)

    INTEGER,OPTIONAL,INTENT(IN)    :: doy

    ! Local variables
    REAL(wp) ::        &
      & k_effective ,  &  ! total heat conductivity of ice/snow
      & deltaT      ,  &  ! temperature increment 
      & F_A         ,  &  ! atmospheric net flux, positive=upward
      & F_S         ,  &  ! conductive flux, positive=upward
      & deltaTdenominator     ! prefactor of deltaT in sfc. flux
                              ! balance
    
!    REAL(wp) :: one_minus_I_0 ! 1.0 - I_0 for use with SWin

    INTEGER :: k, jc ! loop indices

    

    ! initialization of output variables
    Qbot(:,:) = 0._wp
    Qtop(:,:) = 0._wp

    ! --- initialization
!    one_minus_I_0 = 1.0_wp

    DO k=1,kice
      DO jc = i_startidx_c,i_endidx_c
        IF (hi(jc,k) > 0._wp) THEN
          
          ! --- total heat conductivity for the ice-snow system
          ! #slo/Josiane# 2014-11: rho_ref/rhos missing - TODO: check
          k_effective = ki*ks/(ks*hi(jc,k) + ki*hs(jc,k))


! --- calculate (1-I_0)
!          IF (hs(jc,k) > 0.0_wp ) THEN
!            one_minus_I_0=1.0_wp
!          ELSE
!            one_minus_I_0=1.0_wp-I_0
!          END IF
                  
          ! --- F_A, F_S : pos=upward flux

          ! F_A: flux ice-atmosphere
          IF (i_therm_model == 1) THEN
            F_A = - nonsolar(jc,k) - SWnet(jc,k) !* one_minus_I_0
          ELSE IF (i_therm_model ==3) THEN
            ! #achim: first draft: hard-coding simpler form of
            ! atmospheric fluxes (from Dirk's thesis, p.193)
            !
            ! atm. flx = LWout - (LWin, sens, lat) - SWin
            F_A =   zemiss_def * StBo * (Tsurf(jc,k) +tmelt)**4   &
              &            - ( 118.0_wp * EXP(-0.5_wp*((doy-206)/53.1_wp)**2) + 179.0_wp ) &
              &            - 314.0_wp * EXP(-0.5_wp*((doy-164)/47.9_wp)**2)  & !SW, NO 1-I_0 factor
              &                * ( 0.431_wp / (1.0_wp+((doy-207)/44.5_wp)**2) - 0.086_wp)!1-albedo
          END IF

          ! F_S conductive heat flux through the ice
          F_S = k_effective * (Tfw(jc) - Tsurf(jc,k))


          IF (i_therm_model == 1 ) THEN
          ! We add constant heat capacity to deltaTdenominator to stabilize the atmosphere
            deltaTdenominator = k_effective  - dnonsolardT(jc,k) + rhoi*hci_layer*ci/pdtime
          ELSE IF (i_therm_model == 3) THEN
            ! dLWdT is missing!
            deltaTdenominator = k_effective + 4.0_wp*zemiss_def*StBo*(Tsurf(jc,k)+tmelt)**3
          END IF

          ! --- temperature increment
          deltaT = (F_S - F_A) / deltaTdenominator


          ! --- ice temperatures over 0 deg C impossible:
          IF (Tsurf(jc,k) + deltaT > 0.0_wp) THEN  ! if new temperature would be over 0 deg C
            deltaT = -Tsurf(jc,k) 
             Tsurf(jc,k) = 0.0_wp
            
            ! pos. flux means into uppermost ice layer
            Qtop(jc,k) = - F_A + F_S - deltaT * deltaTdenominator

            ! pos. flux means into lowest ice layer
            ! correction r20136 by Josiane/Dirk
            Qbot(jc,k) = - F_S + deltaT * k_effective
            !Qbot(jc,k) = - F_S -  deltaT * k_effective
            ! NB: flux from ocean to ice still missing, this is done in
            !                       ice_growth_zerolayer


          ELSE   ! if new temperature is less than 0 deg C, then we can achieve F_A=F_S just as we wanted
            ! new surface temperature
            Tsurf(jc,k) = Tsurf(jc,k) + deltaT
            
            ! surface flux balanced:
            Qtop(jc,k) = 0.0_wp
            
            ! pos. flux means into lowest ice layer
            Qbot(jc,k) =  k_effective * (Tsurf(jc,k) - Tfw(jc))
            !!! ice%Qbot = - new F_S
            ! NB: flux from ocean to ice still missing, this is done in
            !                       ice_growth_zerolayer
            
          END IF

        ELSE
          Tsurf(jc,k) = Tfw(jc)

   !  check whether correct?
   !      Qtop(jc,k) = 0.0_wp
   !      Qbot(jc,k) = 0.0_wp
        END IF

      END DO
    END DO

! ----------------------------------------

  END SUBROUTINE set_ice_temp_zerolayer


  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! set_ice_temp_zerolayer:: calculate new ice + snow temperatures based on:
  !!    Semtner, Albert J., 1976: A Model for the Thermodynamic Growth of Sea
  !!    Ice in Numerical Investigations of Climate. J. Phys. Oceanogr., 6,
  !!    379--389.  doi:
  !!    http://dx.doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2
  !!   (Appendix)
  !!
  !! This function changes:
  !! ice % Tsurf    the new surface temperature   for each ice category     [deg C]
  !! ice % Qbot     Heat flux available for freezing/melting at ice bottom  [W/m^2]
  !! ice % Qtop     Heat flux available for melting at ice surface          [W/m^2]
  !!
  !!           all "dtime" in this function are atmospheric time step
  !! @par Revision History
  !! Initial release by Achim Randelhoff
  !!
  !! New routine set_ice_temp_zero_nogradients for test of simplified scheme
  !! Initial release by Stephan Lorenz, MPI,  2015-07

  SUBROUTINE set_ice_temp_zero_nogradients(i_startidx_c, i_endidx_c, nbdim, kice, i_therm_model, pdtime, &
            &   Tsurf,          & 
            &   hi,             & 
            &   hs,             & 
            &   Qtop,           & 
            &   Qbot,           & 
            &   SWnet,          & 
            &   nonsolar,       & 
            &   dnonsolardT,    &
            &   Tfw,            &
            &   doy)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice, i_therm_model
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)

    INTEGER,OPTIONAL,INTENT(IN)    :: doy

    ! Local variables
    REAL(wp) ::        &
      & k_effective ,  &  ! total heat conductivity of ice/snow
      & deltaT      ,  &  ! temperature increment 
      & F_A         ,  &  ! atmospheric net flux, positive=upward
      & F_S         ,  &  ! conductive flux, positive=upward
      & deltaTdenominator     ! prefactor of deltaT in sfc. flux
                              ! balance
    REAL(wp) :: c_icelayer
    
!    REAL(wp) :: one_minus_I_0 ! 1.0 - I_0 for use with SWin

    INTEGER :: k, jc ! loop indices

    

    ! initialization of output variables
    Qbot(:,:) = 0._wp
    Qtop(:,:) = 0._wp

    ! declaration of constant:
    c_icelayer = rhoi*hci_layer*ci/pdtime

    ! --- initialization
!    one_minus_I_0 = 1.0_wp

    DO k=1,kice
      DO jc = i_startidx_c,i_endidx_c
        IF (hi(jc,k) > 0._wp) THEN
          
          ! --- total heat conductivity for the ice-snow system
          ! #slo/Josiane# 2014-11: rho_ref/rhos missing - TODO: check
          k_effective = ki*ks/(ks*hi(jc,k) + ki*hs(jc,k))


! --- calculate (1-I_0)
!          IF (hs(jc,k) > 0.0_wp ) THEN
!            one_minus_I_0=1.0_wp
!          ELSE
!            one_minus_I_0=1.0_wp-I_0
!          END IF
                  
          ! --- F_A, F_S : pos=upward flux

          ! F_A: flux ice-atmosphere
          IF (i_therm_model == 1) THEN
            F_A = - nonsolar(jc,k) - SWnet(jc,k) !* one_minus_I_0
          ELSE IF (i_therm_model ==3) THEN
            ! #achim: first draft: hard-coding simpler form of
            ! atmospheric fluxes (from Dirk's thesis, p.193)
            !
            ! atm. flx = LWout - (LWin, sens, lat) - SWin
            F_A =   zemiss_def * StBo * (Tsurf(jc,k) +tmelt)**4   &
              &            - ( 118.0_wp * EXP(-0.5_wp*((doy-206)/53.1_wp)**2) + 179.0_wp ) &
              &            - 314.0_wp * EXP(-0.5_wp*((doy-164)/47.9_wp)**2)  & !SW, NO 1-I_0 factor
              &                * ( 0.431_wp / (1.0_wp+((doy-207)/44.5_wp)**2) - 0.086_wp)!1-albedo
          END IF

          ! F_S conductive heat flux through the ice
          F_S = k_effective * (Tfw(jc) - Tsurf(jc,k))


          IF (i_therm_model == 1 ) THEN
          ! We add constant heat capacity to deltaTdenominator to stabilize the atmosphere
            deltaTdenominator = k_effective  - dnonsolardT(jc,k) + rhoi*hci_layer*ci/pdtime
          ELSE IF (i_therm_model == 3) THEN
            ! dLWdT is missing!
            deltaTdenominator = k_effective + 4.0_wp*zemiss_def*StBo*(Tsurf(jc,k)+tmelt)**3
          END IF

          ! --- temperature increment
        ! deltaT = (F_S - F_A) / deltaTdenominator
          deltaT = (k_effective * Tfw(jc) - F_A + c_icelayer * Tsurf(jc,k) ) / &
            &      (c_icelayer + k_effective) - Tsurf(jc,k)

          ! --- ice temperatures over 0 deg C impossible:
          IF (Tsurf(jc,k) + deltaT > 0.0_wp) THEN  ! if new temperature would be over 0 deg C
            deltaT = -Tsurf(jc,k) 
             Tsurf(jc,k) = 0.0_wp
            
            ! pos. flux means into uppermost ice layer
          ! Qtop(jc,k) = -F_A + F_S - deltaT * deltaTdenominator
            Qtop(jc,k) = -F_A + F_S - deltaT * c_icelayer

            ! pos. flux means into lowest ice layer
            ! correction r20136 by Josiane/Dirk
            Qbot(jc,k) = - F_S + deltaT * k_effective
            !Qbot(jc,k) = - F_S -  deltaT * k_effective
            ! NB: flux from ocean to ice still missing, this is done in
            !                       ice_growth_zerolayer


          ELSE   ! if new temperature is less than 0 deg C, then we can achieve F_A=F_S just as we wanted
            ! new surface temperature
            Tsurf(jc,k) = Tsurf(jc,k) + deltaT
            
            ! surface flux balanced:
            Qtop(jc,k) = 0.0_wp
            
            ! pos. flux means into lowest ice layer
            Qbot(jc,k) =  k_effective * (Tsurf(jc,k) - Tfw(jc))
            !!! ice%Qbot = - new F_S
            ! NB: flux from ocean to ice still missing, this is done in
            !                       ice_growth_zerolayer
            
          END IF

        ELSE
          Tsurf(jc,k) = Tfw(jc)

   !  check whether correct?
   !      Qtop(jc,k) = 0.0_wp
   !      Qbot(jc,k) = 0.0_wp
        END IF

      END DO
    END DO

! ----------------------------------------

  END SUBROUTINE set_ice_temp_zero_nogradients


 !
 ! The counterpart to the  ice_growth subroutine
 !
 !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! ice_growth_zerolayer - change ice and snow thickness (Semtner 1976, Appendix)
  !! This function changes:
  !! --- currently not fully --- ice % hs       new snow thickness for each ice category                [m]
  !! ice % hi       new ice  thickness for each ice category                [m]
  !! --- not currently --- ice % evapwi   amount of evaporated water from the mixed layer
  !!                in previously ice covered areas if all ice is gone      [kg/m^3]
  !! ice % heatOceI to contain the energy that is available to the mixed layer
  !!                in previously ice covered areas if all ice is gone      [J]
  !!
  !! @par Revision History
  !! Initial release by Achim Randelhoff
  !!
 
 SUBROUTINE ice_growth_zerolayer(p_patch, p_os, ice, rpreci)
   TYPE(t_patch),             INTENT(IN), TARGET    :: p_patch 
   TYPE(t_hydro_ocean_state), INTENT(IN)            :: p_os
   TYPE (t_sea_ice),          INTENT(INOUT)         :: ice
   REAL(wp),                  INTENT(IN)            :: rpreci(:,:) 
                                   ! water equiv. solid precipitation rate [m/s] DIMENSION (ie,je)

   !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch%alloc_cell_blocks) ::         &
      & Tfw,         & ! Ocean freezing temperature [C]
      & Q_surplus   ! energy surplus during ice growth
    
    REAL(wp) ::      &
      & below_water, & ! Thickness of snow layer below water line           [m]
      & draft          ! depth of ice-ocean interface below sea level       [m]

    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c
    
    all_cells            => p_patch%cells%all
    Tfw(:,:,:)           =  0.0_wp
    ice%zHeatOceI(:,:,:) =  0.0_wp
    Q_surplus(:,:,:)     =  0.0_wp

    
    ! Save ice thickness at previous time step for calculation of heat and salt
    ! flux into ocean in subroutine upper_ocean_TS
    ice % hiold (:,:,:) = ice%hi(:,:,:)
    ice % hsold (:,:,:) = ice%hs(:,:,:)

    ! freezing temperature of uppermost sea water
    IF ( no_tracer < 2 .OR. use_constant_tfreez ) THEN
      Tfw(:,:,:) = Tf
    ELSE
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu * p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ENDIF

!---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('GrowZero: heatOceI bef.grow' , ice%heatOceI   , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qtop bef. growth'  , ice%Qtop       , str_module, 5, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qbot bef. growth'  , ice%Qbot       , str_module, 5, in_subset=p_patch%cells%owned)
!---------------------------------------------------------------------
    
    IF (i_ice_therm /= 3 ) THEN
      ! Heat flux from ocean into ice
      CALL oce_ice_heatflx (p_patch, p_os,ice,Tfw,ice%zHeatOceI)
!!$    ELSE IF ( i_ice_therm == 3) THEN
      ! for i_ice_therm == 3, no ocean-ice heatflx is included!
    END IF

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc, draft, below_water) SCHEDULE(dynamic)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%hi(jc,k,jb) > 0._wp) THEN

            ! Add oceanic heat flux to energy available at the bottom of the ice.
            ice%Qbot(jc,k,jb) = ice%Qbot(jc,k,jb) + ice%zHeatOceI(jc,k,jb)

            ! Add snowfall to snow depth
            ! #slo# 2015-01: bugfix: rpreci is rate of snowfall over ice covered area
            ice%hs(jc,k,jb) = ice%hs(jc,k,jb) + rpreci(jc,jb)*dtime*rho_ref/rhos
            ! #slo# 2015-01: bugfix: rpreci is over whole grid-area
            !ice%hs(jc,k,jb) = ice%hs(jc,k,jb) + rpreci(jc,jb)*ice%conc(jc,k,jb)*dtime*rho_ref/rhos
      
            ! for energy flux surplus
            IF ( ice%Qtop(jc,k,jb) > 0.0_wp ) THEN 
              IF  ( ice%hs(jc,k,jb) > 0.0_wp )  THEN ! melt snow where there's snow
                
                ice%hs (jc,k,jb) =  ice%hs(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhos) 
                ! put remaining heat, if any, into melting ice below
                IF (ice%hs(jc,k,jb) < 0.0_wp) THEN
                  ice%hi(jc,k,jb) = ice%hi(jc,k,jb) + ice%hs(jc,k,jb) * (rhos/rhoi) ! snow thickness loss in ice equivalents
                  ice%hs(jc,k,jb) = 0.0_wp
                ENDIF
                
              ELSE   ! where there's no snow
                ice%hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhoi) 
              ENDIF
            ENDIF
            
            ! bottom melt/freeze
            ice%hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qbot(jc,k,jb) * dtime / (alf*rhoi )
            
            ! heat to remove from water
            !  - heatOceI - positive into ocean - positive=downward i.e. same sign convention as HeatFlux_Total into ocean
            !  - zHeatOceI - positive into ice, i.e. positive=upward - melting energy coming from below, from the ocean
            ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - ice%zHeatOceI(jc,k,jb)
            
            ! hi<0 if melting energy (Qbot+zHeatOceI) is larger than needed to melt all ice and snow, see above
            IF (ice%hi (jc,k,jb) <= 0.0_wp) THEN

              ! remove surplus energy of ice thickness from water
              !  - hi<0, if all ice and snow is already melted
              !  - calculate surplus of heatOceI>0 available for heating of ocean after complete melting
              ! #slo# 2014-11: 3. bugfix: sign error in hi for heatOceI
              !ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) + ice%hi(jc,k,jb)*alf*rhoi/dtime
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - ice%hi(jc,k,jb)*alf*rhoi/dtime

              ! remove latent heat of snow from water
              ! #slo# 2014-11: if there is snow on top of melted ice, hs>0, then heatOceI is reduced by latent heat of snow
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - ice%hs(jc,k,jb)*alf*rhos/dtime
              ! #slo# 2014-11: Attention: if heatOceI is not enough to melt whole snow,
              !                then energy budget is not closed! TODO: check and correct (later, little energy)
              ! IF ( ice%heatOceI(jc,k,jb) >0 ) THEN
              !  - snow is set to rest of ice, since no snow without water is possible
              !   ice%hi(jc,k,jb) = ice%heatOceI(jc,k,jb) * dtime / (alf * rhoi )
              !   ice%heatOceI(jc,k,jb) = 0.0_wp
              ! ELSE  ! melting energy was enough to melt all snow and ice
              
              ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
              !ice%conc (jc,k,jb) = 0.0_wp  !  do not change concentration here, but in ice_conc_change only
              ice%hi   (jc,k,jb) = 0.0_wp
              ice%hs   (jc,k,jb) = 0.0_wp

              ! ENDIF ! ( ice%heatOceI(jc,k,jb) >0 )

            ENDIF

            ! #slo# 2014-11: 2. bugfix: calculation moved down to below recalculaton of hi
            ! #slo# 2015-01: could we update ice%draft here or not?
            draft           = ( rhoi*ice%hi(jc,k,jb) + rhos*ice%hs(jc,k,jb) )/rho_ref
            below_water     = draft - ice%hi(jc,k,jb)  !  thickness to be converted to ice
            
            ! snow -> ice conversion for snow below waterlevel
            ! Currently not quite physical: Snow is pushed together to form new ice, hence snow thickness
            ! decreases more than ice thickness by rhoi/rhos ( analogue to old growth.f90 sea-ice model )
            ! Salt content of snow ice is equal to that of normal ice, salt is removed from the ocean
            ! Temperature of new upper ice is calculated as described in the paragraph below 
            ! Eq. 36
            IF ( below_water > 0.0_wp ) THEN
              ice%snow_to_ice(jc,k,jb) = below_water*rhoi/rhos     ! Thickness of snow that is converted into ice
              ice%hs         (jc,k,jb) = ice%hs(jc,k,jb) - ice%snow_to_ice(jc,k,jb)
              ice%hi         (jc,k,jb) = ice%hi(jc,k,jb) + below_water
            END IF

            IF (ice%hs (jc,k,jb) < 0.0_wp) THEN
               ice % hs(jc,k,jb) = 0.0_wp
               ice % hi(jc,k,jb) = ice%hi(jc,k,jb) + ice%hs(jc,k,jb) * (rhos/rhoi) ! snow thickness loss in ice equivalents
            ENDIF
            
            ! check energy conservation
            ! surplus energy = entering - leaving - latent heat
            !!! what's up with the energy that's put into the ocean?
            ! #slo# 2015-01: snowfall changes energy input - not yet considered
            Q_surplus(jc,k,jb) = &!0.0_wp
              &                   ice%Qbot(jc,k,jb) + ice%Qtop(jc,k,jb) &
              &                   + (ice%hi(jc,k,jb)-ice%hiold(jc,k,jb)) *alf*rhoi/dtime&
              &                   + (ice%hs(jc,k,jb)-ice%hsold(jc,k,jb)) *alf*rhos/dtime

          ELSE  !  hi<=0
            ! #slo# 2014-12: check - heatOceI is set in case of no ice - negative ice possible?
            ice%heatOceI(jc,k,jb) = ice%Qtop(jc,k,jb) + ice%Qbot(jc,k,jb)
            ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
            ice%conc (jc,k,jb) = 0.0_wp
            ice%hi   (jc,k,jb) = 0.0_wp
            ice%hs   (jc,k,jb) = 0.0_wp
          ENDIF

          ! #slo# 2014-12: update zUnderIce better here than in ice_slow?
       !  ice%zUnderIce(:,:) = flat(:,:) + p_os%p_prog(nold(1))%h(:,:) &
       !    &                - (rhos * ice%hs(:,1,:) + rhoi * ice%hi(:,1,:)) * ice%conc(:,1,:) / rho_ref

        END DO
      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

!---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('GrowZero: snow_to_ice', ice%snow_to_ice, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: hi'         , ice%hi         , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: hs'         , ice%hs         , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: zHeatOceI'  , ice%zHeatOceI  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: heatOceI '  , ice%heatOceI   , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Q_surplus'  , Q_surplus      , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qtop'       , ice%Qtop       , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Qbot'       , ice%Qbot       , str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero: Tsurf'      , ice%Tsurf      , str_module, 4, in_subset=p_patch%cells%owned)
!---------------------------------------------------------------------
 
  END SUBROUTINE ice_growth_zerolayer
  
END MODULE mo_sea_ice_zerolayer
