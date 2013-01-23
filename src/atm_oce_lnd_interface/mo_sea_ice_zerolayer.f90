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
!! @par Copyright
!! ! 2002-2007 by DWD and MPI-M
!! ! This software is provided for non-commercial use only.
!! ! See the LICENSE and the WARRANTY conditions.
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
MODULE mo_sea_ice_zerolayer

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs,&
    &                               mu,mus,ci, alf, I_0, alv, albedoW, clw,            &
    &                               cpd, zemiss_def,rd, stbo,tmelt   
  USE mo_ocean_nml,           ONLY: no_tracer, i_sea_ice
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_oce_state,           ONLY: t_hydro_ocean_state 
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean
  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 

  ! # achim: for implementation of simple fluxes
  USE mo_datetime,            ONLY: t_datetime

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=12)           :: str_module    = 'SeaIceZeroLy'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

  PUBLIC :: set_ice_temp_zerolayer
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

  SUBROUTINE set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, SWdim, i_therm_model, &
            &   Tsurf,          & 
            &   hi,             & 
            &   hs,             & 
            &   Qtop,           & 
            &   Qbot,           & 
            &   SWin,           & 
            &   alb,            & 
            &   nonsolar,       & 
            &   dnonsolardT,    &
            &   Tfw,            &
            &   doy)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice, SWdim, i_therm_model
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWin       (nbdim,SWdim)
    REAL(wp),INTENT(IN)    :: alb        (nbdim,kice,SWdim)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    INTEGER, INTENT(IN)    :: doy

    ! Local variables
    REAL(wp) ::        &
      & k_effective ,  &  ! total heat conductivity of ice/snow
      & deltaT      ,  &  ! temperature increment 
      & F_A         ,  &  ! atmospheric net flux, positive=upward
      & F_S         ,  &  ! conductive flux, positive=upward
      & deltaTdenominator     ! prefactor of deltaT in sfc. flux
                              ! balance
    
    REAL(wp) :: one_minus_I_0 ! 1.0 - I_0 for use with SWin

    INTEGER :: k, jc ! loop indices

    

    ! initialization of output variables
    Qbot(:,:) = 0._wp
    Qtop(:,:) = 0._wp

    ! --- initialization
    one_minus_I_0 = 1.0_wp

    DO k=1,kice
      DO jc = i_startidx_c,i_endidx_c
        IF (hi(jc,k) > 0._wp) THEN
          
          ! --- total heat conductivity for the ice-snow system
          k_effective = ki*ks/(ks*hi(jc,k) + ki*hs(jc,k))

! --- calculate (1-I_0)
          IF (hs(jc,k) > 0.0_wp ) THEN
            one_minus_I_0=1.0_wp
          ELSE
            one_minus_I_0=1.0_wp-I_0
          END IF
                  
          ! --- F_A, F_S : pos=upward flux

          ! F_A: flux ice-atmosphere
          IF (i_therm_model == 2) THEN
            F_A = - nonsolar(jc,k) - SUM( (1.0_wp - alb(jc,k,:)) * SWin(jc,:) )* one_minus_I_0
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


          IF (i_therm_model == 2 ) THEN
          ! We add constant heat capacity to deltaTdenominator to stabilize the atmosphere
            deltaTdenominator = k_effective  - dnonsolardT(jc,k) + rhoi*0.05_wp*ci/dtime
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
            Qbot(jc,k) = - F_S -  deltaT * k_effective
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

        END IF
      END DO
    END DO


! ----------------------------------------

  END SUBROUTINE set_ice_temp_zerolayer








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
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch%nblks_c) ::         &
      & zHeatOceI,   & ! Oceanic heat flux                                  [W/m^2]
      & Tfw,         & ! Ocean freezing temperature [C]
      & Q_surplus   ! energy surplus during ice growth
    
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c
    
    all_cells => p_patch%cells%all 
    Tfw(:,:,:) = 0.0_wp
    zHeatOceI(:,:,:) = 0.0_wp
    Q_surplus(:,:,:) = 0.0_wp

    
    ! Save ice thickness at previous time step for calculation of heat and salt
    ! flux into ocean in subroutine upper_ocean_TS
    ice % hiold (:,:,:) = ice%hi(:,:,:)
    ice % hsold (:,:,:) = ice%hs(:,:,:)

    ! freezing temperature of uppermost sea water
    IF ( no_tracer >= 2 ) then
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ELSE
      Tfw(:,:,:) = Tf
    ENDIF
    
    IF (i_sea_ice == 2 ) THEN
      ! Heat flux from ocean into ice
      CALL oce_ice_heatflx (p_os,ice,Tfw,zHeatOceI)
!!$    ELSE IF ( i_sea_ice == 3) THEN
      ! for i_sea_ice == 3, no ocean-ice heatflx is included!
    END IF
    
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%hi(jc,k,jb) > 0._wp) THEN

            ! Add oceanic heat flux to energy available at the bottom of the ice.
            ice%Qbot(jc,k,jb) = ice%Qbot(jc,k,jb) + zHeatOceI(jc,k,jb)
      
            ! for energy flux surplus
            IF ( ice%Qtop(jc,k,jb) > 0.0_wp ) THEN 
              IF  ( ice%hs(jc,k,jb) > 0.0_wp )  THEN ! melt snow where there's snow
                
                ice%hs (jc,k,jb) =  ice%hs(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhos) 
                ! put remaining heat, if any, into melting ice below
                IF (ice%hs(jc,k,jb) < 0.0_wp) THEN
                  ice % hi(jc,k,jb) = ice%hi(jc,k,jb) + ice%hs(jc,k,jb) * (rhos/rhoi) ! snow thickness loss
                  !                                                               in ice equivalents
                  ice % hs(jc,k,jb) = 0.0_wp
                ENDIF
                
              ELSE   ! where there's no snow
                ice % hi(jc,k,jb) = ice%hi(jc,k,jb) -  ice%Qtop(jc,k,jb) * dtime / (alf*rhoi) 
              ENDIF
            ENDIF
            
            ! bottom melt/freeze
            ice % hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qbot(jc,k,jb) * dtime / (alf * rhoi )
            
            
            IF (ice%hi (jc,k,jb) <= 0.0_wp) THEN
              ! #achim: check units -- heatocei in J as opposed to W/m2?
              ! remove surplus energy of ice thickness from water
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) &
                &                     + ice%hi(jc,k,jb)*alf*rhoi/dtime
              ! remove latent heat of snow from water
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) &
                &                     - ice%hs(jc,k,jb)*alf*rhos/dtime
              
              ! 
              ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
              ice%conc (jc,k,jb) = 0.0_wp
              ice%hi   (jc,k,jb) = 0.0_wp
              ice%hs   (jc,k,jb) = 0.0_wp
            ELSE
              ! heat to remove from water
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - zHeatOceI(jc,k,jb)
            ENDIF
            
            ! check energy conservation
            ! surplus energy = entering - leaving - latent heat
            !!! what's up with the energy that's put into the ocean?
            Q_surplus(jc,k,jb) = &!0.0_wp
              &                   ice%Qbot(jc,k,jb) + ice%Qtop(jc,k,jb) &
              &                   + (ice%hi(jc,k,jb)-ice%hiold(jc,k,jb)) *alf*rhoi/dtime&
              &                   + (ice%hs(jc,k,jb)-ice%hsold(jc,k,jb)) *alf*rhos/dtime


          ENDIF
        END DO
      END DO
    END DO
    
!!$    CALL print_maxmin_si(Q_surplus(:,1,:),ice,p_patch,'Q_surplus')
!!$    CALL print_maxmin_si(ice%hi(:,1,:),ice,p_patch,'ice%hi')
!!$    CALL print_cells(Q_surplus(:,1,:),'Q_surplus')
!!$    CALL print_cells(ice%hi(:,1,:),'ice%hi')
!!$    CALL print_cells(ice%Qtop(:,1,:),'ice%Qtop')
!!$    CALL print_cells(ice%Qbot(:,1,:),'ice%Qbot')
!!$    CALL print_cells(ice%hi(:,1,:)-ice%hiold(:,1,:),'new ice')
!!$    CALL print_cells(zHeatOceI(:,1,:),'zHeatOceI')

!!$!#slo# !---------DEBUG DIAGNOSTICS-------------------------------------------
!!$    idt_src=1 !3  ! output print level (1-5, fix)
!!$    CALL dbg_print('GrowZero: Q_surplus'       ,Q_surplus                ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%hi'          ,ice%hi                   ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%Qtop'        ,ice%Qtop                 ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%Qbot'        ,ice%Qbot                 ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%Tsurf'       ,ice%Tsurf                ,str_module,idt_src)
!!$!#slo# !---------------------------------------------------------------------
 
  END SUBROUTINE ice_growth_zerolayer
  
END MODULE mo_sea_ice_zerolayer
