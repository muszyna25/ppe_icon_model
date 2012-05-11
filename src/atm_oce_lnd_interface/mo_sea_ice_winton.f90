MODULE mo_sea_ice_winton
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_exception,           ONLY: finish, message
  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary 
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_math_utilities,      ONLY: t_cartesian_coordinates
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs,&
    &                               mu,mus,ci, alf, I_0, alv, albedoW, clw,            &
    &                               cpd, zemiss_def,rd, stbo,tmelt   
  USE mo_math_constants,      ONLY: rad2deg
  USE mo_ocean_nml,           ONLY: no_tracer, init_oce_prog, iforc_oce, &
    &                               FORCING_FROM_FILE_FLUX, i_sea_ice
  USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base, ocean_var_list
  USE mo_oce_index,           ONLY: print_mxmn, ipl_src
  USE mo_var_list,            ONLY: add_var
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  ! # achim
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean
  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx

  IMPLICIT NONE

  PRIVATE

  ! #achim
  PUBLIC :: ice_growth_winton
  PUBLIC :: set_ice_temp_winton
  ! #

CONTAINS

  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! set_ice_temp_winton:: calculate new ice + snow temperatures according to sec.2a from
  !!           Winton, M., 2000: A Reformulated Three-Layer Sea Ice Model,   
  !!           J. Atmos. Oce. Tech., 17, 525-531. 
  !!
  !!           doi: 10.1175/1520-0426(2000)017<0525:ARTLSI> (put into google)
  !!
  !! This function changes:
  !! ice % Ts       the new surface temperature   for each ice category     [�C]
  !! ice % T1       the new upper ice+snow temp.  for each ice category     [�C]
  !! ice % T2       the new lower ice temperature for each ice category     [�C]
  !! ice % Qbot     Heat flux available for freezing/melting at ice bottom  [W/m�]
  !! ice % Qtop     Heat flux available for melting at ice surface          [W/m�]
  !!
  !!           all "dtime" in this function are atmospheric time step
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE set_ice_temp_winton(ppatch,ice, Tfw, Qatm) 
    TYPE(t_patch),        INTENT(IN)    :: ppatch 
    TYPE(t_sea_ice),      INTENT(INOUT) :: ice
    REAL(wp),             INTENT(IN)    :: Tfw(:,:,:)
    TYPE(t_atmos_fluxes), INTENT(IN)    :: Qatm

    !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, ppatch%nblks_c) ::           &
      & A,           & ! Eq. 7
      & A1,          & ! Eq. 16
      & A1a,         & ! First two terms of Eq. 16 and 19
      & B,           & ! Eq. 8
      & B1,          & ! Eq. 17
      & B1a,         & ! First three terms of Eq. 17 and 20
      & C1,          & ! Eq. 18
      & D,           & ! 1./(6*dT*K2) + rhoi*hi*C for Eq. 16, 17, 19, 20
      & iK1B,        & ! 1./(K1 + B) (used in eq. 16 and 17)
      & K1,          & ! Winton's K 1/2 (eq. 5)
      & K2,          & ! Winton's K 3/2 (eq. 10)
      & SWin3D,      & ! Short-wave radiation field splitted into ice categories
      & Tsurfm         ! Surface melting temperature
    
    REAL(wp) :: idt2 ! 1 / (2*dt)
    
    INTEGER :: i,j,k ! counter for loops
   !-------------------------------------------------------------------------------

    idt2   =  1.0_wp / (2.0_wp*dtime)

    ! Create array of shortwave radiation split up into ice categories
    ! (purely for computational reasons)
    FORALL(i=1:nproma, j=1:ppatch%nblks_c, k=1:ice%kice, ice % isice (i,k,j)) &
      & SWin3d(i,k,j) = Qatm% SWin(i,j)

    ! Calculate new ice temperature wherever there is ice 
    ! lat > 0, sens > 0, LWnet >0 , SWin > 0  for downward flux 
    ! dlatdT, dsensdT, dLWdT >0 for downward flux increasing with increasing Tsurf
    !

    WHERE (ice % isice (:,:,:) )
      B   (:,:,:) = -Qatm% dlatdT(:,:,:) - Qatm% dsensdT(:,:,:) - Qatm% dLWdT(:,:,:)        ! Eq.  8
      A   (:,:,:) = -Qatm% lat(:,:,:) - Qatm% sens(:,:,:) - Qatm% LWnet(:,:,:) -  &
        &           (1.0_wp - ice% alb(:,:,:)) * I_0 *  SWin3d  - ice%Tsurf(:,:,:)* B       ! Eq.  7
      K1  (:,:,:)  =  4.0_wp * ki * ks / (ks * ice%hi(:,:,:) + 4.0_wp * ki * ice%hs(:,:,:)) ! Eq.  5
      K2  (:,:,:)  =  2.0_wp * ki / ice%hi(:,:,:)                                           ! Eq. 10
      D   (:,:,:)  =  1.0_wp / (6.0_wp * dtime * K2 + rhoi*ice%hi(:,:,:)*ci)                 
      iK1B(:,:,:)  =  1.0_wp / (K1 + B)

     ! Set temperature at which surface is fully liquid
      WHERE (ice%hs(:,:,:) > 1e-6_wp) 
        Tsurfm(:,:,:)  =  0.0_wp
      ELSEWHERE
        Tsurfm(:,:,:)  =  - muS
      END WHERE

      
      A1a       (:,:,:)  =  rhoi*ice%hi(:,:,:) * idt2 * ci &
        &                    + K2* (4.0_wp * dtime * K2 + rhoi*ice%hi(:,:,:)*ci)*D 
      A1        (:,:,:)  =  A1a + K1*B * iK1B                                              ! Eq. 16
      B1a       (:,:,:)  =  -rhoi*ice%hi(:,:,:)* (ci*ice%T1(:,:,:) &
        &                    - alf*muS/ice%T1(:,:,:)) * idt2 - I_0 & 
        &                    - K2*(4.0_wp*dtime*K2*Tfw(:,:,:)&
        &                    + rhoi*ice%hi(:,:,:)*ci*ice%T2(:,:,:))*D
      B1        (:,:,:)  =  B1a + A*K1*iK1B                                                ! Eq. 17
      C1        (:,:,:)  =  - rhoi*ice%hi(:,:,:) * alf * muS * idt2                        ! Eq. 18
      ice%T1    (:,:,:)  =  -(B1 + SQRT(B1*B1-4.0_wp*A1*C1)) / (2.0_wp*A1)                 ! Eq. 21
      ice%Tsurf (:,:,:)  =  (K1*ice%T1(:,:,:)-A) * iK1B                                    ! Eq.  6


      WHERE ( ice%Tsurf(:,:,:) > Tsurfm(:,:,:) ) 
        A1           (:,:,:)  =  A1a + K1                                                  ! Eq. 19
        B1           (:,:,:)  =  B1a - K1*Tsurfm                                           ! Eq. 20
        ice%T1       (:,:,:)  =  -(B1 + SQRT(B1*B1-4.0_wp*A1*C1)) / (2.0_wp*A1)            ! Eq. 21
        ice%Tsurf    (:,:,:)  =  Tsurfm                               
        ! Sum up heatfluxes available for melting at ice surface for each atmopheric time step.
        ! ice%Qtop will be averaged in ave_fluxes
        ice%Qtop     (:,:,:)  =  ice% Qtop(:,:,:) + K1*(ice%T1(:,:,:)-ice%Tsurf(:,:,:)) &
          &                       - (A + B*ice%Tsurf(:,:,:))                               ! Eq. 22
      END WHERE
     
     
      ! Eq. 15
      ice%T2     (:,:,:)  =  ( 2.0_wp*dtime*K2*(ice%T1(:,:,:)+2.0_wp*Tfw(:,:,:)) &
        &                     + rhoi*ice%hi(:,:,:)*ci*ice%T2(:,:,:)) * D
      ! Sum up conductive heatflux at ice-ocean interface for each atmospheric time step. ice%Qtop
      ! will be averaged in ave_fluxes The ocean heat flux is calculated at the beginning of
      ! ice_growth
      ice% Qbot  (:,:,:)  =  ice% Qbot(:,:,:) &
        &                     - 4.0_wp*Ki*(Tfw(:,:,:)-ice%T2(:,:,:))/ice%hi(:,:,:)         ! Eq. 23
    END WHERE

    ipl_src=1  ! output print level (1-5, fix)
    CALL print_mxmn('ice%Tsurf',1,ice%Tsurf(:,1,:),1,ppatch%nblks_c,'ice',ipl_src)

  END SUBROUTINE set_ice_temp_winton



  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! set_ice_temp_winton:: ice_growth_winton - change ice and snow thickness (Winton 2000, section 2b)
  !! This function changes:
  !! ice % hs       new snow thickness for each ice category                [m]
  !! ice % hi       new ice  thickness for each ice category                [m]
  !! ice % hsold    old snow thickness for each ice category                [m]
  !! ice % hiold    old ice  thickness for each ice category                [m]
  !! ice % T1       the new upper ice+snow temp.  for each ice category     [�C]
  !! ice % T2       the new lower ice temperature for each ice category     [�C]
  !! ice % evapwi   amount of evaporated water from the mixed layer
  !!                in previously ice covered areas if all ice is gone      [kg/m�]
  !! ice % heatOceI to contain the energy that is available to the mixed layer
  !!                in previously ice covered areas if all ice is gone      [J]
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!
  SUBROUTINE ice_growth_winton(ppatch, p_os, ice, rpreci)!, lat)
    TYPE(t_patch),             INTENT(IN)    :: ppatch 
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os
    TYPE (t_sea_ice),          INTENT(INOUT) :: ice
    REAL(wp),                  INTENT(IN)    :: rpreci(:,:) 
                                   ! water equiv. solid precipitation rate [m/s] DIMENSION (ie,je)
    !REAL(wp),                  INTENT(IN)    :: lat(:,:,:) 
                                   !! lat. heat flux  [W/m�] DIMENSION (ie,je,kice)

    !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, ppatch%nblks_c) ::         &
      & below_water, & ! Thickness of snow layer below water line           [m]
      & C1,          & ! L  * rhos * hs                                     [J/m�]
      & C2,          & ! E1 * rhoi * h1                                     [J/m�]
      & C3,          & ! E2 * rhoi * h2                                     [J/m�]
      & delh2,       & ! increase of bottom layer thickness (Eq. 24)        [m]
      & draft,       & ! depth of ice-ocean interface below sea level       [m]
      & E1,          & ! Energy content of upper ice+snow layer             [J/kg]
      & E2,          & ! Energy content of lower ice      layer             [J/kg]
      & f1,          & ! Fraction of upper ice in new ice layer (Eq. 37)
      & h1,          & ! Thickness of upper ice layer                       [m]
      & h2,          & ! Thickness of lower ice layer                       [m] 
      & new_snow3d,  & ! New snow fallen onto each ice category             [m]
      !& subli,       & ! Amount of ice+snow that is sublimated away         [kg/m�]
      & Tbar,        & ! Dummy temperature for new temperature calculation  [�C]
      & surfmeltsn,  & ! Surface melt water from snow melt with T=0�C       [m]
      & surfmelti1,  & ! Surface melt water from upper ice with T=-muS      [m]
      & surfmelti2,  & ! Surface melt water from lower ice with T=-muS      [m]
      & zHeatOceI,   & ! Oceanic heat flux                                  [W/m^2]
      & Tfw            ! Ocean freezing temperature [°C]

    INTEGER k

    ! Necessary initialisation
    delh2=0._wp
    surfmelti1 = 0.0_wp
    surfmelti2 = 0.0_wp


    !-------------------------------------------------------------------------------
    ! Calculate snow fall and create array split into ice categories
    new_snow3d (:,1,:)   = rpreci (:,:) * dtime * rho_ref / rhos 
    FORALL(k=2:ice%kice)  new_snow3d(:,k,:)  = new_snow3d(:,1,:)


    ! freezing temperature of uppermost sea water
    IF ( no_tracer >= 2 ) THEN
      DO k=1,ice%kice
        Tfw(:,k,:) = -mu*p_os%p_prog(nold(1))%tracer(:,1,:,2)
      ENDDO
    ELSE
      Tfw(:,:,:) = Tf
    ENDIF
    
    ! Heat flux from ocean into ice
    CALL oce_ice_heatflx (p_os,ice,Tfw,zHeatOceI)
    ! Add oceanic heat flux to energy available at the bottom of the ice.
    ice%Qbot(:,:,:) = ice%Qbot(:,:,:) + zHeatOceI(:,:,:)
   
    ! Do the following wherever there is ice
    !isice: &  
    WHERE (ice% isice(:,:,:))
     
      ! Save ice thickness at previous time step for calculation of heat and salt
      ! flux into ocean in subroutine upper_ocean_TS
      ice % hiold (:,:,:) = ice%hi(:,:,:)
      ice % hsold (:,:,:) = ice%hs(:,:,:)

      h1(:,:,:) = ice%hi(:,:,:) / 2.0_wp
      h2(:,:,:) = h1(:,:,:)


      ! Apply mass increasing changes first. 
      ! 1. Snow fall

      ice%hs(:,:,:) = ice%hs(:,:,:) + new_snow3d(:,:,:)
      
      ! 2. Bottom ice-growth  (maybe add frazil ice?)

      ! #eoo# Eqns. 24, 27--29 and 31--36 appear to be missing rhoi or rhos to get the proper units
      ! for Delta h. But these are included in this program
      WHERE (ice%Qbot < 0.0_wp) 
        delh2  (:,:,:) = ice%Qbot(:,:,:) * dtime / (rhoi * (ci * (Tfw(:,:,:) + muS) - alf)) ! Eq. 24 & 25
        ice%T2 (:,:,:) = (delh2*Tfw(:,:,:) + h2 * ice%T2(:,:,:)) / (delh2 + h2)             ! Eq. 26
        h2     (:,:,:) = h2 + delh2
      END WHERE

      ! Now mass decreasing changes. 
      ! 1. Evaporation
      ! #eoo# Not in Winton - does this count the latent fluxes twice?

      !subli(:,:,:) = lat  / als * dtime;    ![kg/m�]
      !WHERE     (subli <= ice%hs*rhos )         
      !  ice%hs(:,:,:) = ice%hs - subli / rhos
      !ELSEWHERE (subli <= ice%hs*rhos + h1*rhoi )              ! if all snow is gone
      !  ice%hs(:,:,:) = 0.0_wp
      !  h1(:,:,:) = h1 - (subli - ice%hs*rhos) / rhoi
      !ELSEWHERE (subli <= ice%hs*rhos + (h1+h2)*rhoi )         ! if upper ice is gone
      !  ice%hs(:,:,:) = 0.0_wp
      !  h1(:,:,:) = 0.0_wp
      !  h2(:,:,:) = h2 - (subli - ice%hs*rhos - h1*rhoi) / rhoi
      !ELSEWHERE                                                ! if all ice is gone
      !  ice%hs(:,:,:) = 0.0_wp
      !  h1(:,:,:) = 0.0_wp
      !  h2(:,:,:) = 0.0_wp
      !  ice% evapwi(:,:,:) = (subli - ice%hs*rhos - (h1+h2)*rhoi) * als / alv
      !END WHERE
     
   
      ! 2. surface ablation (if any) 

      E1(:,:,:) = ci * ( ice%T1(:,:,:)+muS ) - alf*(1.0_wp+muS/ice%T1(:,:,:))    ! Eq.  1 (energy upper layer)
      E2(:,:,:) = ci * ( ice%T2(:,:,:)+muS ) - alf                        ! Eq. 25 (energy lower layer)
      C1(:,:,:) = alf  * rhos * ice%hs(:,:,:)
      C2(:,:,:) = E1(:,:,:) * rhoi * h1
      C3(:,:,:) = E2(:,:,:) * rhoi * h2
    
      WHERE ( ice%Qtop(:,:,:) > 0.0_wp ) 
        surfmeltsn   (:,:,:) = MIN(ice%Qtop(:,:,:)*dtime / (alf * rhos), ice%hs(:,:,:))
        ice%hs       (:,:,:) = ice%hs(:,:,:) - surfmeltsn                               ! Eq. 27
        ice%surfmelt (:,:,:) = surfmeltsn * rhos/rho_ref
        WHERE (ice%hs(:,:,:) <= 0.0_wp) 
          surfmelti1   (:,:,:) = MIN((ice%Qtop(:,:,:)*dtime-C1) / (-E1(:,:,:)*rhoi), h1)
          h1           (:,:,:) = h1 - surfmelti1                                        ! Eq. 28
          ice%surfmelt (:,:,:) = ice%surfmelt(:,:,:) + surfmelti1 * rhoi/rho_ref
          WHERE (h1(:,:,:) <= 0.0_wp) 
            surfmelti2   (:,:,:) = MIN((ice%Qtop(:,:,:)*dtime-C1+C2) / (-E2(:,:,:)*rhoi), h2)
            h2           (:,:,:) = h2 - surfmelti2                                      ! Eq. 29
            ice%surfmelt (:,:,:) = ice%surfmelt(:,:,:) + surfmelti2 * rhoi/rho_ref
            WHERE (h2(:,:,:) <= 0.0_wp) 
              ice% heatOceI(:,:,:) = ice%Qtop(:,:,:) + (-C1 + C2 + C3)/dtime            ! Eq. 30
              !Flux - not heat
              !ice% heatOceI(:,:,:) = ice%Qtop*dtime - C1 + C2 + C3                     ! Eq. 30
            END WHERE
          END WHERE
        END WHERE
        ! Calculate average temperature of surface melt water 
        ! T(snow) = 0�C, T(ice) = -muS �C
        ice%surfmeltT(:,:,:) = (surfmelti1+surfmelti2) * (-muS) /  ice%surfmelt(:,:,:)
      END WHERE
     
      C1(:,:,:) = alf    * rhos * ice%hs(:,:,:)
      C2(:,:,:) = E1(:,:,:) * rhoi * h1
      C3(:,:,:) = E2(:,:,:) * rhoi * h2
   
     ! 3. bottom ablation (if any)

      WHERE ( ice%Qbot(:,:,:) > 0.0_wp ) 
        h2 (:,:,:) = h2 - MIN(ice%Qbot(:,:,:) * dtime/ (-E2(:,:,:)*rhoi), h2)             ! Eq. 31
        WHERE (h2(:,:,:) <= 0.0_wp) 
          h1 (:,:,:) = h1 - MIN((ice%Qbot(:,:,:) * dtime  + C3) / (-E1(:,:,:)*rhoi), h1)  ! Eq. 32
          WHERE (h1(:,:,:) <= 0.0_wp) 
            ice%hs (:,:,:) = ice%hs(:,:,:) - MIN((ice%Qbot(:,:,:) * dtime+C3+C2)&
              & /(alf*rhos), ice%hs(:,:,:))                                             ! Eq. 33
            WHERE (ice%hs (:,:,:) <= 0.0_wp) 
              ice% heatOceI(:,:,:) = ice% heatocei(:,:,:) + ice%Qbot(:,:,:) &
                &                     + (-C1 + C2 + C3)/dtime                           ! Eq. 34
              ! Flux - not heat
              !ice% heatOceI(:,:,:) = ice% heatocei(:,:,:) + ice%Qbot(:,:,:) * dtime - C1 + C2 + C3    ! Eq. 34
            END WHERE
          END WHERE
        END WHERE
      END WHERE

      ! Calculate ice thickness and draft (ice+snow depth below water line)
      ice%hi      (:,:,:) = h1 + h2
      draft       (:,:,:) = (rhoi*ice%hi(:,:,:)+rhos*ice%hs(:,:,:)) / rho_ref
      below_water (:,:,:) = draft(:,:,:) - ice%hi(:,:,:)

      
      ! snow -> ice conversion for snow below waterlevel
      ! Currently not quite physical: Snow is pushed together to form new ice, hence snow thickness
      ! decreases more than ice thickness by rhoi/rhos ( analogue to old growth.f90 sea-ice model )
      ! Salt content of snow ice is equal to that of normal ice, salt is removed from the ocean
      ! Temperature of new upper ice is calculated as described in the paragraph below 
      ! Eq. 36
      WHERE (below_water (:,:,:) > 0.0_wp) 
        ice% snow_to_ice(:,:,:) = below_water * rhoi / rhos
        ice%hs          (:,:,:) = ice%hs - ice% snow_to_ice(:,:,:)
        f1              (:,:,:) = h1 / (h1+below_water)
        Tbar            (:,:,:) = f1  * ( ice%T1(:,:,:) - alf* muS/(ci*ice%T1(:,:,:)) )&
          &                        - (1.0_wp-f1)*muS 
        ice%T1          (:,:,:) = 0.5_wp * ( Tbar - SQRT(Tbar*Tbar + 4.0_wp*muS*alf/ci) )
        h1              (:,:,:) = h1 + below_water
        ice%hi          (:,:,:) = h1 + h2
      END WHERE

      ! Even up upper and lower layer
      WHERE ( h1(:,:,:) < h2(:,:,:)  ) 
        f1    (:,:,:) =  h1 / (0.5_wp*ice%hi(:,:,:))                                
        Tbar  (:,:,:) =  f1 * ( ice%T1(:,:,:) - alf*muS/(ci*ice%T1(:,:,:)) ) &
          &              + (1.0_wp-f1)*ice%T2(:,:,:)  ! Eq. 39
        ice%T1(:,:,:) =  0.5_wp * ( Tbar - SQRT(Tbar*Tbar + 4.0_wp*muS*alf/ci) )     ! Eq. 38
      ELSEWHERE ( h1(:,:,:) > h2(:,:,:) ) 
        f1    (:,:,:) =  h1 / (0.5_wp*ice%hi(:,:,:)) - 1.0_wp
        ice%T2(:,:,:) =  f1 * ( ice%T1(:,:,:) - alf*muS/(ci*ice%T1(:,:,:)) ) &
          &              + (1.0_wp-f1)*ice%T2(:,:,:)  ! Eq. 40
      END WHERE
    
      ! ice%T2 can get above bulk melting temperature. If this happens, use additional energy to
      ! melt equal thickness of upper and lower layer (last para.  section 2)
      ! Energy available for melting: -h2 * ci * (ice%T2+muS)
      ! Energy needed for melting lower layer: L
      ! Energy needed for melting upper layer: -(ci*(ice%T1+muS)-L*(1+muS/ice%T1)) (Eq. 1)
      WHERE (ice%t2 (:,:,:) > -muS)                  
        ice%hi (:,:,:) = ice%hi(:,:,:) - h2*ci*(ice%T2(:,:,:)+muS) / &
          &            ( 0.5_wp*alf - 0.5_wp*(ci*(ice%T1(:,:,:)+muS) &
          &              - alf*(1.0_wp+muS/ice%T1(:,:,:))) )
        ice%T2 (:,:,:) = -muS
      END WHERE

      ! Is this necessary?
      WHERE (ice%hi(:,:,:) <= 0.0_wp) 
        ice%Tsurf(:,:,:) =  Tfw(:,:,:)
        ice%T1   (:,:,:) =  Tfw(:,:,:)
        ice%T2   (:,:,:) =  Tfw(:,:,:)
        ice%isice(:,:,:) =  .FALSE.
        ice%conc (:,:,:) = 0.0_wp
        ice%hi   (:,:,:) = 0.0_wp
      ELSEWHERE
        ice%heatOceI(:,:,:) = ice%heatOceI(:,:,:) - zHeatOceI(:,:,:)
      END WHERE
    
    END WHERE !isice

    ipl_src=1  ! output print level (1-5, fix)
    CALL print_mxmn('ice%hi',1,ice%hi(:,1,:),1,ppatch%nblks_c,'ice',ipl_src)
     
  END SUBROUTINE ice_growth_winton

END MODULE mo_sea_ice_winton
