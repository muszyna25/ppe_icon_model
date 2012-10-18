!>
!! Provide an implementation of the sea-ice model.
!!
!! Provide an implementation of the parameters of the surface module (sea ice)
!! used between the atmopshere and the hydrostatic ocean model.
!!
!! @author 
!! 
!! @par Revision History
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_sea_ice_winton

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,&
    &                               mu,mus,ci, alf, I_0
  USE mo_ocean_nml,           ONLY: no_tracer 
!  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_oce_state,           ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 

  IMPLICIT NONE

  PRIVATE

!  CHARACTER(len=12)           :: str_module    = 'SeaIceWinton'  ! Output of module for 1 line debug
!  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

  PUBLIC :: ice_growth_winton
  PUBLIC :: set_ice_temp_winton

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
  !! ice % Tsurf    the new surface temperature   for each ice category     [C]
  !! ice % T1       the new upper ice+snow temp.  for each ice category     [C]
  !! ice % T2       the new lower ice temperature for each ice category     [C]
  !! ice % Qbot     Heat flux available for freezing/melting at ice bottom  [W/m^2]
  !! ice % Qtop     Heat flux available for melting at ice surface          [W/m^2]
  !!
  !!           all "dtime" in this function are atmospheric time step
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!

  SUBROUTINE set_ice_temp_winton(p_patch,ice, Tfw, Qatm) 
    TYPE(t_patch),        INTENT(IN), TARGET    :: p_patch 
    TYPE(t_sea_ice),      INTENT(INOUT)         :: ice
    REAL(wp),             INTENT(IN)            :: Tfw(:,:,:)
    TYPE(t_atmos_fluxes), INTENT(IN)            :: Qatm

    !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch%nblks_c) ::           &
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
      & Tsurfm,      & ! Surface melting temperature
      & I              ! Penetrating shortwave radiation (taking snow cover into account)
    
    REAL(wp) :: idt2 ! 1 / (2*dt)

    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c     ! loop indices
    
    
   !-------------------------------------------------------------------------------

    ! initialization
    idt2   =  1.0_wp / (2.0_wp*dtime)
    
    all_cells => p_patch%cells%all 
    
!!$    FORALL(jc=1:nproma, jb=1:p_patch%nblks_c, k=1:ice%kice, ice % isice (i,k,j)) &
!!$      & SWin3d(jc,k,jb) = Qatm% SWin(jc,jb)
    

    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%isice(jc,k,jb)) THEN

            ! No surface penetrating shortwave if there's snow on top
            IF (ice%hs(jc,k,jb) > 0.0_wp ) THEN
              I(jc,k,jb)=0._wp
            ELSE
              I(jc,k,jb)=I_0
            END IF
            
            ! Create array of shortwave radiation split up into ice categories
            ! (purely for computational reasons)
            SWin3d(jc,k,jb) = Qatm% SWin(jc,jb)

            ! Calculate new ice temperature wherever there is ice 
            ! lat > 0, sens > 0, LWnet >0 , SWin > 0  for downward flux 
            ! dlatdT, dsensdT, dLWdT >0 for downward flux increasing with increasing Tsurf
            
            B   (jc,k,jb) = - Qatm% dlatdT(jc,k,jb) - Qatm% dsensdT(jc,k,jb) &
              &              - Qatm% dLWdT(jc,k,jb)                                     ! Eq.  8

            A   (jc,k,jb) = - Qatm% lat(jc,k,jb) - Qatm% sens(jc,k,jb) &
              &             - Qatm% LWnet(jc,k,jb) &
              &             - (1.0_wp - ice% alb(jc,k,jb)) * (1._wp-I(jc,k,jb)) *  SWin3d(jc,k,jb)  &
              &             - ice%Tsurf(jc,k,jb)* B(jc,k,jb)                            ! Eq.  7

            K1  (jc,k,jb)  =  4.0_wp * ki * ks &
              &                / (ks * ice%hi(jc,k,jb) + 4.0_wp * ki * ice%hs(jc,k,jb)) ! Eq.  5

            K2  (jc,k,jb)  =  2.0_wp * ki / ice%hi(jc,k,jb)                             ! Eq. 10

            D   (jc,k,jb)  =  1.0_wp / (6.0_wp * dtime * K2(jc,k,jb) &
              &                          + rhoi*ice%hi(jc,k,jb)*ci)                 

            iK1B(jc,k,jb)  =  1.0_wp / ( K1(jc,k,jb) + B(jc,k,jb) )

            ! Set temperature at which surface is fully liquid
            IF(ice%hs(jc,k,jb) > 1e-6_wp) THEN
              Tsurfm(jc,k,jb)  =  0.0_wp
            ELSE
              Tsurfm(jc,k,jb)  =  - muS
            END IF

      
            !
            A1a       (jc,k,jb)  =  rhoi*ice%hi(jc,k,jb) * idt2 * ci &
              &                     + K2(jc,k,jb)* (4.0_wp * dtime * K2(jc,k,jb) &
              &                                      + rhoi*ice%hi(jc,k,jb)*ci)*D(jc,k,jb)
            ! Eq. 16
            A1        (jc,k,jb)  =  A1a(jc,k,jb) + K1(jc,k,jb)*B(jc,k,jb) * iK1B(jc,k,jb)
            !
            B1a       (jc,k,jb)  = - rhoi * ice%hi(jc,k,jb) &
              &                       * (ci * ice%T1(jc,k,jb) &
              &                          - alf * muS/ice%T1(jc,k,jb)) * idt2 &
              &                    - I(jc,k,jb) & 
              &                    - K2(jc,k,jb) * (4.0_wp * dtime &
              &                                     * K2(jc,k,jb) * Tfw(jc,k,jb)&
              &                    + rhoi * ice%hi(jc,k,jb) * ci &
              &                       * ice%T2(jc,k,jb))*D(jc,k,jb)
            ! Eq. 17
            B1        (jc,k,jb)  =  B1a(jc,k,jb) +   A   (jc,k,jb) &
              &                                    * K1  (jc,k,jb) &
              &                                    * iK1B(jc,k,jb)
            
            ! Eq. 18
            C1        (jc,k,jb)  =  - rhoi*ice%hi(jc,k,jb) * alf * muS * idt2
            ! Eq. 21
            ice%T1    (jc,k,jb)  =  - ( B1(jc,k,jb) + SQRT(B1(jc,k,jb)*B1(jc,k,jb) &
              &                                           - 4.0_wp*A1(jc,k,jb)     &
              &                                             * C1(jc,k,jb)) )       &
              &                     / (2.0_wp*A1(jc,k,jb))
            ! Eq.  6
            ice%Tsurf (jc,k,jb)  =   (K1(jc,k,jb) * ice%T1(jc,k,jb)-A(jc,k,jb)) &
              &                    * iK1B(jc,k,jb)     
            
            
            IF ( ice%Tsurf(jc,k,jb) > Tsurfm(jc,k,jb) ) THEN
              ! Eq. 19
              A1           (jc,k,jb)  =  A1a(jc,k,jb) + K1(jc,k,jb)
              ! Eq. 20
              B1           (jc,k,jb)  =  B1a(jc,k,jb) - K1(jc,k,jb) * Tsurfm(jc,k,jb)
              ! Eq. 21
              ice%T1       (jc,k,jb)  =  -( B1(jc,k,jb) + SQRT(B1(jc,k,jb) * B1(jc,k,jb) &
                &                                               - 4.0_wp * A1(jc,k,jb) &
                &                                                 * C1(jc,k,jb))) &
                &                        / (2.0_wp*A1(jc,k,jb))
              ice%Tsurf    (jc,k,jb)  =  Tsurfm(jc,k,jb)                               
              ! Sum up heatfluxes available for melting at ice surface for each atmopheric time step.
              ! ice%Qtop will be averaged in ave_fluxes
              ! Eq. 22
              ice%Qtop     (jc,k,jb)  =  ice% Qtop(jc,k,jb) &
                &                        + K1(jc,k,jb) * (ice%T1(jc,k,jb) - ice%Tsurf(jc,k,jb)) &
                &                        - ( A(jc,k,jb) + B(jc,k,jb)*ice%Tsurf(jc,k,jb) )
            END IF
     
     
            ! Eq. 15
            ice%T2 (jc,k,jb)  =  ( 2.0_wp * dtime &
              &                     * K2(jc,k,jb) * (ice%T1(jc,k,jb) &
              &                                      + 2.0_wp*Tfw(jc,k,jb)) &
              &                    + rhoi*ice%hi(jc,k,jb) * ci &
              &                       * ice%T2(jc,k,jb)) * D(jc,k,jb)
            ! Sum up conductive heatflux at ice-ocean interface for each atmospheric time step. ice%Qtop
            ! will be averaged in ave_fluxes The ocean heat flux is calculated at the beginning of
            ! ice_growth
             ! Eq. 23
            ice% Qbot  (jc,k,jb)  =  ice% Qbot(jc,k,jb) &
              &                     - 4.0_wp*Ki*(Tfw(jc,k,jb)-ice%T2(jc,k,jb))/ice%hi(jc,k,jb)
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE set_ice_temp_winton

  !-------------------------------------------------------------------------------
  !
  !  
  !>
  !! ! ice_growth_winton - change ice and snow thickness (Winton 2000, section 2b)
  !! This function changes:
  !! ice % hs       new snow thickness for each ice category                [m]
  !! ice % hi       new ice  thickness for each ice category                [m]
  !! ice % hsold    old snow thickness for each ice category                [m]
  !! ice % hiold    old ice  thickness for each ice category                [m]
  !! ice % T1       the new upper ice+snow temp.  for each ice category     [C]
  !! ice % T2       the new lower ice temperature for each ice category     [C]
  !! ice % evapwi   amount of evaporated water from the mixed layer
  !!                in previously ice covered areas if all ice is gone      [kg/m^3]
  !! ice % heatOceI to contain the energy that is available to the mixed layer
  !!                in previously ice covered areas if all ice is gone      [J]
  !!
  !! @par Revision History
  !! Initial release by Peter Korn, MPI-M (2010-07). Originally code written by
  !! Dirk Notz, following MPI-OM. Code transfered to ICON.
  !!

  SUBROUTINE ice_growth_winton(p_patch, p_os, ice, rpreci)!, lat)
    TYPE(t_patch),             INTENT(IN), TARGET    :: p_patch 
    TYPE(t_hydro_ocean_state), INTENT(IN)            :: p_os
    TYPE (t_sea_ice),          INTENT(INOUT)         :: ice
    REAL(wp),                  INTENT(IN)            :: rpreci(:,:) 
                                   ! water equiv. solid precipitation rate [m/s] DIMENSION (ie,je)
    !REAL(wp),                  INTENT(IN)    :: lat(:,:,:) 
                                   !! lat. heat flux  [W/m^2] DIMENSION (ie,je,kice)

    !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch%nblks_c) ::         &
      & below_water, & ! Thickness of snow layer below water line           [m]
      & C1,          & ! L  * rhos * hs                                     [J/m^2]
      & C2,          & ! E1 * rhoi * h1                                     [J/m^2]
      & C3,          & ! E2 * rhoi * h2                                     [J/m^2]
      & delh2,       & ! increase of bottom layer thickness (Eq. 24)        [m]
      & draft,       & ! depth of ice-ocean interface below sea level       [m]
      & E1,          & ! Energy content of upper ice+snow layer             [J/kg]
      & E2,          & ! Energy content of lower ice      layer             [J/kg]
      ! for energy conservation
      & Q_surplus,   &
      !
      & f1,          & ! Fraction of upper ice in new ice layer (Eq. 37)
      & h1,          & ! Thickness of upper ice layer                       [m]
      & h2,          & ! Thickness of lower ice layer                       [m] 
      & new_snow3d,  & ! New snow fallen onto each ice category             [m]
      !& subli,       & ! Amount of ice+snow that is sublimated away         [kg/m^3]
      & Tbar,        & ! Dummy temperature for new temperature calculation  [C]
      & surfmeltsn,  & ! Surface melt water from snow melt with T=0C        [m]
      & surfmelti1,  & ! Surface melt water from upper ice with T=-muS      [m]
      & surfmelti2,  & ! Surface melt water from lower ice with T=-muS      [m]
      & zHeatOceI,   & ! Oceanic heat flux                                  [W/m^2]
      & Tfw            ! Ocean freezing temperature [C]

    TYPE(t_subset_range), POINTER :: all_cells
    
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c     ! loop indices

    
    ! Necessary initialisation
    delh2=0._wp
    surfmelti1 = 0.0_wp
    surfmelti2 = 0.0_wp

    !
    ice%E1(:,:,:) =     0.0_wp
    ice%E2(:,:,:) =     0.0_wp
    Q_surplus(:,:,:) = 0.0_wp

    !
    all_cells => p_patch%cells%all 






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
    


    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          ! Do the following wherever there is ice
          !
          IF (ice%isice(jc,k,jb)) THEN
            
            ! Add oceanic heat flux to energy available at the bottom of the ice.
            ice%Qbot(jc,k,jb) = ice%Qbot(jc,k,jb) + zHeatOceI(jc,k,jb)
            
            
              
            ! Save ice thickness at previous time step for calculation of heat and salt
            ! flux into ocean in subroutine upper_ocean_TS
            ice % hiold (jc,k,jb) = ice%hi(jc,k,jb)
            ice % hsold (jc,k,jb) = ice%hs(jc,k,jb)
            
            h1(jc,k,jb) = ice%hi(jc,k,jb) / 2.0_wp
            h2(jc,k,jb) = h1(jc,k,jb)
            
            
            ! Apply mass increasing changes first. 
            ! 1. Snow fall
            
            ice%hs(jc,k,jb) = ice%hs(jc,k,jb) + new_snow3d(jc,k,jb)
            
            ! 2. Bottom ice-growth  (maybe add frazil ice?)
            
            ! #eoo# Eqns. 24, 27--29 and 31--36 appear to be missing rhoi or rhos to get the proper units
            ! for Delta h. But these are included in this program
            IF (ice%Qbot(jc,k,jb) < 0.0_wp) THEN
              ! Eq. 24 & 25
              delh2  (jc,k,jb) = ice%Qbot(jc,k,jb) * dtime &
                &                / (rhoi * (ci * (Tfw(jc,k,jb) + muS) - alf)) 
              
              ! Eq. 26
              ice%T2 (jc,k,jb) =   (delh2(jc,k,jb)*Tfw(jc,k,jb) + h2(jc,k,jb) * ice%T2(jc,k,jb)) &
                &                / (delh2(jc,k,jb) + h2(jc,k,jb))  
              
              h2     (jc,k,jb) = h2(jc,k,jb) + delh2(jc,k,jb)
            END IF
            
            ! Now mass decreasing changes. 
            ! 1. Evaporation
            ! #eoo# Not in Winton - does this count the latent fluxes twice?
            
            !subli(jc,k,jb) = lat  / als * dtime;    ![kg/m^3]
            !WHERE     (subli <= ice%hs*rhos )         
            !  ice%hs(jc,k,jb) = ice%hs - subli / rhos
            !ELSEWHERE (subli <= ice%hs*rhos + h1*rhoi )              ! if all snow is gone
            !  ice%hs(jc,k,jb) = 0.0_wp
            !  h1(jc,k,jb) = h1 - (subli - ice%hs*rhos) / rhoi
            !ELSEWHERE (subli <= ice%hs*rhos + (h1+h2)*rhoi )         ! if upper ice is gone
            !  ice%hs(jc,k,jb) = 0.0_wp
            !  h1(jc,k,jb) = 0.0_wp
            !  h2(jc,k,jb) = h2 - (subli - ice%hs*rhos - h1*rhoi) / rhoi
            !ELSEWHERE                                                ! if all ice is gone
            !  ice%hs(jc,k,jb) = 0.0_wp
            !  h1(jc,k,jb) = 0.0_wp
            !  h2(jc,k,jb) = 0.0_wp
            !  ice% evapwi(jc,k,jb) = (subli - ice%hs*rhos - (h1+h2)*rhoi) * als / alv
            !END WHERE
            
            
            ! 2. surface ablation (if any) 
            
            ! Eq.  1 (energy upper layer)
            E1(jc,k,jb) = ci * ( ice%T1(jc,k,jb)+muS ) - alf*(1.0_wp+muS/ice%T1(jc,k,jb))
            ! Eq. 25 (energy lower layer)
            E2(jc,k,jb) = ci * ( ice%T2(jc,k,jb)+muS ) - alf                        
            
            C1(jc,k,jb) = alf  * rhos * ice%hs(jc,k,jb)
            C2(jc,k,jb) = E1(jc,k,jb) * rhoi * h1(jc,k,jb)
            C3(jc,k,jb) = E2(jc,k,jb) * rhoi * h2(jc,k,jb)
            
            IF ( ice%Qtop(jc,k,jb) > 0.0_wp ) THEN
              surfmeltsn   (jc,k,jb) = MIN(ice%Qtop(jc,k,jb)*dtime / (alf * rhos),&
                &                          ice%hs(jc,k,jb)                          )
              
              ! Eq. 27
              ice%hs       (jc,k,jb) = ice%hs(jc,k,jb) - surfmeltsn(jc,k,jb) 
              ice%surfmelt (jc,k,jb) = surfmeltsn(jc,k,jb) * rhos/rho_ref
              
              IF (ice%hs(jc,k,jb) <= 0.0_wp) THEN
                surfmelti1   (jc,k,jb) = MIN((ice%Qtop(jc,k,jb)*dtime - C1(jc,k,jb)) &
                  &                             / (-E1(jc,k,jb)*rhoi),&
                  &                           h1(jc,k,jb)                            )
                ! Eq. 28
                h1           (jc,k,jb) = h1(jc,k,jb) - surfmelti1(jc,k,jb)
                ice%surfmelt (jc,k,jb) =   ice%surfmelt(jc,k,jb) &
                  &                      + surfmelti1(jc,k,jb) * rhoi/rho_ref
                
                ! Eq. 29
                IF (h1(jc,k,jb) <= 0.0_wp) THEN
                  surfmelti2   (jc,k,jb) = MIN(  ( ice%Qtop(jc,k,jb)*dtime     &
                    &                            - C1(jc,k,jb) + C2(jc,k,jb) ) &
                    &                           / (-E2(jc,k,jb)*rhoi),         &
                    &                          h2(jc,k,jb)                      )
                  h2           (jc,k,jb) = h2(jc,k,jb) - surfmelti2(jc,k,jb)                   
                  ice%surfmelt (jc,k,jb) =   ice%surfmelt(jc,k,jb) &
                    &                      + surfmelti2(jc,k,jb) * rhoi/rho_ref
                  
                  IF (h2(jc,k,jb) <= 0.0_wp) THEN
                    ! Eq. 30
                    ice% heatOceI(jc,k,jb) =   ice%Qtop(jc,k,jb) &
                      &                      + (-C1(jc,k,jb) + C2(jc,k,jb) + C3(jc,k,jb))/dtime
                    !Flux - not heat
                    ! Eq. 30
                    !ice% heatOceI(jc,k,jb) = ice%Qtop*dtime - C1 + C2 + C3
                  END IF
                END IF
              END IF
              ! Calculate average temperature of surface melt water 
              ! T(snow) = 0C, T(ice) = -muS C
              ice%surfmeltT(jc,k,jb) =   (surfmelti1(jc,k,jb) + surfmelti2(jc,k,jb)) &
                &                      * (-muS) /  ice%surfmelt(jc,k,jb)
            END IF
            
            C1(jc,k,jb) = alf    * rhos * ice%hs(jc,k,jb)
            C2(jc,k,jb) = E1(jc,k,jb) * rhoi * h1(jc,k,jb)
            C3(jc,k,jb) = E2(jc,k,jb) * rhoi * h2(jc,k,jb)
            
            ! 3. bottom ablation (if any)
            
            IF ( ice%Qbot(jc,k,jb) > 0.0_wp ) THEN
              ! Eq. 31
              h2 (jc,k,jb) = h2(jc,k,jb) - MIN( ice%Qbot(jc,k,jb) * dtime&
                &                                  / (-E2(jc,k,jb) * rhoi) ,&
                &                               h2(jc,k,jb)                   ) 
              
              IF (h2(jc,k,jb) <= 0.0_wp) THEN
                
                ! Eq. 32
                h1 (jc,k,jb) = h1(jc,k,jb) - MIN((ice%Qbot(jc,k,jb) * dtime + C3(jc,k,jb)) &
                  &                               / (-E1(jc,k,jb)*rhoi),                 &
                  &                              h1(jc,k,jb)                        )
                IF (h1(jc,k,jb) <= 0.0_wp) THEN
                  ! Eq. 33
                  ice%hs (jc,k,jb) =  ice%hs(jc,k,jb) &
                    &               - MIN(  (ice%Qbot(jc,k,jb) * dtime &
                    &                        + C3(jc,k,jb) + C2(jc,k,jb))&
                    &                      / (alf*rhos),&
                    &                     ice%hs(jc,k,jb)                           )
                  IF (ice%hs (jc,k,jb) <= 0.0_wp) THEN
                    ! Eq. 34
                    ice% heatOceI(jc,k,jb) = ice% heatocei(jc,k,jb) + ice%Qbot(jc,k,jb) &
                      &                     + (-C1(jc,k,jb) + C2(jc,k,jb) + C3(jc,k,jb)) / dtime                          
                    ! Flux - not heat
                    !ice% heatOceI(jc,k,jb) = ice% heatocei(jc,k,jb) + ice%Qbot(jc,k,jb) * dtime - C1 + C2 + C3 ! Eq. 34
                  END IF
                END IF
              END IF
            END IF
            
            ! Calculate ice thickness and draft (ice+snow depth below water line)
            ice%hi      (jc,k,jb) = h1(jc,k,jb) + h2(jc,k,jb)
            draft       (jc,k,jb) = (rhoi*ice%hi(jc,k,jb)+rhos*ice%hs(jc,k,jb)) / rho_ref
            below_water (jc,k,jb) = draft(jc,k,jb) - ice%hi(jc,k,jb)
            
            
            ! snow -> ice conversion for snow below waterlevel
            ! Currently not quite physical: Snow is pushed together to form new ice, hence snow thickness
            ! decreases more than ice thickness by rhoi/rhos ( analogue to old growth.f90 sea-ice model )
            ! Salt content of snow ice is equal to that of normal ice, salt is removed from the ocean
            ! Temperature of new upper ice is calculated as described in the paragraph below 
            ! Eq. 36
            IF (below_water (jc,k,jb) > 0.0_wp) THEN
              ice% snow_to_ice(jc,k,jb) = below_water(jc,k,jb) * rhoi / rhos
              ice%hs          (jc,k,jb) = ice%hs(jc,k,jb) - ice% snow_to_ice(jc,k,jb)
              f1              (jc,k,jb) =  h1(jc,k,jb) &
                &                         / ( h1(jc,k,jb)+below_water(jc,k,jb) )
              Tbar            (jc,k,jb) =   f1(jc,k,jb)  * ( ice%T1(jc,k,jb) &
                &                                    - alf* muS / (ci*ice%T1(jc,k,jb)) ) &
                &                         - (1.0_wp-f1(jc,k,jb))*muS 
              ice%T1          (jc,k,jb) = 0.5_wp * ( Tbar(jc,k,jb) &
                &                                  - SQRT(  Tbar(jc,k,jb)*Tbar(jc,k,jb) &
                &                                           + 4.0_wp*muS*alf/ci) )
              h1              (jc,k,jb) = h1(jc,k,jb) + below_water(jc,k,jb)
              ice%hi          (jc,k,jb) = h1(jc,k,jb) + h2(jc,k,jb)
            END IF
            
            ! Even up upper and lower layer
            IF( h1(jc,k,jb) < h2(jc,k,jb)  ) THEN
              f1    (jc,k,jb) =  h1(jc,k,jb) / (0.5_wp*ice%hi(jc,k,jb)) 
              ! Eq. 39
              Tbar  (jc,k,jb) =  f1(jc,k,jb) * ( ice%T1(jc,k,jb) - alf*muS/(ci*ice%T1(jc,k,jb)) ) &
                &              + (1.0_wp-f1(jc,k,jb)) * ice%T2(jc,k,jb) 
              ! Eq. 38
              ice%T1(jc,k,jb) =  0.5_wp * ( Tbar(jc,k,jb) - SQRT(Tbar(jc,k,jb)*Tbar(jc,k,jb) &
                &                                       + 4.0_wp*muS*alf/ci) )
            ELSE IF ( h1(jc,k,jb) > h2(jc,k,jb) ) THEN
              f1    (jc,k,jb) =  h1(jc,k,jb) / (0.5_wp*ice%hi(jc,k,jb)) - 1.0_wp
              ice%T2(jc,k,jb) =  f1(jc,k,jb) &
                &                * ( ice%T1(jc,k,jb) - alf*muS/(ci*ice%T1(jc,k,jb)) ) &
                &              + (1.0_wp-f1(jc,k,jb))*ice%T2(jc,k,jb)  ! Eq. 40
            END IF
            
            ! ice%T2 can get above bulk melting temperature. If this happens, use additional energy to
            ! melt equal thickness of upper and lower layer (last para.  section 2)
            ! Energy available for melting: -h2 * ci * (ice%T2+muS)
            ! Energy needed for melting lower layer: L
            ! Energy needed for melting upper layer: -(ci*(ice%T1+muS)-L*(1+muS/ice%T1)) (Eq. 1)
            IF (ice%t2 (jc,k,jb) > -muS) THEN
              ice%hi (jc,k,jb) = ice%hi(jc,k,jb) - h2(jc,k,jb) * ci * (ice%T2(jc,k,jb)+muS) &
                &               / ( 0.5_wp*alf - 0.5_wp*(ci*(ice%T1(jc,k,jb)+muS) &
                &              - alf*(1.0_wp+muS/ice%T1(jc,k,jb))) )
              ice%T2 (jc,k,jb) = -muS
            END IF
            
            ! Is this necessary?
            IF (ice%hi(jc,k,jb) <= 0.0_wp) THEN
              ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
              ice%T1   (jc,k,jb) =  Tfw(jc,k,jb)
              ice%T2   (jc,k,jb) =  Tfw(jc,k,jb)
              ice%isice(jc,k,jb) =  .FALSE.
              ice%conc (jc,k,jb) = 0.0_wp
              ice%hi   (jc,k,jb) = 0.0_wp
            ELSE
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - zHeatOceI(jc,k,jb)
            END IF
            
            ! Eq.  1 (energy upper layer)
            ice%E1(jc,k,jb) =   ci * ( ice%T1(jc,k,jb)+muS ) &
              &              - alf * (1.0_wp+muS/ice%T1(jc,k,jb))
            ! Eq. 25 (energy lower layer)
            ice%E2(jc,k,jb) = ci * ( ice%T2(jc,k,jb)+muS ) - alf
            
            ! Energy balance: discrepancy in eqn. (2)
            
!            Q_surplus(jc,k,jb) =   (   ice%E1(jc,k,jb) - E1(jc,k,jb)            &
!              &                      + ice%E2(jc,k,jb) - E2(jc,k,jb) ) / dtime  &
!              &                  - (ice%Qtop(jc,k,jb) + ice%Qbot(jc,k,jb) )    &
!              &                / (ice%hs(jc,k,jb)*rhos + ice%hi(jc,k,jb)*rhoi)
            
          END IF  !isice
        END DO
      END DO
    END DO

!!$    CALL print_cells(Q_surplus(:,1,:),'Q_surplus')
!!$    CALL print_cells(ice%hi(:,1,:),'ice%hi')
!!$    CALL print_cells(ice%Qtop(:,1,:),'ice%Qtop')
!!$    CALL print_cells(ice%Qbot(:,1,:),'ice%Qbot')
!!$    CALL print_cells(ice%hi(:,1,:)-ice%hiold(:,1,:),'new ice')
!!$    CALL print_cells(zHeatOceI(:,1,:),'zHeatOceI')
!!$
!!$    CALL print_cells(ice%Tsurf(:,1,:),'ice%Tsurf')

!!$! !---------DEBUG DIAGNOSTICS-------------------------------------------
!!$    idt_src=1 !3  ! output print level (1-5, fix)
!!$    CALL dbg_print('GrowZero: Q_surplus'       ,Q_surplus                ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%hi'          ,ice%hi                   ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%Qtop'        ,ice%Qtop                 ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%Qbot'        ,ice%Qbot                 ,str_module,idt_src)
!!$    CALL dbg_print('GrowZero: ice%Tsurf'       ,ice%Tsurf                ,str_module,idt_src)
!!$! !---------------------------------------------------------------------
 

    
    
  END SUBROUTINE ice_growth_winton

END MODULE mo_sea_ice_winton
