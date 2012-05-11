MODULE mo_sea_ice_zerolayer

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

  PUBLIC :: set_ice_temp_zerolayer
  PUBLIC :: ice_growth_zerolayer
  PUBLIC :: ice_flux_balance_zerolayer

CONTAINS
  
  SUBROUTINE ice_flux_balance_zerolayer (ppatch,ice,Qatm,SWin3d,LWin3d,TsurfK,k_effective,&
    &                                    F_A,F_S)
    TYPE(t_patch),        INTENT(IN)    :: ppatch 
    TYPE(t_sea_ice),      INTENT(IN)    :: ice
    TYPE(t_atmos_fluxes), INTENT(IN)    :: Qatm
    REAL(wp),             INTENT(OUT)   :: F_A(:,:,:), F_S(:,:,:) 
    REAL(wp),             INTENT(IN)    :: SWin3d(:,:,:), LWin3d(:,:,:),&
      &                                    k_effective(:,:,:), TsurfK(:,:,:)

    
    ! Local
    REAL(wp), DIMENSION (nproma,ice%kice, ppatch%nblks_c) :: I_0_eff

    WHERE ( ice%isice(:,:,:) )
      ! fraction of LW radiation that penetrates further into the ice
      WHERE ( ice%hs(:,:,:) <= 0.0_wp )
        ! in case of no snow:
        I_0_eff(:,:,:) = I_0 ! what value to put here really?
      ELSEWHERE ( ice%hs(:,:,:) > 0.0_wp )
        ! in case of snow
        I_0_eff(:,:,:) = 0.0_wp
      END WHERE
      
      ! net flux between ice and atmosphere F_A, positive upward
      F_A(:,:,:) = &
        & StBo * TsurfK(:,:,:) ** 4 &
        & - (1.0_wp - ice%alb(:,:,:)) *  (1.0_wp - I_0_eff) * SWin3d(:,:,:) & 
        & - Qatm%sens(:,:,:) - Qatm%lat(:,:,:) - LWin3d(:,:,:) 
      ! conductive heat flux coming through the ice. F_S > 0 for surface colder than bottom
      F_S(:,:,:) = k_effective(:,:,:) * (Tf - ice%Tsurf(:,:,:))
    END WHERE
  END SUBROUTINE ice_flux_balance_zerolayer



  !
  ! Achim: Start implementing Semtner's 0-Layer Model (1976)
  !


  !!    Semtner, Albert J., 1976: A Model for the Thermodynamic Growth of Sea
  !!    Ice in Numerical Investigations of Climate. J. Phys. Oceanogr., 6,
  !!    379–389.  doi:
  !!    http://dx.doi.org/10.1175/1520-0485(1976)006<0379:AMFTTG>2.0.CO;2
  ! (Appendix)

  SUBROUTINE set_ice_temp_zerolayer(ppatch,ice, Tfw, Qatm) 
    TYPE(t_patch),        INTENT(IN)    :: ppatch 
    TYPE(t_sea_ice),      INTENT(INOUT) :: ice
    REAL(wp),             INTENT(IN)    :: Tfw(:,:,:)
    TYPE(t_atmos_fluxes), INTENT(IN)    :: Qatm

    ! Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, ppatch%nblks_c) ::          &
      & k_effective ,  &  ! total heat conductivity of ice/snow
      & sigma_T_p3  ,  &  ! previous temperature ^3 times Stefan-Boltzmann
      & deltaT      ,  &  ! Temperature increment 
      & TsurfK      ,  &  ! Tsurf in [K]
      & SWin3D      ,  &  ! Short-wave radiation field split into ice categories
      & LWin3D      ,  &  ! Long-wave radiation field split into ice categories
      & F_A_p       ,  &  ! F_A with outg. LW from previous time step
      & F_S_p       ,  &  ! conductive flux with T_S from previous time step
      & F_A         ,  &  ! F_A with
      & F_S             ! conductive flux
!!$      & I_0_eff           ! fraction of the net incoming solar radiation which
!!$                         ! is absorbed in the ice

      
      
!!$    REAL(wp), PARAMETER :: &
!!$      & one  = 1.0_wp , & ! Bottom temperature [K]
!!$      & four = 4.0_wp   ! 4

    INTEGER :: i,j,k ! loop indices

    ! initialize surf. temp. in Kelvin
    TsurfK(:,:,:) = ice%Tsurf(:,:,:)+tmelt
 
    ! Create array of incoming SW/LW radiation split up into ice categories
    ! (purely for computational reasons)
    FORALL(i=1:nproma, j=1:ppatch%nblks_c, k=1:ice%kice, ice % isice (i,k,j)) &
      & SWin3d(i,k,j) = Qatm% SWin(i,j)
    FORALL(i=1:nproma, j=1:ppatch%nblks_c, k=1:ice%kice, ice % isice (i,k,j)) &
      & LWin3d(i,k,j) = Qatm% LWin(i,j)

    WHERE ( ice%isice(:,:,:) )
      ! combined heat conductivity for snow and ice
      k_effective(:,:,:) = ki * ks / (ks * ice%hi(:,:,:) + ki * ice%hs(:,:,:))
      ! one factor in the linearization of Stefan-Boltzmann around the previous Tsurf
      sigma_T_p3 (:,:,:) = StBo *TsurfK(:,:,:) ** 3
      
      !
      ! calculate temperature increment
      !

      ! calculate net upwelling het flux and conductive heat flux through the ice
    END WHERE
    CALL ice_flux_balance_zerolayer (ppatch,ice,Qatm,SWin3d,LWin3d,TsurfK,k_effective,&
      &                                    F_A,F_S)
    WHERE ( ice%isice(:,:,:) )
      deltaT(:,:,:) = (F_S - F_A)&
        & / ( k_effective(:,:,:) + 4.0_wp*sigma_T_p3(:,:,:) )
      
      ! update surface temperature
      ice%Tsurf(:,:,:) = ice%Tsurf(:,:,:) + deltaT(:,:,:)
      ! in case of positive surface temperature, reset to zero
      WHERE ( ice%Tsurf(:,:,:) > 0.0_wp )
        ice%Tsurf(:,:,:) = 0.0_wp
      END WHERE
      
      ! update surface temp. in Kelvin
      TsurfK(:,:,:) = ice%Tsurf(:,:,:)+tmelt
    END WHERE
    !
    CALL ice_flux_balance_zerolayer (ppatch,ice,Qatm,SWin3d,LWin3d,TsurfK,k_effective,&
      &                                    F_A,F_S)
    
    WHERE ( ice%isice(:,:,:) )
      ! flux balance 
      ice%Qtop(:,:,:) = - F_A(:,:,:) + F_S(:,:,:)
      ice%Qbot(:,:,:) = - F_S ! flux from ocean to ice still missing, this is done in
      !                       ice_growth_zerolayer
    END WHERE
    
    ipl_src=1  ! output print level (1-5, fix)
    CALL print_mxmn('ice%Tsurf',1,ice%Tsurf(:,1,:),1,ppatch%nblks_c,'ice',ipl_src)
    
  END SUBROUTINE set_ice_temp_zerolayer



 !
 ! The counterpart to the  ice_growth subroutnine
 !
 
 SUBROUTINE ice_growth_zerolayer(ppatch, p_os, ice, rpreci)
   TYPE(t_patch),             INTENT(IN)    :: ppatch 
   TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os
   TYPE (t_sea_ice),          INTENT(INOUT) :: ice
   REAL(wp),                  INTENT(IN)    :: rpreci(:,:) 
                                   ! water equiv. solid precipitation rate [m/s] DIMENSION (ie,je)

   !!Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, ppatch%nblks_c) ::         &
      & zHeatOceI,   & ! Oceanic heat flux                                  [W/m^2]
      & Tfw            ! Ocean freezing temperature [°C]
    
    INTEGER :: k ! Loop index for thickness categories


    ! freezing temperature of uppermost sea water
    IF ( no_tracer >= 2 ) then
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
   
    WHERE ( ice%isice(:,:,:) )
      ! surface ablation
      WHERE ( ice%Qtop(:,:,:) > 0.0_wp ) 
        WHERE ( ice%hs(:,:,:) > 0.0_wp )  ! melt snow where there's snow
          ice%hs (:,:,:) =  ice%hs(:,:,:) - ice%Qtop(:,:,:) * dtime / (alf*rhos) 
          ! put remaining heat, if any, into melting ice below
          WHERE (ice%hs(:,:,:) < 0.0_wp) 
            ice % hi(:,:,:) = ice%hi(:,:,:) + ice%hs(:,:,:) * (rhos/rhoi) ! snow thickness loss
            !                                                               in ice equivalents
            ice % hs(:,:,:) = 0.0_wp
          END WHERE
          
        ELSEWHERE   ! where there's no snow
          ice % hi(:,:,:) = ice%hi(:,:,:) -  ice%Qtop(:,:,:) * dtime / (alf*rhoi) 
        END WHERE
      END WHERE
      
      ! bottom melt/freeze
      ice % hi(:,:,:) = ice%hi(:,:,:) + ice%Qbot(:,:,:) * dtime / (alf * rhoi )
      
      
      WHERE (ice%hi(:,:,:) <= 0.0_wp) 
        ice%Tsurf(:,:,:) =  Tfw(:,:,:)
        ice%isice(:,:,:) =  .FALSE.
        ice%conc (:,:,:) = 0.0_wp
        ice%hi   (:,:,:) = 0.0_wp
      ELSEWHERE
        ice%heatOceI(:,:,:) = ice%heatOceI(:,:,:) - zHeatOceI(:,:,:)
      END WHERE
      
    END WHERE

    ipl_src=1  ! output print level (1-5, fix)
    CALL print_mxmn('ice%hi',1,ice%hi(:,1,:),1,ppatch%nblks_c,'ice',ipl_src)
    
  END SUBROUTINE ice_growth_zerolayer

END MODULE mo_sea_ice_zerolayer
