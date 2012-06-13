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

!!$  USE mo_kind,                ONLY: wp
!!$  USE mo_parallel_config,     ONLY: nproma
!!$  USE mo_run_config,          ONLY: dtime
!!$  USE mo_dynamics_config,     ONLY: nold
!!$  USE mo_model_domain,        ONLY: t_patch
!!$  !
!!$  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 
!!$
!!$!  USE mo_exception,           ONLY: finish, message
!!$!  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell, sea_boundary 
!!$  USE mo_loopindices,         ONLY: get_indices_c
!!$!  USE mo_math_utilities,      ONLY: t_cartesian_coordinates
!!$  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs,&
!!$    &                               mu,mus,ci, alf, I_0, alv, albedoW, clw,            &
!!$    &                               cpd, zemiss_def,rd, stbo,tmelt   
!!$!  USE mo_math_constants,      ONLY: rad2deg
!!$  USE mo_ocean_nml,           ONLY: no_tracer, init_oce_prog, iforc_oce, &
!!$    &                               FORCING_FROM_FILE_FLUX, i_sea_ice
!!$  USE mo_util_dbg_prnt,       ONLY: dbg_print
!!$  USE mo_oce_state,           ONLY: t_hydro_ocean_state, v_base, ocean_var_list
!!$  USE mo_var_list,            ONLY: add_var
!!$!  USE mo_master_control,      ONLY: is_restart_run
!!$!  USE mo_cf_convention
!!$!  USE mo_grib2
!!$!  USE mo_cdi_constants
!!$  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
!!$    &                               t_atmos_for_ocean
!!$  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch
!  USE mo_exception,           ONLY: finish, message
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref,ki,ks,Tf,albi,albim,albsm,albs,&
    &                               mu,mus,ci, alf, I_0, alv, albedoW, clw,            &
    &                               cpd, zemiss_def,rd, stbo,tmelt   
!  USE mo_math_constants,      ONLY: rad2deg
  USE mo_ocean_nml,           ONLY: no_tracer !, &
!    &                               init_oce_prog, iforc_oce, &
!    &                               FORCING_FROM_FILE_FLUX, i_sea_ice
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_oce_state,           ONLY: t_hydro_ocean_state !, v_base, ocean_var_list
!  USE mo_var_list,            ONLY: add_var
!  USE mo_master_control,      ONLY: is_restart_run
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_sfc_flx, t_atmos_fluxes, &
    &                               t_atmos_for_ocean
  USE mo_sea_ice_shared_sr,   ONLY: oce_ice_heatflx
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range 


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
  !!    379–389.  doi:
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

  SUBROUTINE set_ice_temp_zerolayer(p_patch,ice, Tfw, Qatm) 
    TYPE(t_patch),        INTENT(IN),TARGET    :: p_patch 
    TYPE(t_sea_ice),      INTENT(INOUT)        :: ice
    REAL(wp),             INTENT(IN)           :: Tfw(:,:,:)
    TYPE(t_atmos_fluxes), INTENT(IN)           :: Qatm

    ! Local variables
    REAL(wp), DIMENSION (nproma,ice%kice, p_patch%nblks_c) ::          &
      & k_effective ,  &  ! total heat conductivity of ice/snow
      & deltaT      ,  &  ! temperature increment 
      & F_A         ,  &  ! atmospheric net flux, positive=upward
      & F_S             ! conductive flux, positive=upward
    
    TYPE(t_subset_range), POINTER :: all_cells

    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c     ! loop indices
    
    ! initialization
    all_cells => p_patch%cells%all 
    k_effective(:,:,:) = 0.0_wp
    deltaT     (:,:,:) = 0.0_wp
    F_A        (:,:,:) = 0.0_wp
    F_S        (:,:,:) = 0.0_wp   

    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%isice(jc,k,jb)) THEN
            
            ! total heat conductivity for the ice-snow system
            k_effective(jc,k,jb) = ki*ks/(ks*ice%hi(jc,k,jb) + ki*ice%hs(jc,k,jb))
                    
            ! F_A, F_S : pos=upward flux
            ! flux ice-atmosphere
            !            F_A(jc,k,jb) = zemiss_def * StBo * (ice%Tsurf(jc,k,jb) + tmelt)**4&
            F_A(jc,k,jb) = Qatm%LWout(jc,k,jb) &
              &            - (1.0_wp - ice%alb(jc,k,jb)) * Qatm%SWin(jc,jb) &
              &            - Qatm%LWin(jc,jb) - Qatm%sens(jc,k,jb) - Qatm%lat(jc,k,jb)
            
            ! conductive heat flux through the ice
            F_S(jc,k,jb) = k_effective(jc,k,jb) * (Tf - ice%Tsurf(jc,k,jb))  

            ! temperature increment
            deltaT(jc,k,jb) = (F_S(jc,k,jb) - F_A(jc,k,jb))&
              &               / (k_effective(jc,k,jb) &
!              &                  + 4.0_wp*zemiss_def * StBo * &
!              &                    (ice%Tsurf(jc,k,jb) + tmelt)**3 &
              &                  + Qatm%  dLWdT(jc,k,jb) &
              &                  - Qatm%dsensdT(jc,k,jb))


            ! ice temperatures over 0 deg C impossible:
            IF (ice%Tsurf(jc,k,jb) + deltaT(jc,k,jb) > 0.0_wp) THEN  ! if new temperature would be over 0 deg C
              deltaT(jc,k,jb) = -ice%Tsurf(jc,k,jb) 
              ice% Tsurf(jc,k,jb) = 0.0_wp
              
              ! pos. flux means into uppermost ice layer
              ice%Qtop(jc,k,jb) = - F_A(jc,k,jb) + F_S(jc,k,jb) &
                &               - deltaT(jc,k,jb) * (k_effective(jc,k,jb) &
                !              &                                    + 4.0_wp*zemiss_def * StBo *       &
                !              &                                      (ice%Tsurf(jc,k,jb) + tmelt)**3  &
                &                                    + Qatm%  dLWdT(jc,k,jb) &
                &                                    - Qatm%dsensdT(jc,k,jb))

              ! pos. flux means into lowest ice layer
              ice%Qbot(jc,k,jb) = - F_S(jc,k,jb) -  deltaT(jc,k,jb) * k_effective(jc,k,jb)
              ! NB: flux from ocean to ice still missing, this is done in
              !                       ice_growth_zerolayer


            ELSE   ! if new temperature would be under 0 deg C, then we can achieve F_A=F_S just as we wanted
              ! new surface temperature
              ice%Tsurf(jc,k,jb) = ice%Tsurf(jc,k,jb) + deltaT(jc,k,jb)
              
              ! surface flux balanced:
              ice%Qtop(jc,k,jb) = 0.0_wp
              
              ! pos. flux means into lowest ice layer
              ice%Qbot(jc,k,jb) =  k_effective(jc,k,jb) * (ice%Tsurf(jc,k,jb) - Tf)
              !!! ice%Qbot = - new F_S
              ! NB: flux from ocean to ice still missing, this is done in
              !                       ice_growth_zerolayer
              
            END IF

          END IF
        END DO
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
  !!                in previously ice covered areas if all ice is gone      [kg/m�]
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
      & Tfw,         & ! Ocean freezing temperature [°C]
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
    
    ! Heat flux from ocean into ice
    CALL oce_ice_heatflx (p_os,ice,Tfw,zHeatOceI)
    ! Add oceanic heat flux to energy available at the bottom of the ice.
    ice%Qbot(:,:,:) = ice%Qbot(:,:,:) + zHeatOceI(:,:,:)
    
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c) 
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
          IF (ice%isice(jc,k,jb)) THEN

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
              ! put surplus energy from ice thickness into water
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) &
                &                     - ice%hi(jc,k,jb)*alf*rhoi/dtime
              ! put snow energy into water
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) &
                &                     + ice%hs(jc,k,jb)*alf*rhos/dtime
              
              ! 
              ice%Tsurf(jc,k,jb) =  Tfw(jc,k,jb)
              ice%isice(jc,k,jb) =  .FALSE.
              ice%conc (jc,k,jb) = 0.0_wp
              ice%hi   (jc,k,jb) = 0.0_wp
              ice%hs   (jc,k,jb) = 0.0_wp
            ELSE
              ! heat to remove from water
              ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb) - zHeatOceI(jc,k,jb)
            ENDIF
            
            ! check energy conservation
            ! surplus energy = entering - leaving - latent heat
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
