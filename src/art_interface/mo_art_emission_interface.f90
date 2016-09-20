!>
!! Provides interface to ART-routines dealing with emissions
!!
!! This module provides an interface to the ART emission routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!! @author Daniel Reinert, DWD
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2012-01-27)
!! Modification by Kristina Lundgren, KIT (2012-01-30)
!! - Modification for dealing with the ART-routine emission_volc.
!! Rewritten by Daniel Rieger, KIT (2013-09-30)
!! - Usage of the generalized ART infrastructure
!! Modification by Daniel Rieger, KIT (2014-11-12)
!! - First adaptions for unified use of parameterizations within
!!   ICON-ART and COSMO-ART. Changes performed for: Mineral dust, Sea salt
!! Modification by Roland Ruhnke, Jennifer Schroeter, KIT (2014-11-27)
!! - Chemical tracer section
!! Modification by Michael Weimer, KIT (2015-08-06)
!! - Emission of chemical tracer
!! Modification by Jonas Straub, Daniel Rieger KIT (2015-09-15)
!! - Included OpenMP statements
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_emission_interface

  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int,dzsoil
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_parallel_config,               ONLY: nproma
  USE mo_exception,                     ONLY: finish
  USE mo_nonhydro_types,                ONLY: t_nh_state
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_nwp_lnd_types,                 ONLY: t_lnd_diag
  USE mo_run_config,                    ONLY: lart
  USE mo_datetime,                      ONLY: t_datetime
#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc,t_fields_1mom
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_config,                    ONLY: art_config,                   &
                                          &   iCS137,iI131,iTE132,          &
                                          &   iZR95,iXE133,iI131g,          &
                                          &   iI131o,iBA140,iRU103,         &
                                          &   iasha, iashb, iashc
  USE mo_art_integration,               ONLY: art_integrate_explicit
! Emission Routines
  USE mo_art_emission_volc_1mom,        ONLY: art_organize_emission_volc
  USE mo_art_emission_volc_2mom,        ONLY: art_prepare_emission_volc,    &
                                          &   art_calculate_emission_volc
  USE mo_art_emission_radioact,         ONLY: art_emiss_radioact
  USE mo_art_emission_seas,             ONLY: art_seas_emiss_martensson,    &
                                          &   art_seas_emiss_monahan,       &
                                          &   art_seas_emiss_smith
  USE mo_art_emission_dust,             ONLY: art_emission_dust,art_prepare_emission_dust
  USE mo_art_emission_dust_simple,      ONLY: art_prepare_emission_dust_simple
  USE mo_art_emission_chemtracer,       ONLY: art_emiss_chemtracer
  USE mo_art_emission_gasphase,         ONLY: art_emiss_gasphase
  USE omp_lib 
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_emission_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_interface(ext_data,p_patch,dtime,p_nh_state,prm_diag,p_diag_lnd,rho,datetime,tracer)
  !! Interface for ART: Emissions
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-01-27)
  !! Modification by Kristina Lundgren, KIT (2012-01-30)
  !! Rewritten by Daniel Rieger, KIT (2013-09-30)
  TYPE(t_external_data), INTENT(in) ::  &
    &  ext_data                !< Atmosphere external data
  TYPE(t_patch), TARGET, INTENT(in) ::  & 
    &  p_patch                 !< Patch on which computation is performed
  REAL(wp), INTENT(in)    :: &
    &  dtime                   !< Time step (advection)
  TYPE(t_nh_state),INTENT(in)       :: &
    &  p_nh_state              !< State variables (prognostic, diagnostic, metrics)
  TYPE(t_nwp_phy_diag), INTENT(in)  :: &
    &  prm_diag                !< List of diagnostic fields (physics)
  TYPE(t_lnd_diag), INTENT(in)      :: &
    &  p_diag_lnd              !< List of diagnostic fields (land)
  REAL(wp), INTENT(inout) :: &
    &  rho(:,:,:)              !< Density of air [kg/m3]
  TYPE(t_datetime), INTENT(IN) :: &
    &  datetime                !< Date and time information
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:,:)         !< Tracer mixing ratios [kg kg-1]
  ! Local variables
  INTEGER                 :: & 
    &  jg, jb, ijsp, jk, jc, & !< Patch id, counter for block loop, jsp loop, vertical loop
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,   & !< Relaxation start and end
    &  i_nchdom,             & !< Number of child domains
    &  nlev                    !< Number of levels (equals index of lowest full level)
  REAL(wp),ALLOCATABLE    :: &
    &  emiss_rate(:,:),      & !< Emission rates [UNIT m-3 s-1], UNIT might be mug, kg or just a number
    &  dz(:,:)                 !< Height of model layer
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: &
    &  this_mode               !< pointer to current aerosol mode

  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

  IF (lart) THEN
    
    ALLOCATE(emiss_rate(nproma,nlev))
    ALLOCATE(dz(nproma,nlev))
  
    CALL art_air_properties(p_patch,p_art_data(jg))
       
    IF (art_config(jg)%lart_aerosol) THEN

!$omp parallel do default (shared) private(jb, istart, iend, dz)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        
        ! Get model layer heights
        DO jk = 1, nlev
          DO jc = istart, iend
            dz(jc,jk) = p_nh_state%metrics%ddqz_z_full(jc,jk,jb)
          ENDDO
        ENDDO
        
        ! ----------------------------------
        ! --- Preparations for emission routines
        ! ----------------------------------
        SELECT CASE(art_config(jg)%iart_dust)
          CASE(0)
            ! Nothing to do, no dust emissions
          CASE(1)
            CALL art_prepare_emission_dust(p_nh_state%diag%u(:,nlev,jb),p_nh_state%diag%v(:,nlev,jb),           &
              &          rho(:,nlev,jb),prm_diag%tcm(:,jb),p_diag_lnd%w_so(:,1,jb),p_diag_lnd%w_so_ice(:,1,jb), &
              &          dzsoil(1),p_diag_lnd%h_snow(:,jb),jb,istart,iend,p_art_data(jg)%soil_prop)
          CASE(2)
            ! not available yet
          CASE default
            CALL finish('mo_art_emission_interface:art_emission_interface', &
              &         'ART: Unknown dust emissions configuration')
        END SELECT
        SELECT CASE(art_config(jg)%iart_volcano)
          CASE(0)
            ! Nothing to do, no volcano emissions
          CASE(1)
            ! bulk emissions see below
          CASE(2)
            CALL art_prepare_emission_volc(jb,nlev,p_nh_state%metrics%z_ifc(:,:,jb),p_art_data(jg)%volc_data)
          CASE default
            CALL finish('mo_art_emission_interface:art_emission_interface', &
                 &      'ART: Unknown volc emissions configuration')
        END SELECT
      ENDDO !jb   
!$omp end parallel do


        ! ----------------------------------
        ! --- Call the emission routines
        ! ----------------------------------
        
        this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
      
        DO WHILE(ASSOCIATED(this_mode))
          ! Check how many moments the mode has
          SELECT TYPE (fields=>this_mode%fields)
            TYPE is (t_fields_2mom)
              ! Now the according emission routine has to be found
!$omp parallel do default (shared) private (jb, istart, iend, emiss_rate, dz)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                  &                istart, iend, i_rlstart, i_rlend)
                  
                emiss_rate(:,:) = 0._wp
                
                ! Get model layer heights
                DO jk = 1, nlev
                  DO jc = istart, iend
                    dz(jc,jk) = p_nh_state%metrics%ddqz_z_full(jc,jk,jb)
                  ENDDO
                ENDDO

                SELECT CASE(TRIM(fields%name))
                  CASE ('seasa')
                    CALL art_seas_emiss_martensson(prm_diag%u_10m(:,jb), prm_diag%v_10m(:,jb),                        &
                      &             dz(:,nlev), p_diag_lnd%t_s(:,jb),                                                 &
                      &             ext_data%atm%fr_land(:,jb),p_diag_lnd%fr_seaice(:,jb),ext_data%atm%fr_lake(:,jb), &
                      &             istart,iend,emiss_rate(:,nlev))
                  CASE ('seasb')
                    CALL art_seas_emiss_monahan(prm_diag%u_10m(:,jb), prm_diag%v_10m(:,jb),                           &
                      &             dz(:,nlev), ext_data%atm%fr_land(:,jb),                                           &
                      &             p_diag_lnd%fr_seaice(:,jb),ext_data%atm%fr_lake(:,jb), istart,iend,emiss_rate(:,nlev))
                  CASE ('seasc')
                    CALL art_seas_emiss_smith(prm_diag%u_10m(:,jb), prm_diag%v_10m(:,jb),                             &
                      &             dz(:,nlev), ext_data%atm%fr_land(:,jb),                                           &
                      &             p_diag_lnd%fr_seaice(:,jb),ext_data%atm%fr_lake(:,jb), istart,iend,emiss_rate(:,nlev))
                  CASE ('dusta')
                    CALL art_emission_dust(dz(:,nlev),                                                 &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),   &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil),  &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),     &
                      &             jb,istart,iend,'dusta',p_art_data(jg)%soil_prop,emiss_rate(:,nlev))
                  CASE ('dustb')
                    CALL art_emission_dust(dz(:,nlev),                                                 &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),   &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil),  &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),     &
                      &             jb,istart,iend,'dustb',p_art_data(jg)%soil_prop,emiss_rate(:,nlev))
                  CASE ('dustc')
                    CALL art_emission_dust(dz(:,nlev),                                                 &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),   &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil),  &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),     &
                      &             jb,istart,iend,'dustc',p_art_data(jg)%soil_prop,emiss_rate(:,nlev))
                  CASE ('asha')
                    CALL art_calculate_emission_volc( jb, p_nh_state%metrics%ddqz_z_full(:,:,jb),      &
                      &             p_patch%cells%area(:,jb), nlev, p_art_data(jg)%volc_data,          &
                      &             iasha, emiss_rate(:,:) )
                  CASE ('ashb')
                    CALL art_calculate_emission_volc( jb, p_nh_state%metrics%ddqz_z_full(:,:,jb),      &
                      &             p_patch%cells%area(:,jb), nlev, p_art_data(jg)%volc_data,          &
                      &             iashb, emiss_rate(:,:) )
                  CASE ('ashc')
                    CALL art_calculate_emission_volc( jb, p_nh_state%metrics%ddqz_z_full(:,:,jb),      &
                      &             p_patch%cells%area(:,jb), nlev, p_art_data(jg)%volc_data,          &
                      &             iashc, emiss_rate(:,:) )
                END SELECT
                
                ! Update mass mixing ratios
                DO ijsp = 1, fields%ntr-1
                  CALL art_integrate_explicit(tracer(:,:,jb,fields%itr3(ijsp)),  emiss_rate(:,:),      &
                    &                         dtime, istart, iend, nlev, opt_rho = rho(:,:,jb))
                ENDDO
                ! Update mass-specific number
                CALL art_integrate_explicit(tracer(:,:,jb,fields%itr0), emiss_rate(:,:), dtime,        &
                  &                         istart, iend, nlev, opt_rho = rho(:,:,jb),                 &
                  &                         opt_fac=(fields%info%mode_fac * fields%info%factnum))
              ENDDO !jb
!$omp end parallel do
            CLASS is (t_fields_1mom)
              ! drieg: This needs to be done here instead of the version outside the jb loop below in the future
            CLASS default
              CALL finish('mo_art_emission_interface:art_emission_interface', &
                &         'ART: Unknown mode field type')
          END SELECT
          this_mode => this_mode%next_mode
        ENDDO !while(associated)
    
      DEALLOCATE(emiss_rate)
      DEALLOCATE(dz)
      
      ! START OLD BLOCK: Needs to be realized within jb loop above in the future
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE is (t_fields_2mom)
            ! Nothing to do here
          CLASS is (t_fields_radio)
            SELECT CASE(TRIM(fields%name))
              CASE ('Cs-137')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iCS137),373, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('I-131')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iI131),340, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('Te-132')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iTE132),325, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('Zr-95')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iZR95),184, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('Xe-133')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iXE133),355, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('I-131g')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iI131g),870, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('I-131o')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iI131o),880, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('Ba-140')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iBA140),384, & 
                  &                     p_art_data(jg)%radioact_data)
              CASE ('Ru-103')
                CALL art_emiss_radioact(p_patch,dtime,rho,tracer(:,:,:,iRU103),220, & 
                  &                     p_art_data(jg)%radioact_data)
              ! And Default...
              CASE default
                CALL finish('mo_art_emission_interface:art_emission_interface', &
                  &         'No according emission routine to mode')
            END SELECT
          CLASS is (t_fields_volc)
            ! nothing to do here, see below
          CLASS default
            CALL finish('mo_art_emission_interface:art_emission_interface', &
              &         'ART: Unknown mode field type')
        END SELECT
        this_mode => this_mode%next_mode
      ENDDO
  
      ! ----------------------------------
      ! --- volcano emissions
      ! ----------------------------------
    
      IF (art_config(jg)%iart_volcano == 1) THEN
        CALL art_organize_emission_volc(p_patch,dtime,rho,p_art_data(jg)%volc_data,tracer) 
      ENDIF
      ! END OLD BLOCK
    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- emissions of chemical tracer
    ! ----------------------------------
  
    IF (art_config(jg)%lart_chem) THEN
      SELECT CASE(art_config(jg)%iart_chem_mechanism)
        CASE(0)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              &                istart, iend, i_rlstart, i_rlend)
            
            CALL art_emiss_chemtracer(datetime,                       &
              &                       dtime,                          &
              &                       tracer,                         &
              &                       p_nh_state%diag%pres,           &
              &                       p_nh_state%diag%temp,           &
              &                       p_nh_state%metrics,             &
              &                       ext_data%atm%llsm_atm_c,        &
              &                       p_patch,                        &
              &                       jb,istart,iend,nlev,nproma,     &
              &                       prm_diag%swflx_par_sfc)
          ENDDO
        CASE(1)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              &                istart, iend, i_rlstart, i_rlend)
            
            CALL art_emiss_chemtracer(datetime,                       &
              &                       dtime,                          &
              &                       tracer,                         &
              &                       p_nh_state%diag%pres,           &
              &                       p_nh_state%diag%temp,           &
              &                       p_nh_state%metrics,             &
              &                       ext_data%atm%llsm_atm_c,        &
              &                       p_patch,                        &
              &                       jb,istart,iend,nlev,nproma,     &
              &                       prm_diag%swflx_par_sfc)
          ENDDO
        CASE(2)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              &                istart, iend, i_rlstart, i_rlend)
            
            CALL art_emiss_gasphase(datetime,                       &
              &                     dtime,                          &
              &                     tracer,                         &
              &                     p_nh_state%diag%pres,           &
              &                     p_nh_state%diag%temp,           &
              &                     p_nh_state%metrics,             &
              &                     ext_data%atm%llsm_atm_c,        &
              &                     p_patch,                        &
              &                     jb,istart,iend,nlev,nproma,     &
              &                     prm_diag%swflx_par_sfc)
          ENDDO
        CASE DEFAULT
          CALL finish('mo_art_emission_interface:art_emission_interface', &
               &      'ART: Unknown iart_chem_mechanism')
      END SELECT !iart_chem_mechanism
    ENDIF !lart_chem
  ENDIF !lart
       
#endif
END SUBROUTINE art_emission_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_interface


