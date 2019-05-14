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
!! Modification by Jonas Straub, KIT (2017-02-08)
!! - Emission of pollen
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
  USE mo_run_config,                    ONLY: lart,ntracer,iforcing 
  USE mo_time_config,                   ONLY: time_config
  USE mo_impl_constants,                ONLY: iecham, inwp
  USE mo_echam_phy_memory,              ONLY: prm_field
  USE mtime,                            ONLY: datetime
  USE mo_util_mtime,                    ONLY: getElapsedSimTimeInSeconds
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_emissInt
#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_pollen, t_fields_volc
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_integration,               ONLY: art_integrate_explicit
! Emission Routines
  USE mo_art_emission_volc_1mom,        ONLY: art_organize_emission_volc
  USE mo_art_emission_volc_2mom,        ONLY: art_prepare_emission_volc,    &
                                          &   art_calculate_emission_volc
  USE mo_art_emission_seas,             ONLY: art_seas_emiss_martensson,    &
                                          &   art_seas_emiss_monahan,       &
                                          &   art_seas_emiss_smith
  USE mo_art_emission_dust,             ONLY: art_emission_dust,art_prepare_emission_dust
  USE mo_art_emission_dust_simple,      ONLY: art_prepare_emission_dust_simple
  USE mo_art_emission_chemtracer,       ONLY: art_emiss_chemtracer
  USE mo_art_emission_gasphase,         ONLY: art_emiss_gasphase
  USE mo_art_emission_pollen,           ONLY: art_emiss_pollen
  USE mo_art_emission_pntSrc,           ONLY: art_emission_pntSrc
  USE mo_art_read_emissions,            ONLY: art_add_emission_to_tracers,  &
                                          &   art_init_emission_struct,     &
                                          &   art_read_emissions
  USE mo_art_prescribed_state,          ONLY: art_prescribe_tracers
  USE omp_lib 
  USE mo_sync,                          ONLY: sync_patch_array_mult, SYNC_C

#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_emission_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_interface(ext_data,p_patch,dtime,p_nh_state,prm_diag,p_diag_lnd,rho,current_date,nnow,tracer)
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
  TYPE(datetime), INTENT(in), POINTER :: &
    &  current_date            !< Date and time information
  INTEGER, INTENT(in)               :: &
    &  nnow                    !< Time level
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
    &  dz_3d(:,:,:),         & !< Height of model layer (3dimensional)
    &  dz(:,:)                 !< Height of model layer
  REAL(wp), POINTER :: &
    &  land_sea(:,:),  &
    &  u_10m(:,:),     &
    &  v_10m(:,:)
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: &
    &  this_mode               !< pointer to current aerosol mode
  REAL(wp)                :: &
    &  p_sim_time              !< elapsed simulation time on this grid level
  
  ! calculate elapsed simulation time in seconds (local time for
  ! this domain!)
  p_sim_time = getElapsedSimTimeInSeconds(current_date) 


  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_emissInt)


    IF (art_config(jg)%lart_pntSrc) THEN
      ! Point sources
      CALL art_emission_pntSrc(jg, p_art_data(jg)%pntSrc, current_date, dtime, rho,        &
        &                      p_patch%cells%area, p_nh_state%metrics%ddqz_z_full, tracer)
    ENDIF

    ALLOCATE(emiss_rate(nproma,nlev))
    ALLOCATE(dz(nproma,nlev))
    ALLOCATE(dz_3d(nproma,nlev,p_patch%nblks_c))

    IF (art_config(jg)%lart_aerosol .OR. art_config(jg)%lart_chem &
        .OR. art_config(jg)%lart_passive) THEN

      CALL art_prescribe_tracers(tracer, p_art_data(jg)%prescr_list,     &
               &                 p_patch, current_date, p_nh_state%diag, &
               &                 p_nh_state%metrics%z_mc,                &
               &                 i_startblk, i_endblk, i_rlstart, i_rlend)

      IF (p_art_data(jg)%emiss%is_init) THEN
        IF (iforcing == inwp) THEN
          CALL art_add_emission_to_tracers(tracer,p_art_data(jg)%emiss,p_patch,p_nh_state%metrics, &
                                      &  p_nh_state%diag%temp,p_nh_state%diag%pres,dtime,        &
                                      &  current_date,prm_diag%swflx_par_sfc(:,:))
        ELSE IF (iforcing == iecham) THEN
          CALL art_add_emission_to_tracers(tracer,p_art_data(jg)%emiss,p_patch,p_nh_state%metrics, &
                                      &  p_nh_state%diag%temp,p_nh_state%diag%pres,dtime,        &
                                      &  current_date)
        ENDIF
      ENDIF
    ENDIF

    IF (art_config(jg)%lart_aerosol) THEN
!$omp parallel do default (shared) private(jb, istart, iend, dz)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        CALL art_air_properties(p_nh_state%diag%pres(:,:,jb),p_nh_state%diag%temp(:,:,jb), &
          &                     istart,iend,1,nlev,jb,p_art_data(jg))
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
              &          dzsoil(1),p_diag_lnd%h_snow(:,jb),jb,istart,iend,p_art_data(jg)%ext%soil_prop,         &
              &          p_art_data(jg)%diag%ustar_threshold(:,jb),p_art_data(jg)%diag%ustar(:,jb))
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
            ! No preparations necessary
          CASE(2)
            CALL art_prepare_emission_volc(current_date,jb,nlev,p_nh_state%metrics%z_ifc(:,:,jb),  &
              &                            p_art_data(jg)%dict_tracer, p_art_data(jg)%ext%volc_data)
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
                      &             jb,istart,iend,'dusta',p_art_data(jg)%ext%soil_prop,emiss_rate(:,nlev))
                    ! Compute emission rate [mug m-2 s-1] and collect for output
                    p_art_data(jg)%diag%emiss_rate_dusta(:,jb) = emiss_rate(:,nlev)*dz(:,nlev)
                    ! Accumulate emission rate 
                    p_art_data(jg)%diag%acc_emiss_dusta(:,jb) = p_art_data(jg)%diag%acc_emiss_dusta(:,jb) &
                                                              + dtime * p_art_data(jg)%diag%emiss_rate_dusta(:,jb)
                  CASE ('dustb')
                    CALL art_emission_dust(dz(:,nlev),                                                 &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),   &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil),  &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),     &
                      &             jb,istart,iend,'dustb',p_art_data(jg)%ext%soil_prop,emiss_rate(:,nlev))
                    ! Compute emission rate [mug m-2 s-1] and collect for output
                    p_art_data(jg)%diag%emiss_rate_dustb(:,jb) = emiss_rate(:,nlev)*dz(:,nlev)
                    ! Accumulate emission rate 
                    p_art_data(jg)%diag%acc_emiss_dustb(:,jb) = p_art_data(jg)%diag%acc_emiss_dustb(:,jb) &
                                                              + dtime * p_art_data(jg)%diag%emiss_rate_dustb(:,jb)
                  CASE ('dustc')
                    CALL art_emission_dust(dz(:,nlev),                                                 &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),   &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),      &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil),  &
                      &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),     &
                      &             jb,istart,iend,'dustc',p_art_data(jg)%ext%soil_prop,emiss_rate(:,nlev))
                    ! Compute emission rate [mug m-2 s-1] and collect for output
                    p_art_data(jg)%diag%emiss_rate_dustc(:,jb) = emiss_rate(:,nlev)*dz(:,nlev)
                    ! Accumulate emission rate 
                    p_art_data(jg)%diag%acc_emiss_dustc(:,jb) = p_art_data(jg)%diag%acc_emiss_dustc(:,jb) &
                                                              + dtime * p_art_data(jg)%diag%emiss_rate_dustc(:,jb)
                  CASE ('asha')
                    CALL art_calculate_emission_volc( jb, p_nh_state%metrics%ddqz_z_full(:,:,jb),      &
                      &             p_patch%cells%area(:,jb), nlev, p_art_data(jg)%ext%volc_data,      &
                      &             fields%itr3(1), emiss_rate(:,:) ) !< itr3(1) assumes only 1 mass component of mode
                  CASE ('ashb')
                    CALL art_calculate_emission_volc( jb, p_nh_state%metrics%ddqz_z_full(:,:,jb),      &
                      &             p_patch%cells%area(:,jb), nlev, p_art_data(jg)%ext%volc_data,      &
                      &             fields%itr3(1), emiss_rate(:,:) ) !< itr3(1) assumes only 1 mass component of mode
                  CASE ('ashc')
                    CALL art_calculate_emission_volc( jb, p_nh_state%metrics%ddqz_z_full(:,:,jb),      &
                      &             p_patch%cells%area(:,jb), nlev, p_art_data(jg)%ext%volc_data,      &
                      &             fields%itr3(1), emiss_rate(:,:) ) !< itr3(1) assumes only 1 mass component of mode
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
            TYPE is (t_fields_pollen)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                  &                istart, iend, i_rlstart, i_rlend)

                DO jc = istart, iend
                  dz(jc,nlev) = p_nh_state%metrics%ddqz_z_full(jc,nlev,jb)
                ENDDO

                SELECT CASE(TRIM(fields%name))
                  CASE('pollbetu') ! Betula --> birch
                    CALL art_emiss_pollen(dtime,current_date,                             &
                      &                   rho(:,nlev,jb),                                 &
                      &                   p_art_data(jg)%ext%pollen_prop%pollen_type(1),  & !< 1=pollbetu
                      &                   tracer(:,nlev,jb,:),                            &
                      &                   p_art_data(jg)%dict_tracer,                     &
                      &                   p_nh_state%diag%temp(:,nlev,jb),                &
                      &                   p_nh_state%diag%pres_sfc(:,jb),                 &
                      &                   p_nh_state%prog(nnow)%tke(:,nlev,jb),           &
                      &                   prm_diag%rain_gsp_rate(:,jb),                   &
                      &                   prm_diag%rain_con_rate(:,jb),                   &
                      &                   prm_diag%rh_2m(:,jb),                           &
                      &                   dz(:,nlev),                                     &
                      &                   ext_data%atm%llsm_atm_c(:,jb),                  &
                      &                   jb, istart, iend )
                  CASE DEFAULT
                    CALL finish('mo_art_emission_interface:art_emission_interface', &
                      &         'ART: Pollen emissions for '//TRIM(fields%name)//' not yet implemented')
                END SELECT
              ENDDO
            CLASS is (t_fields_radio)
              ! handled via pntSrc structure
            CLASS is (t_fields_volc)
              ! handled below
            CLASS default
              CALL finish('mo_art_emission_interface:art_emission_interface', &
                &         'ART: Unknown mode field type')
          END SELECT
          this_mode => this_mode%next_mode
        ENDDO !while(associated)

      DEALLOCATE(emiss_rate)

      ! volcano emissions
      IF (art_config(jg)%iart_volcano == 1) THEN
        CALL art_organize_emission_volc(p_patch, current_date, dtime,rho,p_art_data(jg)%dict_tracer, &
          &                             p_art_data(jg)%ext%volc_data,tracer)
      ENDIF


    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- emissions of chemical tracer
    ! ----------------------------------
  
    IF ((art_config(jg)%lart_chem) .OR. (art_config(jg)%lart_passive)) THEN
      ! Get model layer heights
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        DO jk = 1, nlev
          DO jc = istart, iend
            dz_3d(jc,jk,jb) = p_nh_state%metrics%z_ifc(jc,jk,jb) - p_nh_state%metrics%z_ifc(jc,jk+1,jb)
          ENDDO
        ENDDO
      END DO

      IF (iforcing == iecham) THEN
        land_sea => prm_field(jg)%lsmask
        u_10m => prm_field(jg)%uas
        v_10m => prm_field(jg)%vas
      ELSE IF (iforcing == inwp) THEN
        land_sea => ext_data%atm%fr_land
        u_10m => prm_diag%u_10m
        v_10m => prm_diag%v_10m
      ELSE
         CALL finish('mo_art_emission_interface:art_emission_interface', &
              &         'ART: no land sea mask for this iforcing available')
      ENDIF

        
      SELECT CASE(art_config(jg)%iart_chem_mechanism)
        CASE(0)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              &                istart, iend, i_rlstart, i_rlend)
            CALL art_emiss_chemtracer(current_date,                   &
              &                       dtime,                          &
              &                       tracer,                         &
              &                       p_nh_state%diag%pres,           &
              &                       p_nh_state%diag%temp,           &
              &                       land_sea,                       &
              &                       p_patch,                        &
              &                       p_art_data(jg)%dict_tracer,     &
              &                       jb,istart,iend,nlev,nproma,     &
              &                       u_10m, v_10m, dz_3d(:,:,jb))
          ENDDO
        CASE(1)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              &                istart, iend, i_rlstart, i_rlend)
            CALL art_emiss_chemtracer(current_date,                   &
              &                       dtime,                          &
              &                       tracer,                         &
              &                       p_nh_state%diag%pres,           &
              &                       p_nh_state%diag%temp,           &
              &                       land_sea,                       &
              &                       p_patch,                        &
              &                       p_art_data(jg)%dict_tracer,     &
              &                       jb,istart,iend,nlev,nproma,     &
              &                       u_10m, v_10m, dz_3d(:,:,jb))
          ENDDO
        CASE(2)
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
              &                istart, iend, i_rlstart, i_rlend)
            CALL art_emiss_gasphase(current_date,                   &
              &                     dtime,                          &
              &                     tracer,                         &
              &                     p_nh_state%diag%pres,           &
              &                     p_nh_state%diag%temp,           &
              &                     land_sea,                       &
              &                     p_patch,                        &
              &                     jb,istart,iend,nlev,nproma,     &
              &                     u_10m, v_10m, dz_3d(:,:,jb))
          ENDDO
        CASE DEFAULT
          CALL finish('mo_art_emission_interface:art_emission_interface', &
               &      'ART: Unknown iart_chem_mechanism')
      END SELECT !iart_chem_mechanism
    ENDIF !lart_chem
    
    DEALLOCATE(dz)
    DEALLOCATE(dz_3d)

    CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer,  f4din=tracer(:,:,:,:))

    IF (timers_level > 3) CALL timer_stop(timer_art_emissInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)

  ENDIF !lart
       
#endif
END SUBROUTINE art_emission_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_interface


