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
  USE mo_run_config,                    ONLY: lart,                         &
                                          &   iCS137,iI131,iTE132,          &
                                          &   iZR95,iXE133,iI131g,          &
                                          &   iI131o,iBA140,iRU103
    
#ifdef __ICON_ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc,t_fields_1mom
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties
  USE mo_art_diagnostics_interface,     ONLY: art_diagnostics_interface
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_integration,               ONLY: art_integrate_explicit
! Emission Routines
  USE mo_art_emission_volc,             ONLY: art_organize_emission_volc
  USE mo_art_radioactive,               ONLY: art_emiss_radioact
  USE mo_art_emission_seas,             ONLY: art_seas_emiss_martensson, &
                                          &   art_seas_emiss_monahan, &
                                          &   art_seas_emiss_smith
  USE mo_art_emission_dust,             ONLY: art_emission_dust,art_prepare_emission_dust
  USE mo_art_emission_dust_simple,      ONLY: art_prepare_emission_dust_simple
  USE mo_art_chemtracer,                ONLY: art_emiss_chemtracer
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_emission_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_interface(ext_data,p_patch,dtime,p_nh_state,prm_diag,p_diag_lnd,rho,p_trac)
  !! Interface for ART: Emissions
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-01-27)
  !! Modification by Kristina Lundgren, KIT (2012-01-30)
  !! Rewritten by Daniel Rieger, KIT (2013-09-30)
  TYPE(t_external_data), INTENT(in) ::  &
    &  ext_data                !< atmosphere external data
  TYPE(t_patch), TARGET, INTENT(in) ::  & 
    &  p_patch                 !< patch on which computation is performed
  REAL(wp), INTENT(in)    :: &
    &  dtime                   !< time step
  TYPE(t_nh_state),INTENT(in)       :: &
    &  p_nh_state              !< State variables (prognostic, diagnostic, metrics)
  TYPE(t_nwp_phy_diag), INTENT(in)  :: &
    &  prm_diag                !< list of diagnostic fields (physics)
  TYPE(t_lnd_diag), INTENT(in)      :: &
    &  p_diag_lnd              !< list of diagnostic fields (land)
  REAL(wp), INTENT(inout) :: &
    &  rho(:,:,:)              !< density of air [kg/m3]
  REAL(wp), INTENT(inout) :: &
    &  p_trac(:,:,:,:)         !< tracer mixing ratios [kg kg-1]
  ! Local variables
  INTEGER                 :: & 
    &  jg, jb, ijsp, jk,     & !< patch id, counter for block loop, conuter for jsp loop, counter for vertical loop
    &  i_startblk, i_endblk, & !< Start and end of block loop
    &  istart, iend,         & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,   & !< 
    &  i_nchdom,             & !< 
    &  nlev                    !< Number of levels (equals index of lowest full level)
  REAL(wp),ALLOCATABLE    :: &
    &  emiss_rate(:),        & !< Emission rates [UNIT kg-1 s-1] or [UNIT m-3 s-1], UNIT might be mug, kg or just a number
    &  dz(:,:)                 !< Height of lowest layer
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode     !< pointer to current aerosol mode
  
  
  
  ! --- Get the loop indizes
  i_nchdom   = MAX(1,p_patch%n_childdom)
  jg         = p_patch%id
  nlev       = p_patch%nlev
  i_rlstart  = grf_bdywidth_c+1
  i_rlend    = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)


  IF (lart) THEN
  
    ALLOCATE(emiss_rate(nproma),dz(nproma,nlev))
  
    CALL art_air_properties(p_patch,p_art_data(jg))
       
    IF (art_config(jg)%lart_aerosol) THEN
    
!DEVSTART new emission interface structure for unified use of emission routines in ICON/COSMO
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                istart, iend, i_rlstart, i_rlend)
        
        DO jk = 1, nlev
          dz(:,jk) = ABS(p_nh_state%metrics%z_ifc(:,jk,jb) - p_nh_state%metrics%z_ifc(:,jk+1,jb))
        ENDDO
        
        ! Call the ART diagnostics
        CALL art_diagnostics_interface(rho(:,:,jb),p_trac(:,:,jb,:), dz(:,:), istart, iend, nlev, jb, &
          &                            art_config(jg), p_art_data(jg))
        
        ! ----------------------------------
        ! --- Preparations for emission routines
        ! ----------------------------------
        select case(art_config(jg)%iart_dust)
          case(0)
            ! Nothing to do, no dust emissions
          case(1)
            CALL art_prepare_emission_dust(p_nh_state%diag%u(:,nlev,jb),p_nh_state%diag%v(:,nlev,jb), &
              &          rho(:,nlev,jb),prm_diag%tcm(:,jb),p_diag_lnd%w_so(:,1,jb),                   &
              &          dzsoil(1),p_diag_lnd%h_snow(:,jb),jb,istart,iend,p_art_data(jg)%soil_prop)
          case(2)
            ! not available yet
!            call art_prepare_emission_dust_simple(...)
          case default
            CALL finish('mo_art_emission_interface:art_emission_interface', &
                 &      'ART: Unknown dust emissions configuration')
        end select
        
        ! ----------------------------------
        ! --- Call the emission routines
        ! ----------------------------------
        
        this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
      
        DO WHILE(ASSOCIATED(this_mode))
          emiss_rate(:) = 0._wp
          
          ! Check how many moments the mode has
          select type (fields=>this_mode%fields)
            type is (t_fields_2mom)
              ! Now the according emission routine has to be found
              select case(TRIM(fields%info%name))
                case ('seasa')
                  CALL art_seas_emiss_martensson(prm_diag%u_10m(:,jb),prm_diag%v_10m(:,jb),dz(:,nlev),p_diag_lnd%t_s(:,jb),&
                    &             ext_data%atm%fr_land(:,jb),p_diag_lnd%fr_seaice(:,jb),ext_data%atm%fr_lake(:,jb),   &
                    &             istart,iend,emiss_rate(:))
                case ('seasb')
                  CALL art_seas_emiss_monahan(prm_diag%u_10m(:,jb),prm_diag%v_10m(:,jb),dz(:,nlev),                   &
                    &             ext_data%atm%fr_land(:,jb),p_diag_lnd%fr_seaice(:,jb),ext_data%atm%fr_lake(:,jb),   &
                    &             istart,iend,emiss_rate(:))
                case ('seasc')
                  CALL art_seas_emiss_smith(prm_diag%u_10m(:,jb),prm_diag%v_10m(:,jb),dz(:,nlev),                     &
                    &             ext_data%atm%fr_land(:,jb),p_diag_lnd%fr_seaice(:,jb),ext_data%atm%fr_lake(:,jb),   &
                    &             istart,iend,emiss_rate(:))
                case ('dusta')
                  CALL art_emission_dust(dz(:,nlev),ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil), &
                    &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),                      &
                    &             jb,istart,iend,'dusta',p_art_data(jg)%soil_prop,emiss_rate(:))
                case ('dustb')
                  CALL art_emission_dust(dz(:,nlev),ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil), &
                    &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),                      &
                    &             jb,istart,iend,'dustb',p_art_data(jg)%soil_prop,emiss_rate(:))
                case ('dustc')
                  CALL art_emission_dust(dz(:,nlev),ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil), &
                    &             ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),                      &
                    &             jb,istart,iend,'dustc',p_art_data(jg)%soil_prop,emiss_rate(:))
              end select
              
              ! Update mass mixing ratios
              DO ijsp = 1, fields%info%njsp
                CALL art_integrate_explicit(p_trac(:,nlev,jb,fields%info%jsp(ijsp)),  emiss_rate(:), dtime,          &
                  &                         istart,iend, opt_rho = rho(:,nlev,jb))
              ENDDO
              ! Update mass-specific number
              CALL art_integrate_explicit(p_trac(:,nlev,jb,fields%info%i_number_conc), emiss_rate(:), dtime,         &
                &                         istart,iend, opt_rho = rho(:,nlev,jb),                                     &
                &                         opt_fac=(fields%info%mode_fac * fields%info%factnum))
                
            class is (t_fields_1mom)
              ! ...
            class default
              call finish('mo_art_emission_interface:art_emission_interface', &
                   &      'ART: Unknown mode field type')
          end select
          this_mode => this_mode%next_mode
        ENDDO !while(associated)
      ENDDO !jb
    
      DEALLOCATE(emiss_rate,dz)
! DEV END -> Leave radioact,volc and chemtracer alone so far
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        select type (fields=>this_mode%fields)
          type is (t_fields_2mom)
            ! Nothing to do here
          class is (t_fields_radio)
            select case(TRIM(fields%info%name))
              case ('CS137')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iCS137),373)
              case ('I131')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iI131),340)
              case ('TE132')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iTE132),325)
              case ('ZR95')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iZR95),184)
              case ('XE133')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iXE133),355)
              case ('I131g')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iI131g),870)
              case ('I131o')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iI131o),880)
              case ('BA140')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iBA140),384)
              case ('RU103')
                call art_emiss_radioact(p_patch,dtime,rho,p_trac(:,:,:,iRU103),220)
              ! And Default...
              case default
                call finish('mo_art_emission_interface:art_emission_interface', &
                     &      'No according emission routine to mode')
            end select
          class is (t_fields_volc)
            ! nothing to do here, see below
          class default
            call finish('mo_art_emission_interface:art_emission_interface', &
                 &      'ART: Unknown mode field type')
        end select
        this_mode => this_mode%next_mode
      ENDDO
  
      ! ----------------------------------
      ! --- volcano emissions
      ! ----------------------------------
    
      IF (art_config(jg)%iart_volcano == 1) THEN
        CALL art_organize_emission_volc(p_patch,dtime,rho,p_trac) 
      ENDIF
    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- emissions of chemical tracer
    ! ----------------------------------
  
    IF (art_config(jg)%lart_chem) THEN
    
      IF (art_config(jg)%iart_chem_mechanism == 0) THEN
        CALL art_emiss_chemtracer(ext_data,p_patch,dtime,p_nh_state%diag,p_trac)
      ENDIF
      
    ENDIF
    
  ENDIF !lart
       
#endif

END SUBROUTINE art_emission_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_interface

