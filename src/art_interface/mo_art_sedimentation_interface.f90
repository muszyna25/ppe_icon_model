!>
!! Provides interface to ART-routines dealing with sedimentation
!!
!! This module provides an interface to the ART-routine sedi_volc.
!! The interface is written in such a way, that ICON will compile and run
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-01)
!! Modifications by Daniel Rieger, KIT (2014-05-22)
!! - Adaption to changes in ART data structure
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_sedi_interface

  USE mo_kind,                          ONLY: wp
  USE mo_parallel_config,               ONLY: nproma
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_nonhydro_types,                ONLY: t_nh_prog, t_nh_metrics,t_nh_diag
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mo_run_config,                    ONLY: lart,iqr,iqc,iqv
  USE mo_exception,                     ONLY: finish
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_advection_vflux,               ONLY: upwind_vflux_ppm_cfl
  USE mo_loopindices,                   ONLY: get_indices_c
#ifdef __ICON_ART
! infrastructure routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_volc
  USE mo_art_data,                      ONLY: p_art_data, UNDEF_INT_ART
  USE mo_art_clipping,                  ONLY: art_clip_lt
  USE mo_art_config,                    ONLY: art_config
! sedimentation and deposition routines
  USE mo_art_sedi_volc,                 ONLY: art_sedi_volc
  USE mo_art_sedi_2mom,                 ONLY: art_calc_v_sed, art_calc_sed_flx
  USE mo_art_depo_2mom,                 ONLY: art_calc_v_dep, art_store_v_dep
  USE mo_art_radioactive,               ONLY: art_drydepo_radioact
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_sedi_interface


CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_sedi_interface( p_patch, &
           &                   p_dtime, &
           &                   p_prog,  &
           &                   p_metrics, rho, p_diag, &
           &                   prm_diag, &
           &                   p_trac, &
           &                   dz, &
           &                   lprint_cfl )
!! Interface for ART routines for calculation of
!! sedimentation and deposition velocities
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-01)
!! Modifications by Daniel Rieger, KIT (2014-05-22)
!! - Adaption to changes in ART data structure
!! - Calculation of deposition velocities
  TYPE(t_patch), TARGET, INTENT(IN) :: &
    &  p_patch                            !< patch on which computation is performed
  REAL(wp), INTENT(IN)              :: &
    &  p_dtime
  TYPE(t_nh_prog), INTENT(IN)       :: &  
    &  p_prog                             !< current prognostic state
  TYPE(t_nh_metrics), INTENT(IN)    :: &
    &   p_metrics                         !< metrical fields
  REAL(wp), INTENT(IN)              :: & 
    &  rho(:,:,:)                       !< density of air at full levels
  TYPE(t_nh_diag), INTENT(IN)       :: &
    &  p_diag                             !< diagnostic variables
  TYPE(t_nwp_phy_diag), INTENT(in)  :: &
    &  prm_diag                !< list of diagnostic fields (physics)
  REAL(wp), INTENT(INOUT)           :: &
    &  p_trac(:,:,:,:)              !< tracer mixing ratio [kg/kg]
  REAL(wp), INTENT(IN)              :: &      
    &  dz(:,:,:)                        !< cell height defined at full levels
  LOGICAL, INTENT(IN)               :: &               
    &  lprint_cfl                         !< determines if vertical CFL number shall be printed in upwind_vflux_ppm_cfl
  REAL(wp), ALLOCATABLE             :: &  
    &  p_upflux_sed(:,:,:),            &  !< upwind flux at half levels due to sedimentation
    &  rhodz_new(:,:,:)                   !< density * height of full layer
    
  REAL(wp),POINTER                  :: &
    &  mflx_contra_vsed(:,:,:),        &
    &  nflx_contra_vsed(:,:,:),        &
    &  flx_contra_vsed(:,:,:)
      
  REAL(wp) ::       &
    &  sedim_update    !< tracer tendency due to sedimentation
      
  INTEGER  ::       &  
    &  jc,          &  !< loop index for: index in cell
    &  jk,          &  !< loop index for: index in level full
    &  ikp1,        &  !< loop index for: index in level half
    &  jb,          &  !< loop index for: index in block
    &  nlev,        &  !< number of full levels
    &  nlevp1,      &  !< number of half levels
    &  nblks,       &  !< number of blocks
    &  i_nchdom,    &
    &  i_rlstart,   &
    &  i_rlend,     &
    &  i_startblk,  &
    &  i_endblk,    &
    &  istart,iend, &
    &  jg,          & !< patch id
    &  jsp, i,      & !< counter
    &  iubc=0,      & !< upper boundary condition 0 = none
    &  itype_vlimit=2 !< Monotone limiter
  LOGICAL  :: lcompute_gt
  LOGICAL  :: lcleanup_gt

#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
  REAL(wp), ALLOCATABLE   :: &
    &  vsed0(:,:),           & !< Sedimentation velocities 0th moment [m s-1]
    &  vsed3(:,:),           & !< Sedimentation velocities 3th moment [m s-1]
    &  vdep0(:),             & !< Deposition velocities 0th moment [m s-1]
    &  vdep3(:)                !< Deposition velocities 3th moment [m s-1]
  
  
  lcompute_gt=.TRUE. ! compute geometrical terms
  lcleanup_gt=.TRUE. ! clean up geometrical terms. obs. this i currently done for all components. improvement:
                     ! compute values for first component, cleanup after last component.

  NULLIFY(mflx_contra_vsed)
  NULLIFY(nflx_contra_vsed)
  NULLIFY(flx_contra_vsed)
  
  jg     = p_patch%id
  nlevp1 = p_patch%nlevp1      !< Number of vertical half levels
  nlev   = p_patch%nlev    !< Number of vertical full levels
  nblks  = p_patch%nblks_c
  
  ALLOCATE(vsed0(nproma,nlev),vsed3(nproma,nlev))
  ALLOCATE(vdep0(nproma),vdep3(nproma))
  
  !Get all cell enitities, except halos
  i_nchdom  = MAX(1,p_patch%n_childdom)
  i_rlstart = grf_bdywidth_c+1
  i_rlend   = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
  
  IF(lart) THEN
    IF (art_config(jg)%lart_aerosol) THEN
    
      ALLOCATE(p_upflux_sed(nproma,nlevp1,nblks))
      ALLOCATE(rhodz_new(nproma,nlev,nblks))
      
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
          &                istart, iend, i_rlstart, i_rlend)
        DO jk = 1, nlev
          DO jc = istart, iend
            rhodz_new(jc,jk,jb) = rho(jc,jk,jb) * dz(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO
      
      
  !drieg: Necessary depending on parameterizations that will be used
  !    CALL art_air_properties(p_patch,p_art_data(jg))
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
     
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        select type (fields=>this_mode%fields)
        
          class is (t_fields_2mom)
            ! Before sedimentation/deposition velocity calculation, the modal parameters have to be calculated
            call fields%modal_param(p_art_data(jg),p_patch,p_trac)
            
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                &                istart, iend, i_rlstart, i_rlend)
              ! Calculate sedimentation velocities for 0th and 3rd moment
              CALL art_calc_v_sed(rho(:,:,jb),dz(:,:,jb),p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb), &
                &     fields%density(:,:,jb),fields%diameter(:,:,jb), fields%info%exp_aero,         &
                &     fields%knudsen_nr(:,:,jb), istart, iend, nlev, vsed0(:,:), vsed3(:,:))
              ! Calculate massflux due to sedimentation for 0th and 3rd moment
              CALL art_calc_sed_flx(vsed0(:,:),vsed3(:,:),p_metrics%wgtfac_c(:,:,jb),               &
                &     p_metrics%wgtfacq_c(:,:,jb),rho(:,:,jb),                                      &
                &     p_diag%rho_ic(:,:,jb), istart, iend, nlev,                                    &
                &     fields%flx_contra_vsed0(:,:,jb),fields%flx_contra_vsed3(:,:,jb))
              ! Calculate deposition velocities for 0th and 3rd moment
              CALL art_calc_v_dep(p_diag%temp(:,nlev,jb), p_diag%temp_ifc(:,nlev+1,jb),             &
                &     p_diag%u(:,nlev,jb), p_diag%v(:,nlev,jb), rho(:,nlev,jb),prm_diag%tcm(:,jb),  &
                &     prm_diag%tch(:,jb),p_trac(:,nlev,jb,iqv),p_trac(:,nlev,jb,iqc),               &
                &     p_trac(:,nlev,jb,iqr), p_prog%theta_v(:,nlev,jb), prm_diag%gz0(:,jb),         &
                &     dz(:,nlev,jb), p_art_data(jg)%air_prop%art_dyn_visc(:,nlev,jb),               &
                &     fields%diameter(:,nlev,jb),fields%info%exp_aero,fields%knudsen_nr(:,nlev,jb), &
                &     vsed0(:,nlev), vsed3(:,nlev), istart, iend, nlev, vdep0(:), vdep3(:))
              ! Store deposition velocities for the use in turbulence scheme
              CALL art_store_v_dep(vdep0(:), vdep3(:), fields%info%njsp, fields%info%jsp(:),        &
                &     fields%info%i_number_conc, art_config(jg)%nturb_tracer,                       &
                &     p_prog%turb_tracer(jb,:),istart,iend,p_art_data(jg)%turb_fields%vdep(:,jb,:))
            ENDDO
            mflx_contra_vsed => fields%flx_contra_vsed3
            nflx_contra_vsed => fields%flx_contra_vsed0
            
          class is (t_fields_volc)
            call art_sedi_volc( p_patch,p_metrics,p_prog,           &
              &              rho, p_diag,                         & 
              &              fields%diam,fields%rho,                &
              &              fields%flx_contra_vsed3,               &
              &              fields%itracer) 
            mflx_contra_vsed => fields%flx_contra_vsed3
          class is (t_fields_radio)
            ! Sedimentation velocity is zero for radioact. tracers
            fields%flx_contra_vsed3(:,:,:) = 0.0_wp
            mflx_contra_vsed => fields%flx_contra_vsed3
            ! However, a deposition velocity is required
            call art_drydepo_radioact( p_patch,p_prog,           &
              &              fields%itracer) 
          class default
            call finish('mo_art_sedimentation_interface:art_sedimentation_interface', &
              &         'ART: Unknown mode field type')
        end select
        
        ! here the number needs to be accounted for, too
        do i=0, this_mode%fields%info%njsp !< loop through the tracer mass mixing ratios contained in the mode
          if (i .NE. 0) then !< get index of mass mixing ratio of the species contained in this_mode
            jsp = this_mode%fields%info%jsp(i)
            flx_contra_vsed => mflx_contra_vsed
          else               !< get index of number mixing ratio of this_mode
            jsp = this_mode%fields%info%i_number_conc
            flx_contra_vsed => nflx_contra_vsed
          endif
          
          if (jsp .NE. UNDEF_INT_ART) then !< for monodisperse tracer, i_number_conc = UNDEF_INT_ART
          
            ! ----------------------------------
            ! --- calculate vertical flux term due to sedimentation
            ! ----------------------------------
            ! upwind_vflux_ppm_cfl is internally OpenMP parallelized
            CALL upwind_vflux_ppm_cfl(p_patch, p_trac(:,:,:,jsp),           & !< in
              &                       iubc, flx_contra_vsed, p_dtime,             & !< in
              &                       lcompute_gt, lcleanup_gt, itype_vlimit,     & !< in
              &                       dz, rhodz_new, lprint_cfl,                & !< in
              &                       p_upflux_sed(:,:,:), opt_elev=nlevp1 )        !< out
  
            ! ----------------------------------
            ! --- update mixing ratio after sedimentation
            ! ----------------------------------        
          
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                &                istart, iend, i_rlstart, i_rlend)
              DO jk = 1, nlev
                ! index of bottom half level
                ikp1 = jk + 1
                DO jc = istart, iend
    
                  sedim_update = (p_dtime * (  p_upflux_sed(jc,jk,  jb)   &
                     &                      - p_upflux_sed(jc,ikp1,jb) )  &
                     &         / rhodz_new(jc,jk,jb))
              
                  p_trac(jc,jk,jb,jsp) =   p_trac(jc,jk,jb,jsp) - sedim_update
              
                END DO!jc
              END DO !jk
            END DO !jb
          endif !jsp .ne. -999
        end do !i
        
        this_mode => this_mode%next_mode
      END DO
    
    ! ----------------------------------
    ! --- Clip the tracers
    ! ----------------------------------
  
      CALL art_clip_lt(p_trac,0.0_wp)
      
      DEALLOCATE(p_upflux_sed)
    ENDIF !lart_aerosol
  ENDIF !lart
  
  DEALLOCATE(vsed0,vsed3,vdep0,vdep3)
  
#endif

END SUBROUTINE art_sedi_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_sedi_interface

