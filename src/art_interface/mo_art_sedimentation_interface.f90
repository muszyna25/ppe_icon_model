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
!!
!! @par Copyright
!! 2002-2010 by DWD, MPI-M, and KIT.
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
!!    an according license agreement with DWD, MPI-M, and KIT.
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
MODULE mo_art_sedi_interface

    USE mo_kind,                   ONLY: wp
    USE mo_parallel_config,        ONLY: nproma
    USE mo_model_domain,           ONLY: t_patch
    USE mo_impl_constants,         ONLY: min_rlcell
    USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_metrics,t_nh_diag
    USE mo_art_config,             ONLY: art_config
    USE mo_exception,              ONLY: finish
    USE mo_linked_list,            ONLY: t_var_list, t_list_element
    USE mo_var_metadata_types,     ONLY: t_var_metadata, t_tracer_meta
    USE mo_advection_vflux,        ONLY: upwind_vflux_ppm_cfl
    USE mo_run_config,             ONLY: ntracer
    USE mo_loopindices,            ONLY: get_indices_c
#ifdef __ICON_ART
! infrastructure routines
    USE mo_art_modes_linked_list,  ONLY: p_mode_state,t_mode
    USE mo_art_modes,              ONLY: t_fields_2mom,t_fields_radio, &
        &                                t_fields_volc
    USE mo_art_data,               ONLY: p_art_data, UNDEF_INT_ART
    USE mo_art_clipping,           ONLY: art_clip_tracers_zero
! sedimentation and deposition routines
    USE mo_art_sedi_volc,          ONLY: art_sedi_volc
!    USE mo_art_aerosol,            ONLY: p_mflx_contra_vsed, vdep_ash
    USE mo_art_sedi_depo,          ONLY: art_calc_v_sed_dep
!    USE mo_art_aerosol_utilities,  ONLY: art_modal_parameters,art_air_properties
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_sedi_interface


CONTAINS


  !>
  !! Interface for ART-routine sedi_volc
  !!
  !! This interface calls the ART-routine sedi_volc, if ICON has been
  !! built including the ART-package. Otherwise, this is simply a dummy
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Kristina Lundgren, KIT (2012-06-01)
  !!

  SUBROUTINE art_sedi_interface( p_patch, &
             &                   p_dtime, &
             &                   p_prog_list, p_prog, &
             &                   p_metrics, p_rho, p_diag, &
             &                   p_tracer_new, &
             &                   p_cellhgt_mc_now, p_rhodz_new, &
             &                   lprint_cfl, &
             &                   opt_topflx_tra )


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN)              :: &
      &  p_dtime

    TYPE(t_var_list), INTENT(IN)      :: &   !< current prognostic state list
     &  p_prog_list                          ! drieg: I think we do not need this

    TYPE(t_nh_prog), INTENT(IN)       :: &   !< current prognostic state
     &  p_prog

    TYPE(t_nh_metrics), INTENT(IN)    :: &
     &   p_metrics


    REAL(wp), INTENT(IN)              :: &   !<density of air at full levels
      &  p_rho(:,:,:)                        !< [kg/m3]
                                             !< dim: (nproma,nlev,nblks_c)

    TYPE(t_nh_diag), INTENT(IN)       :: &   !<diagnostic variables
      &  p_diag


    REAL(wp), INTENT(INOUT) ::  &            !< tracer mixing ratios (specific concentrations)
      &  p_tracer_new(:,:,:,:)               !< at current time level n+1 (after transport)
                                             !< [kg/kg]
                                             !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(IN) ::          &       !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)             !<
                                             !< NH: \Delta z       [m]
                                             !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &               !< NH: density weighted cell depth (\rho*(z_half(top)-z_half(bottom))
      &  p_rhodz_new(:,:,:)                  !< at full levels and time step n+1 [kg/m**2]
                                             !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &               !< determines if vertical CFL number shall be printed
      &  lprint_cfl                          !< in routine upwind_vflux_ppm_cfl

    REAL(wp), INTENT(IN), OPTIONAL:: &       !< vertical tracer flux at upper boundary
      &  opt_topflx_tra(:,:,:)               !< NH: [kg/m**2/s]
                                             !< dim: (nproma,nblks_c,ntracer)
                                             
    REAL(wp), ALLOCATABLE :: &      !< upwind flux at half levels due to sedimentation
      &  p_upflux_sed(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)
      
    REAL(wp),POINTER :: &
      &  mflx_contra_vsed(:,:,:),  &
      &  nflx_contra_vsed(:,:,:),  &
      &  flx_contra_vsed(:,:,:)
      
    REAL(wp) :: sedim_update            !< tracer tendency due to sedimentation
      
    INTEGER  :: jc,jk,ikp1,jb           !< loop index for: index in block,full and half levels,block
    INTEGER  :: nlev,nlevp1,nblks, &
    &           i_nchdom, i_rlstart, i_rlend,  i_startblk, i_endblk,i_startidx, i_endidx
      
    INTEGER  :: jg                      !< patch id
    INTEGER  :: jsp
    INTEGER  :: i
    INTEGER  :: iubc=0                  !< upper boundary condition 0 = none
    INTEGER  :: itype_vlimit=2          !< Monotone limiter
    LOGICAL  :: lcompute_gt
    LOGICAL  :: lcleanup_gt

#ifdef __ICON_ART
    TYPE(t_mode), POINTER   :: this_mode
#endif


#ifdef __ICON_ART

  lcompute_gt=.TRUE. ! compute geometrical terms
  lcleanup_gt=.TRUE. ! clean up geometrical terms. obs. this i currently done for all components. improvement:
                     ! compute values for first component, cleanup after last component.

  NULLIFY(mflx_contra_vsed)
  NULLIFY(nflx_contra_vsed)
  NULLIFY(flx_contra_vsed)
  
  jg     = p_patch%id
  nlevp1 = p_patch%nlevp1      !< Number of vertical half levels
  nlev       = p_patch%nlev    !< Number of vertical full levels
  nblks  = p_patch%nblks_c
  
  !Get all cell enitities, except halos
  i_nchdom  = MAX(1,p_patch%n_childdom)
  i_rlstart = 1  !is always one
  i_rlend   = min_rlcell
  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)
  
  IF(art_config(jg)%lart) THEN 
  
    ALLOCATE(p_upflux_sed(nproma,nlevp1,nblks))
    
!drieg: i think we do not need the air properties at this point, but i am not sure  yet
!    CALL art_air_properties(p_patch,p_art_data(jg))
    this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
   
    DO WHILE(ASSOCIATED(this_mode))
      ! Select type of mode
      select type (fields=>this_mode%fields)
      
        class is (t_fields_2mom)
          ! Before sedimentation/deposition velocity calculation, the modal parameters have to be calculated
          call fields%modal_param(p_art_data(jg),p_patch,p_tracer_new)
          
          call art_calc_v_sed_dep(p_patch,p_metrics,p_prog,p_diag,fields,p_rho,p_tracer_new)
          mflx_contra_vsed => fields%flx_contra_vsed3
          nflx_contra_vsed => fields%flx_contra_vsed0
          
        class is (t_fields_volc)
          call art_sedi_volc( p_patch,p_metrics,p_rho,            &
            &              p_diag,                                & 
            &              fields%diam,fields%rho,                &
            &              fields%flx_contra_vsed3) 
          mflx_contra_vsed => fields%flx_contra_vsed3
        class is (t_fields_radio)
         
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
        
          CALL upwind_vflux_ppm_cfl(p_patch, p_tracer_new(:,:,:,jsp),           & !< in
            &                       iubc, flx_contra_vsed, p_dtime,            & !< in! we need here nflx if i = 0
            &                       lcompute_gt, lcleanup_gt, itype_vlimit,     & !< in
            &                       p_cellhgt_mc_now, p_rhodz_new, lprint_cfl,  & !< in
            &                       p_upflux_sed(:,:,:), opt_elev=nlevp1 )        !< out

          ! ----------------------------------
          ! --- update mixing ratio after sedimentation
          ! ----------------------------------        
        
          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
              &                i_startidx, i_endidx, i_rlstart, i_rlend)
            DO jk = 1, nlev
              ! index of bottom half level
              ikp1 = jk + 1
              DO jc = i_startidx, i_endidx
  
                sedim_update = (p_dtime * (  p_upflux_sed(jc,jk,  jb)   &
                   &                      - p_upflux_sed(jc,ikp1,jb) )  &
                   &         / p_rhodz_new(jc,jk,jb))
            
                p_tracer_new(jc,jk,jb,jsp) =   p_tracer_new(jc,jk,jb,jsp) - sedim_update
            
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

    CALL art_clip_tracers_zero(p_tracer_new)
    
    DEALLOCATE(p_upflux_sed)
  ENDIF

#endif

  END SUBROUTINE art_sedi_interface


END MODULE mo_art_sedi_interface

