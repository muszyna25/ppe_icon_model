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

    USE mo_kind,                ONLY: wp
    USE mo_parallel_config,      ONLY: nproma
    USE mo_model_domain,        ONLY: t_patch
    USE mo_impl_constants,       ONLY: min_rlcell
    USE mo_nonhydro_types,       ONLY: t_nh_metrics,t_nh_diag
    USE mo_art_config,          ONLY: art_config
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_linked_list,         ONLY: t_var_list, t_list_element
    USE mo_var_metadata,        ONLY: t_var_metadata, t_tracer_meta
    USE mo_advection_vflux,     ONLY: upwind_vflux_ppm_cfl 
    USE mo_run_config,          ONLY: ntracer
    USE mo_loopindices,          ONLY: get_indices_c
#ifdef __ICON_ART
    USE mo_art_sedi_volc,       ONLY: art_sedi_volc
    USE mo_art_aerosol,         ONLY: p_mflx_contra_vsed
    USE mo_art_aerosol,         ONLY: p_art_mode,nmodes,imode_seasa,imode_seasb,imode_seasc
    USE mo_art_sedi_depo,       ONLY: art_calc_v_sed_dep
    USE mo_art_aerosol_utilities,  ONLY: art_modal_parameters,art_air_properties
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

  SUBROUTINE art_sedi_interface( p_patch,&
             &                   p_dtime,&
             &                   p_prog_list,&
             &                   p_metrics, &
             &                   p_rho,&
             &                   p_diag,     &
             &                   p_tracer_new,&
             &                   p_cellhgt_mc_now,&
             &                   p_rhodz_new,&
             &                   opt_topflx_tra)


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN)              :: &
      &  p_dtime

    TYPE(t_var_list), INTENT(IN)      :: &   !< current prognostic state list
     &  p_prog_list

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
      &  p_rhodz_new(:,:,:)               !< at full levels and time step n+1 [kg/m**2]
                                             !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN), OPTIONAL:: &       !< vertical tracer flux at upper boundary
      &  opt_topflx_tra(:,:,:)               !< NH: [kg/m**2/s]
                                             !< dim: (nproma,nblks_c,ntracer)

    ! local variables:

    REAL(wp), ALLOCATABLE :: &      !< upwind flux at half levels due to sedimentation 
      &  p_upflux_sed(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    TYPE(t_list_element), POINTER :: current_element !< returns the reference to 
                                                     !< current element in list
    TYPE(t_var_metadata), POINTER :: info            !< returns reference to tracer 
                                                     !< metadata of current element

    INTEGER, POINTER :: jsp                          !< returns index of element
    
    INTEGER          :: n                            !<loop variable

    CHARACTER(len=32), POINTER :: var_name            !< returns a character containing the name
                                                     !< of current ash component without the time level
                                                     !< suffix at the end. e.g. qash1(.TL1)

    REAL(wp), POINTER  :: diameter_ash, &
    &                     rho_ash                !<  resturns diameter and density of volcanic ash particles

    CHARACTER(*), PARAMETER :: art_routine = TRIM("mo_art_sedimentation_interface:art_sedi_interface")
    
    INTEGER  :: jg,jc,jk,ikp1,jb           !< loop index for: patch,index in block,full and half levels,block
    INTEGER  :: nlev,nlevp1,nblks,istat, &
    &           i_nchdom, i_rlstart, i_rlend,  i_startblk, i_endblk,i_startidx, i_endidx
    INTEGER  :: p_iubc, &                !< Upper boundary condition. Default value=0, no upper bc cond.
    &           p_itype_vlimit           !< Type of limiter for vertical transport. Default val. =1, semi-monotone slope limiter.
    LOGICAL  :: lcompute_gt, lcleanup_gt !Compute and clean up geometrical terms in connection to flux calculation.  

    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART

jg  = p_patch%id

IF(art_config(jg)%lart) THEN

    nlev      = p_patch%nlev        !< Number of vertical full levels
    nlevp1    = p_patch%nlevp1      !< Number of vertical half levels
    nblks     = p_patch%nblks_c

    !Get all cell enitities, except halos
    i_nchdom  = MAX(1,p_patch%n_childdom)
    i_rlstart = 1  !is always one
    i_rlend   = min_rlcell
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    ! First calculation of sedimentation and deposition velocity for the modal aerosol
    CALL art_air_properties(p_patch)
    DO n=1, nmodes
      CALL art_modal_parameters(p_patch,p_art_mode(n),p_tracer_new,'SEDIMENTATION')
      WRITE(*,*) 'Calculating sedimentation velocity for ', p_art_mode(n)%zname
      CALL art_calc_v_sed_dep(p_patch,p_metrics,p_diag,p_art_mode(n),p_rho,p_tracer_new)
    ENDDO
    
    p_iubc =0        ! No upper boundary condition
    p_itype_vlimit=2 ! Monotone limiter


    ALLOCATE(p_mflx_contra_vsed(nproma,nlevp1,nblks),p_upflux_sed(nproma,nlevp1,nblks),stat=istat)


    lcompute_gt=.TRUE. ! compute geometrical terms
    lcleanup_gt=.TRUE. ! clean up geometrical terms. obs. this i currently done for all components. improvement:
                       ! compute values for first component, cleanup after last component.


    current_element=>p_prog_list%p%first_list_element 

      !start DO-loop over elements in list:
    DO WHILE (ASSOCIATED(current_element))

      !get meta data of current element:
      info=>current_element%field%info

      ! assure that current element is tracer
      IF (info%tracer%lis_tracer) THEN
        IF (info%tracer%lsed_tracer) THEN

        ! ----------------------------------
        ! --- retrieve  running index
        ! ----------------------------------

          jsp           =>  info%ncontained 
          var_name      =>  info%name
          diameter_ash  =>  info%tracer%rdiameter_tracer
          rho_ash       =>  info%tracer%rrho_tracer
 
          WRITE (message_text,*) 'Sedimentation of ',var_name,' with idx= ',jsp,info%tracer%lsed_tracer
          CALL message(TRIM(art_routine),message_text)
          
        ! ----------------------------------
        ! --- calculate sedimentation velocities
        ! ----------------------------------

          SELECT CASE(info%tracer%tracer_class)

            CASE('volcash')

              CALL art_sedi_volc(p_patch,p_dtime,p_metrics,     &
                &                p_rho,                         &
                &                p_diag,                        &
                &                diameter_ash,rho_ash,          &
                &                p_mflx_contra_vsed) 

            CASE('radioact')
 
              WRITE (message_text,*) 'Sedimentation of ',var_name,' currently not possible'
              CALL message(TRIM(art_routine),message_text)
              p_mflx_contra_vsed = 0.0_wp 
              
            CASE('mode_seasa')
              p_mflx_contra_vsed = p_art_mode(imode_seasa)%mflx_contra_vsed3
            CASE('mode_seasa_number')
              p_mflx_contra_vsed = p_art_mode(imode_seasa)%mflx_contra_vsed0
            CASE('mode_seasb')
              p_mflx_contra_vsed = p_art_mode(imode_seasb)%mflx_contra_vsed3
            CASE('mode_seasb_number')
              p_mflx_contra_vsed = p_art_mode(imode_seasb)%mflx_contra_vsed0
            CASE('mode_seasc')
              p_mflx_contra_vsed = p_art_mode(imode_seasc)%mflx_contra_vsed3
            CASE('mode_seasc_number')
              p_mflx_contra_vsed = p_art_mode(imode_seasc)%mflx_contra_vsed0
              
          END SELECT

        ! ----------------------------------
        ! --- calculate vertical flux term due to sedimentation
        ! ----------------------------------
          
          CALL upwind_vflux_ppm_cfl( p_patch, p_tracer_new(:,:,:,jsp), p_iubc,    &! in
            &                  p_mflx_contra_vsed, p_dtime, lcompute_gt,          &! in
            &                  lcleanup_gt, p_itype_vlimit,                       &! in
            &                  p_cellhgt_mc_now, p_rhodz_new,                     &! in
            &                  p_upflux_sed(:,:,:))                                ! out
            
        ! ----------------------------------
        ! --- update mixing ratio after sedimentation
        ! ----------------------------------

            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)
              DO jk =1, nlev-1
              ! index of top half level
                ikp1 = jk + 1
                DO jc = i_startidx, i_endidx
                
                  p_tracer_new(jc,jk,jb,jsp) =                       &
                  &    p_tracer_new(jc,jk,jb,jsp)                    &
                  &  - p_dtime * (  p_upflux_sed(jc,jk,jb)           &
                  &               - p_upflux_sed(jc,ikp1  ,jb) )     &
                  &  / p_rhodz_new(jc,jk,jb)
                  
                  IF (p_tracer_new(jc,jk,jb,jsp) .LT. 0.0_wp) THEN
                    WRITE(*,*) 'After Sedi: Tracer ',var_name,' below 0: ',p_tracer_new(jc,jk,jb,jsp)
                    p_tracer_new(jc,jk,jb,jsp) = 0.0_wp
                    WRITE(*,*) 'p_upflux_sed_jk',p_upflux_sed(jc,jk,jb),'p_upflux_sed_ikp1',p_upflux_sed(jc,ikp1  ,jb)
                    WRITE(*,*) 'Diameter:',p_art_mode(imode_seasc)%diameter(jc,jk,jb)
                  ENDIF
                    
                END DO!jc
              END DO !jk
            END DO !jb

        ENDIF !lsed_tracer
      ENDIF !lis_tracer

      ! ----------------------------------
      ! --- select the next element in the list
      ! ----------------------------------

      current_element => current_element%next_list_element

    ENDDO !loop elements

    DEALLOCATE(p_mflx_contra_vsed,p_upflux_sed)

ENDIF !lart
#endif

  END SUBROUTINE art_sedi_interface


END MODULE mo_art_sedi_interface

