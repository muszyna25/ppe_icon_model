!>
!! Provides interface to ART-routines dealing with sedimentation of volcanic ash particles
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
!! 2002-2010 by DWD and MPI-M
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
    USE mo_art_data_volc,       ONLY: p_mflx_contra_vsed
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


    CHARACTER(len=32), POINTER :: var_name            !< returns a character containing the name
                                                     !< of current ash component without the time level
                                                     !< suffix at the end. e.g. qash1(.TL1)

    REAL(wp), POINTER  :: diameter_ash, &
    &                     rho_ash                !<  resturns diameter and density of volcanic ash particles

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_art_sedimentation_interface:art_sedi_interface")
    
    INTEGER  :: jg,jc,jk,ikp1,jb           !< loop index for: patch,index in block,full and half levels,block
    INTEGER  :: nlev,nlevp1,nblks,istat, &
    &           i_nchdom, i_rlstart, i_rlend,  i_startblk, i_endblk,i_startidx, i_endidx
    INTEGER  :: p_iubc, &                !< Upper boundary condition. Default value=0, no upper bc cond.
    &           p_itype_vlimit           !< Type of limiter for vertical transport. Default val. =1, semi-monotone slope limiter.
    LOGICAL  :: lcompute_gt, lcleanup_gt !Compute and clean up geometrical terms in connection to flux calculation.  
    LOGICAL,SAVE  :: ltest=.TRUE.  
    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
    
     jg  = p_patch%id
     p_iubc =0        ! No upper boundary condition
!     p_itype_vlimit=4 ! Positive definite limiter1 
!     p_itype_vlimit=1 ! Semi monotone limiter
     p_itype_vlimit=2 ! Monotone limiter


     nlev      = p_patch%nlev        !< Number of vertical full levels
     nlevp1    = p_patch%nlevp1      !< Number of vertical half levels
     nblks     = p_patch%nblks_c

     !Get all cell enitities, except halos
     i_nchdom  = MAX(1,p_patch%n_childdom)
     i_rlstart = 1  !is always one
     i_rlend   = min_rlcell
     i_startblk = p_patch%cells%start_blk(i_rlstart,1)
     i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

     ALLOCATE(p_mflx_contra_vsed(nproma,nlevp1,nblks),p_upflux_sed(nproma,nlevp1,nblks),stat=istat)


    lcompute_gt=.TRUE. ! compute geometrical terms
    lcleanup_gt=.TRUE. ! clean up geometrical terms. obs. this i currently done for all components. improvement:
                       ! compute values for first component, cleanup after last component.


!      WRITE(0,*) 'K.L. in sed_ifc'
      current_element=>p_prog_list%p%first_list_element 

      !start DO-loop over elements in list:
      DO WHILE (ASSOCIATED(current_element))

      !get meta data of current element:
      info=>current_element%field%info

      ! assure that current element is tracer
      IF (info%tracer%lis_tracer) THEN
        IF (info%tracer%lsed_tracer) THEN

         !
         ! retrieve  running index:
         !
          jsp           =>  info%ncontained 
          var_name      =>  info%name
          diameter_ash  =>  info%tracer%rdiameter_tracer
          rho_ash       =>  info%tracer%rrho_tracer
 
          WRITE (message_text,*) 'Sedimentation of ',var_name,' with idx= ',jsp,info%tracer%lsed_tracer
          CALL message(TRIM(routine),message_text)

         !
         ! calculate sedimentation velocities of volcanic ash particles:
         !
            CALL art_sedi_volc(p_patch,p_dtime,p_metrics,     &
              &                p_rho,&
              &                p_diag,    &
              &                diameter_ash,rho_ash,&
              &                p_mflx_contra_vsed) 
!*******************************************
! MaBa: SET CONC IN LAYERS FOR SEDI TEST
!      DO jb = i_startblk, i_endblk
!
!        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
!                       i_startidx, i_endidx, i_rlstart, i_rlend)
!
!!          DO  jk =1, nlev
!          DO jk =1, nlev-1
!            ! index of top half level
!            DO jc = i_startidx, i_endidx
!           IF(ltest) THEN
!              p_tracer_new(jc,jk+1,jb,jsp) =                                         &
!                 &    1.0E-30    !kg m-3
!            IF(jk.eq.45) THEN
!              p_tracer_new(jc,jk,jb,jsp) =                                         &
!                 &    0.1    !kg m-3
!            ENDIF
!            IF(jsp.eq.11) THEN
!              p_tracer_new(jc,jk,jb,jsp) =                                         &
!                 &    1.1    !kg m-3
!            ENDIF
!            ELSE 
!            p_tracer_new(jc,jk,jb,jsp) = 1.0E-30 
!            ENDIF !layer
!           ENDIF !ltest

!              IF(jb.eq.7.and.jc.eq.724.and.jk.ge.45) THEN
!                WRITE(0,*) 'K.L. tracer INI 1 mass:',p_tracer_new(jc,jk,jb,jsp), p_cellhgt_mc_now(jc,jk,jb)
!              ENDIF !jb..jk..
!            p_tracer_new(jc,jk,jb,jsp) =                                         &
!                & p_tracer_new(jc,jk,jb,jsp)/p_rho(jc,jk,jb) !kg kg-1
!             IF(jb.eq.7.and.jc.eq.724.and.jk.eq.45) THEN
!               WRITE(0,*) 'K.L. tracer INI 2 q:',p_tracer_new(jc,jk,jb,jsp)
!             ENDIF  ! jb..jk..
!          ENDDO  !jc
!          ENDDO  !jk
!          ENDDO  !jb

!*******************************

         !
         ! calculate vertical flux term due to sedimentation
         !
          
          CALL upwind_vflux_ppm_cfl( p_patch, p_tracer_new(:,:,:,jsp), p_iubc,    &! in
            &                  p_mflx_contra_vsed, p_dtime, lcompute_gt, &! in
            &                  lcleanup_gt, p_itype_vlimit,       &! in
            &                  p_cellhgt_mc_now, p_rhodz_new,             &! in
            &                  p_upflux_sed(:,:,:))!,                           &! out
!            &                  opt_topflx_tra=opt_topflx_tra(:,:,jt),        &! in
!            &                  opt_slev=p_iadv_slev(jt),                     &! in
!            &                  opt_rlstart=opt_rlstart,                      &! in
!            &                  opt_rlend=opt_rlend                           )! in




         !
         ! update mixing ratio after sedimentation
         ! 
WRITE(0,*) 'K.L. sedi update bef',jsp, MAXVAL(p_tracer_new(:,:,:,jsp)),MAXLOC(p_tracer_new(:,:,:,jsp))

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

!          DO jk =1, nlev
          DO jk =1, nlev-1

            ! index of top half level
            ikp1 = jk + 1

            DO jc = i_startidx, i_endidx

!***           IF(.NOT. ltest) THEN

            p_tracer_new(jc,jk,jb,jsp) =                                         &
                &    p_tracer_new(jc,jk,jb,jsp)   &
                &  - p_dtime * (  p_upflux_sed(jc,jk,jb)                &
                &               - p_upflux_sed(jc,ikp1  ,jb) )             &
                &  / p_rhodz_new(jc,jk,jb)

 
!***           ENDIF !ltest
!***             p_tracer_new(jc,jk,jb,jsp) =                                         &
!***                & p_tracer_new(jc,jk,jb,jsp)*p_rho(jc,jk,jb) !kg m-3

          
            END DO!jc
          END DO !jk
       END DO !jb
WRITE(0,*) 'K.L. sedi update after',jsp, MAXVAL(p_tracer_new(:,:,:,jsp)),MAXLOC(p_tracer_new(:,:,:,jsp))

      ENDIF !lsed_tracer

      ENDIF !lis_tracer

      ! select the next element in the list
      current_element => current_element%next_list_element

     ENDDO !loop elements
           ltest=.FALSE. 

       DEALLOCATE(p_mflx_contra_vsed,p_upflux_sed)
#endif

  END SUBROUTINE art_sedi_interface


END MODULE mo_art_sedi_interface

