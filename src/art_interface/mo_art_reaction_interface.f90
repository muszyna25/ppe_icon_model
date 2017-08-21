!>
!! Provides interface to ART-routines dealing with chemical reactions and radioactive decay
!!
!! This module provides an interface to the ART-routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Max Bangert, KIT
!!
!! @par Revision History
!! Initial revision by Max Bangert, KIT (2013-02-15)
!! Modifications by Daniel Rieger, KIT (2014-05-22)
!! - Adaption to changes in ART data structure
!! Modifications by Roland Ruhnke, Jennifer Schröter, KIT (2015-08-06)
!! - Splitting between simple chemtracer and full gas phase chemistry
!!   called by iart_chem_mechanism == 0 or 1
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_reaction_interface

  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: min_rlcell_int, iecham, inwp
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_nonhydro_types,                ONLY: t_nh_diag
  USE mo_run_config,                    ONLY: lart, iforcing
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_nonhydro_types,                ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag
  USE mtime,                            ONLY: datetime
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_reacInt,            &
                                          &   timer_art_losschem, timer_art_photo
  USE mo_echam_phy_memory,              ONLY: t_echam_phy_tend
#ifdef __ICON_ART
  USE mo_art_decay_radioact,            ONLY: art_decay_radioact
  USE mo_art_chemtracer,                ONLY: art_loss_chemtracer
  USE mo_art_gasphase,                  ONLY: art_loss_gasphase
  USE mo_art_photolysis,                ONLY: art_photolysis
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio
  USE mo_art_config,                    ONLY: art_config
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_reaction_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_reaction_interface(ext_data, p_patch,current_date,p_dtime,p_prog_list,p_prog, &
  &                               p_metrics,p_diag,tracer, prm_diag, tend)

  !>
  !! Interface for ART-routines treating reactions of any kind (chemistry, radioactive decay)
  !!
  !! This interface calls the ART-routines, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Max Bangert, KIT (2013-02-25)
  !!
  ! atmosphere external data                                
  TYPE(t_external_data), INTENT(INOUT) :: &
    &  ext_data
  TYPE(t_patch), TARGET, INTENT(IN)   :: & 
    &  p_patch                           !< patch on which computation is performed
  TYPE(datetime), POINTER, INTENT(IN) :: &
    &  current_date                      !< Actual time and date
  REAL(wp), INTENT(IN)                :: &
    &  p_dtime                           !< time step
  TYPE(t_var_list), INTENT(IN)        :: &
    &  p_prog_list                       !< current prognostic state list
  TYPE(t_nh_prog), INTENT(IN)         :: &
    &  p_prog
  TYPE(t_nh_metrics), INTENT(IN)      :: &
    &  p_metrics                         !< NH metrics state
  TYPE(t_nh_diag), INTENT(IN)       :: &
    &  p_diag                            !< list of diagnostic fields
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:,:,:)                   !< tracer mixing ratios (specific concentrations)
  TYPE(t_nwp_phy_diag),OPTIONAL, INTENT(IN)  :: &
    &  prm_diag                          !< NH metrics state
  TYPE(t_echam_phy_tend) , OPTIONAL,  POINTER  :: tend
! Local variables
  INTEGER                           :: &
    &  jb,                             & !< loop index
    &  jg,                             & !< domain index
    &  i_startblk, i_endblk,           & !< Start and end of block loop
    &  istart, iend,                   & !< Start and end of nproma loop
    &  i_rlstart, i_rlend,             & !< Relaxation start and end
    &  i_nchdom,                       & !< Number of child domains
    &  nlev                              !< Number of levels (equals index of lowest full level)
#ifdef __ICON_ART
  TYPE(t_mode), POINTER   :: this_mode
#endif
  
  !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
  IF(lart) THEN

    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_reacInt)

    ! --- Get the loop indizes
    i_nchdom   = MAX(1,p_patch%n_childdom)
    jg         = p_patch%id
    nlev       = p_patch%nlev
    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

    IF (art_config(jg)%lart_aerosol) THEN
      ! ----------------------------------
      ! --- Radioactive particles
      ! ----------------------------------
      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
    
      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_radio)
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                &                istart, iend, i_rlstart, i_rlend)
              CALL art_decay_radioact(p_dtime,istart,iend,nlev,     &
                &                     fields%halflife,tracer(:,:,jb,fields%itr))
            ENDDO
        END SELECT
        this_mode => this_mode%next_mode
      END DO
    
      NULLIFY(this_mode)
      
    ENDIF !lart_aerosol
    
    ! ----------------------------------
    ! --- chemical tracer reactions
    ! ----------------------------------

    IF (art_config(jg)%lart_chem) THEN
      SELECT CASE(art_config(jg)%iart_chem_mechanism)
        CASE(0)
          IF (timers_level > 3) CALL timer_start(timer_art_losschem)
          CALL art_loss_chemtracer(ext_data, p_patch, &
                 & current_date,                      &
                 & p_dtime,                           &
                 & p_prog,                            &
                 & p_prog_list,                       &
                 & p_diag,                            &
                 & p_metrics,                         &
                 & tracer)
          IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
        CASE(1)
          IF (timers_level > 3) CALL timer_start(timer_art_photo)
            IF (iforcing == inwp) THEN
          CALL art_photolysis(ext_data,               &
                 & p_patch,                           &
                 & current_date,                      &
                 & p_dtime,                           &
                 & p_prog_list,                       &
                 & p_prog,                            &
                 & p_diag,                            &
                 & p_prog%rho,                        &
                 & p_metrics,                         &
                 & tracer,                            &
                 & prm_diag = prm_diag                )
         ELSEIF (iforcing == iecham) THEN
          CALL art_photolysis(ext_data,               &
                 & p_patch,                           &
                 & current_date,                      &
                 & p_dtime,                           &
                 & p_prog_list,                       &
                 & p_prog,                            &
                 & p_diag,                            &
                 & p_prog%rho,                        &
                 & p_metrics,                         &
                 & tracer                             )
          ENDIF

          IF (timers_level > 3) CALL timer_stop(timer_art_photo)

          IF (timers_level > 3) CALL timer_start(timer_art_losschem)

          CALL art_loss_chemtracer(ext_data, p_patch, &
                 & current_date,                      &
                 & p_dtime,                           &
                 & p_prog,                            &
                 & p_prog_list,                       &
                 & p_diag,                            &
                 & p_metrics,                         &
                 & tracer)
          IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
        CASE(2)
          IF (iforcing == inwp) THEN
            IF (timers_level > 3) CALL timer_start(timer_art_photo)

            CALL art_photolysis(ext_data,               &
                   & p_patch,                           &
                   & current_date,                      &
                   & p_dtime,                           &
                   & p_prog_list,                       &
                   & p_prog,                            &
                   & p_diag,                            &
                   & p_prog%rho,                        &
                   & p_metrics,                         &
                   & tracer,                            &
                   & prm_diag = prm_diag)

            IF (timers_level > 3) CALL timer_stop(timer_art_photo)
            IF (timers_level > 3) CALL timer_start(timer_art_losschem)

            CALL art_loss_gasphase(current_date,        &
                   & ext_data,                          &
                   & p_patch,                           &
                   & p_dtime,                           &
                   & p_prog_list,                       &
                   & p_diag,                            &
                   & p_metrics,                         &
                   & tracer)

            IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
          ELSEIF (iforcing == iecham) THEN
            IF (timers_level > 3) CALL timer_start(timer_art_photo)

            CALL art_photolysis(ext_data,               &
                   & p_patch,                           &
                   & current_date,                      &
                   & p_dtime,                           &
                   & p_prog_list,                       &
                   & p_prog,                            &
                   & p_diag,                            &
                   & p_prog%rho,                        &
                   & p_metrics,                         &
                   & tracer)

            IF (timers_level > 3) CALL timer_stop(timer_art_photo)
            IF (timers_level > 3) CALL timer_start(timer_art_losschem)

            CALL art_loss_gasphase(current_date,        &
                   & ext_data,                          &
                   & p_patch,                           &
                   & p_dtime,                           &
                   & p_prog_list,                       &
                   & p_diag,                            &
                   & p_metrics,                         &
                   & tracer)

            IF (timers_level > 3) CALL timer_stop(timer_art_losschem)
          ENDIF

        CASE DEFAULT
          CALL finish('mo_art_reaction_interface:art_reaction_interface', &
               &      'ART: Unknown iart_chem_mechanism')
      END SELECT
    ENDIF !lart_chem

    IF (timers_level > 3) CALL timer_stop(timer_art_reacInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF ! lart
#endif

END SUBROUTINE art_reaction_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_reaction_interface
