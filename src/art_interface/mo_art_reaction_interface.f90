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
!!
!! @par Copyright
!! 2002-2010 by DWD, MPI-M, and KIT
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
MODULE mo_art_reaction_interface

    USE mo_kind,                   ONLY: wp
    USE mo_model_domain,           ONLY: t_patch
    USE mo_art_config,             ONLY: art_config
    USE mo_exception,              ONLY: message, message_text, finish
    USE mo_linked_list,            ONLY: t_var_list,t_list_element
    USE mo_var_metadata_types,     ONLY: t_var_metadata
    USE mo_nonhydro_types,         ONLY: t_nh_diag
#ifdef __ICON_ART
    USE mo_art_radioactive,        ONLY: art_decay_radioact
!    USE mo_art_chemtracer,      ONLY: art_loss_chemtracer
    USE mo_art_modes_linked_list,  ONLY: p_mode_state,t_mode
    USE mo_art_modes,              ONLY: t_fields_radio
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_reaction_interface


CONTAINS

  !>
  !! Interface for ART-routines treating reactions of any kind (chemistry, radioactive decay)
  !!
  !! This interface calls the ART-routines, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Max Bangert, KIT (2013-02-25)
  SUBROUTINE art_reaction_interface( p_patch,p_dtime,p_prog_list,p_diag,p_tracer_now)


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN) ::p_dtime           !< time step

    TYPE(t_nh_diag), INTENT(IN) ::p_diag

    TYPE(t_var_list), INTENT(IN) :: &        !< current prognostic state list
      &  p_prog_list

    REAL(wp), INTENT(INOUT) ::  &
      &  p_tracer_now(:,:,:,:)               !< tracer mixing ratios (specific concentrations)

#ifdef __ICON_ART
    TYPE(t_mode), POINTER   :: this_mode
#endif
    INTEGER  :: jg                           !< domain index

    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART

jg  = p_patch%id

IF(art_config(jg)%lart) THEN
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
  
  DO WHILE(ASSOCIATED(this_mode))
    ! Select type of mode
    select type (fields=>this_mode%fields)
        type is (t_fields_radio)
          CALL  art_decay_radioact(p_patch,p_dtime,p_tracer_now(:,:,:,fields%itracer),fields%halflife) 
    end select                  
    this_mode => this_mode%next_mode
  END DO
  
ENDIF
#endif

  END SUBROUTINE art_reaction_interface


END MODULE mo_art_reaction_interface

