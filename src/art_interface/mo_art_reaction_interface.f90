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

    USE mo_kind,                ONLY: wp
    USE mo_model_domain,        ONLY: t_patch
    USE mo_art_config,          ONLY: art_config
    USE mo_exception,           ONLY: message, message_text, finish
    USE mo_linked_list,         ONLY: t_var_list,t_list_element
    USE mo_var_metadata_types,  ONLY: t_var_metadata
#ifdef __ICON_ART
    USE mo_art_radioactive,       ONLY: art_decay_radioact
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
  SUBROUTINE art_reaction_interface( p_patch,p_dtime,p_prog_list,p_tracer_now)


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN) ::p_dtime           !< time step

    TYPE(t_var_list), INTENT(IN) :: &        !< current prognostic state list
      &  p_prog_list

    REAL(wp), INTENT(INOUT) ::  &  !< tracer mixing ratios (specific concentrations)
      &  p_tracer_now(:,:,:,:)     !< at current time level n (before transport)      ! MaBa Timelevel has to fit!!!
                                   !< [kg/kg]
                                   !< dim: (nproma,nlev,nblks_c,ntracer)



    TYPE(t_list_element), POINTER :: current_element !< returns the reference to
                                                     !< current element in list
    TYPE(t_var_metadata), POINTER :: info            !< returns reference to tracer
                                                     !< metadata of current element
    
    INTEGER  :: jg                !< loop index

    INTEGER, POINTER :: jsp                          !< returns index of element
    CHARACTER(len=32), POINTER :: var_name           !< returns a character containing the name
                                                     !< of current ash component without the time level
                                                     !< suffix at the end. e.g. qash1(.TL1)

    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART

jg  = p_patch%id

IF(art_config(jg)%lart) THEN

      current_element=>p_prog_list%p%first_list_element

      ! ----------------------------------
      ! --- start DO-loop over elements in list:
      ! ----------------------------------

      DO WHILE (ASSOCIATED(current_element))

      ! ----------------------------------
      ! --- get meta data of current element:
      ! ----------------------------------

      info=>current_element%field%info

      ! ----------------------------------
      ! ---  assure that current element is tracer
      ! ----------------------------------

      IF (info%tracer%lis_tracer) THEN

        IF (art_config(jg)%lart_radioact .AND. art_config(jg)%lart_decay_radioact) THEN

          SELECT CASE(info%tracer%tracer_class)

          CASE('radioact')

          ! ----------------------------------
          ! --- retrieve  running index:
          ! ----------------------------------

          jsp=>info%ncontained
          var_name=>info%name

          CALL  art_decay_radioact(p_patch,p_dtime,p_tracer_now(:,:,:,jsp),info%tracer%halflife_tracer) 

          END SELECT

        ENDIF
      ENDIF !lis_tracer

      ! ----------------------------------
      ! --- select the next element in the list
      ! ----------------------------------

      current_element => current_element%next_list_element

     ENDDO !loop elements

ENDIF
#endif

  END SUBROUTINE art_reaction_interface


END MODULE mo_art_reaction_interface

