!>
!! Provides interface to ART-routines dealing with emissions
!!
!! This module provides an interface to the ART-routine emission_volc.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Reinert, DWD
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2012-01-27)
!! Modification by Kristina Lundgren, KIT (2012-01-30)
!! - Modification for dealing with the ART-routine emission_volc.
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
MODULE mo_art_emission_interface

    USE mo_kind,                  ONLY: wp
    USE mo_model_domain,          ONLY: t_patch
    USE mo_art_config,            ONLY: art_config
    USE mo_exception,             ONLY: message, message_text, finish
    USE mo_datetime,              ONLY: t_datetime
    USE mo_linked_list,           ONLY: t_var_list,t_list_element
    USE mo_var_metadata,          ONLY: t_var_metadata
#ifdef __ICON_ART
    USE mo_art_emission_volc,     ONLY: art_organize_emission_volc
    USE mo_art_radioactive,       ONLY: art_emiss_radioact
    USE mo_art_emission_seas,     ONLY: art_emission_seas
    USE mo_art_emission_dust,     ONLY: art_emission_dust
    USE mo_art_aerosol,           ONLY: p_art_mode,imode_seasa,imode_seasb,imode_seasc, &
        &                               imode_dusta,imode_dustb,imode_dustc
    USE mo_art_aerosol_utilities, ONLY: art_modal_parameters,art_air_properties
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_emission_interface


CONTAINS


  !>
  !! Interface for ART-routine emission_volc
  !!
  !! This interface calls the ART-routine emission_volc, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-01-27)
  !! Modification by Kristina Lundgren, KIT (2012-01-30)
  !! - Call modified for emission of volcanic ash
  SUBROUTINE art_emission_interface( p_patch,p_dtime,datetime,p_prog_list,p_rho,p_tracer_now)


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed
    REAL(wp), INTENT(IN) ::p_dtime
     TYPE(t_datetime), INTENT(IN)::datetime

    TYPE(t_var_list), INTENT(IN) :: &        !< current prognostic state list
      &  p_prog_list

    REAL(wp), INTENT(INOUT) ::  &  !< density of air 
      &  p_rho(:,:,:)              !< [kg/m3] 
                                   !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &  !< tracer mixing ratios (specific concentrations)
      &  p_tracer_now(:,:,:,:)     !< at current time level n (before transport)
                                   !< [kg/kg]
                                   !< dim: (nproma,nlev,nblks_c,ntracer)
    
    INTEGER  :: jg                !< patch id

    TYPE(t_list_element), POINTER :: current_element !< returns the reference to
                                                     !< current element in list
    TYPE(t_var_metadata), POINTER :: info            !< returns reference to tracer
                                                     !< metadata of current element
    INTEGER, POINTER :: jsp                          !< returns index of element
    CHARACTER(len=32), POINTER :: var_name           !< returns a character containing the name
                                                     !< of current ash component without the time level
                                                     !< suffix at the end. e.g. qash1(.TL1)
    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
    
   jg  = p_patch%id
     
   IF (art_config(jg)%lart .AND. art_config(jg)%lart_emiss) THEN

     ! First: Modal Parameters 
       CALL art_air_properties(p_patch)
       
   ! ----------------------------------
   ! --- sea salt emissions
   ! ----------------------------------
       
       IF (art_config(jg)%lart_seasalt) THEN
         CALL art_modal_parameters(p_patch,p_art_mode(imode_seasa),p_tracer_now,'EMISSION')
         CALL art_modal_parameters(p_patch,p_art_mode(imode_seasb),p_tracer_now,'EMISSION')
         CALL art_modal_parameters(p_patch,p_art_mode(imode_seasc),p_tracer_now,'EMISSION')
         CALL art_emission_seas(p_patch,p_dtime,p_rho,p_tracer_now) 
       ENDIF
   
   ! ----------------------------------
   ! --- mineral dust emissions
   ! ----------------------------------
       
       IF (art_config(jg)%lart_dust) THEN
         CALL art_modal_parameters(p_patch,p_art_mode(imode_dusta),p_tracer_now,'EMISSION')
         CALL art_modal_parameters(p_patch,p_art_mode(imode_dustb),p_tracer_now,'EMISSION')
         CALL art_modal_parameters(p_patch,p_art_mode(imode_dustc),p_tracer_now,'EMISSION')
         CALL art_emission_dust(p_patch,p_dtime,p_rho,p_tracer_now) 
       ENDIF
   
   ! ----------------------------------
   ! --- volcano emissions
   ! ----------------------------------

       IF (art_config(jg)%lart_volcano) THEN
         CALL art_organize_emission_volc(p_patch,p_dtime,p_rho,p_tracer_now) 
       ENDIF

   ! ----------------------------------
   ! --- radioactive nuclide emissions
   ! ----------------------------------

       IF (art_config(jg)%lart_radioact) THEN

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
           ! ----------------------------------
           ! --- retrieve  running index:
           ! ----------------------------------

             jsp=>info%ncontained
             var_name=>info%name

             IF(info%tracer%tracer_class=='radioact') THEN
               CALL art_emiss_radioact(p_patch,p_dtime,p_tracer_now(:,:,:,jsp),p_rho,info%tracer%imis_tracer)
             ENDIF

           ENDIF
           ! ----------------------------------
           ! --- select the next element in the list
           ! ----------------------------------

           current_element => current_element%next_list_element

         ENDDO ! loop elements

       ENDIF
       
   ENDIF

#endif

  END SUBROUTINE art_emission_interface


END MODULE mo_art_emission_interface

