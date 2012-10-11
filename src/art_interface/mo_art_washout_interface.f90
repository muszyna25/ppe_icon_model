!>
!! Provides interface to ART-routines dealing with washout of volcanic ash particles
!!
!! This module provides an interface to the ART-routine emission_volc.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-15)
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
MODULE mo_art_washout_interface

    USE mo_kind,                ONLY: wp
    USE mo_model_domain,        ONLY: t_patch
    USE mo_art_config,          ONLY: art_config
    USE mo_exception,             ONLY: message, message_text, finish
    USE mo_linked_list,         ONLY: t_var_list,t_list_element
    USE mo_var_metadata,        ONLY: t_var_metadata
#ifdef __ICON_ART
    USE mo_art_washout_volc,       ONLY:art_washout_volc
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_washout_interface


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
  SUBROUTINE art_washout_interface(            & !>in
      &          dt_phy_jg,                    & !>in
      &          p_patch,                      & !>in
      &          p_prog_list,                  & !>in
      &          p_rain_gsp_rate,              & !>in
      &          p_snow_gsp_rate,              & !>in
      &          p_rain_con_rate,              & !>in
      &          p_snow_con_rate,              & !>in
      &          p_rho,                        & !>in
      &          p_tracer_new)                  !>inout

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN) ::dt_phy_jg         !< time interval, fast physics

    TYPE(t_var_list), INTENT(IN) :: & !current prognostic state list
      &  p_prog_list

    REAL(wp), INTENT(IN)    ::  &       !< grid-scale surface rain rate  [kg/m2/s]
      &  p_rain_gsp_rate(:,:)           !< dim: (nproma,nblks_c)

    REAL(wp), INTENT(IN)    ::  &       !< grid-scale surface snow rate  [kg/m2/s]
      &  p_snow_gsp_rate(:,:)           !< dim: (nproma,nblks_c)

    REAL(wp), INTENT(IN)    ::  &       !< convective surface rain rate  [kg/m2/s]
      &  p_rain_con_rate(:,:)           !< dim: (nproma,nblks_c)

    REAL(wp), INTENT(IN)    ::  &       !< convective surface snow rate  [kg/m2/s]
      &  p_snow_con_rate(:,:)           !< dim: (nproma,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &       !< density of air 
      &  p_rho(:,:,:)               !< 
                                        !< [kg/m3]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &       !< tracer mixing ratios (specific concentrations)
      &  p_tracer_new(:,:,:,:)          !< after transport 
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    TYPE(t_list_element), POINTER :: current_element !< returns the reference to
                                                     !< current element in list
    TYPE(t_var_metadata), POINTER :: info            !< returns reference to tracer
                                                     !< metadata of current element
    INTEGER  :: jg                !< loop index
    INTEGER, POINTER :: jsp                          !< returns index of element
    CHARACTER(len=5), POINTER :: var_name            !< returns a character containing the name
                                                     !< of current ash component without the time level
                                                     !< suffix at the end. e.g. qash1(.TL1)

    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
    
    jg  = p_patch%id

      current_element=>p_prog_list%p%first_list_element

      !start DO-loop over elements in list:
      DO WHILE (ASSOCIATED(current_element))

      !get meta data of current element:
      info=>current_element%field%info

      ! assure that current element is tracer
      IF (info%tracer%lis_tracer) THEN
        IF (art_config(jg)%lart_wash_volcano .AND. info%tracer%lwash_tracer) THEN

         !
         ! retrieve  running index:
         !
          jsp=>info%ncontained
          var_name=>info%name

          WRITE(0,*) 'WASHOUT of ', var_name,' with idx= ',jsp

!!$         CALL art_washout_volc(dt_phy_jg,       & !>in
!!$         &          p_patch,                    & !>in
!!$         &          p_rain_gsp_rate,              & !>in
!!$         &          p_snow_gsp_rate,              & !>in
!!$         &          p_rain_con_rate,              & !>in
!!$         &          p_snow_con_rate,              & !>in
!!$         &          p_rho,                   & !>in
!!$         &          p_tracer_new,jsp)                !>inout 
        ENDIF
      ENDIF !lis_tracer

      ! select the next element in the list
      current_element => current_element%next_list_element

     ENDDO !loop elements

#endif

  END SUBROUTINE art_washout_interface


END MODULE mo_art_washout_interface

