!>
!! Provides interface to ART-routines dealing with washout
!!
!! This module provides an interface to the ART-routines.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Daniel Rieger, KIT
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-06-15)
!! Rewritten by Daniel Rieger, KIT (2013-30-09)
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
MODULE mo_art_washout_interface

    USE mo_kind,                ONLY: wp
    USE mo_model_domain,        ONLY: t_patch
    USE mo_art_config,          ONLY: art_config
    USE mo_exception,           ONLY: finish
    USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
    USE mo_run_config,          ONLY: iqr,iqs
    USE mo_nonhydro_state,      ONLY: p_nh_state

#ifdef __ICON_ART
! Infrastructure Routines
    USE mo_art_modes_linked_list,  ONLY: p_mode_state,t_mode, &
        &                               t_fields_2mom,t_fields_1mom
    USE mo_art_data,            ONLY: p_art_data
    USE mo_art_aerosol_utilities,  ONLY: art_air_properties
    USE mo_art_clipping,           ONLY: art_clip_tracers_zero
! Washout Routines
    USE mo_art_washout_volc,       ONLY:art_washout_volc
    USE mo_art_radioactive,        ONLY:art_washout_radioact
    USE mo_art_washout_aerosol,    ONLY:art_aerosol_washout
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_washout_interface


CONTAINS


  !>
  !! Interface for ART-routines dealing with washout
  !!
  !! This interface calls the ART-routines, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!

  SUBROUTINE art_washout_interface(            & !>in
      &          dt_phy_jg,                    & !>in
      &          p_patch,                      & !>in
      &          prm_diag,                     & !>in
      &          p_rho,                        & !>in
      &          p_tracer_new)                  !>inout

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN) ::dt_phy_jg         !< time interval, fast physics

    TYPE(t_nwp_phy_diag),       INTENT(IN) ::& !< diagnostic variables
     &  prm_diag 

    REAL(wp), INTENT(IN) ::  &          !< density of air 
      &  p_rho(:,:,:)                   !< 
                                        !< [kg/m3]
                                        !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(INOUT) ::  &       !< tracer mixing ratios (specific concentrations)
      &  p_tracer_new(:,:,:,:)          !< after transport 
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)

    INTEGER  :: jg                      !< patch id

    TYPE(t_mode), POINTER   :: this_mode
    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART

  jg  = p_patch%id

  IF(art_config(jg)%lart .AND. art_config(jg)%lart_wash) THEN 
      
    CALL art_air_properties(p_patch,p_art_data(jg))
    
    this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
   
    DO WHILE(ASSOCIATED(this_mode))
      ! Select type of mode
      select type (fields=>this_mode%fields)
        type is (t_fields_2mom)
          ! Before washout, the modal parameters have to be calculated
          call fields%modal_param(p_art_data(jg),p_patch,p_tracer_new)
          ! calc washout
          call art_aerosol_washout(fields,p_art_data(jg),     &
                 &                 p_patch,dt_phy_jg,p_rho,   &
                 &                 p_tracer_new,prm_diag,p_nh_state(jg))
        type is (t_fields_1mom)
          ! do the stuff
        class default
          CALL finish('mo_art_washout_interface:art_washout_interface', &
               &      'ART: Unknown mode field type')
      end select
      this_mode => this_mode%next_mode
    END DO
   
   
  ! ----------------------------------
  ! --- Clip the tracers as the washout
  ! --- can produce negative values when
  ! --- timestep is too high
  ! ----------------------------------
  
  CALL art_clip_tracers_zero(p_tracer_new)
  
  
  
  
! OLD STUFF  
  
  
      ! ----------------------------------
      ! --- second for monodisperse aerosol
      ! ----------------------------------
    
!      current_element=>p_prog_list%p%first_list_element

      ! ----------------------------------
      ! --- start DO-loop over elements in list:
      ! ----------------------------------

!      DO WHILE (ASSOCIATED(current_element))

      ! ----------------------------------
      ! --- get meta data of current element:
      ! ----------------------------------

!      info=>current_element%field%info

      ! ----------------------------------
      ! ---  assure that current element is tracer
      ! ----------------------------------

!      IF (info%tracer%lis_tracer) THEN
!        IF (info%tracer%lwash_tracer) THEN

          ! ----------------------------------
          ! --- retrieve  running index:
          ! ----------------------------------

!          jsp=>info%ncontained
!          var_name=>info%name

          ! ----------------------------------
          ! --- Choose parameterisation
          ! ----------------------------------

!          SELECT CASE(info%tracer%tracer_class)

!          CASE('volcash')

!             CALL art_washout_volc(dt_phy_jg,       & !>in
!             &          p_patch,                    & !>in
!             &          prm_diag,                   & !>in
!             &          p_rho,                   & !>in
!             &          p_tracer_new,jsp)                !>inout 

!          CASE('radioact')

!             CALL art_washout_radioact(p_patch,    & !>in
!             &          dt_phy_jg,                 & !>in
!             &          prm_diag,                  & !>in
!             &          p_tracer_new(:,:,:,jsp),   & !>inout
!             &          info%tracer%imis_tracer,       & !>in  needed for nuclide specific factors
!             &          p_tracer_new(:,:,:,iqr),   & !>in
!             &          p_tracer_new(:,:,:,iqs)    ) !>in

!          END SELECT
          
!        ENDIF
!      ENDIF !lis_tracer

      ! ----------------------------------
      ! --- select the next element in the list
      ! ----------------------------------

!      current_element => current_element%next_list_element

!     ENDDO !loop elements

ENDIF
#endif

  END SUBROUTINE art_washout_interface


END MODULE mo_art_washout_interface

