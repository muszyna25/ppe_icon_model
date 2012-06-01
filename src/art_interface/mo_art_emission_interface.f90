!>
!! Provides interface to ART-routines dealing with emissions of volcanic ash particles
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
MODULE mo_art_emission_interface

    USE mo_kind,                ONLY: wp
    USE mo_model_domain,        ONLY: t_patch
    USE mo_art_config,          ONLY: art_config
    USE mo_exception,             ONLY: message, message_text, finish
    USE mo_datetime,             ONLY:t_datetime
#ifdef __ICON_ART
    USE mo_art_emission_volc,       ONLY:art_organize_emission_volc
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
  SUBROUTINE art_emission_interface( p_patch,p_dtime,datetime,p_rho_now,p_tracer_now)


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed
    REAL(wp), INTENT(IN) ::p_dtime
     TYPE(t_datetime), INTENT(IN)::datetime

    REAL(wp), INTENT(INOUT) ::  &  !<density of air 
      &  p_rho_now(:,:,:)          !< at current time level n (before transport)
                                   !< [kg/m3]
                                   !< dim: (nproma,nlev,nblks_c)
    REAL(wp), INTENT(INOUT) ::  &  !< tracer mixing ratios (specific concentrations)
      &  p_tracer_now(:,:,:,:)          !< at current time level n (before transport)
                                        !< [kg/kg]
                                        !< dim: (nproma,nlev,nblks_c,ntracer)
!    REAL(wp), INTENT(IN) ::          &  !< cell height defined at full levels for
!      &  p_cellhgt_mc_now(:,:,:)        !< time step n
!                                        !< NH: \Delta z       [m]
!                                        !< HA: \Delta p       [Pa]
!                                        !< dim: (nproma,nlev,nblks_c)
!

    INTEGER  :: jg                !< loop index

    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
    
    jg  = p_patch%id
     
     IF (art_config(jg)%lart_volc) THEN
      CALL art_organize_emission_volc(p_patch,p_dtime,datetime,p_rho_now,p_tracer_now) 
     ENDIF


#endif

  END SUBROUTINE art_emission_interface


END MODULE mo_art_emission_interface

