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


 
#ifdef __ICON_ART

#endif

  END SUBROUTINE art_sedi_interface


END MODULE mo_art_sedi_interface

