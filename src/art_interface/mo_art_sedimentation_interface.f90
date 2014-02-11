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

    USE mo_kind,                   ONLY: wp
    USE mo_parallel_config,        ONLY: nproma
    USE mo_model_domain,           ONLY: t_patch
    USE mo_impl_constants,         ONLY: min_rlcell
    USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_metrics,t_nh_diag
    USE mo_art_config,             ONLY: art_config
    USE mo_exception,              ONLY: message, message_text, finish
    USE mo_linked_list,            ONLY: t_var_list, t_list_element
    USE mo_var_metadata_types,     ONLY: t_var_metadata, t_tracer_meta
    USE mo_advection_vflux,        ONLY: upwind_vflux_ppm_cfl
    USE mo_run_config,             ONLY: ntracer
    USE mo_loopindices,            ONLY: get_indices_c
#ifdef __ICON_ART
!    USE mo_art_sedi_volc,          ONLY: art_sedi_volc
!    USE mo_art_aerosol,            ONLY: p_mflx_contra_vsed, vdep_ash
!    USE mo_art_aerosol,            ONLY: p_art_mode,nmodes,imode_seasa,imode_seasb,imode_seasc
!    USE mo_art_sedi_depo,          ONLY: art_calc_v_sed_dep
!    USE mo_art_aerosol_utilities,  ONLY: art_modal_parameters,art_air_properties
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

  SUBROUTINE art_sedi_interface( p_patch, &
             &                   p_dtime, &
             &                   p_prog_list, p_prog, &
             &                   p_metrics, p_rho, p_diag, &
             &                   p_tracer_new, &
             &                   p_cellhgt_mc_now, p_rhodz_new, &
             &                   lprint_cfl, &
             &                   opt_topflx_tra )


    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation
      &  p_patch                             !< is performed

    REAL(wp), INTENT(IN)              :: &
      &  p_dtime

    TYPE(t_var_list), INTENT(IN)      :: &   !< current prognostic state list
     &  p_prog_list

    TYPE(t_nh_prog), INTENT(IN)       :: &   !< current prognostic state
     &  p_prog

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
      &  p_rhodz_new(:,:,:)                  !< at full levels and time step n+1 [kg/m**2]
                                             !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &               !< determines if vertical CFL number shall be printed
      &  lprint_cfl                          !< in routine upwind_vflux_ppm_cfl

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

    INTEGER, ALLOCATABLE  :: idx_trac_arr(:)     !< Array to map jsp of tracer list element to idx_trac
    INTEGER  :: idx_trac

    !-----------------------------------------------------------------------

#ifdef __ICON_ART


#endif

  END SUBROUTINE art_sedi_interface


END MODULE mo_art_sedi_interface

