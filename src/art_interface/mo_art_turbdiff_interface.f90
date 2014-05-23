!>
!! Provides interface to ART-routines dealing with turbulent diffusion
!! using the COSMO turbulence scheme of Matthias Raschendorfer.
!!
!! This module provides an interface to the ART-routines.
!! The interface is written in such a way, that ICON will compile and run
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Jochen Foerstner, DWD
!!
!! @par Revision History
!! Initial revision by Jochen Foerstner, DWD (2013-09-02)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_turbdiff_interface


  USE mo_kind,                          ONLY: wp
  USE mo_parallel_config,               ONLY: nproma
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_art_config,                    ONLY: art_config
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_nonhydro_types,                ONLY: t_nh_metrics, t_nh_diag, t_nh_prog
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE src_turbdiff,                     ONLY: modvar
#ifdef __ICON_ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_surface_value,             ONLY: art_surface_value
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_turbdiff_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_turbdiff_interface( defcase,  & !>in
    &          p_patch,                      & !>in
    &          p_prog_rcf,                   & !>in
    &          prm_nwp_tend,                 & !>in
    &          ptr,                          & !>out
    &          p_metrics, p_diag, prm_diag,  & !>in, optional
    &          jb,                           & !>in, optional
    &          opt_sv, opt_fc,               & !>in, optional
    &          i_st, i_en, dt )                !>in, optional
!! Interface for ART-routines dealing with turbdiff
!!
!! This interface calls the ART-routines, if ICON has been
!! built including the ART-package. Otherwise, this is simply a dummy
!! routine.
!!
  CHARACTER(len=*), INTENT(in)              :: &
    &  defcase                           !< definition of case
  TYPE(t_patch), TARGET, INTENT(IN)         :: &
    &  p_patch                           !< patch on which computation is performed
  TYPE(t_nh_prog),  INTENT(in)              :: &
    &  p_prog_rcf                        !< the prog vars
  TYPE(t_nwp_phy_tend), INTENT(in)          :: &
    &  prm_nwp_tend                      !< atm phys tendencies
  TYPE(modvar), DIMENSION(:), INTENT(inout) :: &
    &  ptr                               !< passive tracer pointer type structure for diffusion
  TYPE(t_nh_metrics), INTENT(in), OPTIONAL  :: &
    &  p_metrics                         !< metrical fields
  TYPE(t_nh_diag), INTENT(in), OPTIONAL     :: &
    &  p_diag                            !< the diag vars
  TYPE(t_nwp_phy_diag), INTENT(in), OPTIONAL:: &
    &  prm_diag                          !< atm phys vars
  INTEGER, INTENT(in), OPTIONAL             :: &
    &  jb, i_st, i_en
  REAL(wp), DIMENSION(:,:,:), INTENT(out), TARGET, OPTIONAL :: &
    &  opt_sv                            !< surface value according to transfer coeff.
  LOGICAL, INTENT(in), OPTIONAL             :: &
    &  opt_fc
  REAL(wp), INTENT(in), OPTIONAL            :: &
    &  dt
  REAL(wp),POINTER                 ::  &
    &  sv(:,:,:),                      & !< surface value of tracer
    &  vdep(:,:,:)                       !< deposition velocity of tracer
    
  !Local variables
  ! ---------------------------------------

  INTEGER  :: jg, idx_trac, jk, jc       !< loop indices
  INTEGER  :: nblks, istat, nlev, i_startidx, i_endidx

  !-----------------------------------------------------------------------

#ifdef __ICON_ART

  jg  = p_patch%id
  IF ( art_config(jg)%lart ) THEN
    SELECT CASE(TRIM(defcase))
    
    CASE('setup_ptr')
      sv => p_art_data(jg)%turb_fields%sv
      vdep => p_art_data(jg)%turb_fields%vdep

      IF ( .NOT. PRESENT(opt_sv) ) THEN
        CALL art_surface_value( p_patch, p_prog_rcf, p_metrics, p_diag, prm_diag, &
          &                     jb, vdep, sv )
      END IF

      DO idx_trac = 1, art_config(jg)%nturb_tracer

        ! set up pointer to tracer type structure for diffusion
        ptr(idx_trac)%av => p_prog_rcf%turb_tracer(jb,idx_trac)%ptr
        ptr(idx_trac)%at => prm_nwp_tend%turb_tracer_tend(jb,idx_trac)%ptr
        ptr(idx_trac)%at =  0._wp
  
        IF ( PRESENT(opt_sv) ) THEN
          ptr(idx_trac)%sv => opt_sv(:,jb,idx_trac)
          IF ( PRESENT(opt_fc) ) THEN
            ptr(idx_trac)%fc = opt_fc
          ELSE
            ptr(idx_trac)%fc = .FALSE.
          END IF
        ELSE
          ptr(idx_trac)%sv => sv(:,jb,idx_trac)
          ptr(idx_trac)%fc = .FALSE.
        END IF
      END DO

    CASE('update_ptr')
      nlev = p_patch%nlev !< Number of vertical full levels
      ! update tracers due to diffusion
      DO idx_trac = 1, art_config(jg)%nturb_tracer
        DO jk = 1, nlev
  !DIR$ IVDEP
          DO jc = i_st, i_en
            ptr(idx_trac)%av(jc,jk) = MAX( 0._wp, ptr(idx_trac)%av(jc,jk)     &
              &                     + dt * ptr(idx_trac)%at(jc,jk) )
          END DO
        END DO
      END DO
    
    END SELECT
  END IF

#endif

END SUBROUTINE art_turbdiff_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_turbdiff_interface

