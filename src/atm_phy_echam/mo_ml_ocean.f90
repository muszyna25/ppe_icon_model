!>
!! @brief Subroutine ml_ocean computes the surface temperature of a mixed
!! layer slab ocean based only on the vertical energy budget.  currently
!! there is no land/lakes, no ice, and only a very simple Qflux correction
!!
!! original purpose: use for a Radiative Convective Equilibrium APE exp
!! 
!! @author Levi G Silvers, MPI-M; Andrew I L Williams, Oxford
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_ml_ocean

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rho_ref  
  USE mo_nh_testcases_nml,   ONLY: qflux_case, dmixsea
  USE mo_model_domain,       ONLY: p_patch
  USE mo_ml_qflux,           ONLY: ml_qflux

  IMPLICIT NONE
  PUBLIC :: ml_ocean

CONTAINS

  SUBROUTINE ml_ocean ( jg, kproma, jcs, jce, jb, pdtime,pahflw,pahfsw,ptrflw,psoflw,ptsw )

    !
    ! Arguments
    !
    INTEGER, INTENT(IN)    :: jg
    INTEGER, INTENT(IN)    :: kproma
    INTEGER, INTENT(IN)    :: jcs, jce
    INTEGER, INTENT(IN)    :: jb
    REAL(wp),INTENT(IN)    :: pdtime                   !< time step
    REAL(wp),INTENT(IN)    :: pahflw(kproma)           ! latent heat flux over water
    REAL(wp),INTENT(IN)    :: pahfsw(kproma)           ! sensible heat flux over water
    REAL(wp),INTENT(IN)    :: ptrflw(kproma)           ! LW flux over water
    REAL(wp),INTENT(IN)    :: psoflw(kproma)           ! SW flux over water
    REAL(wp),INTENT(INOUT) :: ptsw(kproma)             ! surface temperature of water
    !
    ! Local variables
    !
    INTEGER  :: jl
    REAL(wp) :: cpsea,zmixcap,zfluxw !dmixsea
    REAL(wp) :: qflux_val, zlat, zlon
   
    !$ACC DATA PRESENT( pahflw, pahfsw, ptrflw, psoflw, ptsw )

    !
    ! Parameters/constants set in ml_flux(ECHAM) originally:
    !  zdmix=50._dp
    !  zcpwater=3994._dp
    !  zrho_sea=1025._dp
    !  zmixcap=zrho_sea*zcpwater*zdmix
    !  zmonlen=2592000._dp
    !
    !dmixsea   = 50._wp   ! [m] depth of mixed layer (zdmix)
    
    cpsea     = 3994._wp ! [J/(kg K)] specific heat capacity of water (what temp? looks low to me)

    zmixcap   = rho_ref*cpsea*dmixsea ! mixed layer heat capacity? 

    !$ACC PARALLEL DEFAULT(PRESENT)
    !$ACC LOOP GANG VECTOR PRIVATE( zfluxw )
    DO jl = jcs, jce

      zlat = p_patch(jg)%cells%center(jl,jb)%lat 
      zlon = p_patch(jg)%cells%center(jl,jb)%lon

      qflux_val          = ml_qflux(qflux_case, zlat, zlon)

      zfluxw             = pahflw(jl) + pahfsw(jl) + ptrflw(jl) + psoflw(jl) - qflux_val
      ptsw(jl)           = ptsw(jl) + pdtime * zfluxw/zmixcap

    END DO
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE ml_ocean

END MODULE mo_ml_ocean
