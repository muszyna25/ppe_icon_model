!>
!! @brief Subroutine ml_ocean computes the surface temperature of a mixed
!! layer slab ocean based only on the vertical energy budget.  currently
!! there is no land/lakes, no ice, and no flux heat correction
!!
!! original purpose: use for a Radiative Convective Equilibrium APE exp
!!
!! @author Levi G Silvers, MPI-M
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
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

MODULE mo_ml_ocean

  USE mo_kind,               ONLY: wp
  USE mo_physical_constants, ONLY: rho_ref  

  IMPLICIT NONE
  PUBLIC :: ml_ocean

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

  SUBROUTINE ml_ocean ( kproma, start_column, end_column, pdtime,pahflw,pahfsw,ptrflw,psoflw,ptsw )

    !
    ! Arguments
    !
    INTEGER, INTENT(IN)    :: kproma
    INTEGER, INTENT(IN)    :: start_column, end_column
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
    REAL(wp) :: dmixsea,cpsea,zmixcap,zfluxw

    !
    ! Parameters/constants set in ml_flux(ECHAM) originally:
    !  zdmix=50._dp
    !  zcpwater=3994._dp
    !  zrho_sea=1025._dp
    !  zmixcap=zrho_sea*zcpwater*zdmix
    !  zmonlen=2592000._dp
    !
    dmixsea   = 50._wp   ! [m] depth of mixed layer (zdmix)
    cpsea     = 3994._wp ! [J/(kg K)] specific heat capacity of water (what temp? looks low to me)

    zmixcap   = rho_ref*cpsea*dmixsea ! mixed layer heat capacity? 

    DO jl = start_column,  end_column

      zfluxw             = pahflw(jl) + pahfsw(jl) + ptrflw(jl) + psoflw(jl)
      ptsw(jl)           = ptsw(jl) + pdtime * zfluxw/zmixcap

    END DO

  END SUBROUTINE ml_ocean

END MODULE mo_ml_ocean
