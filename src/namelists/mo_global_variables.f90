!>
!!                 Contains global variables used by the shallow water model.
!!
!!                 More comments on many of the variable defined here
!!                 to be found in the <i>user_introduction</i>.
!!
MODULE mo_global_variables

!  USE mo_kind,               ONLY: wp
!  USE mo_io_units,           ONLY: nnml, nnml_output
!  USE mo_namelist,           ONLY: position_nml, positioned
!  USE mo_mpi,                ONLY: p_pe, p_io
!  USE mo_ocean_nml,          ONLY: lmpiom_radiation, lmpiom_convection,           &
!    &                              lmpiom_gentmcwill

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC
!--------------------------------------------------------------------

  CONTAINS


 SUBROUTINE setup_physics
!!-----------------------------------------------------------------------
!! read namelist for physics
!! (chose the physical package and subsequent parameters)
!!-----------------------------------------------------------------------
!!
!  INTEGER :: i_status
!
!  SELECT CASE (iforcing)
!  CASE (impiom)
!    !
!    lmpiom_radiation  = .FALSE.
!    lmpiom_convection = .FALSE.
!    lmpiom_gentmcwill = .FALSE.
!    !
!    CALL position_nml ('mpiom_phy_ctl', status=i_status)
!    !
!    SELECT CASE (i_status)
!    CASE (positioned)
!      READ (nnml, mpiom_phy_ctl)
!    END SELECT
!    !  write the contents of the namelist to an ASCII file
!    IF(p_pe == p_io) WRITE(nnml_output,nml=mpiom_phy_ctl)
!    !
!   CASE DEFAULT
!    !
!  END SELECT
!!
END SUBROUTINE setup_physics

!-------------------------------------------------------------------------
!
END MODULE mo_global_variables
