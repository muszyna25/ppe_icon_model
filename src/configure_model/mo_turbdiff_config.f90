!>
!! @brief configuration setup for turbulent diffusion (turbdiff)
!!
!! configuration setup for turbulent diffusion
!!
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-09-19)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_turbdiff_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for turbulent diffusion (turbdiff)
  !!--------------------------------------------------------------------------
  TYPE :: t_turbdiff_config

    ! namelist variables

    INTEGER :: &   ! type of surface-atmosphere transfer
      &  itype_tran
    INTEGER :: &   ! mode of surface-atmosphere transfer
      &  imode_tran 
    INTEGER :: &   ! mode of cloud representation in transfer parametr.
      &  icldm_tran
    INTEGER :: &   ! mode of turbulent diffusion parametrization
      &  imode_turb
    INTEGER :: &   ! mode of cloud representation in turbulence parametr.
      &  icldm_turb
    INTEGER :: &   ! type of shear production for TKE
      &  itype_sher

    LOGICAL :: &   ! calculation SSO-wake turbulence production for TKE
      &  ltkesso
    LOGICAL :: &   ! consider convective buoyancy production for TKE
      &  ltkecon
    LOGICAL :: &   ! explicit corrections of the implicit calculated turbul. diff.
      &  lexpcor
    LOGICAL :: &   ! consideration of thermal TKE-sources in the enthalpy budget
      &  ltmpcor
    LOGICAL :: &   ! using the profile values of the lowest main level instead of
      &  lprfcor   ! the mean value of the lowest layer for surface flux calulations

    LOGICAL :: &   ! nonlocal calculation of vertical gradients used for turbul. diff.
      &  lnonloc 
    LOGICAL :: &   ! consideration of fluctuations of the heat capacity of air
      &  lcpfluc
    LOGICAL :: &   ! use semi-implicit TKE diffusion
      &  limpltkediff

    LOGICAL :: &   ! TRUE: horizontally homogeneous roughness length 
      &  lconst_z0 ! (for idealized testcases)

    REAL(wp):: &   ! horizontally homogeneous roughness length 
      &  const_z0  ! (for idealized testcases)


    !
    ! Switches controlling other physical parameterizations:
    !
    INTEGER :: &  ! type of water cloud diagnosis
      &  itype_wcld
    INTEGER :: &  ! type of diagnostics of synoptical near surface variables
      &  itype_synd


    !
    ! derived variables
    !

                   
  END TYPE t_turbdiff_config

  !>
  !!
  TYPE(t_turbdiff_config), TARGET :: turbdiff_config(0:max_dom)


!CONTAINS


END MODULE mo_turbdiff_config
