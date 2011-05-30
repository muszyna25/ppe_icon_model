!>
!!  Namelist for surface physics
!!
!!  these Subroutines are called by control model and construct the
!!  surface scheme composition
!!
!! @author <Kristina Froehlich, DWD>
!!
!!
!! @par Revision History
!! First implementation by Kristina Froehlich, DWD (2010-06-20>)
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
MODULE mo_lnd_nwp_nml

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_run_nml,            ONLY: iforcing, inwp
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_io_units,           ONLY: nnml

  IMPLICIT NONE

  PRIVATE

!> Action Variables for physical schemes
! --------------------------------------
  INTEGER ::  nlev_soil, nztlev    !! number of soil layers, time integration scheme
  INTEGER ::  nlev_snow        !! number of snow layers
  INTEGER ::  nsfc_subs        !! number of TILES

  REAL(wp), DIMENSION(8):: zml_soil=(/ 0.005_wp,0.02_wp,0.06_wp,0.18_wp,0.54_wp,1.62_wp,&
    & 4.86_wp,14.58_wp /)

  LOGICAL ::       &
       lseaice,    & !> forecast with sea ice model
       llake,      & !! forecst with lake model FLake
       lmelt     , & !! soil model with melting process
       lmelt_var , & !! freezing temperature dependent on water content
       lmulti_snow   !! run the multi-layer snow model

!--------------------------------------------------------------------
! nwp forcing (right hand side)
!--------------------------------------------------------------------

  NAMELIST/lnd_ctl/ nlev_soil, nztlev, nlev_snow, nsfc_subs, &
 &                  lseaice, llake, lmulti_snow  
   
   PUBLIC :: setup_nwp_lnd
   PUBLIC :: zml_soil
   PUBLIC :: nlev_soil, nztlev
   PUBLIC :: nlev_snow, nsfc_subs
   PUBLIC :: lseaice, llake, lmulti_snow


 CONTAINS



  !-------------------------------------------------------------------------
  !
  !>
  !! Setup NWP physics
  !!
  !! Read namelist for physics. Choose the physical package and subsequent
  !! parameters.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-10-06)
  !!
  SUBROUTINE setup_nwp_lnd

    INTEGER :: i_status

    !-------------------------------------------------------------------------

    nlev_soil       = 7     !> 7 = default value for number of soil layers
    nztlev          = 2     !> 2 = default value for time integration scheme
    nlev_snow       = 1     !> 0 = default value for number of snow layers
    nsfc_subs       = 2     !> 1 = default value for number of TILES



  !> KF  current settings to get NWP turbulence running
        lseaice    = .FALSE.
        llake      = .FALSE.
        lmulti_snow= .FALSE.


        CALL position_nml ('lnd_ctl', status=i_status)
           IF (i_status == POSITIONED) THEN
              READ (nnml, lnd_ctl)
           ENDIF


  END SUBROUTINE setup_nwp_lnd



!==============================================================================

END MODULE mo_lnd_nwp_nml

