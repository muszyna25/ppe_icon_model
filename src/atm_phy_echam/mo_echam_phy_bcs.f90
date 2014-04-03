!>
!! @brief Prepares boundary conditions needed for ECHAM physics
!!
!! @author Marco Giorgetta (MPI-M)
!!
!! @par Revision History
!!
!! @par Copyright
!! 2002-2014 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!      violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!      copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_echam_phy_bcs

  USE mo_kind                       ,ONLY: wp
  USE mo_math_constants             ,ONLY: pi
!!$  USE mo_datetime                   ,ONLY: t_datetime, add_time , OPERATOR(==) ! for t_datetime typed variables
  USE mo_datetime                   ,ONLY: t_datetime, add_time
  USE mo_model_domain               ,ONLY: t_patch

  USE mo_master_nml                 ,ONLY: lrestart
!!$  USE mo_time_config                ,ONLY: time_config
  USE mo_echam_phy_config           ,ONLY: echam_phy_config
  USE mo_radiation_config           ,ONLY: ighg, isolrad, tsi, tsi_radt, ssi_radt, irad_o3, irad_aero

  USE mo_greenhouse_gases           ,ONLY: ghg_time_interpolation

  USE mo_amip_bc                    ,ONLY: get_current_amip_bc_year, read_amip_bc, amip_time_weights, amip_time_interpolation
  USE mo_icoham_sfc_indices         ,ONLY: iwtr

  USE mo_time_interpolation         ,ONLY: time_weights_limm
  USE mo_time_interpolation_weights ,ONLY: wi_limm

  USE mo_solar_irradiance           ,ONLY: read_ssi_bc, ssi_time_interpolation

  USE mo_impl_constants             ,ONLY: io3_amip
  USE mo_o3                         ,ONLY: read_amip_o3

  USE mo_aero_kinne                 ,ONLY: read_aero_kinne
  USE mo_aero_stenchikov            ,ONLY: read_aero_stenchikov

  USE mo_echam_phy_memory           ,ONLY: prm_field

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_bcs_global

  CHARACTER(len=*), PARAMETER :: version    = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_echam_phy_bcs'

CONTAINS
  !>
  !! SUBROUTINE echam_phy_bcs_global
  !!
  !! This subroutine is called in the time loop of the model and serves
  !! to prepare the boundary conditions for the ECHAM physics.
  !!
  !! This routine contains all time dependent preparations that are valid
  !! for the whole globe and need to be done outside of the loop over
  !! rows, in which the echam_phy_main routine is called for each row
  !! of the block.
  !!
  !! Note that each call of this subroutine deals with a single grid
  !! with index jg rather than the entire grid tree.

  SUBROUTINE echam_phy_bcs_global( datetime     ,&! in
    &                              jg           ,&! in
    &                              patch        ,&! in
    &                              dtadv_loc    ,&! in
    &                              ltrig_rad    ,&! out
    &                              time_radtran ) ! out

    ! Arguments

    TYPE(t_datetime)         ,INTENT(in)    :: datetime      !< date and time fo this timestep
    INTEGER                  ,INTENT(in)    :: jg            !< grid index
    TYPE(t_patch)    ,TARGET ,INTENT(in)    :: patch         !< description of grid jg
    REAL(wp)                 ,INTENT(in)    :: dtadv_loc     !< timestep of advection and physics on grid jg
    LOGICAL                  ,INTENT(out)   :: ltrig_rad     !< trigger for radiation transfer computation
    REAL(wp)                 ,INTENT(out)   :: time_radtran  !< time of day (in radian) at which radiative transfer is computed

    ! Local variables

    TYPE(t_datetime) :: datetime_radtran  !< date and time of zenith angle for radiative transfer comp.
    REAL(wp)         :: dsec              !< [s] time increment of datetime_radtran wrt. datetime

!!$    LOGICAL          :: is_initial_datetime
!!$    LOGICAL          :: is_radtran_datetime
    LOGICAL          :: is_1st_call = .TRUE.

    ! Local parameters

!!$    CHARACTER(*), PARAMETER :: method_name = "echam_phy_bcs_global"

    !-------------------------------------------------------------------------
    ! Prepare some global parameters or parameter arrays
    !-------------------------------------------------------------------------

    ! Check whether the radiative transfer needs to be calculated in this
    ! timestep. Then boundary conditions must be prepared for this purpose.
    !
    IF (echam_phy_config%lrad) THEN

!!$      is_initial_datetime = (datetime == time_config%ini_datetime) ! here the overloaded == from mo_datetime is used
!!$      is_radtran_datetime = (MOD(NINT(datetime%daysec),NINT(echam_phy_config%dt_rad)) == 0)
!!$      ltrig_rad = ( is_initial_datetime .OR. is_radtran_datetime )
      ltrig_rad   = ( is_1st_call .AND. (.NOT.lrestart)                             ) .OR. &
        &           ( MOD(NINT(datetime%daysec),NINT(echam_phy_config%dt_rad)) == 0 )

    ELSE
      ltrig_rad = .FALSE.
    END IF

    ! Set the time instance datetime_radtran for the zenith angle to be used
    ! in the radiative transfer. All other input for the radiative transfer
    ! is for datatime, i.e. the current timestep.
    !
    datetime_radtran = datetime                         ! copy current date and time
    IF (ltrig_rad) THEN
      dsec = 0.5_wp*(echam_phy_config%dt_rad - dtadv_loc) ! [s] time increment for zenith angle
      CALL add_time(dsec,0,0,0,datetime_radtran)          ! datetime_radtran = datetime_radtran + dsec
    END IF
    time_radtran  = 2._wp*pi * datetime_radtran%daytim  ! time of day in radian

    ! Read and interpolate in time monthly mean SST for AMIP simulations
    ! SST is needed for turbulent vertical fluxes and for radiation.
    !
    IF (echam_phy_config%lamip) THEN
      IF (datetime%year /= get_current_amip_bc_year()) THEN
        CALL read_amip_bc(datetime%year, patch)
      END IF
      CALL amip_time_weights(datetime)
      CALL amip_time_interpolation( prm_field(jg)%seaice(:,:) ,&
        &                           prm_field(jg)%tsurfw(:,:) ,&
        &                           prm_field(jg)%siced(:,:)  ,&
        &                           prm_field(jg)%lsmask(:,:) )

      prm_field(jg)%tsfc_tile(:,:,iwtr) = prm_field(jg)%tsurfw(:,:)

      ! The ice model should be able to handle different thickness classes, 
      ! but for AMIP we ONLY USE one ice class.
      prm_field(jg)%conc(:,1,:) = prm_field(jg)%seaice(:,:)
      prm_field(jg)%hi  (:,1,:) = prm_field(jg)%siced (:,:)
    END IF

    ! interpolation weights for linear interpolation
    ! of monthly means onto the actual integration time step
    CALL time_weights_limm(datetime, wi_limm)


    ! total solar irradiation at the mean sun earth distance
    IF (isolrad==1) THEN
      CALL read_ssi_bc(datetime%year,.FALSE.)
      CALL ssi_time_interpolation(wi_limm,.FALSE.,tsi)
    END IF

    ! quantities needed for the radiative transfer only
    !
    IF (ltrig_rad) THEN
      !
      ! total and spectral solar irradiation at the mean sun earth distance
      IF (isolrad==1) THEN
        CALL read_ssi_bc(datetime%year,.TRUE.)
        CALL ssi_time_interpolation(wi_limm,.TRUE.,tsi_radt,ssi_radt)
      END IF
      !
      ! greenhouse gas concentrations, assumed constant in horizontal dimensions
      IF (ighg > 0) THEN
        CALL ghg_time_interpolation(datetime)
      END IF
      !
      ! ozone concentration
      IF (irad_o3 == io3_amip) THEN
        CALL read_amip_o3(datetime%year, patch)
      END IF
      !
      ! tropospheric aerosol optical properties
      IF (irad_aero == 13) THEN
        CALL read_aero_kinne(datetime%year, patch)
      END IF
      !
      ! stratospheric aerosol optical properties
      IF (irad_aero == 15) THEN
        CALL read_aero_kinne     (datetime%year, patch)
        CALL read_aero_stenchikov(datetime%year, patch)
      END IF
      !
    END IF ! ltrig_rad

    is_1st_call = .FALSE.

  END SUBROUTINE echam_phy_bcs_global

END MODULE mo_echam_phy_bcs
