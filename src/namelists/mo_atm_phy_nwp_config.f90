!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_atm_phy_nwp_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: max_dom,MAX_CHAR_LENGTH,itconv,itccov,&
    &                               itrad,itradheat, itsso,itgscp,itsatad,itupdate,&
    &                               itturb, itsfc,  itgwd, iphysproc
  USE mo_exception,        ONLY: message, message_text, finish

  USE mo_data_turbdiff,       ONLY: imode_turb,                              &
    &                               limpltkediff, ltkesso, lexpcor,          &
    &                               tur_len, pat_len, a_stab,                &
    &                               tkhmin, tkmmin, c_diff,                  &
    &                               itype_wcld, icldm_turb,                  &
    &                               itype_tran, rlam_heat, rlam_mom, rat_sea,&
    &                               llake, lseaice 

  USE mo_icoham_sfc_indices,  ONLY: init_sfc_indices, nsfc_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_phy_nwp_config, t_atm_phy_nwp_config, tcall_phy
  PUBLIC :: configure_atm_phy_nwp

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for atm dynamics
  !!--------------------------------------------------------------------------
  TYPE :: t_atm_phy_nwp_config

    ! namelist variables

  INTEGER ::  inwp_gscp        !> microphysics
  INTEGER ::  inwp_satad       !! saturation adjustment
  INTEGER ::  inwp_convection  !! convection
  INTEGER ::  inwp_radiation   !! radiation
  INTEGER ::  inwp_sso         !! sso
  INTEGER ::  inwp_gwd         !! non-orographic gravity wave drag
  INTEGER ::  inwp_cldcover    !! cloud cover
  INTEGER ::  inwp_turb        !! turbulence
  INTEGER ::  inwp_surface     !! surface including soil, ocean, ice,lake

  REAL(wp) :: dt_conv   !> field element for convection
  REAL(wp) :: dt_ccov   !! field element for subscale cloud cover
  REAL(wp) :: dt_rad    !! "-"                     radiation
  REAL(wp) :: dt_radheat!! "-" rad. heating from radiative fluxes with updated cosmu0 
  REAL(wp) :: dt_sso    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gwd    !! "-"  for subscale gravity waves
  REAL(wp) :: dt_gscp   !! field element for microphysics
  REAL(wp) :: dt_turb   !! field element for turbulence
  REAL(wp) :: dt_sfc    !! field element for surface
  REAL(wp) :: dt_satad  !! field element for sat. adjustment
  REAL(wp) :: dt_update !! field element for tracer phys update

  LOGICAL ::  lseaice     !> forecast with sea ice model
  LOGICAL ::  llake       !! forecst with lake model FLake

  INTEGER :: imode_turb, itype_wcld, icldm_turb, itype_tran
  LOGICAL :: limpltkediff, ltkesso, lexpcor
  REAL(wp):: tur_len, pat_len, a_stab, tkhmin, tkmmin, c_diff, &
    &        rlam_heat, rlam_mom, rat_sea 

  REAL(wp):: qi0, qc0

  END TYPE t_atm_phy_nwp_config

  !>
  !!
  TYPE(t_atm_phy_nwp_config) :: atm_phy_nwp_config(max_dom) !< shape: (n_dom)


  REAL(wp) ::  &                    !> Field of calling-time interval (seconds) for
    &  tcall_phy(max_dom,iphysproc) !! each domain and phys. process


CONTAINS

SUBROUTINE configure_atm_phy_nwp(n_dom,ltestcase, iadv_rcf, dtime)
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
  !! revision for restructurring by Kristina Froehlich MPI-M (2011-07-12)


  REAL(wp),INTENT(IN) :: dtime
  INTEGER, INTENT(IN) :: n_dom
  INTEGER, INTENT(IN) :: iadv_rcf
  LOGICAL, INTENT(IN) :: ltestcase

  INTEGER :: jg
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
    &      routine = 'mo_atm_phy_nwp_config:configure_atm_phy_nwp'

    tcall_phy(:,:) = 0._wp

    DO jg = 1,n_dom
      ! Slow physics:
      ! currently for each domain the same time intervals are set
      !
      ! Fast physics:
      ! time interval are set equal to the time interval for advection
      !
      ! note that time intervalls are set to 0 if the physical process
      ! is not considered at all.
      IF (  atm_phy_nwp_config(jg)%inwp_convection == 0 ) THEN   ! 0 = no convection
        tcall_phy(jg,itconv) = 0._wp
      ELSE
        tcall_phy(jg,itconv) = atm_phy_nwp_config(jg)% dt_conv    ! seconds
      ENDIF

       !> KF always call clouds after convection
       !! to ensure the proper output of Qx_tot

      IF (  atm_phy_nwp_config(jg)% inwp_convection > 0 ) THEN   ! 0 = no convection
        tcall_phy(jg,itccov) =  atm_phy_nwp_config(jg)% dt_conv    ! seconds
      ELSE 
        tcall_phy(jg,itccov) =  atm_phy_nwp_config(jg)% dt_gscp    ! seconds
      ENDIF
      !really switch off the clouds
      IF (  atm_phy_nwp_config(jg)% inwp_cldcover == 0 ) THEN     ! 0 = no cloud cover
        tcall_phy(jg,itccov) = 0._wp
      ENDIF

!      ELSE
!        tcall_phy(jg,itccov) = dt_rad(jg)           ! cloud cover
!                                                    ! should be coupled
!      ENDIF                                         ! to convection
                                                    ! and radiation!

      IF (  atm_phy_nwp_config(jg)% inwp_radiation == 0 ) THEN    ! 0 = no radiation
        tcall_phy(jg,itrad)     =  0._wp
        tcall_phy(jg,itradheat) =  0._wp
      ELSE
        tcall_phy(jg,itrad)     =   atm_phy_nwp_config(jg)% dt_rad       ! seconds
        tcall_phy(jg,itradheat) =   atm_phy_nwp_config(jg)% dt_radheat   ! seconds       
      ENDIF

      IF (  atm_phy_nwp_config(jg)% inwp_sso == 0 ) THEN          ! 0 = no sso
        tcall_phy(jg,itsso) =  0._wp
      ELSE
        tcall_phy(jg,itsso) =   atm_phy_nwp_config(jg)% dt_sso           ! seconds
      ENDIF

      IF (  atm_phy_nwp_config(jg)% inwp_gwd == 0 ) THEN          ! 0 = no sso
        tcall_phy(jg,itgwd) =  0._wp
      ELSE
        tcall_phy(jg,itgwd) =   atm_phy_nwp_config(jg)% dt_gwd    ! seconds
      ENDIF

      IF (  atm_phy_nwp_config(jg)% inwp_gscp == 0 ) THEN         ! 0 = no microphysics
        tcall_phy(jg,itgscp) =  0._wp
      ELSE
        tcall_phy(jg,itgscp) =  atm_phy_nwp_config(jg)% dt_gscp    ! seconds
      ENDIF

      IF (  atm_phy_nwp_config(jg)%inwp_satad == 0 ) THEN         ! 0 = no satad
        tcall_phy(jg,itsatad)  =  0._wp
        tcall_phy(jg,itupdate) =  0._wp   ! no moist update needed if no satad
      ELSE
        tcall_phy(jg,itsatad)  =   atm_phy_nwp_config(jg)%dt_satad  !seconds
        tcall_phy(jg,itupdate) =   atm_phy_nwp_config(jg)%dt_update !seconds
      ENDIF

      IF (  atm_phy_nwp_config(jg)%inwp_turb == 0 ) THEN         ! 0 = no turbulence
        tcall_phy(jg,itturb) =  0._wp 
      ELSE
        tcall_phy(jg,itturb) =   atm_phy_nwp_config(jg)%dt_turb   !seconds
      ENDIF

      IF (  atm_phy_nwp_config(jg)%inwp_surface == 0 ) THEN       ! 0 = no soil
        tcall_phy(jg,itsfc) =  0._wp 
      ELSE
        tcall_phy(jg,itsfc) =   atm_phy_nwp_config(jg)%dt_turb  !seconds
      ENDIF

    ENDDO

    IF( atm_phy_nwp_config(1)%inwp_turb == 2) THEN
       CALL init_sfc_indices( ltestcase, 'APE' ) !call of a hydrostatic testcase
                                             ! to obtain the demanded parameters
    ENDIF

    IF( MOD( REAL(  iadv_rcf,wp)*dtime, &
      &         atm_phy_nwp_config(jg)%dt_conv) /= 0._wp )  THEN
      WRITE(message_text,'(a,I4,2F10.2)') &
        &'advective and convective timesteps are not- but will be synchronized ', &
        &     1, REAL(  iadv_rcf,wp)*dtime,tcall_phy(1,itconv)
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

 
END SUBROUTINE configure_atm_phy_nwp


END MODULE mo_atm_phy_nwp_config
