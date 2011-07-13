!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Kristina Froehlich, MPI-M (2011-07-12)
!! @author Hui Wan, MPI-M (2011-07-12)
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_atm_nml_crosscheck

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_impl_constants,      ONLY: max_char_length, max_dom,itconv,itccov,&
    &                               itrad,itradheat, itsso,itgscp,itsatad,itupdate,&
    &                               itturb, itsfc,  itgwd, iphysproc,iecham, ildf_echam,&
    &                               inwp
  USE mo_time_config,         ONLY: time_config
  USE mo_run_config
  USE mo_gridref_config,      ONLY: gridref_config
  USE mo_interpol_config       ! all, ONLY: interpol_config
  USE mo_grid_configuration    !,     ONLY: global_cell_type
  USE mo_sleve_config          ! all, ONLY: sleve_config

  USE mo_dynamics_config,     ONLY: dynamics_config
  USE mo_advection_config,    ONLY: advection_config

  USE mo_nh_dyn_config       !For now: all, later  ONLY: nh_dyn_config
  USE mo_ha_dyn_config,     ONLY: ha_dyn_config
  USE mo_diffusion_config,  ONLY: diffusion_config

  USE mo_io_config           !all,         ONLY: io_config

  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, tcall_phy
  USE mo_lnd_nwp_config,     ONLY: nwp_lnd_config

  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_radiation_config,   ONLY: radiation_config  
  USE mo_echam_conv_config,  ONLY: echam_conv_config
  USE mo_gw_hines_config,    ONLY: gw_hines_config
  USE mo_vdiff_config,       ONLY: vdiff_config


  IMPLICIT NONE

!  PRIVATE

!  PUBLIC  

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

SUBROUTINE atm_crosscheck

  INTEGER :: jg

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                              'atm_crosscheck'


  IF (lplane .AND. global_cell_type==3) THEN
    CALL finish( TRIM(routine),&
      'currently, only the hexagon version runs on a plane')
  ENDIF


  !--------------------------------------------------------------------
  ! checking of interpolation  parameters
  !--------------------------------------------------------------------

  IF (global_cell_type == 6) THEN
    ! ... check i_cori_method
    IF (i_cori_method <1 .OR. i_cori_method>4) THEN
      CALL finish( TRIM(routine),'value of i_cori_method out of range [1,2,3,4]')
    ENDIF
  ENDIF
  !--------------------------------------------------------------------
  ! checking of nh-dynamics parameters
  !--------------------------------------------------------------------

  ! reset l_nest_rcf to false if iadv_rcf = 1
  IF (iadv_rcf == 1) l_nest_rcf = .FALSE.

  IF (upstr_beta > 1.0_wp .OR. upstr_beta < 0.0_wp) THEN
    CALL finish(TRIM(routine), 'upstr_beta out of range 0..1')
  ENDIF

  ! for reduced calling frequency of tracer advection / fast physics:
  ! odd values of iadv_rcf are allowed only if nest calls are synchronized with advection
  IF ( .NOT. l_nest_rcf .AND. MOD(iadv_rcf,2) /= 0 .AND. iadv_rcf /= 1 .OR. iadv_rcf == 0) THEN
      CALL finish( TRIM(routine), 'Invalid reduced-calling-frequency parameter& '//&
        &'Value must be even or 1 if l_nest_rcf=.FALSE.')
  ENDIF


    !--------------------------------------------------------------------
    ! checking the meanings of the nwp physics namelist
    !--------------------------------------------------------------------

  DO jg =1,max_dom

    IF( (atm_phy_nwp_config(jg)%inwp_convection >0 ) .OR. (atm_phy_nwp_config(jg)%inwp_gscp > 0) &
      .AND. atm_phy_nwp_config(jg)%inwp_satad == 0)& 
     & CALL finish( TRIM(routine),'satad has to be switched on')


     IF( MOD( REAL(  iadv_rcf,wp)*dtime, &
       &         atm_phy_nwp_config(jg)%dt_conv) /= 0._wp )  THEN
       WRITE(message_text,'(a,I4,2F10.2)') &
      &'advective and convective timesteps are not- but will be synchronized ', &
      &     1, REAL(  iadv_rcf,wp)*dtime,tcall_phy(1,itconv)
      CALL message(TRIM(routine), TRIM(message_text))
     ENDIF

   IF( (atm_phy_nwp_config(jg)%inwp_gscp==0) .AND. &
     & (atm_phy_nwp_config(jg)%inwp_convection==0) .AND.&
     & (atm_phy_nwp_config(jg)%inwp_radiation==0) .AND.&
     & (atm_phy_nwp_config(jg)%inwp_sso==0)  .AND. &
     & (atm_phy_nwp_config(jg)%inwp_surface == 0) .AND.&
     & (atm_phy_nwp_config(jg)%inwp_turb> 0) )   &
     CALL message(TRIM(routine),' WARNING! NWP forcing set but only turbulence selected!')


  IF (iforcing == iecham     .OR. &
&     iforcing == ildf_echam .OR. &
&     iforcing == inwp          ) &
&     dynamics_config(jg)%ldry_dycore     = .FALSE.

! check radiation scheme in relation to chosen ozone


    IF (  atm_phy_nwp_config(jg)%inwp_radiation > 0 )  THEN

      SELECT CASE (radiation_config(jg)%irad_o3)
      CASE (0,6)
        ! ok
      CASE default
        CALL finish(TRIM(routine),'irad_o3 currently has to be 0 or 6.')
      END SELECT
    ENDIF

  ENDDO

    !--------------------------------------------------------------------
    ! checking the meanings of the io settings
    !--------------------------------------------------------------------

  SELECT CASE(iforcing)
  CASE ( inwp )
    ! Do nothing. Keep the initial values, if not specified in namelist.
    ! consider special idealized testcase with turbulence only
    IF( .NOT. ltransport  )   THEN
      lwrite_precip    = .FALSE.
      lwrite_cloud     = .FALSE.
      lwrite_radiation = .FALSE.
      lwrite_tke       = .TRUE.
      lwrite_surface   = .FALSE.
      CALL message('io_nml_setup',' ATTENTION! Only TKE output for TURBULENCE ONLY test')
    ENDIF
  CASE ( iecham,ildf_echam )
    ! Do nothing. Keep the initial values, if not specified in namelist.
  CASE DEFAULT
    ! Do nothing. Keep the initial values, if not specified in namelist.
  END SELECT
  
  IF (( inextra_2D > 0) .OR. (inextra_3D > 0) ) THEN 
    lwrite_extra = .TRUE.
    WRITE(message_text,'(a,2I4,a,L4)') &
      &'inextra is',inextra_2d,inextra_3d ,' lwrite_extra has been set', lwrite_extra
    CALL message('io_namelist', TRIM(message_text))
  ENDIF
  
  IF (inextra_2D == 0 .AND. inextra_3D == 0 .AND. lwrite_extra) &
    CALL finish('io_namelist','need to specify extra fields for extra output')




END  SUBROUTINE atm_crosscheck




END MODULE mo_atm_nml_crosscheck

