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
  USE mo_exception,           ONLY: message, message_text, finish, print_value
  USE mo_master_nml,          ONLY: lrestart
  USE mo_impl_constants,      ONLY: max_char_length, max_dom,itconv,itccov,&
    &                               itrad,itradheat, itsso,itgscp,itsatad,itupdate,&
    &                               itturb, itsfc,  itgwd, iphysproc,iecham, ildf_echam,&
    &                               inwp, iheldsuarez, ildf_dry,  &
    &                               IHS_ATM_TEMP,IHS_ATM_THETA, &
    &                               tracer_only, inh_atmosphere, ishallow_water
  USE mo_parallel_configuration, ONLY: check_parallel_configuration
  USE mo_run_config,          ONLY: lrestore_states, dtime, iforcing, ltransport, &
                                    &ntracer, nlev, io3, inextra_2D, inextra_3D,&
                                    & configure_run, ltestcase
  USE mo_time_config,         ONLY: time_config, configure_time
  USE mo_gridref_config
  USE mo_interpol_config      
  USE mo_grid_configuration   
  USE mo_sleve_config         

  USE mo_dynamics_config,     ONLY: configure_dynamics,&
    &                            iequations, itime_scheme, idiv_method, divavg_cntrwgt,&
    &                            sw_ref_height, ldry_dycore, lcoriolis, lshallow_water, ltwotime
  USE mo_advection_config,    ONLY: advection_config !, configure_advection

  USE mo_nonhydrostatic_config      
  USE mo_ha_dyn_config,     ONLY: ha_dyn_config
  USE mo_diffusion_config,  ONLY: diffusion_config, configure_diffusion

  USE mo_io_config          

  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, tcall_phy, configure_atm_phy_nwp
  USE mo_lnd_nwp_config,     ONLY: nlev_soil, nztlev ,nlev_snow ,nsfc_subs,&
    &                              lseaice,  llake, lmelt , lmelt_var, lmulti_snow
  USE mo_echam_phy_config,   ONLY: echam_phy_config, configure_echam_phy
  USE mo_radiation_config
  USE mo_echam_conv_config,  ONLY: echam_conv_config, configure_echam_convection
  USE mo_gw_hines_config,    ONLY: gw_hines_config
  USE mo_vdiff_config,       ONLY: vdiff_config
  USE mo_nh_testcases,       ONLY: linit_tracer_fv,nh_test_name
  USE mo_hydro_testcases,    ONLY: ctest_name

  IMPLICIT NONE

!  PRIVATE

  PUBLIC :: atm_crosscheck !, atmospheric_configuration

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  SUBROUTINE atm_crosscheck

    INTEGER :: jg
    CHARACTER(len=*), PARAMETER :: routine =  'atm_crosscheck'

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration(lrestore_states)

    !--------------------------------------------------------------------
    ! Horizontal interpolation
    !--------------------------------------------------------------------
    IF (global_cell_type == 6) THEN
    ! ... check i_cori_method
      IF (i_cori_method <1 .OR. i_cori_method>4) THEN
        CALL finish( TRIM(routine),'value of i_cori_method out of range [1,2,3,4]')
      ENDIF
    ENDIF

    !--------------------------------------------------------------------
    ! Tracer transport
    !--------------------------------------------------------------------
    IF((itime_scheme==tracer_only).AND.(.NOT.ltransport)) THEN
      WRITE(message_text,'(A,i2,A)') &
      'itime_scheme set to ', tracer_only, 'but ltransport to .FALSE.'
      CALL finish( TRIM(routine),TRIM(message_text))
    END IF

    IF(ltransport .AND. ntracer <= 0) THEN
      CALL finish( TRIM(routine),'Tracer transport switched on but ntracer <= 0')
      ! [for nwp forcing ntracer setting is treated in setup_transport]
    ENDIF

    SELECT CASE (iforcing)
    CASE(IECHAM,ILDF_ECHAM)
      IF (ntracer < 3) &
      CALL finish(TRIM(routine),'ECHAM physics needs at least 3 tracers')

    CASE(INWP)
      IF (ntracer < 5) &
      CALL finish(TRIM(routine),'NWP physics needs at least 3 tracers')
    END SELECT

    !--------------------------------------------------------------------
    ! Grid and dynamics 
    !--------------------------------------------------------------------
    IF (lplane .AND. global_cell_type==3) CALL finish( TRIM(routine),&
      'Currently only the hexagon model can run on a plane')

    IF (global_cell_type==6.AND.idiv_method==2) THEN
      CALL finish( TRIM(ROUTINE),'idiv_method =2 not valid for the hexagonal model') 
    ENDIF

    IF ((iforcing==IHELDSUAREZ.OR.iforcing==ILDF_DRY).AND.(.NOT.ldry_dycore)) &
    CALL finish( TRIM(ROUTINE),'ldry_dycore should be .TRUE. for the '//&
               'Held-Suarez test and the dry local diabatic forcing test.')

    !--------------------------------------------------------------------
    ! testcases 
    !--------------------------------------------------------------------

    IF(global_cell_type==3) THEN
      linit_tracer_fv  = .TRUE. ! like default
    ELSE
      linit_tracer_fv  = .FALSE.
    ENDIF



    !--------------------------------------------------------------------
    ! Shallow water 
    !--------------------------------------------------------------------
    IF (iequations==ISHALLOW_WATER.AND.ha_dyn_config%lsi_3d) THEN
      CALL message( TRIM(routine), 'lsi_3d = .TRUE. not appicable to shallow water model')
    ENDIF

    IF ((iequations==ISHALLOW_WATER).AND.(nlev/=1)) &
    CALL finish(TRIM(routine),'Multiple vertical level specified for shallow water model')

    !--------------------------------------------------------------------
    ! Hydrostatic atm
    !--------------------------------------------------------------------
    IF (iequations==IHS_ATM_THETA) ha_dyn_config%ltheta_dyn = .TRUE.

    !--------------------------------------------------------------------
    ! Nonhydrostatic atm
    !--------------------------------------------------------------------
    ! reset l_nest_rcf to false if iadv_rcf = 1
    IF (iadv_rcf == 1) l_nest_rcf = .FALSE.

    IF (upstr_beta > 1.0_wp .OR. upstr_beta < 0.0_wp) THEN
      CALL finish(TRIM(routine), 'upstr_beta out of range 0..1')
    ENDIF

    ! for reduced calling frequency of tracer advection / fast physics:
    ! odd values of iadv_rcf are allowed only if nest calls are synchronized 
    ! with advection
    IF ( .NOT. l_nest_rcf .AND. MOD(iadv_rcf,2) /= 0 &
         .AND. iadv_rcf /= 1 .OR. iadv_rcf == 0) THEN
      CALL finish( TRIM(routine), 'Invalid reduced-calling-frequency parameter. '//&
        &'Value must be even or 1 if l_nest_rcf=.FALSE.')
    ENDIF

    !--------------------------------------------------------------------
    ! Atmospheric physics, general
    !--------------------------------------------------------------------
    IF ((iforcing==INWP).AND.(iequations/=INH_ATMOSPHERE)) &
    CALL finish( TRIM(routine), 'NWP physics only implemented in the '//&
               'nonhydrostatic atm model')

    IF ((iforcing==IECHAM).AND.(iequations==INH_ATMOSPHERE)) &
    CALL finish( TRIM(routine), 'ECHAM physics not implemented in the '//&
               'nonhydrostatic atm model')

    !--------------------------------------------------------------------
    ! NWP physics
    !--------------------------------------------------------------------
    IF (iforcing==inwp) THEN
      DO jg =1,n_dom

        IF( (atm_phy_nwp_config(jg)%inwp_convection >0 ) .OR. &
            (atm_phy_nwp_config(jg)%inwp_gscp > 0)       .AND.&
             atm_phy_nwp_config(jg)%inwp_satad == 0)& 
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
        CALL message(TRIM(routine),' WARNING! NWP forcing set but '//&
                    'only turbulence selected!')

        IF ((iequations==INH_ATMOSPHERE).AND.(iforcing==inwp).AND.ldry_dycore) &
        CALL finish(TRIM(routine),'ldry_dycore = .TRUE. not allowed for '//&
                   'the nonhydrostaic atm model with NWP physics.')

        ! check radiation scheme in relation to chosen ozone

        IF (  atm_phy_nwp_config(jg)%inwp_radiation > 0 )  THEN

          SELECT CASE (irad_o3)
          CASE (0,6) ! ok
          CASE default
            CALL finish(TRIM(routine),'irad_o3 currently has to be 0 or 6.')
          END SELECT
        ENDIF

      ENDDO
    END IF

    !--------------------------------------------------------------------
    ! ECHAM physics
    !--------------------------------------------------------------------
    IF ((iforcing==IECHAM).AND.(echam_phy_config%lrad).AND. &
        (irad_o3 > 0) .AND. (io3 > ntracer) ) THEN

      CALL print_value('irad_o3' ,irad_o3)
      CALL print_value('io3    ' ,io3)
      CALL print_value('ntracer' ,ntracer)
      CALL finish(TRIM(routine), 'Not enough tracers for ECHAM physics with RRTM.')
    END IF

    !--------------------------------------------------------------------
    ! checking the meanings of the diffusion settings
    !--------------------------------------------------------------------

 DO jg =1,n_dom

   SELECT CASE( diffusion_config(jg)%hdiff_order)
   CASE(-1)
     CALL message(TRIM(routine),'Horizontal diffusion switched off.')
   CASE(2,4)
     CONTINUE

   CASE(3)
     IF (global_cell_type==3) CALL finish(TRIM(routine), &
     ' hdiff_order = 3 invalid for triangular model.')

   CASE(5)
     IF (global_cell_type==6) CALL finish(TRIM(routine), &
     ' hdiff_order = 5 invalid for hexagonal model.')

   CASE(24,42)
     IF (.NOT.( iequations==IHS_ATM_TEMP)) CALL finish(TRIM(routine), &
     ' hdiff_order = 24 or 42 only implemented for the hydrostatic atm model')

   CASE DEFAULT
     CALL finish(TRIM(routine),                     &
       & 'Error: Invalid choice of  hdiff_order. '// &                
       & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
   END SELECT

   IF (  diffusion_config(jg)%hdiff_efdt_ratio<=0._wp) THEN
     CALL message(TRIM(routine),'No horizontal background diffusion is used')
  ENDIF

  IF ( lshallow_water )  diffusion_config(jg)%lhdiff_temp=.FALSE.

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

!  SUBROUTINE atmospheric_configuration
!  INTEGER :: jg
!  CHARACTER(len=*), PARAMETER :: routine =  'atm_setup'
!  END SUBROUTINE atmospheric_configuration

END MODULE mo_atm_nml_crosscheck

