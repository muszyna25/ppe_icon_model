!>
!!  Namelist for NWP physics
!!
!!  these Subroutines are called by control model and construct the
!!  physics composition
!!
!! @author <Kristina Froehlich, DWD>
!!
!!
!! @par Revision History
!! First implementation by Kristina Froehlich, DWD (2010-06-20>)
!!  namelist varaibles for calling frequency of fast physics schemes
!!   have been removed to ensure the high frequent calls
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
MODULE mo_nwp_phy_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_dom !,MAX_CHAR_LENGTH
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run

  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_nwp_phy_namelist 

   !
   ! user defined calling intervals
   !
  REAL(wp) :: dt_conv(max_dom)   !> field element for convection
  REAL(wp) :: dt_ccov(max_dom)   !! field element for subscale cloud cover
  REAL(wp) :: dt_rad(max_dom)    !! "-"                     radiation
  REAL(wp) :: dt_sso(max_dom)    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gwd(max_dom)    !! "-"  for subscale gravity waves

  ! switches defining physics packages
  INTEGER  :: inwp_convection(max_dom)    !! convection
  INTEGER  :: inwp_cldcover(max_dom)      !! cloud cover
  INTEGER  :: inwp_radiation(max_dom)     !! radiation
  INTEGER  :: inwp_sso(max_dom)           !! sso
  INTEGER  :: inwp_gwd(max_dom)           !! non-orographic gravity wave drag
  INTEGER  :: inwp_gscp(max_dom)          !! microphysics
  INTEGER  :: inwp_satad(max_dom)         !! saturation adjustment
  INTEGER  :: inwp_turb(max_dom)          !! turbulence
  INTEGER  :: inwp_surface(max_dom)       !! surface including soil, ocean, ice,lake

  INTEGER  :: itype_z0           !! type of roughness length data
  REAL(wp) :: qi0, qc0           !! variables for hydci_pp
  REAL(wp) :: ustart_raylfric    !! velocity at which extra Rayleigh friction starts
  REAL(wp) :: efdt_min_raylfric  !! e-folding time corresponding to maximum relaxation coefficient
  LOGICAL  :: latm_above_top(max_dom) !! use extra layer above model top for radiation (reduced grid only)
  ! parameter for cloud microphysics
  real(wp) :: mu_rain            !! shape parameter in gamma distribution for rain
  real(wp) :: mu_snow            !! ...for snow
  


  NAMELIST /nwp_phy_nml/ inwp_convection, inwp_cldcover,           &
    &                    inwp_radiation, inwp_sso, inwp_gwd,       &
    &                    inwp_gscp, inwp_satad,                    &
    &                    inwp_turb, inwp_surface,                  &
    &                    dt_conv, dt_ccov, dt_rad, dt_sso, dt_gwd, &
    &                    qi0, qc0,                                 &
    &                    ustart_raylfric, efdt_min_raylfric,       &
    &                    latm_above_top, itype_z0, mu_rain, mu_snow



CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read physics Namelist
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
 SUBROUTINE read_inwp_nml

  INTEGER :: i_status

  !0!CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
  !0!  &  routine = 'mo_atm_phy_nwp_nml:read_inwp_nml'

    CALL position_nml ('nwp_phy_nml', status=i_status)
    !
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_nml)
    END SELECT
  !  write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_phy_nml)

 END SUBROUTINE read_inwp_nml



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_nwp_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, jg
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = 'mo_atm_phy_nwp_nml:read_nwp_phy_namelist'

    !-----------------------
    ! 0a. dummy settings; will be replaced with defaults after reading the namelist
    !     wherever nothing is specified explicitly   
    !-----------------------
    inwp_gscp(:)       = -1
    inwp_satad(:)      = -1
    inwp_convection(:) = -1
    inwp_radiation(:)  = -1
    inwp_sso(:)        = -1
    inwp_gwd(:)        = -1
    inwp_cldcover(:)   = -1
    inwp_turb(:)       = -1
    inwp_surface(:)    = -1

    dt_conv (:) = -999._wp
    dt_ccov (:) = -999._wp
    dt_rad  (:) = -999._wp
    dt_sso  (:) = -999._wp
    dt_gwd  (:) = -999._wp

   



    !------------------------------------------------------------------
    ! 0b. Real defaults for some other variables
    !------------------------------------------------------------------

    itype_z0 = 1  !  1 = include orographic contribution to roughness length
    qi0      = 0.0_wp 
    qc0      = 0.0_wp 

    ! shape paramter for gamma distribution for rain and snow
    mu_rain = 0.0_wp
    mu_snow = 0.0_wp

    ustart_raylfric    = 160._wp
    efdt_min_raylfric  = 10800._wp

    latm_above_top(:)  = .FALSE.  ! no extra layer above model top for radiation computation

    !------------------------------------------------------------------
    ! 1. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('nwp_phy_nml')
      READ(funit,NML=nwp_phy_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_phy_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_nml)
    END SELECT
    CALL close_nml


    !-----------------------
    ! 3. apply default settings where nothing is specified explicitly (except for restart)
    !-----------------------

    IF (.NOT. is_restart_run()) THEN

      ! 3a. Set default values for global domain where nothing at all has been specified

      ! Physics packages
      IF (inwp_gscp(1)       < 0) inwp_gscp(1)       = 1  !> 1 = hydci (COSMO-EU microphysics)
      IF (inwp_satad(1)      < 0) inwp_satad(1)      = 1  !> 1 = saturation adjustment on
      IF (inwp_convection(1) < 0) inwp_convection(1) = 1  !> 1 = Tiedtke/Bechthold convection
      IF (inwp_radiation(1)  < 0) inwp_radiation(1)  = 1  !> 1 = RRTM radiation
      IF (inwp_sso(1)        < 0) inwp_sso(1)        = 1  !> 1 = Lott and Miller scheme (COSMO)
      IF (inwp_gwd(1)        < 0) inwp_gwd(1)        = 1  !> 1 = Orr-Ern-Bechthold scheme (IFS)
      IF (inwp_cldcover(1)   < 0) inwp_cldcover(1)   = 3  !> 3 = clouds from COSMO SGS cloud scheme
      IF (inwp_turb(1)       < 0) inwp_turb(1)       = 1  !> 1 = turbdiff (COSMO diffusion oand transfer)
      IF (inwp_surface(1)    < 0) inwp_surface(1)    = 1  !> 1 = TERRA

      ! Time steps
      IF (dt_conv (1) < 0._wp) dt_conv (1) = 600._wp    !seconds
      IF (dt_ccov (1) < 0._wp) dt_ccov (1) = dt_conv(1) !cloud cover is synchronized with convection by default
      IF (dt_sso  (1) < 0._wp) dt_sso  (1) = 1200._wp   !seconds
      IF (dt_gwd  (1) < 0._wp) dt_gwd  (1) = 1200._wp   !seconds
      IF (dt_rad  (1) < 0._wp) dt_rad  (1) = 1800._wp   !seconds

      ! 3b. Copy values of parent domain (in case of linear nesting) to nested domains where nothing has been specified

      DO jg = 2, max_dom

        ! Physics packages
        IF (inwp_gscp(jg)       < 0) inwp_gscp(jg)       = inwp_gscp(jg-1)
        IF (inwp_satad(jg)      < 0) inwp_satad(jg)      = inwp_satad(jg-1)
        IF (inwp_convection(jg) < 0) inwp_convection(jg) = inwp_convection(jg-1)
        IF (inwp_radiation(jg)  < 0) inwp_radiation(jg)  = inwp_radiation(jg-1)
        IF (inwp_sso(jg)        < 0) inwp_sso(jg)        = inwp_sso(jg-1)
        IF (inwp_gwd(jg)        < 0) inwp_gwd(jg)        = inwp_gwd(jg-1)
        IF (inwp_cldcover(jg)   < 0) inwp_cldcover(jg)   = inwp_cldcover(jg-1)
        IF (inwp_turb(jg)       < 0) inwp_turb(jg)       = inwp_turb(jg-1)
        IF (inwp_surface(jg)    < 0) inwp_surface(jg)    = inwp_surface(jg-1)

        ! Time steps
        IF (dt_conv (jg) < 0._wp) dt_conv (jg) = dt_conv (jg-1) 
        IF (dt_ccov (jg) < 0._wp) dt_ccov (jg) = dt_ccov (jg-1)
        IF (dt_sso  (jg) < 0._wp) dt_sso  (jg) = dt_sso  (jg-1)
        IF (dt_gwd  (jg) < 0._wp) dt_gwd  (jg) = dt_gwd  (jg-1)
        IF (dt_rad  (jg) < 0._wp) dt_rad  (jg) = dt_rad  (jg-1)

      ENDDO

    ENDIF

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    
    ! check for valid parameters in namelists

    DO jg = 1, max_dom

      IF (inwp_gscp(jg) /= 0 .AND. inwp_gscp(jg) /= 1 .AND. inwp_gscp(jg) /= 2 .AND. &
        & inwp_gscp(jg) /= 3 .AND. inwp_gscp(jg) /= 4 .AND. inwp_gscp(jg) /= 9 .AND. &
        & inwp_gscp(jg) /= 10 ) THEN
   
        CALL finish( TRIM(routine), 'Incorrect setting for inwp_gscp. Must be 0,1,2,3,4,9 or 10.')
      END IF

    ENDDO


    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg=1,max_dom
      atm_phy_nwp_config(jg)%inwp_convection = inwp_convection(jg)
      atm_phy_nwp_config(jg)%inwp_cldcover   = inwp_cldcover(jg)
      atm_phy_nwp_config(jg)%inwp_radiation  = inwp_radiation(jg)
      atm_phy_nwp_config(jg)%inwp_sso        = inwp_sso(jg)
      atm_phy_nwp_config(jg)%inwp_gwd        = inwp_gwd(jg) 
      atm_phy_nwp_config(jg)%inwp_gscp       = inwp_gscp(jg)
      atm_phy_nwp_config(jg)%inwp_satad      = inwp_satad(jg)
      atm_phy_nwp_config(jg)%inwp_turb       = inwp_turb(jg)
      atm_phy_nwp_config(jg)%inwp_surface    = inwp_surface(jg)
      atm_phy_nwp_config(jg)%itype_z0        = itype_z0

      atm_phy_nwp_config(jg)%dt_conv         = dt_conv (jg) 
      atm_phy_nwp_config(jg)%dt_ccov         = dt_ccov (jg)
      atm_phy_nwp_config(jg)%dt_rad          = dt_rad  (jg)
      atm_phy_nwp_config(jg)%dt_sso          = dt_sso  (jg)
      atm_phy_nwp_config(jg)%dt_gwd          = dt_gwd  (jg)
      atm_phy_nwp_config(jg)%qi0             = qi0 
      atm_phy_nwp_config(jg)%qc0             = qc0 
      atm_phy_nwp_config(jg)%ustart_raylfric = ustart_raylfric 
      atm_phy_nwp_config(jg)%efdt_min_raylfric = efdt_min_raylfric
      atm_phy_nwp_config(jg)%latm_above_top  = latm_above_top(jg)
      atm_phy_nwp_config(jg)%mu_rain         = mu_rain
      atm_phy_nwp_config(jg)%mu_snow         = mu_snow

    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nwp_phy_nml)                    
      CALL store_and_close_namelist(funit, 'nwp_phy_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_phy_nml)

  END SUBROUTINE read_nwp_phy_namelist


END MODULE mo_nwp_phy_nml

