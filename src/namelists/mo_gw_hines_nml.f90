!>
!! Namelist for the configuration of the vertical diffusion
!!
!!
!! @par Revision History
!! Revision history in mo_echam_vdiff_params.f90 (r4300)
!! Modification by Constantin Junk, MPI-M (2011-05-05)
!! - moved echam_vdiff namelist variables and subroutine setup_vdiff
!!   from mo_echam_vdiff_params to namelists/mo_echam_vdiff_nml
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gw_hines_nml

  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_exception,           ONLY: finish

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom

  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist, &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_gw_hines_config,     ONLY: gw_hines_config
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_gw_hines_namelist

  !-----------------------------------!
  ! gw_hines_nml namelist variables   !
  !-----------------------------------!

  LOGICAL  :: lheatcal      !< true : compute momentum flux dep., heating and diffusion coefficient
                            !< false: compute only momentum flux deposition

  INTEGER  :: emiss_lev     !< number of levels above the ground at which gw are emitted
  REAL(wp) :: rmscon        !< [m/s] root mean square gravity wave wind at emission level
  REAL(wp) :: kstar         !< [1/m] typical gravity wave horizontal wavenumber
  REAL(wp) :: m_min         !< [1/m] minimum bound in  vertical wavenumber

!!$  LOGICAL  :: lfront        !< true: compute gw sources emerging from fronts and background
!!$                            !< (Charron and Manzini, 2002)
!!$  REAL(wp) :: rms_front     !< [m/s] rms frontal gw wind at source level
!!$  REAL(wp) :: front_thres   !< [(K/m)^2/hr] minimum value of the frontogenesis function,
!!$                            !< for which gravity waves are emitted from fronts
!!$
!!$  LOGICAL  :: lozpr         !< true: for background enhancement associated with precipitation
!!$                            !< (Manzini et al., 1997)
!!$  REAL(wp) :: pcrit         !< [mm/d] critical precipitation value, above which 
!!$                            !< gravity wave rms wind enhancement is applied
!!$  REAL(wp) :: pcons         !< [] adimensional factor for background enhancement 
!!$                            !< associated with precipitation

  LOGICAL  :: lrmscon_lat   !< true:  use latitude dependent rmscon
                            !< - |latitude| >= lat_rmscon:
                            !<      use rmscon
                            !< - |latitude| <= lat_rmscon_eq:
                            !<      use rmscon_eq
                            !< - lat_rmscon_eq < |latitude| < lat_rmscon: 
                            !<      use linear interpolation between rmscon_eq and rmscon
                            !< false: use rmscon for all latitudes
                            !< attention: may be overwritten if lfront or lozpr is true
  REAL(wp) :: lat_rmscon_eq !< [degN] rmscon_eq is used equatorward of this latitude
  REAL(wp) :: lat_rmscon    !< [degN] rmscon is used poleward of this latitude
  REAL(wp) :: rmscon_eq     !< [m/s]  rms constant used equatorward of lat_rmscon_eq


NAMELIST /gw_hines_nml/ lheatcal, rmscon, emiss_lev, kstar, m_min,         &
!!$  &                     lfront, rms_front, front_thres,                    &
!!$  &                     lozpr, pcrit, pcons,                               &
  &                     lrmscon_lat, lat_rmscon_eq, lat_rmscon, rmscon_eq

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for Hines gw parameterization. 
  !!
  !! This subroutine 
  !! - reads the Namelist for Hines gw parameterization
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
  SUBROUTINE read_gw_hines_namelist( filename )

    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg 
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER :: routine = 'mo_gw_hines_nml:read_gw_hines_namelist'

    !------------------------------------------------------------------
    ! 1. Set default values
    !------------------------------------------------------------------
    lheatcal      = .FALSE.

    emiss_lev     = 10         ! is correct for L31 and L47
    rmscon        = 0.87_wp    ! default value used in ECHAM5
    kstar         = 5.0e-5_wp  ! = 2*pi/(126000 m)
    m_min         = 0.0_wp

    lrmscon_lat   = .FALSE.    ! ICON default: no latitude dependence
    lat_rmscon_eq =  5.0_wp    ! as in ECHAM6 T63
    lat_rmscon    = 10.0_wp    ! as in ECHAM6 T63
    rmscon_eq     = 1.2_wp     ! as in ECHAM6 T63

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('gw_hines_nml')
      READ(funit,NML=gw_hines_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------
    ! 3. Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gw_hines_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, gw_hines_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, gw_hines_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, gw_hines_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! 4. Sanity Check
    !------------------------------------------------------------------
    IF ( emiss_lev < 0 )   CALL finish(TRIM(routine),'emiss_lev < 0 is not allowed')
    IF ( rmscon < 0.0_wp ) CALL finish(TRIM(routine),'rmscon < 0. is not allowed')
    IF ( kstar  < 0.0_wp ) CALL finish(TRIM(routine),'kstar  < 0. is not allowed')
    IF ( m_min  < 0.0_wp ) CALL finish(TRIM(routine),'m_min  < 0. is not allowed')

    IF ( lrmscon_lat ) THEN
      !
      IF ( lat_rmscon_eq <  0.0_wp  ) &
        & CALL finish(TRIM(routine),'lat_rmscon_eq <  0. degN is not valid')
      IF ( lat_rmscon_eq > 90.0_wp ) &
        & CALL finish(TRIM(routine),'lat_rmscon_eq > 90. degN is not valid')
      IF ( lat_rmscon    <  0.0_wp ) &
        & CALL finish(TRIM(routine),'lat_rmscon    <  0. degN is not valid')
      IF ( lat_rmscon    > 90.0_wp ) &
        & CALL finish(TRIM(routine),'lat_rmscon    > 90. degN is not valid')
      IF ( lat_rmscon_eq > lat_rmscon ) &
        & CALL finish(TRIM(routine),'lat_rmscon_eq > lat_rmscon is not allowed')
      IF ( rmscon_eq < 0.0_wp ) &
        & CALL finish(TRIM(routine),'rmscon_eq < 0. is not allowed')
      !
    END IF

    !------------------------------------------------------------------
    ! 5. Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=gw_hines_nml)                    
      CALL store_and_close_namelist(funit, 'gw_hines_nml') 
    ENDIF

    !------------------------------------------------------------------
    ! 6. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=gw_hines_nml)

    !------------------------------------------------------------------
    ! 7. Fill the configuration state
    !------------------------------------------------------------------

    DO jg = 1,max_dom
      
      gw_hines_config(jg) %lheatcal      = lheatcal
      gw_hines_config(jg) %emiss_lev     = emiss_lev
      gw_hines_config(jg) %rmscon        = rmscon
      gw_hines_config(jg) %kstar         = kstar
      gw_hines_config(jg) %m_min         = m_min

      gw_hines_config(jg) %lrmscon_lat   = lrmscon_lat
      gw_hines_config(jg) %lat_rmscon_eq = lat_rmscon_eq
      gw_hines_config(jg) %lat_rmscon    = lat_rmscon
      gw_hines_config(jg) %rmscon_eq     = rmscon_eq

    ENDDO

  END SUBROUTINE read_gw_hines_namelist

END MODULE mo_gw_hines_nml
