!>
!! Namelist specific to the hydro atm dynamical core.
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_ha_dyn_nml

  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: AB2, TRACER_ONLY, TWO_TL_SI, LEAPFROG_EXPL, &
                                      LEAPFROG_SI, RK4, SSPRK54
  USE mo_mpi,                   ONLY: my_process_is_stdio
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_exception,             ONLY: finish
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,        ONLY: use_restart_namelists
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &
                                      open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,          ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_ha_dyn_namelist

  !---------------------
  ! namelist variables
  !---------------------

  INTEGER  :: itime_scheme   ! variable used to select the time stepping scheme
                             ! = 1, preparational computations for advection test
                             ! = 12, semi implicit 2 time level scheme
                             ! = 13, explicit leapfrog
                             ! = 14, leapfrog with semi implicit correction
                             ! = 15, 4-stage Runge-Kutta method
                             ! = 16, SSPRK(5,4) (Runge-Kutta) method

  INTEGER  :: ileapfrog_startup  ! choice of time stepping scheme for 
                                 ! the first step in
                                 ! a leapfrog time stepping scheme
                                 ! 1 = Euler forward
                                 ! 2 = several sub-steps

  REAL(wp) :: asselin_coeff      ! parameter used in Asselin filter

  INTEGER  :: si_expl_scheme     ! scheme for the explicit part of the
                                 ! 2-time-level semi-implicit time integration.
                                 ! See mo_impl_constants for the options.
  REAL(wp) :: si_2tls

  REAL(wp) :: si_rtol     ! relative tolerance

  REAL(wp) :: si_coeff    !  = 0 : explicit scheme(for *d*,*t*,*alps*).
                          !  = 1 : semi implicit scheme.
                          !  in (0,1): a weighted scheme

  REAL(wp) :: si_offctr   ! weighting parameter used in calculating the
                          ! second temporal derivatives in the semi-implicit
                          ! correction scheme. The value read from namelist are
                          ! assumed to be the offcentering (i.e. between 0 and 1).



  REAL(wp) :: si_cmin     ! min. phase speed of the decomposed modes to be
                          ! solved by the semi-implicit correction scheme
  LOGICAL  :: lsi_3d      ! if .true., solve the 3D equation

  LOGICAL :: ldry_dycore  ! if .TRUE., ignore the effact of water vapor,
                          ! cloud liquid and cloud ice on virtual temperature.

  LOGICAL  :: lref_temp   ! if .TRUE., involve the reference temperature profile
                          ! in the calculation of pressure gradient force.

  NAMELIST/ha_dyn_nml/ itime_scheme,                     &
                       ileapfrog_startup, asselin_coeff, &
                       si_expl_scheme, si_2tls, si_rtol, &
                       si_coeff, si_offctr, si_cmin,     &
                       lsi_3d, ldry_dycore, lref_temp

CONTAINS
  !>
  !!
  SUBROUTINE read_ha_dyn_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    CHARACTER(len=*),PARAMETER :: routine = 'mo_ha_dyn_nml:read_ha_dyn_namelist'

    !------------------------------------------------------------
    ! Default values
    !------------------------------------------------------------
    itime_scheme      = LEAPFROG_SI
    ileapfrog_startup = 1
    asselin_coeff     = 0.1_wp

    si_expl_scheme    = AB2
    si_2tls           = 0.6_wp
    si_rtol           = 1.e-3_wp

    si_coeff          = 1.0_wp
    si_offctr         = 0.7_wp
    si_cmin           = 30._wp
    lsi_3d            = .FALSE.

    ldry_dycore       = .FALSE.
    lref_temp         = .FALSE.

    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('ha_dyn_nml')
      READ(funit,NML=ha_dyn_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ha_dyn_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, ha_dyn_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ha_dyn_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ha_dyn_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Sanity Check
    !-----------------------------------------------------
    SELECT CASE (itime_scheme)
    CASE (TRACER_ONLY,TWO_TL_SI,LEAPFROG_EXPL,LEAPFROG_SI,RK4,SSPRK54) !OK
    CASE DEFAULT
      CALL finish( TRIM(routine),'wrong value for ha_dyn_nml:itime_scheme. '//&
                 & 'See mo_impl_constants.f90 for possible options.' )
    END SELECT

    IF (asselin_coeff<0._wp) CALL finish( TRIM(routine), &
      'wrong (negative) coefficient of Asselin filter')

    IF (si_offctr>1._wp.OR.si_offctr<0._wp) CALL finish( TRIM(routine), &
      'Invalid value for parameter si_offctr. Valid range is [0,1].')

    IF (si_2tls<0.5_wp.OR.si_2tls>1._wp) CALL finish( TRIM(routine), &
      'Improper value for parameter si_2tls. Should be in [0.5,1]')


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=ha_dyn_nml)
      CALL store_and_close_namelist(funit, 'ha_dyn_nml')
    ENDIF
    
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=ha_dyn_nml)

    !-----------------------------------------------------
    ! Fill configuration state
    !-----------------------------------------------------
    ha_dyn_config% itime_scheme      = itime_scheme
    ha_dyn_config% ileapfrog_startup = ileapfrog_startup
    ha_dyn_config% asselin_coeff     = asselin_coeff
    ha_dyn_config% si_expl_scheme    = si_expl_scheme
    ha_dyn_config% si_2tls           = si_2tls
    ha_dyn_config% si_rtol           = si_rtol
    ha_dyn_config% si_coeff          = si_coeff
    ha_dyn_config% si_offctr         = si_offctr
    ha_dyn_config% si_cmin           = si_cmin
    ha_dyn_config% lsi_3d            = lsi_3d
    ha_dyn_config% lref_temp         = lref_temp
    ha_dyn_config% ldry_dycore       = ldry_dycore

  END SUBROUTINE read_ha_dyn_namelist

END MODULE mo_ha_dyn_nml
