!>
!!        Contains the variables to set up the coupling.
!!
!!        
!! @par Revision History
!!   Created by Rene Redler (2011-03-22)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

MODULE mo_coupling_nml

  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------

  USE mo_impl_constants,  ONLY: max_char_length
  USE mo_exception,       ONLY: finish
  USE mo_io_units,        ONLY: nnml
  USE mo_namelist,        ONLY: open_nml, close_nml, position_nml, POSITIONED
  USE mo_nml_annotate,    ONLY: temp_defaults, temp_settings
  USE mo_mpi,             ONLY: my_process_is_stdio

#ifdef YAC_coupling
  USE mo_coupling_config, ONLY: config_coupled_mode

#else
  USE mo_coupling_config, ONLY: t_cpl_field_nml, config_cpl_fields, &
    & number_of_coupled_variables, config_debug_coupler_level
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_coupling_namelist

  ! ------------------------------------------------------------------------
  ! Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------

  ! initialization
  ! --------------

  LOGICAL            :: l_time_average
  LOGICAL            :: l_time_accumulation
  LOGICAL            :: l_diagnostic
  LOGICAL            :: l_activated
  INTEGER            :: dt_coupling
  INTEGER            :: dt_model
  INTEGER            :: lag
  CHARACTER(len=132) :: name

CONTAINS

  !>
  !!  Initialization of variables that contain general information.
  !!
  !!               Initialization of variables that contain general information
  !!               about the coupled model run. The configuration is read from
  !!               namelist 'icon_cpl'.
  !!
  !! @par Revision History
  !!
  SUBROUTINE read_coupling_namelist (namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename

    !
    ! Local variables
    !

#ifdef YAC_coupling
    LOGICAL            :: coupled_mode
#else
    TYPE(t_cpl_field_nml), POINTER :: new_cpl_fields(:)
#endif

    INTEGER :: debug_coupler_level
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: new_dim

    INTEGER :: nbr_cpl_fields
    INTEGER :: istat
    LOGICAL :: first
    LOGICAL :: l_redirect_stdout
    INTEGER :: iunit

    CHARACTER(len=max_char_length), PARAMETER :: &
         &   routine = 'mo_coupling_nml:read_coupling_namelist'

#ifdef YAC_coupling
    NAMELIST /coupling_mode_nml/ coupled_mode
#else
    NAMELIST /coupling_mode_nml/ debug_coupler_level
    NAMELIST /coupling_nml/ name,              &
                          dt_coupling,         &
                          dt_model,            &
                          lag,                 &
                          l_time_average,      &
                          l_time_accumulation, &
                          l_diagnostic,        &
                          l_activated
#endif

    !--------------------------------------------------------------------
    ! 1. Set default values
    !--------------------------------------------------------------------

    name                = " "

    dt_coupling         = 0
    dt_model            = 0

    l_time_average      = .FALSE.
    l_time_accumulation = .FALSE.
    l_redirect_stdout   = .FALSE.

    l_diagnostic        = .FALSE.
    l_activated         = .FALSE.

    debug_coupler_level = 0

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml (TRIM(namelist_filename))
    
    CALL position_nml('coupling_mode_nml',STATUS=istat)
    IF (istat==POSITIONED) THEN
      READ (nnml, coupling_mode_nml)
    ENDIF

#ifdef YAC_coupling

    config_coupled_mode = coupled_mode

#else

    !--------------------------------------------------------------------
    ! 3. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !--------------------------------------------------------------------

!rr    IF (is_restart_run()) THEN
!rr      funit = open_and_restore_namelist('coupling_nml')
!rr      READ(funit,NML=coupling_nml)
!rr      CALL close_tmpfile(funit)
!rr    END IF


    first = .TRUE.
    i     = 0

    !--------------------------------------------------------------------
    ! 3.a loop over occurences of namelist group /coupling_nml/
    !--------------------------------------------------------------------

    DO

       CALL position_nml('coupling_nml', lrewind=first, status=istat)

       dt_coupling = 0
       dt_model    = 0
       lag         = 0

       l_time_average      = .FALSE.
       l_time_accumulation = .FALSE.

       l_diagnostic        = .FALSE.
       l_activated         = .FALSE.

       first = .FALSE.

       !-----------------------------------------------------------------
       ! 3.b if namelist group is present ...
       !-----------------------------------------------------------------

       IF (my_process_is_stdio()) THEN
         iunit = temp_defaults()
         WRITE(iunit, coupling_nml)  ! write defaults to temporary text file
       END IF

       SELECT CASE (istat)

       CASE (POSITIONED)

         READ (nnml, coupling_nml, iostat=istat)                          ! overwrite default settings
         IF (my_process_is_stdio()) THEN
           iunit = temp_settings()
           WRITE(iunit, coupling_nml)  ! write settings to temporary text file
         END IF

          i = i + 1

          !--------------------------------------------------------------
          ! 3.c allocate memory
          !--------------------------------------------------------------

          IF ( .NOT. ASSOCIATED(config_cpl_fields) ) THEN

            nbr_cpl_fields = 8

            ALLOCATE(config_cpl_fields(nbr_cpl_fields), STAT = ierr )

            IF ( ierr > 0 ) &
               CALL finish ( TRIM(routine), ' Error allocating cpl_fields ' )

          ENDIF

          !--------------------------------------------------------------
          ! 3.d increase memory if needed
          !--------------------------------------------------------------

          IF ( i > nbr_cpl_fields ) THEN

             ! we need to allocated more memory

             new_dim = nbr_cpl_fields + 8

             ALLOCATE ( new_cpl_fields(new_dim), STAT = ierr )
             IF ( ierr > 0 ) &
                  CALL finish ( TRIM(routine), ' Error allocating cpl_fields ' )

             new_cpl_fields (1:nbr_cpl_fields) = config_cpl_fields (1:nbr_cpl_fields)

             DEALLOCATE ( config_cpl_fields, STAT = ierr )
             IF ( ierr > 0 ) &
                  CALL finish ( TRIM(routine), ' Error deallocating cpl_fields ' )

             config_cpl_fields => new_cpl_fields

             ! update size

             nbr_cpl_fields = new_dim

          ENDIF

          !--------------------------------------------------------------
          ! 4. Consistency check
          !--------------------------------------------------------------

          if ( dt_coupling < dt_model ) &
               CALL finish (TRIM(routine), &
               'Coupling interval must be larger orequal to model time step' )

          !--------------------------------------------------------------
          ! 5. Fill the configuration state
          !--------------------------------------------------------------

          config_cpl_fields(i)%name                = name
          config_cpl_fields(i)%dt_coupling         = dt_coupling
          config_cpl_fields(i)%dt_model            = dt_model
          config_cpl_fields(i)%lag                 = lag
          config_cpl_fields(i)%l_time_average      = l_time_average
          config_cpl_fields(i)%l_time_accumulation = l_time_accumulation
          config_cpl_fields(i)%l_diagnostic        = l_diagnostic
          config_cpl_fields(i)%l_activated         = l_activated

       END SELECT

       IF ( istat /= POSITIONED ) EXIT

    END DO

    number_of_coupled_variables = i
    config_debug_coupler_level  = debug_coupler_level

#endif

    CALL close_nml

  END SUBROUTINE read_coupling_namelist

END MODULE mo_coupling_nml
