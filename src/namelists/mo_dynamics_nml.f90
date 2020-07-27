!>
!! Namelist variables shared by the hydrostatic and nonhydrostatic
!! dynamical cores
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
!!
MODULE mo_dynamics_nml

  USE mo_dynamics_config,     ONLY: config_iequations     => iequations,     &
                                  & config_idiv_method    => idiv_method,    &
                                  & config_divavg_cntrwgt => divavg_cntrwgt, &
                                  & config_sw_ref_height  => sw_ref_height,  &
                                  & config_lcoriolis      => lcoriolis,      &
                                  & config_ldeepatmo      => ldeepatmo

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: INH_ATMOSPHERE
  USE mo_physical_constants,  ONLY: grav
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 

  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,   &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_dynamics_namelist

  !---------------------------------------------------------------
  ! Namelist variables 
  !---------------------------------------------------------------
  ! time stepping scheme 

  INTEGER  :: iequations

  ! way of computing the divergence operator in the triangular model -------------

  INTEGER  :: idiv_method    ! 1: Hydrostatic atmospheric model: 
                             !    Gauss integral with original normal 
                             !    velocity components
                             ! 1: Non-Hydrostatic atmospheric model: 
                             !    Gauss integral with averged normal 
                             !    velocity components
                             ! Thus, in a linear equilateral grid, methods 1 and 2 for
                             ! the non-hydrostatic model are the same.
                             ! 2: divergence averaging with bilinear averaging

  REAL(wp) :: divavg_cntrwgt ! weight of central cell for divergence averaging

  REAL(wp) :: sw_ref_height  ! reference height to linearize around if using
                             ! lshallow_water and semi-implicit correction

  LOGICAL  :: lcoriolis      ! if .TRUE.,  the Coriolis force is switched on

  LOGICAL  :: ldeepatmo      ! if .TRUE., deep-atmosphere modification is applied 
                             ! to the non-hydrostatic model equations 

  NAMELIST/dynamics_nml/ iequations,                  &
                         idiv_method, divavg_cntrwgt, &
                         sw_ref_height,  lcoriolis,   &
                         ldeepatmo

CONTAINS
  !>
  !!
  SUBROUTINE read_dynamics_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_nml:read_dynamics_namelist'

    !------------------------------------------------------------
    ! Set up the default values
    !------------------------------------------------------------
    iequations     = INH_ATMOSPHERE
    idiv_method    = 1
    divavg_cntrwgt = 0.5_wp
    sw_ref_height  = 0.9_wp*2.94e4_wp/grav
    lcoriolis      = .TRUE.
    ldeepatmo      = .FALSE.
 
    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('dynamics_nml')
      READ(funit,NML=dynamics_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('dynamics_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, dynamics_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, dynamics_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, dynamics_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Sanity check
    !-----------------------------------------------------

    IF (idiv_method > 2 .OR. idiv_method < 1 )THEN
      CALL finish(TRIM(routine),'Error: idiv_method must be 1 or 2 !')
    ENDIF

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=dynamics_nml)
      CALL store_and_close_namelist(funit, 'dynamics_nml')
    ENDIF
    
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=dynamics_nml)

    !-----------------------------------------------------
    ! 5. Fill configuration state
    !-----------------------------------------------------

    config_iequations     = iequations
    config_idiv_method    = idiv_method
    config_divavg_cntrwgt = divavg_cntrwgt
    config_sw_ref_height  = sw_ref_height
    config_lcoriolis      = lcoriolis
    config_ldeepatmo      = ldeepatmo

  END SUBROUTINE read_dynamics_namelist
  !-------------

END MODULE mo_dynamics_nml
