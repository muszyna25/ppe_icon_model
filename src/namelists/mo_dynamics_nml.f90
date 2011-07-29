!>
!! Namelist variables shared by the hydrostatic and nonhydrostatic
!! dynamical cores
!!        
!! @par Revision History
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
!!
MODULE mo_dynamics_nml

  USE mo_dynamics_config,     ONLY: config_iequations     => iequations,     &
                                  & config_idiv_method    => idiv_method,    &
                                  & config_divavg_cntrwgt => divavg_cntrwgt, &
                                  & config_sw_ref_height  => sw_ref_height,  &
                                  & config_lcoriolis      => lcoriolis

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: IHS_ATM_TEMP
  USE mo_physical_constants,  ONLY: grav
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 

  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,   &
                                  & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_dynamics_namelist

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

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

  NAMELIST/dynamics_nml/ iequations,                  &
                         idiv_method, divavg_cntrwgt, &
                         sw_ref_height,  lcoriolis

CONTAINS
  !>
  !!
  SUBROUTINE read_dynamics_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_nml:read_dynamics_namelist'

    !------------------------------------------------------------
    ! Set up the default values
    !------------------------------------------------------------
    iequations     = IHS_ATM_TEMP
    idiv_method    = 1
    divavg_cntrwgt = 0.5_wp
    sw_ref_height  = 0.9_wp*2.94e4_wp/grav
    lcoriolis      = .TRUE.
 
    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('dynamics_nml')
      READ(funit,NML=dynamics_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('dynamics_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, dynamics_nml)
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

  END SUBROUTINE read_dynamics_namelist
  !-------------

END MODULE mo_dynamics_nml
