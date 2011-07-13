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

  USE mo_dynamics_config
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, finish
  USE mo_impl_constants,     ONLY: UNKNOWN, LEAPFROG_SI, IHS_ATM_TEMP
  USE mo_physical_constants, ONLY: grav
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io

  USE mo_master_nml,            ONLY: lrestart
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist,   &       
                                    & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_dynamics_namelist
  PUBLIC :: dynamics_nml_setup

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'


  !---------------------------------------------------------------
  ! Namelist variables and auxiliary parameters setting up the
  ! configuration of the dynamical core
  !---------------------------------------------------------------
  ! time stepping scheme 

  INTEGER  :: nml_iequations

  INTEGER  :: nml_itime_scheme   ! parameter used to select the time stepping scheme
                                 ! = 1, explicit 2 time level scheme
                                 ! = 2, semi implicit 2 time level scheme
                                 ! = 3, explicit leapfrog
                                 ! = 4, leapfrog with semi implicit correction
                                 ! = 5, 4-stage Runge-Kutta method
                                 ! = 6, SSPRK(5,4) (Runge-Kutta) method

  ! way of computing the divergence operator in the triangular model -------------

  INTEGER  :: nml_idiv_method  ! 1: Hydrostatic atmospheric model: 
                               !    Gauss integral with original normal 
                               !    velocity components
                               ! 1: Non-Hydrostatic atmospheric model: 
                               !    Gauss integral with averged normal 
                               !    velocity components
                               ! Thus, in a linear equilateral grid, methods 1 and 2 for
                               ! the non-hydrostatic model are the same.
                               ! 2: divergence averaging with bilinear averaging

  REAL(wp) :: nml_divavg_cntrwgt  ! weight of central cell for divergence averaging

  LOGICAL :: nml_ldry_dycore  ! if .TRUE., ignore the effact of water vapor,
                              ! cloud liquid and cloud ice on virtual temperature.

  REAL(wp) :: nml_sw_ref_height  ! reference height to linearize around if using
                                 ! lshallow_water and semi-implicit correction

  LOGICAL  :: nml_lcoriolis      ! if .TRUE.,  the Coriolis force is switched on

  NAMELIST/dynamics_nml/ nml_iequations,  nml_itime_scheme,   &
                         nml_idiv_method, nml_divavg_cntrwgt, &
                         nml_ldry_dycore, nml_sw_ref_height,  &
                         nml_lcoriolis

CONTAINS
  !>
  !!
  SUBROUTINE read_dynamics_namelist()

    INTEGER :: istat, funit
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_nml:read_dynamics_namelist'

    !------------------------------------------------------------
    ! Set up the default values
    !------------------------------------------------------------
    nml_iequations     = IHS_ATM_TEMP
    nml_itime_scheme   = LEAPFROG_SI
    nml_idiv_method    = 1
    nml_divavg_cntrwgt = 0.5_wp
    nml_ldry_dycore    = .FALSE.
    nml_sw_ref_height  = 0.9_wp*2.94e4_wp/grav
    nml_lcoriolis      = .TRUE.
 
    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('dynamics_nml')
      READ(funit,NML=dynamics_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL position_nml ('dynamics_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, dynamics_nml)
    END SELECT

    !-----------------------------------------------------
    ! Sanity check
    !-----------------------------------------------------
    IF((nml_itime_scheme<=0).OR.(nml_itime_scheme>=unknown)) THEN
      WRITE(message_text,'(A,i2)') &
      'wrong value of itime_scheme, must be 1 ...', unknown -1
      CALL finish( TRIM(routine),TRIM(message_text))
    ENDIF

    IF (nml_idiv_method > 2 .OR. nml_idiv_method < 1 )THEN
      CALL finish(TRIM(routine),'Error: idiv_method must be 1 or 2 !')
    ENDIF

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=dynamics_nml)
    CALL store_and_close_namelist(funit, 'dynamics_nml')

    ! write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=dynamics_nml)

    !-----------------------------------------------------
    ! 5. Fill configuration state
    !-----------------------------------------------------

    iequations     = nml_iequations
    itime_scheme   = nml_itime_scheme
    idiv_method    = nml_idiv_method
    divavg_cntrwgt = nml_divavg_cntrwgt
    ldry_dycore    = nml_ldry_dycore
    sw_ref_height  = nml_sw_ref_height
    lcoriolis      = nml_lcoriolis

  END SUBROUTINE read_dynamics_namelist

  !>
  !!
  SUBROUTINE dynamics_nml_setup(i_ndom)

    INTEGER, INTENT(IN) :: i_ndom
 
    INTEGER :: jdom

    CHARACTER(len=*),PARAMETER ::             &
             & routine = 'mo_dynamics_nml:dynamics_setup'
 
   !CHARACTER(len=MAX_CHAR_LENGTH) :: string

    !------------------------------------------------------------
    ! 4. Check the consistency of the parameters
    !------------------------------------------------------------
    ! for the time stepping scheme
  
  ! IF((itime_scheme==tracer_only).AND.(.NOT.ltransport)) THEN
  !   WRITE(string,'(A,i2,A)') &
  !   'itime_scheme set to ', tracer_only, 'but ltransport to .FALSE.'
  !   CALL finish( TRIM(routine),TRIM(string))
  ! ENDIF
  
  ! IF(ltransport .AND. ntracer <= 0 .AND. iforcing/=INWP) THEN
  !   CALL finish( TRIM(routine),'tracer advection not possible for ntracer=0')
  !   ! [for nwp forcing ntracer setting is treated in setup_transport]
  ! ENDIF
  
  ! IF (global_cell_type/=3) THEN
  !   idiv_method = 1
  ! ENDIF
  
  END SUBROUTINE dynamics_nml_setup
  !-------------
END MODULE mo_dynamics_nml
