!>
!! Contains the setup of the variables for io.
!!
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3592)
!!   Modification by Constantin Junk (2011-02-24)
!!     - added new module mo_io_nml
!!     - separated declaration of namelist io_ctl from 
!!       mo_global_variables and moved it mo_io_nml
!!     - minor changes to variable declaration section
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_io_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_impl_constants,     ONLY: max_char_length, max_ntracer, max_dom
  USE mo_datetime,           ONLY: t_datetime, proleptic_gregorian,          &
    &                              date_to_time, add_time, print_datetime_all
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_config,          ONLY: io_config
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist,   &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC


  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters
  ! ------------------------------------------------------------------------
  !
  CHARACTER(len=max_char_length) :: out_expname
  INTEGER :: out_filetype               ! 1 - GRIB1, 2 - netCDF
  LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
  REAL(wp):: dt_data                    ! output timestep [seconds]
  REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
  REAL(wp):: dt_file                    ! timestep [seconds] for triggering new output file
  REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file
  !
  !
  !
  LOGICAL :: lwrite_vorticity           ! if .true., write out vorticity
  LOGICAL :: lwrite_divergence          ! if .true., write out divergence
  LOGICAL :: lwrite_pres                ! if .true., write out full level pressure
  LOGICAL :: lwrite_tend_phy            ! if .true., write out physics-induced tendencies
  LOGICAL :: lwrite_radiation           ! if .true., write out fields related to radiation
  LOGICAL :: lwrite_precip              ! if .true., write out precip
  LOGICAL :: lwrite_cloud               ! if .true., write out cloud variables
  LOGICAL :: lwrite_z3                  ! if .true., write out geopotential on full levels
  LOGICAL :: lwrite_omega               ! if .true., write out the vertical velocity
                                        ! in pressure coordinate
  LOGICAL :: lwrite_tke                 ! if .true., write out TKE
  LOGICAL :: lwrite_surface             ! if .true., write out surface related fields
  LOGICAL :: lwrite_tracer(max_ntracer) ! for each tracer, if .true. write out
                                        ! tracer on full levels
  LOGICAL :: lwrite_extra               ! if .true., write out extra fields


  NAMELIST/io_nml/ out_expname, out_filetype, dt_data, dt_file, dt_diag, dt_checkpoint, &
    &              lwrite_vorticity, lwrite_divergence, lwrite_omega, lwrite_pres, lwrite_z3, &
    &              lwrite_tracer, lwrite_tend_phy, lwrite_radiation, lwrite_precip,           &
    &              lwrite_cloud, lkeep_in_sync,lwrite_tke,lwrite_surface,lwrite_extra


  !
  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent control variables 
  ! -----------------------------------------------------------------------
  !
  LOGICAL :: l_outputtime         ! if .true., output is written at the end of the time step.
  LOGICAL :: l_checkpoint_time    ! if .true., restart file is written at the end of the time step.
  LOGICAL :: l_diagtime           ! if .true., diagnostic output is computed and written
                                  ! at the end of the time step.

  LOGICAL, ALLOCATABLE :: lprepare_output(:) ! For each grid level:
                                             ! if .true., save the prognostic
                                             ! variables to p_prog_out and
                                             ! update p_diag_out.

  PUBLIC :: read_io_namelist

  CONTAINS
!
!-------------------------------------------------------------------------
!
!
!>
!!  Initialization of variables that determine io.
!!
!!               Initialization of variables that determine
!!               some settings of the io.
!!               The configuration is read from namelist 'io_nml'.
!!
!! @par Revision History
!!  Initial version by Almut Gassmann, MPI-M (2008-09-30)
!!  Modofied by Constantin Junk, MPI-M (2010-22-02):
!!     - renamed subroutine setup_io to io_nml_setup
!!     - moved subroutine to new module mo_io_nml
!!
SUBROUTINE io_nml_setup

  INTEGER :: istat, funit

  !------------------------------------------------------------
  ! 5.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
 !SELECT CASE(iforcing)
 !CASE ( inwp )
 !  ! Do nothing. Keep the initial values, if not specified in namelist.
 !  ! consider special idealized testcase with turbulence only
 !  IF( .NOT. ltransport  )   THEN
 !    lwrite_precip    = .FALSE.
 !    lwrite_cloud     = .FALSE.
 !    lwrite_radiation = .FALSE.
 !    lwrite_tke       = .TRUE.
 !    lwrite_surface   = .FALSE.
 !    CALL message('io_nml_setup',' ATTENTION! Only TKE output for TURBULENCE ONLY test')
 !  ENDIF
 !CASE ( iecham,ildf_echam )
 !  ! Do nothing. Keep the initial values, if not specified in namelist.
 !CASE DEFAULT
 !  ! Do nothing. Keep the initial values, if not specified in namelist.
 !END SELECT

 !IF (( inextra_2D > 0) .OR. (inextra_3D > 0) ) THEN 
 !   lwrite_extra = .TRUE.
 !     WRITE(message_text,'(a,2I4,a,L4)') &
 !          &'inextra is',inextra_2d,inextra_3d ,' lwrite_extra has been set', lwrite_extra
 !    CALL message('io_namelist', TRIM(message_text))
 !  ENDIF

 ! IF (inextra_2D == 0 .AND. inextra_3D == 0 .AND. lwrite_extra) &
 !      CALL finish('io_namelist','need to specify extra fields for extra output')

  END SUBROUTINE io_nml_setup


  !>
  !! Read Namelist for I/O. 
  !!
  !! This subroutine 
  !! - reads the Namelist for I/O
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
  SUBROUTINE read_io_namelist
    !
    INTEGER :: istat, funit
    INTEGER :: jg           ! loop index

    CHARACTER(len=*), PARAMETER :: routine = 'mo_io_nml: read_io_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    out_expname   = 'IIIEEEETTTT'
    out_filetype  = 2

    dt_data       = 21600.0_wp   !  6 hours
    dt_file       = 2592000._wp  ! 30 days
    dt_checkpoint = 2592000._wp  ! 30 days
    dt_diag       = 86400._wp    !  1 time step
    lkeep_in_sync = .FALSE.

    lwrite_vorticity   = .TRUE.
    lwrite_divergence  = .TRUE.
    lwrite_pres        = .TRUE.
    lwrite_z3          = .TRUE.
    lwrite_omega       = .TRUE.
    lwrite_tracer(:)   = .TRUE.

    lwrite_precip    = .FALSE.
    lwrite_cloud     = .FALSE.
    lwrite_radiation = .FALSE.
    lwrite_tend_phy  = .FALSE.
    lwrite_surface   = .FALSE.
    lwrite_tke       = .FALSE.
    lwrite_extra     = .FALSE. 

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('io_nml')
      READ(funit,NML=io_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processors)
    !-------------------------------------------------------------------------
    CALL position_nml ('io_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, io_nml)
    END SELECT

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg= 1,max_dom
      io_config(jg)%out_expname      = out_expname
      io_config(jg)%out_filetype     = out_filetype
      io_config(jg)%dt_data          = dt_data
      io_config(jg)%dt_file          = dt_file
      io_config(jg)%dt_diag          = dt_diag
      io_config(jg)%dt_checkpoint    = dt_checkpoint
      io_config(jg)%lwrite_vorticity = lwrite_vorticity
      io_config(jg)%lwrite_divergence= lwrite_divergence 
      io_config(jg)%lwrite_omega     = lwrite_omega
      io_config(jg)%lwrite_pres      = lwrite_pres 
      io_config(jg)%lwrite_z3        = lwrite_z3
      io_config(jg)%lwrite_tracer    = lwrite_tracer
      io_config(jg)%lwrite_tend_phy  = lwrite_tend_phy
      io_config(jg)%lwrite_radiation = lwrite_radiation
      io_config(jg)%lwrite_precip    = lwrite_precip
      io_config(jg)%lwrite_cloud     = lwrite_cloud
      io_config(jg)%lkeep_in_sync    = lkeep_in_sync
      io_config(jg)%lwrite_tke       = lwrite_tke
      io_config(jg)%lwrite_surface   = lwrite_surface
      io_config(jg)%lwrite_extra     = lwrite_extra
    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=io_nml)                    
    CALL store_and_close_namelist(funit, 'io_nml') 

    ! 6. write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=io_nml)

  END SUBROUTINE read_io_namelist


END MODULE mo_io_nml
