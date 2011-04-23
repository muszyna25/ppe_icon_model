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
  USE mo_impl_constants,     ONLY: max_char_length, max_ntracer
  USE mo_datetime,           ONLY: t_datetime, proleptic_gregorian,          &
    &                              date_to_time, add_time, print_datetime_all
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_run_nml,            ONLY: inwp,iecham,ltransport,dtime,ntracer,          &
    &                              iforcing,lshallow_water,iequations,            &
    &                              inextra_2d, inextra_3d,ildf_echam

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

  CHARACTER(len=*), PARAMETER :: modelname    = 'icon'
  CHARACTER(len=*), PARAMETER :: modelversion = 'dev'

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


  NAMELIST/io_ctl/ out_expname, out_filetype, dt_data, dt_file, dt_diag,  &
    &              lwrite_vorticity, lwrite_divergence, lwrite_omega, lwrite_pres, lwrite_z3, &
    &              lwrite_tracer, lwrite_tend_phy, lwrite_radiation, lwrite_precip,           &
    &              lwrite_cloud, lkeep_in_sync,lwrite_tke,lwrite_surface,lwrite_extra


  !
  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent control variables 
  ! -----------------------------------------------------------------------
  !
  LOGICAL :: l_outputtime         ! if .true., output is written at the end of the time step.
  LOGICAL :: l_diagtime           ! if .true., diagnostic output is computed and written
                                  ! at the end of the time step.

  LOGICAL, ALLOCATABLE :: lprepare_output(:) ! For each grid level:
                                             ! if .true., save the prognostic
                                             ! variables to p_prog_out and
                                             ! update p_diag_out.
  !
  !
  !
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
!!               The configuration is read from namelist 'io_ctl'.
!!
!! @par Revision History
!!  Initial version by Almut Gassmann, MPI-M (2008-09-30)
!!  Modofied by Constantin Junk, MPI-M (2010-22-02):
!!     - renamed subroutine setup_io to io_nml_setup
!!     - moved subroutine to new module mo_io_nml
!!
SUBROUTINE io_nml_setup
!
! !local variable
  INTEGER :: i_status

!-----------------------------------------------------------------------

  !------------------------------------------------------------
  ! 3.0 set up the default values for run_ctl
  !------------------------------------------------------------
  
  out_expname   = 'IIIEEEETTTT'
  out_filetype  = 2

  dt_data       = 21600.0_wp
  dt_file       = 2592000._wp
  dt_diag       = dtime
  lkeep_in_sync = .FALSE.

  lwrite_vorticity   = .TRUE.
  lwrite_divergence  = .TRUE.
  lwrite_pres        = .TRUE.
  lwrite_z3          = .TRUE.
  lwrite_omega       = .TRUE.
  IF (ntracer > 0) THEN
    lwrite_tracer(:)   = .TRUE.
  ENDIF

  SELECT CASE(iforcing)
  CASE ( inwp )
    lwrite_precip    = .TRUE.
    lwrite_cloud     = .FALSE.
    lwrite_radiation = .TRUE.
    lwrite_tke       = .FALSE.
    lwrite_surface   = .FALSE.
    lwrite_tend_phy  = .FALSE.
    lwrite_extra     = .FALSE.
    lwrite_extra     = .FALSE. 
  CASE ( iecham, ildf_echam )
    lwrite_precip    = .TRUE.
    lwrite_cloud     = .TRUE.
    lwrite_radiation = .TRUE.
    lwrite_tend_phy  = .TRUE.
    lwrite_surface   = .FALSE.
    lwrite_tke       = .FALSE.
    lwrite_extra     = .FALSE. 
  CASE DEFAULT
    lwrite_precip    = .FALSE.
    lwrite_cloud     = .FALSE.
    lwrite_radiation = .FALSE.
    lwrite_tend_phy  = .FALSE.
    lwrite_surface   = .FALSE.
    lwrite_tke       = .FALSE.
    lwrite_extra     = .FALSE. 
  END SELECT

  !------------------------------------------------------------
  ! 4.0 Read the namelist to evaluate if a restart file shall
  !     be used to configure and initialize the model, and if
  !     so take read the run_ctl parameters from the restart
  !     file.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)

  CALL position_nml ('io_ctl', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, io_ctl)
  END SELECT

  IF (lshallow_water) THEN
     lwrite_z3     = .FALSE.
     lwrite_omega  = .FALSE.
     lwrite_pres   = .FALSE.
  ENDIF
  IF (iequations == 3) THEN
     lwrite_omega  = .FALSE.
  ENDIF

  IF(.NOT. ltransport .AND. ntracer >0 ) THEN
    lwrite_tracer(:) =.FALSE.
      CALL message( 'io_nml_setup', 'lwrite_tracer changed: NO TRACER OUTPUT IF ltransport FALSE ')
  ENDIF


  !------------------------------------------------------------
  ! 5.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
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


!KF  This has to be moved to the module io_namelist

  IF (( inextra_2D > 0) .OR. (inextra_3D > 0) ) THEN 
     lwrite_extra = .TRUE.
       WRITE(message_text,'(a,2I4,a,L4)') &
            &'inextra is',inextra_2d,inextra_3d ,' lwrite_extra has been set', lwrite_extra
      CALL message('io_namelist', TRIM(message_text))
    ENDIF

   IF (inextra_2D == 0 .AND. inextra_3D == 0 .AND. lwrite_extra) &
        CALL finish('io_namelist','need to specify extra fields for extra output')

  ! write the contents of the namelist to an ASCII file

  IF(p_pe == p_io) WRITE(nnml_output,nml=io_ctl)

!
END SUBROUTINE io_nml_setup

END MODULE mo_io_nml

!!$  !>
!!$  !! "read_restart_radiation_nml" reads the parameters of the radiation
!!$  !! namelist from the global attributes of a restart file in netcdf format.
!!$  !!
!!$  !! @par Revision History
!!$  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!$  !!
!!$  SUBROUTINE read_restart_radiation_nml
!!$
!!$    CALL finish('mo_radiation_nml/read_restart_radiation_nml', &
!!$      &         'This subroutine is not yet available')
!!$
!!$  END SUBROUTINE read_restart_radiation_nml
!!$
!!$  !>
!!$  !! "write_restart_radiation_nml" writes the parameters of the radiation
!!$  !! namelist as global attributes to a restart file in netcdf format.
!!$  !!
!!$  !! @par Revision History
!!$  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!$  !!
!!$  SUBROUTINE write_restart_radiation_nml
!!$
!!$    CALL finish('mo_radiation_nml/write_restart_radiation_nml', &
!!$      &         'This subroutine is not yet available')
!!$
!!$  END SUBROUTINE write_restart_radiation_nml
!!$
!!$  !>
!!$  !! "write_rawdata_radiation_nml" writes the parameters of the radiation
!!$  !! namelist as global attributes to a raw data file in netcdf format.
!!$  !!
!!$  !! @par Revision History
!!$  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!$  !!
!!$  SUBROUTINE write_rawdata_radiation_nml
!!$
!!$    CALL finish('mo_radiation_nml/write_rawdata_radiation_nml', &
!!$      &         'This subroutine is not yet available')
!!$
!!$  END SUBROUTINE write_rawdata_radiation_nml
