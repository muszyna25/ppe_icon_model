#define NOMPI
!>
!! Contains output routines for CDI output
!! This module acts merely as a wrapper for either calling the direct
!! or the asynchronous output routines.
!!
!!
!! @par Revision History
!! Initial implementation by Rainer Johanni (2010-12-02)
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
!!
MODULE mo_output

  USE mo_exception,           ONLY: finish, message_text, get_filename_noext
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_parallel_configuration, ONLY: num_io_procs
  USE mo_io_units,            ONLY: filename_max
  USE mo_model_domain_import, ONLY: n_dom, nroot, lplane
  USE mo_ocean_nml,           ONLY: n_zlev
  USE mo_io_nml,              ONLY: out_expname
  USE mo_run_nml,             ONLY: iequations,           &
     &                              ihs_atm_temp,         &
     &                              ihs_atm_theta,        &
     &                              inh_atmosphere,       &
     &                              ishallow_water,ihs_ocean
  USE mo_atmo_control,        ONLY: p_patch
  USE mo_io_vlist,            ONLY: setup_vlist, destruct_vlist,           &
     &                              open_output_vlist, close_output_vlist, &
     &                              write_vlist
  USE mo_io_vlist,            ONLY: setup_vlist_oce
  USE mo_io_async,            ONLY: setup_io_procs, shutdown_io_procs, &
    &                               output_async, set_output_file
  USE mo_datetime,            ONLY: t_datetime,iso8601
  USE mo_io_restart,          ONLY: set_restart_time, set_restart_vct,       &
                                  & init_restart, open_writing_restart_files,  &
                                  & write_restart, close_writing_restart_files,&
                                  & finish_restart
  USE mo_io_restart_attributes,ONLY: set_restart_attribute
  USE mo_dynamics_nml,        ONLY: nold, nnow, nnew, nnow_rcf, nnew_rcf
  USE mo_model_domain,        ONLY: t_patch

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_output'

  LOGICAL :: l_omit_dom   ! flag if "'_DOM',jg" should be omitted in the filenames

  !-------------------------------------------------------------------------------------------------

  ! Public routines:

  PUBLIC :: init_output_files, close_output_files, write_output
  PUBLIC :: create_restart_file


CONTAINS

  !-------------------------------------------------------------------------------------------------
  !>
  !! Initialize output file(s)
  !! 
  !! This can be done for the first  or any other time:
  !! * Set the file name
  !! * Add current namelist values as global attributes (netcdf)
  !! * Create data variable definitions
  !! For jfile=1 some additional refinement initializations are done.

  SUBROUTINE init_output_files(jfile, lclose)

    INTEGER, INTENT(IN) :: jfile !> Number of fileset to open
    LOGICAL, INTENT(IN) :: lclose !> lclose old file

    INTEGER :: jg, jlev
    INTEGER :: nlev              !< number of full levels
    CHARACTER(LEN=filename_max) :: gridtype, outputfile
    CHARACTER(LEN=filename_max) :: grid_filename

    IF(.NOT.lclose) THEN

      ! This is the first call - initialize

      ! If n_dom=1, i.e. if no grid refinement is used on the sphere, then
      ! do not USE "'_DOM',jg" in output file name, if it is not present
      ! in the input file name. Use l_omit_dom to check if this is the case.
      l_omit_dom = .FALSE.

      IF (lplane) THEN
        gridtype='plan'
      ELSE
        gridtype='icon'
      END IF

      DO jg = 1, n_dom

        jlev = p_patch(jg)%level

        ! Grid file name(s) for input
        !
        ! Allow file names without "DOM" specifier if n_dom=1.
!         IF (n_dom == 1) THEN
!           ! Check if file name without "DOM" specifier exists.
!           WRITE (gridfile(jg),'(a,a,i0,a,i2.2,a)') &
!             &    TRIM(gridtype),'R',nroot,'B',jlev,'-grid.nc'
!           INQUIRE (FILE=gridfile(jg), EXIST=l_omit_dom)
!           ! Otherwise use file name with "DOM" specifier
!           IF (.NOT. l_omit_dom)                                            &
!             &    WRITE (gridfile(jg),'(a,a,i0,2(a,i2.2),a)')                &
!             &    TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
!         ELSE
!           ! n_dom >1 --> "'_DOM',jg" required in file name
!           WRITE (gridfile(jg),'(a,a,i0,2(a,i2.2),a)') &
!             &    TRIM(gridtype),'R',nroot,'B',jlev,'_DOM',jg,'-grid.nc'
!         ENDIF

        ! Set up vlist for this grid level
        ! Please note: setup_vlist only sets up the vlist, it does not open any output file!
       
        IF(iequations /= ihs_ocean) CALL setup_vlist( TRIM(p_patch(jg)%grid_filename), jg )

      ENDDO

    ELSE

      ! If not called for the first time, close previous output files
      ! (only if we are actually doing output!)
      IF(num_io_procs == 0 .AND. p_pe == p_io) THEN
        DO jg = n_dom, 1, -1
          CALL close_output_vlist(jg)
        ENDDO
      ENDIF

    ENDIF

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level
      nlev = p_patch(jg)%nlev
      grid_filename = get_filename_noext(p_patch(jg)%grid_filename)
      ! Raw data file name(s) for output
      !
      WRITE (outputfile,'(a,a,i4.4,a)')  &
        &  TRIM(out_expname), TRIM(grid_filename),'_', jfile, '.nc'
      
!       SELECT CASE (iequations)
!         !
!       CASE (ishallow_water)
!         IF (l_omit_dom) THEN
!           WRITE (outputfile,'(a,a,i0,a,i2.2,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_R', nroot, 'B', jlev, '_', jfile, '.nc'
!         ELSE
!           WRITE (outputfile,'(a,a,i2.2,a,i0,a,i2.2,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_DOM', jg, '_R', nroot, 'B', jlev, '_', jfile, '.nc'
!         END IF
!         !
!       CASE (ihs_atm_temp, ihs_atm_theta, inh_atmosphere)
!         IF (l_omit_dom) THEN
!           WRITE (outputfile,'(a,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_R', nroot, 'B', jlev, 'L', nlev, '_', jfile, '.nc'
!         ELSE
!           WRITE (outputfile,'(a,a,i2.2,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
!             &  TRIM(out_expname), '_DOM', jg, '_R', nroot, 'B', jlev, 'L', nlev, '_', jfile, '.nc'
!         END IF
!       CASE (ihs_ocean)
!         WRITE (outputfile,'(a,a,i0,a,i2.2,a,i0,a,i4.4,a)')  &
!           &  TRIM(out_expname), '_O.R', nroot, 'B', jlev, 'L', n_zlev, '_', jfile, '.nc'
!         WRITE(*,'(a,a)') ' control_model: Initial output file for setup_vlist_oce is ', &
!           &              TRIM(outputfile)
!         !
!       CASE DEFAULT
!         CALL finish(modname,'Unsupported value of iequations in init_output_files')
!         !
!       END SELECT

      WRITE(0,'(a,a)') ' Initial output file for setup_vlist is ', &
          &              TRIM(outputfile)
      IF (iequations == ihs_ocean) THEN
        ! #slo# must be aligned with general output
        CALL setup_vlist_oce( p_patch(1:), TRIM(p_patch(jg)%grid_filename), TRIM(outputfile), jg )
      ELSE
        IF(num_io_procs == 0) THEN
          IF(p_pe == p_io) CALL open_output_vlist(TRIM(outputfile), jg)
        ELSE
          CALL set_output_file(outputfile, jg)
        ENDIF
      ENDIF

    ENDDO

    ! Setup I/O PEs if this is the initial call and I/O PEs are enabled
    ! Note that this has to be done AFTER the output files are set!

    IF(jfile == 1 .AND. num_io_procs>0) CALL setup_io_procs()


  END SUBROUTINE init_output_files

  !-------------------------------------------------------------------------------------------------
  !>
  !! Closes output files and finalizes I/O setup
  !! Note: This routine must only be called for the final close, not when the output files
  !! are switched during the run!

  SUBROUTINE close_output_files

    INTEGER jg

    DO jg = n_dom, 1, -1
      IF(num_io_procs == 0 .AND. p_pe == p_io) CALL close_output_vlist(jg)
      CALL destruct_vlist( jg )
    ENDDO

    IF(num_io_procs>0) CALL shutdown_io_procs

  END SUBROUTINE close_output_files

  !-------------------------------------------------------------------------------------------------
  !>
  SUBROUTINE write_output(datetime, z_sim_time)

    TYPE(t_datetime),   INTENT(in) :: datetime
    REAL(wp), OPTIONAL, INTENT(in) :: z_sim_time(n_dom)

!    Proposal by Matthias Raschendorfer for correct output
!
!    INTEGER sec
!
!    sec=NINT(datetime%second)-INT(datetime%second)
!    outptime=datetime
!    IF (sec.NE.0) THEN
!       CALL add_time(REAL(sec,wp),0,0,0,outptime)
!    END IF
!
!    IF(num_io_procs == 0) THEN
!      CALL write_vlist(outptime)
!    ELSE
!      CALL output_async(outptime)
!    ENDIF
!

    IF ( PRESENT(z_sim_time) ) THEN  
      IF(num_io_procs == 0) THEN
        CALL write_vlist(datetime, z_sim_time(1))
      ELSE
        CALL output_async(datetime,z_sim_time(1))
      ENDIF
    ELSE
      IF(num_io_procs == 0) THEN
        CALL write_vlist(datetime)
      ELSE
        CALL output_async(datetime)
      ENDIF
    ENDIF
  
  END SUBROUTINE write_output

  !-------------
  !>
  !! 
  !! Hui Wan (MPI-M, 2011-05)
  !!
  SUBROUTINE create_restart_file( patch, datetime, pvct,      &
                                & jfile, l_have_output     )

    TYPE(t_patch),   INTENT(IN) :: patch
    TYPE(t_datetime),INTENT(IN) :: datetime
    REAL(wp),INTENT(IN) :: pvct(:) 
    INTEGER, INTENT(IN) :: jfile  ! current output file index
    LOGICAL, INTENT(IN) :: l_have_output
    
    INTEGER :: klev, jg, kcell, kvert, kedge, icelltype
   
    CHARACTER(LEN=132) :: string

    !----------------
    ! Initialization
    klev = patch%nlev
    jg   = patch%id
    kcell =patch%n_patch_cells_g
    kvert = patch%n_patch_verts_g
    kedge = patch%n_patch_edges_g
    icelltype = patch%cell_type
    
    CALL set_restart_attribute( 'current_year'  , datetime%year   )   
    CALL set_restart_attribute( 'current_month' , datetime%month  )
    CALL set_restart_attribute( 'current_day'   , datetime%day    )   
    CALL set_restart_attribute( 'current_hour'  , datetime%hour   )   
    CALL set_restart_attribute( 'current_minute', datetime%minute )
    CALL set_restart_attribute( 'current_second', datetime%second )

    CALL set_restart_attribute( 'nold'    , nold    (jg))
    CALL set_restart_attribute( 'nnow'    , nnow    (jg))
    CALL set_restart_attribute( 'nnew'    , nnew    (jg))
    CALL set_restart_attribute( 'nnow_rcf', nnow_rcf(jg))
    CALL set_restart_attribute( 'nnew_rcf', nnew_rcf(jg))

    IF (l_have_output) THEN
      CALL set_restart_attribute( 'next_output_file', jfile+1 )
    ELSE
      CALL set_restart_attribute( 'next_output_file', jfile   )
    END IF

    CALL set_restart_vct( pvct )  ! Vertical coordinate (A's and B's)

    CALL init_restart( TRIM(out_expname), &! exp name
                     & '1.2.2',           &! model version
                     & kcell, icelltype,  &! total # of cells, # of vertices per cell
                     & kvert, 9-icelltype,&! total # of vertices, # of vertices per dual cell
                     & kedge, 4,          &! total # of cells, shape of control volume for edge 
                     & klev              ) ! total # of vertical layers

    CALL set_restart_time( iso8601(datetime) )  ! Time tag

    ! Open new file, write data, close and then clean-up.
    message_text = get_filename_noext(patch%grid_filename)
    WRITE(string,'(a,a)') 'restart.',TRIM(message_text)

    CALL open_writing_restart_files( TRIM(string) )

#ifdef NOMPI
    CALL write_restart
#else
    CALL write_restart( p_patch )
#endif

    CALL close_writing_restart_files
    CALL finish_restart

  END SUBROUTINE create_restart_file

END MODULE mo_output
