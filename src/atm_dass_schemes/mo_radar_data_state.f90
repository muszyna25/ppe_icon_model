!>
!! Allocation/deallocation and reading of radar datasets for LHN
!!
!! This module contains routines for setting up the radar data state.
!!
!! @author Klaus Stephan, DWD
!!
!! @par Revision History
!! Initial revision by Klaus Stephan, DWD (2015-01-30)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_radar_data_state

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: inwp, max_char_length, SUCCESS
  USE mo_run_config,         ONLY: iforcing, check_uuid_gracefully
  USE mo_assimilation_config,ONLY: assimilation_config
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_linked_list,        ONLY: t_var_list

  USE mo_radar_data_types,   ONLY: t_radar_fields,t_radar_td_fields, t_radar_ct_fields, t_lhn_diag

  USE mo_var_list,           ONLY: default_var_list_settings,   &
    &                              add_var, new_var_list,                &
    &                              delete_var_list
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                ONLY: DATATYPE_PACK16, DATATYPE_FLT32,                 &
    &                              TSTEP_CONSTANT, TSTEP_INSTANT, cdi_undefid, &
    &                              streamClose, gridInqUUID, GRID_UNSTRUCTURED,    &
    &                              FILETYPE_GRB2, streamOpenRead,   &
    &                              streamInqFileType, streamInqVlist, vlistNtsteps,&
    &                              vlistInqTaxis, streamInqTimestep, taxisInqVdate,&
    &                              taxisInqVtime, vlistInqVarGrid, vlistInqVarName
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_zaxis_type,         ONLY: ZA_HYBRID, ZA_SURFACE
  USE mo_util_cdi,           ONLY: get_cdi_varID, read_cdi_2d, t_inputParameters,  &
    &                              makeInputParameters, deleteInputParameters
  USE mo_util_uuid_types,    ONLY: t_uuid, uuid_string_length
  USE mo_util_uuid,          ONLY: OPERATOR(==), uuid_unparse
  USE mo_dictionary,         ONLY: t_dictionary, dict_init, dict_finalize,         &
    &                              dict_loadfile
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE


  ! required for reading radar data
  INCLUDE 'netcdf.inc'
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_radar_data_state'

  PRIVATE

  ! Number of landcover classes provided by radar parameter data
  ! Needs to be changed into a variable if landcover classifications 
  ! with a different number of classes become available
  PUBLIC :: radar_data
  PUBLIC :: init_radar_data  
  PUBLIC :: lhn_fields
  PUBLIC :: construct_lhn
  PUBLIC :: destruct_lhn
!  PUBLIC :: destruct_radar_data

  TYPE(t_radar_fields),TARGET, ALLOCATABLE :: &
    &  radar_data(:)  ! n_dom

  TYPE(t_lhn_diag),TARGET, ALLOCATABLE :: &
    &  lhn_fields(:)  ! n_dom

!-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Init radar data for atmosphere
  !!
  !! 1. Build data structure, including field lists and 
  !!    memory allocation.
  !! 2. External data are read in from netCDF file or set analytically
  !!
  !! @par Revision History
  !! Initial revision by Klaus Stephan, DWD (2015-01-30) based on mo_radar_data_state.f90
  !!
  SUBROUTINE init_radar_data (p_patch, radar_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_radar_fields), INTENT(INOUT) :: radar_data(:)


    INTEGER              :: jg, ist
    INTEGER, ALLOCATABLE :: cdi_radar_id(:)  !< CDI stream ID (for each domain)
    INTEGER, ALLOCATABLE :: cdi_filetype(:)   !< CDI filetype (for each domain)
    ! dictionary which maps internal variable names onto
    ! GRIB2 shortnames or NetCDF var names.
    TYPE (t_dictionary) :: radar_varnames_dict

!    TYPE(t_datetime) :: datetime
    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':init_radar_data',     &
      radar_varnames_map_file = 'radar_dict_file.dat'
    LOGICAL :: lread_process

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    !-------------------------------------------------------------------------
    IF (ANY (assimilation_config(:)%luse_rad)) THEN
      WRITE(message_text,'(a)') 'start reading radar data'
      CALL message (TRIM(routine),message_text)
    ELSE
      WRITE(message_text,'(a)') 'llhn and llhnverif is set to false, nothing to do'
      CALL message (TRIM(routine),message_text)
      RETURN
    ENDIF
    !-------------------------------------------------------------------------
    !  1.  inquire radar files for their data structure
    !-------------------------------------------------------------------------

    ! Allocate and open CDI stream (files):
    ALLOCATE (cdi_radar_id(n_dom), cdi_filetype(n_dom), stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'ALLOCATE failed!')
    CALL inquire_radar_files(p_patch, cdi_radar_id, cdi_filetype, lread_process)


    ! read the map file (internal -> GRIB2) into dictionary data structure:
    CALL dict_init(radar_varnames_dict, lcase_sensitive=.FALSE.)
    IF (ANY(cdi_filetype(:) == FILETYPE_GRB2)) THEN
      IF(radar_varnames_map_file /= ' ') THEN
        CALL dict_loadfile(radar_varnames_dict, TRIM(radar_varnames_map_file))
      END IF
    END IF

    !------------------------------------------------------------------
    !  2.  construct radar fields for the model
    !------------------------------------------------------------------

    ! top-level procedure for building data structures for 
    ! radar data.
    CALL construct_radar_data(p_patch, radar_data)

    !-------------------------------------------------------------------------
    !  3.  read the data into the fields
    !-------------------------------------------------------------------------

    ! Check, whether radar data should be read from file

      CALL message( TRIM(routine),'Start reading radar data from file' )

      CALL read_radar_data (p_patch, radar_data, cdi_radar_id, &
        &                     radar_varnames_dict)

    CALL message( TRIM(routine),'Finished reading radar data' )

    ! close CDI stream (file):
    IF (lread_process) THEN
      DO jg=1,n_dom
        IF (cdi_radar_id(jg) /= cdi_undefid) CALL streamClose(cdi_radar_id(jg))
      END DO
    END IF
    DEALLOCATE (cdi_radar_id, cdi_filetype, stat=ist)
    IF (ist /= SUCCESS)  CALL finish(TRIM(routine),'DEALLOCATE failed!')

    ! destroy variable name dictionary:
    CALL dict_finalize(radar_varnames_dict)

  END SUBROUTINE init_radar_data


  !-------------------------------------------------------------------------
  !>
  !! Top-level procedure for building radar data structure
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE construct_radar_data (p_patch, radar_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_radar_fields),  INTENT(INOUT) :: radar_data(:)

    INTEGER :: jg

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':construct_radar_data'

!-------------------------------------------------------------------------


    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'radar data started')

    ! Build radar data list for constant-in-time fields for the atm model
    DO jg = 1, n_dom
      IF (assimilation_config(jg)%luse_rad) THEN
        WRITE(listname,'(a,i2.2)') 'radar_data_ct_dom_',jg
        CALL new_radar_data_ct_list(p_patch(jg),radar_data(jg)%radar_ct,radar_data(jg)%radar_ct_list, TRIM(listname), &
                                    assimilation_config(jg)%nradar)
      ENDIF
    END DO

    ! Build radar data list for time-dependent fields
    DO jg = 1, n_dom
      IF (assimilation_config(jg)%luse_rad) THEN
        WRITE(listname,'(a,i2.2)') 'radar_data_td_dom_',jg
        CALL new_radar_data_td_list(p_patch(jg),radar_data(jg)%radar_td,radar_data(jg)%radar_td_list,TRIM(listname),  &
                                    assimilation_config(jg)%nobs_times)
      ENDIF
    END DO

    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'radar data finished')

  END SUBROUTINE construct_radar_data

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric radar data structure
  !!
  !! Allocation of atmospheric radar data structure (constant in time 
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Klaus Stephan (2011-05-03)
  !! Statements that assign initial value added by Hui Wan (MPI-M, 2011-05-30)
  !!
  SUBROUTINE new_radar_data_ct_list ( p_patch, p_radar_ct, p_radar_ct_list, &
    &                                listname, nradar)
!
    TYPE(t_patch), TARGET, INTENT(IN)     :: & !< current patch
      &  p_patch

    TYPE(t_radar_ct_fields), INTENT(INOUT) :: & !< current radar data structure
      &  p_radar_ct 

    TYPE(t_var_list) :: p_radar_ct_list !< current radar data list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER, INTENT(IN) :: nradar

    INTEGER :: nblks_c    !< number of cell blocks to allocate

    INTEGER :: shape2d_c(2)
    INTEGER :: shape3d_c(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice


    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! get patch ID
    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    ! number of vertical levels

    ! predefined array shapes
    shape2d_c  = (/ nproma, nblks_c /)
    shape3d_c  = (/ nproma, nradar, nblks_c       /)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_radar_ct_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_radar_ct_list,            &
                                  & lrestart=.FALSE.  )


    ! radar blacklist at cell center
    !
    ! blacklist  p_radar_ct%blacklist(nproma,nblks_c)
    cf_desc    = t_cf_var('rad_bl', '-', &
      &                   'blacklist information', DATATYPE_FLT32)
    grib2_desc = grib2_var( 0, 15, 197, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_radar_ct_list, 'rad_bl', p_radar_ct%blacklist,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,             &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE.,             &
      &           isteptype=TSTEP_CONSTANT )

    ! radar beam height
    !
    ! dxheight  p_radar_ct%dxheight(nproma,nradar,nblks_c)
    cf_desc    = t_cf_var('rad_height', 'm', &
      &                   'radar beam height above ground', DATATYPE_FLT32)
    grib2_desc = grib2_var( 0, 15, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_radar_ct_list, 'rad_height', p_radar_ct%dxheight,                      &
      &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc,              &
      &           grib2_desc, ldims=shape3d_c, loutput=.FALSE. )

  END SUBROUTINE new_radar_data_ct_list

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric radar data structure (time dependent)
  !!
  !! Allocation of atmospheric radar data structure (time dependent  
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Klaus Stephan (2011-05-03)
  !!
  SUBROUTINE new_radar_data_td_list ( p_patch, p_radar_td, &
    &                               p_radar_td_list, listname, nobs)
!
    TYPE(t_patch), TARGET, INTENT(IN)     :: & !< current patch
      &  p_patch

    TYPE(t_radar_td_fields), INTENT(INOUT) :: & !< current radar data structure
      &  p_radar_td 

    TYPE(t_var_list) :: p_radar_td_list  !< current radar data list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER, INTENT(IN) :: nobs
    INTEGER :: nblks_c      !< number of cell blocks to allocate

    INTEGER :: shape3d_c(3)

    INTEGER :: ibits         !< "entropy" of horizontal slice

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':new_radar_data_td_list'
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ibits  = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_c   = (/ nproma, nblks_c, nobs/)

    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_radar_td_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_radar_td_list,         &
                                  & lrestart=.FALSE.,          &
                                  & loutput=.FALSE.  )


    ! radobs       p_radar_td%obs(nproma,nblks_c,nobs_times)
    cf_desc    = t_cf_var('rad_precip', 'mm/h',   &
      &                   'radar derived precipitation rate', DATATYPE_FLT32)
    grib2_desc = grib2_var( 0, 15, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_radar_td_list, 'rad_precip', p_radar_td%obs, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
      &           grib2_desc, ldims=shape3d_c, loutput=.TRUE.,     &
      &           isteptype=TSTEP_INSTANT )  ! Meta info constituentType missing

    ! radqual       p_radar_td%spqual(nproma,nblks_c,nobs_times)
    cf_desc    = t_cf_var('rad_qual', '-',   &
      &                   'quality index of radar observation', DATATYPE_FLT32)
    grib2_desc = grib2_var( 0, 15, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_radar_td_list, 'rad_qual', p_radar_td%spqual, &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,      &
      &           grib2_desc, ldims=shape3d_c, loutput=.FALSE.,     &
      &           isteptype=TSTEP_CONSTANT )  ! Meta info constituentType missing

  END SUBROUTINE new_radar_data_td_list
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  !>
  !! Destruct radar data data structure and lists
  !!
  !! Destruct radar data data structure and lists
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-05-04)
  !!
  SUBROUTINE destruct_radar_data

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = modname//':destruct_radar_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'radar data started')

    DO jg = 1,n_dom
      ! Delete list of constant in time atmospheric elements
      CALL delete_var_list( radar_data(jg)%radar_ct_list)
    ENDDO

    IF (iforcing > 1 ) THEN
    DO jg = 1,n_dom
      ! Delete list of time-dependent atmospheric elements
      CALL delete_var_list( radar_data(jg)%radar_td_list)
    ENDDO
    END IF

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'radar data finished')

  END SUBROUTINE destruct_radar_data
  !-------------------------------------------------------------------------



  !-------------------------------------------------------------------------
  ! Open ExtPar file and investigate the data structure of the
  ! radar parameters.
  !
  ! Note: This subroutine opens the file and returns a CDI file ID.
  !
  ! @author F. Prill, DWD (2014-01-07)
  !-------------------------------------------------------------------------
  SUBROUTINE inquire_radar_file(p_patch, jg, cdi_radar_id, cdi_filetype, &
       lread_process)

    TYPE(t_patch), INTENT(IN)      :: p_patch(:)
    INTEGER,       INTENT(IN)      :: jg
    INTEGER,       INTENT(INOUT)   :: cdi_radar_id     !< CDI stream ID
    INTEGER,       INTENT(INOUT)   :: cdi_filetype     !< CDI filetype
    LOGICAL,       INTENT(in)      :: lread_process

    ! local variables
    CHARACTER(len=max_char_length), PARAMETER :: routine = modname//'::inquire_radar_file'
    INTEGER                 :: vlist_id, radarobs_id
    LOGICAL                 :: l_exist
    CHARACTER(filename_max) :: radar_file !< file name for reading in

    TYPE(t_uuid)            :: radar_uuidOfHGrid             ! same, but converted to TYPE(t_uuid)

    CHARACTER(len=uuid_string_length) :: grid_uuid_unparsed   ! unparsed grid uuid (human readable)
    CHARACTER(len=uuid_string_length) :: radar_uuid_unparsed ! same for radar-file uuid

    LOGICAL                 :: lmatch                         ! for comparing UUIDs

    INTEGER                 :: cdiGridID, taxisID, nt,ntvars,is
    INTEGER, ALLOCATABLE    :: rdate(:),rtime(:)
    CHARACTER(len=MAX_CHAR_LENGTH)  :: zname

    CHARACTER(len=8) cobsdate
    CHARACTER(len=6) cobstime

    !---------------------------------------------!
    ! Check validity of radar parameter file   !
    !---------------------------------------------!


    ! flag. if true, then this PE reads data from file and broadcasts
    IF (lread_process) THEN
      ! generate file name

      radar_file = TRIM(assimilation_config(jg)%radar_in)//TRIM(assimilation_config(jg)%radardata_file)
      CALL message(routine, "radar_file = "//TRIM(radar_file))

      INQUIRE (FILE=radar_file, EXIST=l_exist)
      IF (.NOT.l_exist)  THEN
         CALL message(routine,'radar data file is not found.')
         RETURN
      ENDIF

!return with error_code

      ! open file
      cdi_radar_id = streamOpenRead(TRIM(radar_file))
      cdi_filetype  = streamInqFileType(cdi_radar_id)

      ! get the id of RAD_PRECIP
      radarobs_id          = get_cdi_varID(cdi_radar_id, "RAD_PRECIP")
!      radarheight_id       = get_cdi_varID(cdi_radar_id, "RAD_HEIGHT")

      vlist_id             = streamInqVlist(cdi_radar_id)
      CALL vlistInqVarName(vlist_id, 0, zname)
      ntvars=vlistNtsteps(vlist_id)

      IF (cdi_filetype > FILETYPE_GRB2) THEN
        assimilation_config(jg)%nobs_times=ntvars
      ELSE
        ntvars=assimilation_config(jg)%nobs_times
      ENDIF

      ALLOCATE (rdate(ntvars),rtime(ntvars))

      ! get time dimension from radar data file
      taxisID   = vlistInqTaxis(vlist_id)
      do nt = 0, ntvars-1
         is=streamInqTimestep(cdi_radar_id,nt)
         rdate(nt+1)=taxisInqVdate(taxisID)
         rtime(nt+1)=taxisInqVtime(taxisID)
!       CALL print_cdi_summary (vlist_id)
      enddo

      ! Compare UUID of radar observation file with UUID of grid.
      !
      ! get horizontal grid UUID contained in radar file
      cdiGridID = vlistInqVarGrid(vlist_id, radarobs_id)
      CALL gridInqUUID(cdiGridID, radar_uuidOfHGrid%DATA)
      !
      ! --- compare UUID of horizontal grid file with UUID from radar file
      lmatch = (p_patch(jg)%grid_uuid == radar_uuidOfHGrid)

      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch(jg)%grid_uuid, grid_uuid_unparsed)
        CALL uuid_unparse(radar_uuidOfHGrid   , radar_uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from gridfile: ', TRIM(grid_uuid_unparsed)
        CALL message(routine,message_text)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from radar file: ', TRIM(radar_uuid_unparsed)
        CALL message(routine,message_text)

        WRITE(message_text,'(a)') 'radar file and horizontal grid file do not match!'
        IF (check_uuid_gracefully) THEN
          CALL message(routine, TRIM(message_text))
        ELSE
          CALL finish(routine, TRIM(message_text))
        END IF
      ENDIF      

    ENDIF ! lread_process

    ! broadcast ntvars from I-Pe to WORK Pes
    CALL p_bcast(assimilation_config(jg)%nobs_times, p_io, mpi_comm)
    IF (.NOT. lread_process) THEN
      ALLOCATE (rdate(assimilation_config(jg)%nobs_times),rtime(assimilation_config(jg)%nobs_times))
    ENDIF

    ! broadcast time and date of obs from I-Pe to WORK Pes
    CALL p_bcast(rdate, p_io, mpi_comm)
    CALL p_bcast(rtime, p_io, mpi_comm)

    ALLOCATE (radar_data(jg)%radar_td%obs_date(assimilation_config(jg)%nobs_times))
!    ALLOCATE (radar_data(jg)%radar_td%obs_time(assimilation_config(jg)%nobs_times))

    DO nt = 1, assimilation_config(jg)%nobs_times

      WRITE(cobsdate,'(I8.8)') rdate(nt)
      WRITE(cobstime,'(I6.6)') rtime(nt)
      READ (cobsdate(1:4),'(I4.4)') radar_data(jg)%radar_td%obs_date(nt)%date%year
      READ (cobsdate(5:6),'(I2.2)') radar_data(jg)%radar_td%obs_date(nt)%date%month
      READ (cobsdate(7:8),'(I2.2)') radar_data(jg)%radar_td%obs_date(nt)%date%day
      READ (cobstime(1:2),'(I2.2)') radar_data(jg)%radar_td%obs_date(nt)%time%hour
      READ (cobstime(3:4),'(I2.2)') radar_data(jg)%radar_td%obs_date(nt)%time%minute
      READ (cobstime(5:6),'(I2.2)') radar_data(jg)%radar_td%obs_date(nt)%time%second
      radar_data(jg)%radar_td%obs_date(nt)%time%ms = 0
    ENDDO

!print *,assimilation_config(jg)%nobs_times,radar_data(jg)%radar_td%obs_date

  END SUBROUTINE inquire_radar_file

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE inquire_radar_files(p_patch, cdi_radar_id, cdi_filetype, &
       lread_process)

    !-------------------------------------------------------
    !
    ! open netcdf files and investigate the data structure  
    ! of the radar observations
    !
    !-------------------------------------------------------

    TYPE(t_patch), INTENT(IN)      :: p_patch(:)
    INTEGER,       INTENT(INOUT)   :: cdi_radar_id(:)  !< CDI stream ID
    INTEGER,       INTENT(INOUT)   :: cdi_filetype(:)   !< CDI filetype
    LOGICAL,       INTENT(in)      :: lread_process
    INTEGER :: jg

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':inquire_radar_files'

!--------------------------------------------------------------------------

    ! set stream IDs to "uninitialized":
    IF (lread_process) THEN
      cdi_radar_id(:) = cdi_undefid
    END IF

    DO jg= 1,n_dom

      IF (assimilation_config(jg)%luse_rad) THEN

      !------------------------------------------------!
      ! 1. Check validity of radar files               !
      !------------------------------------------------!
      ! 1.1. radar observations                        !
        CALL inquire_radar_file(p_patch, jg, cdi_radar_id(jg), &
             &                  cdi_filetype(jg), lread_process)
!      ! 1.2. radar blacklist
!        CALL inquire_radar_file(p_patch, jg, cdi_radar_id(jg), cdi_filetype(jg), &
!          &                     nobs_times)
!      ! 1.3. radar height
!        CALL inquire_radar_file(p_patch, jg, cdi_radar_id(jg), cdi_filetype(jg), &
!          &                     nobs_times)

!evt. ist es sinnvoll die grib files vorher zusammen zu ketten
      ENDIF

    ENDDO ! ndom

  END SUBROUTINE inquire_radar_files

  !-------------------------------------------------------------------------
  !>
  !! Read atmospheric radar data
  !!
  !! Read atmospheric radar data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !!
  SUBROUTINE read_radar_data (p_patch, radar_data, cdi_radar_id, &
    &                           radar_varnames_dict)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_radar_fields),    INTENT(INOUT) :: radar_data(:)

    INTEGER,               INTENT(IN)    :: cdi_radar_id(:)      !< CDI stream ID
    TYPE (t_dictionary),   INTENT(IN)    :: radar_varnames_dict  !< variable names dictionary (for GRIB2)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = modname//':read_radar_data'

    INTEGER :: jg
    INTEGER :: jk
    INTEGER :: ret

    TYPE(t_inputParameters) :: parameters

    !----------------------------------------------------------------------


    !------------------------------------------------!
    ! Read data from ExtPar file                     !
    !------------------------------------------------!

    IF (iforcing == inwp) THEN
        DO jg = 1,n_dom
         IF (assimilation_config(jg)%luse_rad) THEN

          radar_data(jg)%radar_td%obs = -0.1_wp

        !--------------------------------------------------------------------
        !
        ! Read radar data for triangle centers (triangular grid)
        !
        !--------------------------------------------------------------------
          parameters = makeInputParameters(cdi_radar_id(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c) ! &
!        &                                opt_dict=radar_varnames_dict)

        ! time dependent
          CALL read_cdi_2d(parameters, assimilation_config(jg)%nobs_times, 'RAD_PRECIP', radar_data(jg)%radar_td%obs)
!        WRITE(message_text,'(a,f10.5)' ) "Max Obs : ",MAXVAL(radar_data(jg)%radar_td%obs)
!        CALL message(routine,message_text)
!        CALL read_cdi_2d(parameters, assimilation_config(jg)%nobs_times, 'RAD_QUAL', radar_data(jg)%radar_td%spqual)
        ! time indipendent
!        CALL read_cdi_2d(parameters, 'RAD_BL', radar_data(jg)%radar_ct%blacklist)
!        CALL read_cdi_3d(parameters, 'RAD_HEIGHT', assimilation_config(jg)%radar_ct%nradar, radar_data(jg)%dxheight)

          CALL deleteInputParameters(parameters)

 
         ENDIF
        END DO ! jg

   END IF ! inwp


  END SUBROUTINE read_radar_data
  !-------------------------------------------------------------------------

 SUBROUTINE construct_lhn (lhn_fields,p_patch)

  TYPE (t_lhn_diag), INTENT(INOUT) :: lhn_fields(:)
  TYPE (t_patch), INTENT(IN) :: p_patch(:)
  INTEGER :: jg

  DO jg = 1, n_dom
    ALLOCATE (lhn_fields(jg)%ttend_lhn(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_c))
    ALLOCATE (lhn_fields(jg)%qvtend_lhn(nproma,p_patch(jg)%nlev,p_patch(jg)%nblks_c))
    ALLOCATE (lhn_fields(jg)%pr_obs_sum(nproma,p_patch(jg)%nblks_c))
    ALLOCATE (lhn_fields(jg)%pr_mod_sum(nproma,p_patch(jg)%nblks_c))
    ALLOCATE (lhn_fields(jg)%pr_ref_sum(nproma,p_patch(jg)%nblks_c))
  
!$OMP PARALLEL 
    CALL init(lhn_fields(jg)%ttend_lhn(:,:,:))
    CALL init(lhn_fields(jg)%qvtend_lhn(:,:,:))
    CALL init(lhn_fields(jg)%pr_obs_sum(:,:))
    CALL init(lhn_fields(jg)%pr_mod_sum(:,:))
    CALL init(lhn_fields(jg)%pr_ref_sum(:,:))
!$OMP END PARALLEL

  ENDDO
 END SUBROUTINE construct_lhn

!===============================================================================

 SUBROUTINE destruct_lhn (lhn_fields)

  TYPE (t_lhn_diag), INTENT(INOUT) :: lhn_fields(:)
  INTEGER :: jg

  DO jg = 1, n_dom
    DEALLOCATE (lhn_fields(jg)%ttend_lhn)
    DEALLOCATE (lhn_fields(jg)%qvtend_lhn)
    DEALLOCATE (lhn_fields(jg)%pr_obs_sum)
    DEALLOCATE (lhn_fields(jg)%pr_mod_sum)
    DEALLOCATE (lhn_fields(jg)%pr_ref_sum)
  
  ENDDO
 END SUBROUTINE destruct_lhn

!===============================================================================

END MODULE mo_radar_data_state

