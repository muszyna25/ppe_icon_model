!>
!! Allocation/deallocation and reading of HAMOCC boundary conditions
!!
!! @author Irene Stemmler, MPI 
!!
!!
!! @par Revision History
!! Initial version based on mo_ocean_ext_data, Irene Stemmler (2015-10-09)
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
#include "hamocc_omp_definitions.inc"
!----------------------------

MODULE mo_bgc_bcond
USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_hamocc_types,       ONLY: t_hamocc_bcond

  USE mo_var_list,           ONLY: default_var_list_settings,   &
    &                              add_var, add_ref,            &
    &                              new_var_list,                &
    &                              delete_var_list
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_cf_convention,      ONLY: t_cf_var
  USE mo_grib2,              ONLY: t_grib2_var, grib2_var
  USE mo_read_interface,     ONLY: openInputFile, closeFile, on_cells, &
    &                              t_stream_id, nf, read_2D, read_2D_int, &
    &                              read_3D
  USE mo_util_string,        ONLY: t_keyword_list,  &
    &                              associate_keyword, with_keywords
  USE mo_datetime,           ONLY: t_datetime
  USE mo_cdi,                ONLY: DATATYPE_FLT32 
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL,  &
    &                              GRID_REFERENCE,         &
    &                              GRID_CELL,  ZA_SURFACE

  USE mo_hamocc_nml,         ONLY: io_stdo_bgc
  USE mo_ext_data_types,     ONLY: t_external_data, t_external_bgc
  USE mo_ocean_ext_data,     ONLY: ext_data

 IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  PRIVATE

  PUBLIC :: construct_bgc_ext_data
  PUBLIC :: destruct_bgc_ext_data
  PUBLIC :: update_bgc_bcond
  PUBLIC :: ext_data_bgc

  TYPE(t_hamocc_bcond),TARGET :: ext_data_bgc


!-------------------------------------------------------------------------

CONTAINS



  !-------------------------------------------------------------------------
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE construct_bgc_ext_data (p_patch, ext_data,pext_data_bgc)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(t_hamocc_bcond), INTENT(INOUT) :: pext_data_bgc

    INTEGER :: jg
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_bgc_bcond:construct_bgc_data'

    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')


    ! top-level procedure for building data structures for 
    ! external data.
    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data started')


    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'ext_data_bgc_D',jg
      CALL new_ext_data_bgc_list(p_patch(jg), ext_data(jg)%bgc, pext_data_bgc, ext_data(jg)%bgc_list, TRIM(listname))
    END DO ! jg = 1,n_dom


    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data finished')


    !
    CALL read_ext_data_bgc (p_patch, ext_data)

  END SUBROUTINE construct_bgc_ext_data
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !!
!<Optimize:inUse>
  SUBROUTINE new_ext_data_bgc_list ( p_patch, p_ext_bgc, p_ext_data_bgc, p_ext_bgc_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)   :: & !< current patch
      &  p_patch

    TYPE(t_external_bgc), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_bgc 
    TYPE(t_hamocc_bcond), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_data_bgc 

    TYPE(t_var_list) :: p_ext_bgc_list !< current external data list

    CHARACTER(len=*), INTENT(IN)  :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
      &        nblks_e       !< number of edge blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2), shape3d_c(3)
    INTEGER :: idim_omip

    INTEGER :: ibits         !< "entropy" of horizontal slice


    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%alloc_cell_blocks
    nblks_e = p_patch%nblks_e


    ibits = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)

    ! OMIP/NCEP or other flux forcing data on cell centers
    idim_omip = 1 ! or now only dust is read in

    shape3d_c = (/ nproma, 12, nblks_c/)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_bgc_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_bgc_list,            &
                                  & lrestart=.FALSE.,          &
                                 & model_type='oce'  )

    cf_desc    = t_cf_var('Dust cell center', 'kg m-2 yr-1', &
      &                   'DUST', DATATYPE_FLT32)
    grib2_desc = grib2_var( 192, 140, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_bgc_list, 'DUST', p_ext_bgc%dust,      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3d_c )
    CALL add_var( p_ext_bgc_list, 'DUSTY', p_ext_data_bgc%dusty,      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )



  END SUBROUTINE new_ext_data_bgc_list
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Destruct external data data structure and lists
  !!
  !!
  !! @par Revision History
  !!
!<Optimize:inUse>
  SUBROUTINE destruct_bgc_ext_data

    INTEGER :: jg, errstat
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_bgc_bcond:destruct_bgc_ext_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data started')

    DO jg = 1,n_dom
      ! Delete list of constant in time oceanic elements
      CALL delete_var_list( ext_data(jg)%bgc_list )
    ENDDO

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE destruct_bgc_ext_data
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !!
  !! Read HAMOCC external data from netcdf
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE read_ext_data_bgc (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'read_ext_data_bgc'

    CHARACTER(filename_max) :: dust_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, i_lev, no_cells, no_verts, no_tst
    INTEGER :: ncid, dimid
    TYPE(t_stream_id) :: stream_id
    INTEGER :: mpi_comm

    REAL(wp):: z_flux(nproma,12,p_patch(1)%alloc_cell_blocks)
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL message (TRIM(routine), 'start')

    !-------------------------------------------------------------------------
    !  READ DUST
    !-------------------------------------------------------------------------

    jg = 1

    i_lev       = p_patch(jg)%level
    z_flux(:,:,:) = 0.0_wp

    CALL associate_keyword("<path>", TRIM(getModelBaseDir()), keywords)
    
      dust_file='dust.nc'

      CALL message( TRIM(routine),'HAMOCC dust file is: '//TRIM(dust_file) )

      IF(my_process_is_stdio()) THEN
        !
        INQUIRE (FILE=dust_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          write(io_stdo_bgc,*)'FORCING FILE: ',TRIM(dust_file)
          CALL finish(TRIM(routine),'DUST file is not found - ABORT')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(dust_file), NF_NOWRITE, ncid), routine)
        CALL message( TRIM(routine),'HAMOCC dust file opened for read' )

        !
        !
        CALL nf(nf_inq_dimid (ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in HAMOCC dust file do not match - ABORT')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst), routine)
        !
        ! check
        !
        WRITE(message_text,'(A,I6,A)')  'HAMOCC dust file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )
        IF(no_tst /= 12 ) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of forcing timesteps is not equal 12 specified in namelist - ABORT')
        ENDIF

        CALL nf(nf_close(ncid), routine)
      ENDIF

      stream_id = openInputFile(dust_file, p_patch(jg))
      
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF
      no_tst = 12
      !-------------------------------------------------------
      !
      ! Read dust for triangle centers
      !
      !-------------------------------------------------------

        CALL read_3D(stream_id, on_cells, 'DUST', z_flux)
        ext_data(jg)%bgc%dust(:,:,:) = z_flux(:,:,:)
     

      !
      ! close file
      !
      CALL closeFile(stream_id)


      CALL message( TRIM(routine),'HAMOCC dust file read' )

  END SUBROUTINE read_ext_data_bgc
  !--------------------------------------------------
!<Optimize:inUse>

   SUBROUTINE update_bgc_bcond(p_patch_3D, bgc_ext, jstep, datetime)
    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hamocc_bcond)                        :: bgc_ext
    INTEGER, INTENT(IN)                         :: jstep
    TYPE(t_datetime), INTENT(INOUT)             :: datetime

  
 ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_bgc_bcond:update_bgc_bcond'
    CHARACTER(LEN=max_char_length)::check_text 
    INTEGER  :: jmon, jdmon, jmon1, jmon2, ylen, yday
    REAL(wp) :: rday1, rday2
    REAL(wp) ::  z_c2(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)

    TYPE(t_patch), POINTER:: p_patch 
    !TYPE(t_subset_range), POINTER :: all_cells
   !   CALL message( TRIM(routine),'start', io_stdo_bgc )

    !-----------------------------------------------------------------------
    p_patch   => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------


    !  calculate day and month
    jmon  = datetime%month         ! integer current month
    jdmon = datetime%day           ! integer day in month
    yday  = datetime%yeaday        ! integer current day in year
    ylen  = datetime%yealen        ! integer days in year (365 or 366)
    

      jmon1=jmon-1
      jmon2=jmon
      rday1=REAL(15-jdmon,wp)/30.0_wp
      rday2=REAL(15+jdmon,wp)/30.0_wp
      IF (jdmon > 15)  THEN
        jmon1=jmon
        jmon2=jmon+1
        rday1=REAL(45-jdmon,wp)/30.0_wp
        rday2=REAL(jdmon-15,wp)/30.0_wp
      END IF

      IF (jmon1 ==  0) jmon1=12
      IF (jmon1 == 13) jmon1=1
      IF (jmon2 ==  0) jmon2=12
      IF (jmon2 == 13) jmon2=1

   

      bgc_ext%dusty(:,:) = rday1*ext_data(1)%bgc%dust(:,jmon1,:) + &
      &                                   rday2*ext_data(1)%bgc%dust(:,jmon2,:)
    !  bgc_ext%dusty(:,:) = ext_data(1)%bgc%dust(:,jmon,:) 

    !  CALL message( TRIM(routine),'end', io_stdo_bgc )
    END SUBROUTINE update_bgc_bcond
END MODULE

