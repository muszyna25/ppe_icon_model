MODULE mo_io_output
  !
  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message, message_text, get_filename_noext
  USE mo_var_metadata,  ONLY: t_var_metadata
  USE mo_linked_list    ! we need all
  USE mo_var_list,      ONLY: t_var_list, nvar_lists, var_lists
  USE mo_cdi_constants
  USE mo_util_string,   ONLY: separator
  USE mo_util_sysinfo,  ONLY: util_user_name, util_os_system, util_node_name 
  USE mo_io_distribute, ONLY: gather_cells, gather_edges, gather_vertices
  USE mo_mpi,           ONLY: my_process_is_stdio, p_pe, p_barrier
  USE mo_model_domain,  ONLY: t_patch, p_patch

  USE mo_impl_constants,        ONLY: ihs_ocean, zml_soil
  USE mo_vertical_coord_table,  ONLY: vct
  USE mo_dynamics_config,       ONLY: iequations
  USE mo_io_config,             ONLY: out_expname, lwrite_pzlev
  USE mo_grid_config,           ONLY: global_cell_type
  USE mo_run_config,            ONLY: num_lev, num_levp1
  USE mo_io_units,              ONLY: filename_max
  USE mo_nh_pzlev_config,       ONLY: nh_pzlev_config
  USE mo_lnd_nwp_config,        ONLY: nlev_snow

  USE mo_ocean_nml,             ONLY: n_zlev
  USE mo_oce_state,             ONLY: set_zlev

  USE mo_io_vlist,              ONLY: addGlobAtts, addAtmAtts, addOceAtts

  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_output_time
  PUBLIC :: init_output
  PUBLIC :: open_output_files
  PUBLIC :: write_output
  PUBLIC :: close_output_files
  PUBLIC :: finish_output
  !
  CHARACTER(len=32) :: private_initial_time = '' 
  CHARACTER(len=32) :: private_output_time = '' 
  !
  !
  !------------------------------------------------------------------------------------------------
CONTAINS
  !------------------------------------------------------------------------------------------------
  !
  ! YYYYMMDDThhmmssZ (T is a separator and Z means UTC as timezone)
  !
  SUBROUTINE set_output_time(iso8601)
    CHARACTER(len=*), INTENT(in) :: iso8601
    !IF (private_initial_time == '') THEN
    !  private_initial_time = iso8601
    !ELSE
      private_output_time = iso8601
    !ENDIF
  END SUBROUTINE set_output_time
  !------------------------------------------------------------------------------------------------
  SUBROUTINE init_output

    INTEGER :: i

    ! This routine currently sets the output flag to TRUE for every var_list.
    ! This is for testing only and has to be made dependend on the different output flags
    ! like in mo_io_vlist - this is experimental code at the moment!

    DO i = 1, nvar_lists

      var_lists(i)%p%loutput = .TRUE.

      ! Some elements of the ext_data_atm_td_.. var_lists have the wrong shape
      ! (/ nproma, nblks_c, ntimes /), so we must omit these lists since otherwise
      ! the program crashes during communication:
      IF(var_lists(i)%p%name(1:15) == 'ext_data_atm_td') var_lists(i)%p%loutput = .FALSE.

      ! Set filetype to NetCDF
      var_lists(i)%p%output_type = FILETYPE_NC2
    ENDDO

  END SUBROUTINE init_output
  !------------------------------------------------------------------------------------------------
  !
  ! Loop over all the output streams and open the associated files. Set
  ! unit numbers (file IDs) for all streams associated with a file.
  !
  SUBROUTINE open_output_files(jfile)

    INTEGER, INTENT(IN) :: jfile ! Number of file set to open

    CHARACTER(len=filename_max) :: output_filename
    CHARACTER(LEN=filename_max) :: grid_filename
    INTEGER :: i ,j, k, k_jg, nlev, nlevp1, nplev, nzlev, nzlevp1, znlev_soil, astatus
    REAL(wp), ALLOCATABLE :: levels(:), levels_i(:), levels_m(:)
    !
    ! first check for output time
    !
!   RJ: Currently private_initial_time not set (and not needed)
!    IF (private_initial_time == '') THEN
!      CALL finish('open_output_files','initial time not set')
!    ENDIF
    !
    ! print header for output var_lists
    !
    CALL message('',separator)
    CALL message('','')
    CALL message('','Open output files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t50,a,t84,a,t94,a)') 'file', 'var list', 'file ID', 'output'
    CALL message('',message_text)
    CALL message('','')
    !
    ! Loop over all var_lists and open the associated files. Set
    ! file IDs if necessary.
    !
    DO i = 1, nvar_lists
      var_lists(i)%p%first = .FALSE.
      var_lists(i)%p%output_opened = .FALSE.
    ENDDO
    !
    DO i = 1, nvar_lists
      !
      ! skip, if file is already opened
      !
      IF (var_lists(i)%p%output_opened) CYCLE
      !
      ! skip, if var_list is not required for outputing
      !
      IF (.NOT. var_lists(i)%p%loutput) CYCLE
      !
      ! check output file type
      !
      SELECT CASE (var_lists(i)%p%output_type)       
      CASE (FILETYPE_NC)
        CALL finish('open_output_files','netCDF classic not supported')
      CASE (FILETYPE_NC2, FILETYPE_NC4)
        ! this is ok, both formats can write more than 2GB files
      CASE (FILETYPE_GRB)
        CALL finish('open_output_files','GRIB1 not supported')
      CASE (FILETYPE_GRB2)
        CALL message('open_output_files','GRIB2 support experimental')
      CASE default
        CALL finish('open_output_files','unknown output_type')
      END SELECT
      !
      var_lists(i)%p%first = .TRUE.
      !
      k_jg = var_lists(i)%p%patch_id
      IF(k_jg<=0) call finish('open_output_files','patch_id <= 0')

      grid_filename = get_filename_noext(p_patch(k_jg)%grid_filename)
      ! Raw data file name(s) for output
      !
      WRITE (output_filename,'(a,a,a,a,i4.4,a)')  &
        &  TRIM(out_expname),"_",TRIM(grid_filename),'_', jfile, '.nc'
      !output_filename = TRIM(private_base_filename)      &
      !     &           //'_'//TRIM(var_lists(i)%p%model_type)//'.nc'

      IF (my_process_is_stdio()) THEN
        SELECT CASE (var_lists(i)%p%output_type)
        CASE (FILETYPE_NC2)
          var_lists(i)%p%cdiFileID_output = streamOpenWrite(output_filename, FILETYPE_NC2)
        CASE (FILETYPE_NC4)
          var_lists(i)%p%cdiFileID_output = streamOpenWrite(output_filename, FILETYPE_NC4)
        CASE (FILETYPE_GRB2)
          var_lists(i)%p%cdiFileID_output = streamOpenWrite(output_filename, FILETYPE_GRB2)
        END SELECT
        !
        var_lists(i)%p%filename = TRIM(output_filename)
        !
        IF (var_lists(i)%p%cdiFileID_output < 0) THEN
          WRITE(message_text,'(a)') cdiStringError(var_lists(i)%p%cdiFileID_output)
          CALL message('',message_text)
          CALL finish ('open_output_files', 'open failed on '//TRIM(output_filename))
        ELSE
          CALL message ('open_output_files', 'opened '//TRIM(output_filename))
          var_lists(i)%p%output_opened = .TRUE.
        END IF
        !
        ! The following sections add the file global properties collected in init_output
        !
        ! 1. create cdi vlist 
        !
        var_lists(i)%p%cdiVlistID = vlistCreate()
        !
        !    set cdi internal time index to 0 for writing time slices in netCDF
        !
        var_lists(i)%p%cdiTimeIndex = 0
        !
        ! 2. add global attributes for netCDF
        !
        ! 3. add horizontal grid descriptions

        ! Cells

        var_lists(i)%p%cdiCellGridID = gridCreate(GRID_UNSTRUCTURED, p_patch(k_jg)%n_patch_cells_g)
        CALL gridDefNvertex(var_lists(i)%p%cdiCellGridID, global_cell_type)
        !
        CALL gridDefXname(var_lists(i)%p%cdiCellGridID, 'clon')
        CALL gridDefXlongname(var_lists(i)%p%cdiCellGridID, 'center longitude')
        CALL gridDefXunits(var_lists(i)%p%cdiCellGridID, 'radians')
        !
        CALL gridDefYname(var_lists(i)%p%cdiCellGridID, 'clat')
        CALL gridDefYlongname(var_lists(i)%p%cdiCellGridID, 'center latitude')
        CALL gridDefYunits(var_lists(i)%p%cdiCellGridID, 'radians')

        ! Verts

        var_lists(i)%p%cdiVertGridID = gridCreate(GRID_UNSTRUCTURED, p_patch(k_jg)%n_patch_verts_g)
        CALL gridDefNvertex(var_lists(i)%p%cdiVertGridID, 9-global_cell_type)
        !
        CALL gridDefXname(var_lists(i)%p%cdiVertGridID, 'vlon')
        CALL gridDefXlongname(var_lists(i)%p%cdiVertGridID, 'vertex longitude')
        CALL gridDefXunits(var_lists(i)%p%cdiVertGridID, 'radians')
        !
        CALL gridDefYname(var_lists(i)%p%cdiVertGridID, 'vlat')
        CALL gridDefYlongname(var_lists(i)%p%cdiVertGridID, 'vertex latitude')
        CALL gridDefYunits(var_lists(i)%p%cdiVertGridID, 'radians')

        ! Edges

        var_lists(i)%p%cdiEdgeGridID = gridCreate(GRID_UNSTRUCTURED, p_patch(k_jg)%n_patch_edges_g)
        CALL gridDefNvertex(var_lists(i)%p%cdiEdgeGridID, 4)
        !
        CALL gridDefXname(var_lists(i)%p%cdiEdgeGridID, 'elon')
        CALL gridDefXlongname(var_lists(i)%p%cdiEdgeGridID, 'edge longitude')
        CALL gridDefXunits(var_lists(i)%p%cdiEdgeGridID, 'radians')
        !
        CALL gridDefYname(var_lists(i)%p%cdiEdgeGridID, 'elat')
        CALL gridDefYlongname(var_lists(i)%p%cdiEdgeGridID, 'edge latitude')
        CALL gridDefYunits(var_lists(i)%p%cdiEdgeGridID, 'radians')

        !
        ! 4. add vertical grid descriptions
        !    RJ: This is copied from mo_io_vlist

        ! surface level
        var_lists(i)%p%cdiZaxisID(ZA_surface) = zaxisCreate(ZAXIS_SURFACE, 1)
        ALLOCATE(levels(1))
        levels(1) = 0.0_wp
        CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_surface), levels)
        DEALLOCATE(levels)

        ! atm (pressure) height, ocean depth
        IF (iequations/=ihs_ocean) THEN ! atm 

          nlev   = num_lev(k_jg)
          nlevp1 = num_levp1(k_jg)
          ! introduce temporary variable znlev_soil, since global variable nlev_soil 
          ! is unknown to the I/O-Processor. Otherwise receive_patch_configuration in 
          ! mo_io_async complains about mismatch of levels. 
          znlev_soil = SIZE(zml_soil)-1

          ! Hybrid

          var_lists(i)%p%cdiZaxisID(ZA_hybrid)      = zaxisCreate(ZAXIS_HYBRID, nlev)
          var_lists(i)%p%cdiZaxisID(ZA_hybrid_half) = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1)

          ALLOCATE(levels(nlev))
          DO k = 1, nlev
            levels(k) = REAL(k,wp)
          END DO
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_hybrid), levels)
          DEALLOCATE(levels)
          CALL zaxisDefVct(var_lists(i)%p%cdiZaxisID(ZA_hybrid), 2*nlevp1, vct(1:2*nlevp1))
          !
          ALLOCATE(levels(nlevp1))
          DO k = 1, nlevp1
            levels(k) = REAL(k,wp)
          END DO
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_hybrid_half), levels)
          DEALLOCATE(levels)
          CALL zaxisDefVct(var_lists(i)%p%cdiZaxisID(ZA_hybrid_half), 2*nlevp1, vct(1:2*nlevp1))

          ! Define axes for soil model

          var_lists(i)%p%cdiZaxisID(ZA_depth_below_land_p1) = &
            & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+2)
          ALLOCATE(levels(znlev_soil+2))
          levels(1) = 0._wp
          DO k = 1, znlev_soil+1
            levels(k+1) = zml_soil(k)*100._wp
          END DO
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_depth_below_land_p1), levels)
          DEALLOCATE(levels)

          var_lists(i)%p%cdiZaxisID(ZA_depth_below_land) = &
            & zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_depth_below_land), zml_soil*100._wp)
          !
          var_lists(i)%p%cdiZaxisID(ZA_generic_snow_p1) = zaxisCreate(ZAXIS_GENERIC, nlev_snow+1)
          ALLOCATE(levels(nlev_snow+1))
          DO k = 1, nlev_snow+1
            levels(k) = REAL(k,wp)
          END DO
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_generic_snow_p1), levels)
          DEALLOCATE(levels)
          !
          var_lists(i)%p%cdiZaxisID(ZA_generic_snow) = zaxisCreate(ZAXIS_GENERIC, nlev_snow)
          ALLOCATE(levels(nlev_snow))
          DO k = 1, nlev_snow
            levels(k) = REAL(k,wp)
          END DO
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_generic_snow), levels)
          DEALLOCATE(levels)

          ! Define axes for output on p- and z-levels
          !
          IF (lwrite_pzlev) THEN
            nplev = nh_pzlev_config(k_jg)%nplev
            var_lists(i)%p%cdiZaxisID(ZA_pressure) = zaxisCreate(ZAXIS_PRESSURE, nplev)
            ALLOCATE(levels(nplev))
            DO k = 1, nplev
              levels(k) = nh_pzlev_config(k_jg)%plevels(k)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_pressure), levels)
            CALL zaxisDefVct(var_lists(i)%p%cdiZaxisID(ZA_pressure), nplev, levels)
            DEALLOCATE(levels)

            nzlev = nh_pzlev_config(k_jg)%nzlev
            var_lists(i)%p%cdiZaxisID(ZA_height)  = zaxisCreate(ZAXIS_HEIGHT, nzlev)
            ALLOCATE(levels(nzlev))
            DO k = 1, nzlev
              levels(k) = nh_pzlev_config(k_jg)%zlevels(k)
            END DO
            CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_height), levels)
            CALL zaxisDefVct(var_lists(i)%p%cdiZaxisID(ZA_height), nzlev, levels)
            DEALLOCATE(levels)
          ENDIF

        ELSE ! oce
          var_lists(i)%p%cdiZaxisID(ZA_depth)      = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, n_zlev)
          nzlevp1 = n_zlev + 1
          var_lists(i)%p%cdiZaxisID(ZA_depth_half) = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, nzlevp1)

          ALLOCATE(levels_i(nzlevp1))
          ALLOCATE(levels_m(n_zlev))
          CALL set_zlev(levels_i, levels_m)
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_depth)     , levels_m)
          CALL zaxisDefLevels(var_lists(i)%p%cdiZaxisID(ZA_depth_half), levels_i)
          DEALLOCATE(levels_i)
          DEALLOCATE(levels_m)
          var_lists(i)%p%cdiZaxisID(ZA_generic_ice) = zaxisCreate(ZAXIS_GENERIC, 1)
        ENDIF

        !
        ! 5. output does contain absolute time 
        !
        var_lists(i)%p%cdiTaxisID = taxisCreate(TAXIS_ABSOLUTE)
        CALL vlistDefTaxis(var_lists(i)%p%cdiVlistID, var_lists(i)%p%cdiTaxisID)

        !
        ! 6. global attributes
        !
        CALL addGlobAtts(var_lists(i)%p%cdiVlistID,k_jg,astatus)
        IF (iequations/=ihs_ocean) THEN
          CALL addAtmAtts(var_lists(i)%p%cdiVlistID,k_jg,astatus)
        ELSE
          CALL addOceAtts(var_lists(i)%p%cdiVlistID,astatus)
        END IF

      ENDIF
      !
      ! add variables
      !
      IF (my_process_is_stdio()) THEN
        !
        CALL addVarListToVlist(var_lists(i), var_lists(i)%p%cdiVlistID)
        !
        WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
             output_filename, var_lists(i)%p%name,             &
             var_lists(i)%p%cdiFileID_output, var_lists(i)%p%loutput
        CALL message('',message_text)
      ENDIF
      !
      ! loop over all other output var_lists eventually corresponding to the same file
      !
      DO j = 1, nvar_lists
        !
        IF (var_lists(j)%p%output_opened) CYCLE
        IF (.NOT. var_lists(j)%p%loutput) CYCLE
        !
        IF (var_lists(j)%p%output_type /= var_lists(i)%p%output_type &
             .AND. var_lists(j)%p%cdiFileID_output == var_lists(i)%p%cdiFileID_output) THEN
          CALL finish('open_output_streams', 'different file types for the same output file')
        ENDIF
        !
        IF (var_lists(i)%p%model_type == var_lists(j)%p%model_type .AND. &
          & var_lists(i)%p%patch_id   == var_lists(j)%p%patch_id   ) THEN 
          var_lists(j)%p%output_opened = .TRUE.
          var_lists(j)%p%filename = var_lists(i)%p%filename          
          !
          ! set file IDs of all associated output files
          !
          var_lists(j)%p%cdiFileID_output     = var_lists(i)%p%cdiFileID_output
          var_lists(j)%p%cdiVlistID           = var_lists(i)%p%cdiVlistID
          var_lists(j)%p%cdiCellGridID        = var_lists(i)%p%cdiCellGridID
          var_lists(j)%p%cdiVertGridID        = var_lists(i)%p%cdiVertGridID
          var_lists(j)%p%cdiEdgeGridID        = var_lists(i)%p%cdiEdgeGridID
          var_lists(j)%p%cdiZaxisID(:)        = var_lists(i)%p%cdiZaxisID(:)
          var_lists(j)%p%cdiTaxisID           = var_lists(i)%p%cdiTaxisID
          !
          ! add variables to already existing cdi vlists
          !
          IF (my_process_is_stdio()) THEN
            !
            CALL addVarListToVlist(var_lists(j), var_lists(j)%p%cdiVlistID)
            !
            WRITE(message_text,'(t1,a49,t50,a31,t84,i6,t94,l5)')        &
                 output_filename, var_lists(j)%p%name,             &
                 var_lists(j)%p%cdiFileID_output, var_lists(j)%p%loutput
            CALL message('',message_text)
          ENDIF
        ENDIF
      ENDDO
      !      
      IF (my_process_is_stdio() .AND. var_lists(i)%p%first) THEN
        CALL streamDefVlist(var_lists(i)%p%cdiFileID_output, var_lists(i)%p%cdiVlistID)
      ENDIF
      !
    END DO
    !    
    CALL message('','')
    !
  END SUBROUTINE open_output_files
  !------------------------------------------------------------------------------------------------
  !
  ! define variables and attributes
  !
  SUBROUTINE addVarListToVlist(this_list, vlistID)
    TYPE (t_var_list), INTENT(inout) :: this_list
    INTEGER,           INTENT(inout) :: vlistID
    !      
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    INTEGER :: varID, gridID, zaxisID, k_jg, nlev, nlevp1, znlev_soil
    !
    REAL(wp) :: casted_missval
    !
    element => start_with
    element%next_list_element => this_list%p%first_list_element
print *,'addVarListToVlist: ',TRIM(this_list%p%name)

    k_jg   = this_list%p%patch_id
    nlev   = num_lev(k_jg)
    nlevp1 = num_levp1(k_jg)

    ! See above ...
    znlev_soil = SIZE(zml_soil)-1
    !
    for_all_list_elements: DO
      !
      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT
      !
      ! retrieve information from actual linked list element
      !
      info => element%field%info
      !
      ! skip this field ?
      !
      IF (.NOT. info%loutput) CYCLE
      !
      ! set grid ID
      !
      SELECT CASE (info%hgrid)
      CASE(GRID_UNSTRUCTURED_CELL)
        info%cdiGridID = this_list%p%cdiCellGridID
        gridID = info%cdiGridID
      CASE(GRID_UNSTRUCTURED_VERT)
        info%cdiGridID = this_list%p%cdiVertGridID
        gridID = info%cdiGridID
      CASE(GRID_UNSTRUCTURED_EDGE)
        info%cdiGridID = this_list%p%cdiEdgeGridID
        gridID = info%cdiGridID
      CASE DEFAULT
        CALL finish('addVarListToVlist', 'GRID definition missing for '//TRIM(info%name))
      END SELECT

      !
      ! set z axis ID
      !
      SELECT CASE (info%vgrid)

      CASE (ZAXIS_SURFACE)
        info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_surface)

      CASE (ZAXIS_HYBRID)
        info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_hybrid)

      CASE (ZAXIS_HYBRID_HALF)
        info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_hybrid_half)

      CASE (ZAXIS_DEPTH_BELOW_SEA)
        ! RJ: Not sure about this ...
        IF (info%used_dimensions(2) == n_zlev) THEN
          info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_depth)
        ELSE IF (info%used_dimensions(2) == n_zlev+1) THEN
          info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_depth_half)
        ELSE
          PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_DEPTH_BELOW_SEA')
        ENDIF

      CASE (ZAXIS_DEPTH_BELOW_LAND)
        IF (info%used_dimensions(2) == znlev_soil+1) THEN
          info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_depth_below_land)
        ELSE IF (info%used_dimensions(2) == znlev_soil+2) THEN
          info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_depth_below_land_p1)
        ELSE
          PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
          CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_DEPTH_BELOW_LAND')
        ENDIF

      CASE (ZAXIS_HEIGHT)
        IF(info%name(1:8)=='lnd_prog' .AND. INDEX(info%name,'snow') /= 0) THEN
          ! This is a special case - use ZA_generic_snow[_p1]
          IF(info%used_dimensions(2) == nlev_snow) THEN
            info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_generic_snow)
          ELSE IF(info%used_dimensions(2) == nlev_snow+1) THEN
            info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_generic_snow_p1)
          ELSE
            PRINT *,'SNOW variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
            CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_HEIGHT')
          ENDIF
        ELSE
          ! In all other cases, ZAXIS_HEIGHT seems to be equivalent to ZAXIS_HYBRID/ZAXIS_HYBRID_HALF
          ! TODO: Is there a difference to ZAXIS_HYBRID/ZAXIS_HYBRID_HALF ???
          IF (info%used_dimensions(2) == nlevp1) THEN
            info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_hybrid_half)
          ELSE IF (info%used_dimensions(2) == nlev) THEN
            info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_hybrid)
          ELSE
            PRINT *,'Variable: ',TRIM(info%name),' Dimension 2: ',info%used_dimensions(2)
            CALL finish('addVarListToVlist','Dimension mismatch for ZAXIS_HEIGHT')
          ENDIF
        ENDIF

      CASE (ZAXIS_PRESSURE)
        info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_pressure)

      CASE (ZAXIS_ALTITUDE)
        info%cdiZaxisID =  this_list%p%cdiZaxisID(ZA_height)

      CASE DEFAULT
        PRINT *,'Variable: ',TRIM(info%name),' ZAXIS: ',info%vgrid
        CALL finish('addVarListToVlist', 'ZAXIS definition missing for '//TRIM(info%name))

      END SELECT

      zaxisID = info%cdiZaxisID
      !
      info%cdiVarID = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE)
      varID = info%cdiVarID 
      !
      CALL vlistDefVarDatatype(vlistID, varID, DATATYPE_FLT64)
print *,'vlistDefVarName: ',info%name,this_list%p%patch_id
      CALL vlistDefVarName(vlistID, varID, info%name)
      !
      IF (info%cf%long_name /= '') CALL vlistDefVarLongname(vlistID, varID, info%cf%long_name)
      IF (info%cf%units /= '') CALL vlistDefVarUnits(vlistID, varID, info%cf%units)

      IF (info%lmiss) THEN
        IF (ASSOCIATED(element%field%r_ptr)) THEN
          casted_missval = info%missval%rval
        ELSE IF (ASSOCIATED(element%field%i_ptr)) THEN
          casted_missval = REAL(info%missval%ival,wp)
        ELSE
          IF (info%missval%lval) THEN
            casted_missval = 1.0_wp
          ELSE
            casted_missval = 0.0_wp
          ENDIF
        ENDIF
        CALL vlistDefVarMissval(vlistID, varID, casted_missval)
      ENDIF
      !
    ENDDO for_all_list_elements
    !
  END SUBROUTINE addVarListToVlist
  !------------------------------------------------------------------------------------------------
  !
  SUBROUTINE close_output_files
    !
    ! Loop over all the output streams and close the associated files, set
    ! opened to false
    !
    INTEGER :: i, j, fileID, vlistID
    !
    CALL message('',separator)
    CALL message('','')
    CALL message('','Close output files:')
    CALL message('','')
    WRITE(message_text,'(t1,a,t84,a)') 'file', 'file ID'
    CALL message('',message_text)
    CALL message('','')
    !
    close_all_lists: DO i = 1, nvar_lists
      !
      IF (var_lists(i)%p%output_opened) THEN
        IF (my_process_is_stdio() .AND. var_lists(i)%p%first) THEN
          !
          fileID = var_lists(i)%p%cdiFileID_output
          !
          IF (fileID /= CDI_UNDEFID) THEN
            CALL streamClose(var_lists(i)%p%cdiFileID_output)
            DO j = 1, nvar_lists
              IF (fileID == var_lists(j)%p%cdiFileID_output) THEN
                var_lists(j)%p%cdiFileID_output = CDI_UNDEFID
              ENDIF
            ENDDO
          ENDIF
          !
          WRITE(message_text,'(t1,a,t84,i6)') &
               TRIM(var_lists(i)%p%filename), fileID
          CALL message('',message_text)
          !
        ENDIF
        var_lists(i)%p%filename   = ''
      ENDIF
    ENDDO close_all_lists
    !
    IF(my_process_is_stdio()) THEN ! not really necessary but prevents a compiler bug (???)
    for_all_vlists: DO i = 1, nvar_lists
      vlistID = var_lists(i)%p%cdiVlistID
      IF (vlistID /= CDI_UNDEFID) THEN
        CALL vlistDestroy(var_lists(i)%p%cdiVlistID)
        DO j = 1, nvar_lists
          IF (vlistID == var_lists(j)%p%cdiVlistID) THEN
            var_lists(j)%p%cdiVlistID = CDI_UNDEFID
          ENDIF
        ENDDO
      ENDIF
    ENDDO for_all_vlists
    ENDIF
    CALL message('','')
    !
    ! reset all var list properties related to cdi files
    !
    reset_all_lists: DO i = 1, nvar_lists
      var_lists(i)%p%output_opened = .FALSE.
      var_lists(i)%p%first  = .FALSE.
    ENDDO reset_all_lists
    !
    private_output_time = ''
    !
  END SUBROUTINE close_output_files
  !------------------------------------------------------------------------------------------------
  !
  ! loop over all var_lists for output
  !
  SUBROUTINE write_output

    INTEGER :: i,j
    LOGICAL :: write_info
    !
    write_info   = .TRUE.
    !
    ! pick up first stream associated with each file
    !
    DO i = 1, nvar_lists

      IF (.NOT. var_lists(i)%p%loutput) CYCLE

      IF (var_lists(i)%p%first) THEN
        IF (write_info) THEN
          SELECT CASE (var_lists(i)%p%output_type)
          CASE (FILETYPE_NC2)
            CALL message('','Write netCDF2 output for : '//TRIM(private_output_time))
          CASE (FILETYPE_NC4)
            IF (var_lists(i)%p%compression_type == COMPRESS_ZIP) THEN
              CALL message('', &
                   'Write compressed netCDF4 output for : '//TRIM(private_output_time))
            ELSE
              CALL message('','Write netCDF4 output for : '//TRIM(private_output_time))
            END IF
          CASE (FILETYPE_GRB2)
            IF (var_lists(i)%p%compression_type == COMPRESS_SZIP) THEN
              CALL message('','Write compressed GRIB2 output for : '//TRIM(private_output_time))
            ELSE
              CALL message('','Write GRIB2 output for : '//TRIM(private_output_time))
            ENDIF
          END SELECT
        ENDIF
        write_info = .FALSE.
        !
        ! write time information to netCDF file
        !
        IF (my_process_is_stdio()) THEN
          CALL write_time_to_output(var_lists(i))
        ENDIF
        !
        ! loop over all streams associated with the file
        !
        DO j = i, nvar_lists

          IF (.NOT. var_lists(j)%p%loutput) CYCLE

          IF (var_lists(i)%p%model_type == var_lists(j)%p%model_type .AND. &
            & var_lists(i)%p%patch_id   == var_lists(j)%p%patch_id   ) THEN 
            !
            !
            ! write variables
            !
            CALL write_output_var_list(var_lists(j))
            !
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !
  END SUBROUTINE write_output
  !------------------------------------------------------------------------------------------------
  !
  ! set time for output in cdi format
  !
  SUBROUTINE write_time_to_output (this_list)
    TYPE (t_var_list), INTENT(inout) :: this_list
    !
    INTEGER :: fileID, idate, itime, iret
    !
    fileID = this_list%p%cdiFileID_output
    !
    CALL get_date_components(private_output_time, idate, itime)
    !
    CALL taxisDefVdate(this_list%p%cdiTaxisID, idate)
    CALL taxisDefVtime(this_list%p%cdiTaxisID, itime)
    !
    iret = streamDefTimestep(fileID, this_list%p%cdiTimeIndex)
    this_list%p%cdiTimeIndex = this_list%p%cdiTimeIndex + 1
    ! 
  CONTAINS
    !
    SUBROUTINE get_date_components(iso8601, idate, itime)
      CHARACTER(len=*), INTENT(in)  :: iso8601
      INTEGER,          INTENT(out) :: idate, itime
      !
      INTEGER :: it, iz 
      !
      it = INDEX(iso8601, 'T')
      iz = INDEX(iso8601, 'Z')
      READ(iso8601(1:it-1), '(i10)') idate 
      READ(iso8601(it+1:iz-1), '(i10)') itime 
      !
    END SUBROUTINE get_date_components
    !
  END SUBROUTINE write_time_to_output
  !------------------------------------------------------------------------------------------------
  !
  ! write variables of a list for output
  !
  SUBROUTINE write_output_var_list(this_list)
    TYPE (t_var_list) ,INTENT(in) :: this_list
    !
    INTEGER           :: gridtype
    !
    ! variables of derived type used in linked list
    !
    TYPE (t_var_metadata), POINTER :: info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with
    !
    REAL(wp), POINTER :: rptr2d(:,:)   ! 2d field distributed over processors
    REAL(wp), POINTER :: rptr3d(:,:,:) ! 3d field distributed over processors
    !
    REAL(wp), POINTER :: r5d(:,:,:,:,:) ! field gathered on I/O processor
    !
    INTEGER :: gdims(5), nindex, k_jg
    !
    ! Loop over all fields in linked list
    !
    element => start_with
    element%next_list_element => this_list%p%first_list_element    

    k_jg = this_list%p%patch_id
    IF(k_jg<=0) call finish('write_output_var_list','patch_id <= 0') ! Safety only
    !
    for_all_list_elements: DO
      !
      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT
      !
      rptr2d => NULL()
      rptr3d => NULL()
      !
      ! retrieve information from actual linked list element
      !
      info => element%field%info 
      !
      ! skip this field ?
      !
      IF (.NOT. info%loutput) CYCLE
      !
      IF (info%lcontained) THEN 
        nindex = info%ncontained
      ELSE
        nindex = 1
      ENDIF
      !
      SELECT CASE (info%ndims)
      CASE (1)
        CALL finish('write_output_var_list','1d arrays not handled yet.')
      CASE (2)
        rptr2d => element%field%r_ptr(:,:,nindex,1,1)
      CASE (3)
        rptr3d => element%field%r_ptr(:,:,:,nindex,1)
      CASE (4)
        CALL finish('write_output_var_list','4d arrays not handled yet.')
      CASE (5)
        CALL finish('write_output_var_list','5d arrays not handled yet.')
      CASE DEFAULT 
        CALL finish('write_output_var_list','dimension not set.')        
      END SELECT
      !
      gridtype = info%hgrid
      !
      ! allocate temporary global array on output processor
      ! and gather field from other processors
      !
      r5d => NULL()
      !
      SELECT CASE (gridtype)
      CASE (GRID_UNSTRUCTURED_CELL)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ p_patch(k_jg)%n_patch_cells_g, 1, 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_cells(rptr2d, r5d, p_patch=p_patch(k_jg))
        ELSE
          gdims(:) = (/ p_patch(k_jg)%n_patch_cells_g, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_cells(rptr3d, r5d, p_patch=p_patch(k_jg))
        ENDIF
      CASE (GRID_UNSTRUCTURED_VERT)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ p_patch(k_jg)%n_patch_verts_g, 1, 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_vertices(rptr2d, r5d, p_patch=p_patch(k_jg))
        ELSE
          gdims(:) = (/ p_patch(k_jg)%n_patch_verts_g, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_vertices(rptr3d, r5d, p_patch=p_patch(k_jg))
        ENDIF
      CASE (GRID_UNSTRUCTURED_EDGE)
        IF (info%ndims == 2) THEN
          gdims(:) = (/ p_patch(k_jg)%n_patch_edges_g, 1, 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_edges(rptr2d, r5d, p_patch=p_patch(k_jg))
        ELSE
          gdims(:) = (/ p_patch(k_jg)%n_patch_edges_g, info%used_dimensions(2), 1, 1, 1 /)
          ALLOCATE(r5d(gdims(1),gdims(2),gdims(3),gdims(4),gdims(5)) )
          CALL gather_edges(rptr3d, r5d, p_patch=p_patch(k_jg))
        ENDIF
      CASE default
        CALL finish('out_stream','unknown grid type')
      END SELECT
      !
      ! write data
      !
      IF (my_process_is_stdio()) THEN
        CALL write_var (this_list, info, r5d)
      END IF
      !
      ! deallocate temporary global arrays
      !
      IF (ASSOCIATED (r5d)) DEALLOCATE (r5d)
      !
    END DO for_all_list_elements
    !
  END SUBROUTINE write_output_var_list
  !------------------------------------------------------------------------------------------------
  !
  ! finally write data ...
  !
  SUBROUTINE write_var (this_list, info, array)
    TYPE (t_var_list),     INTENT(in) :: this_list
    TYPE (t_var_metadata), INTENT(in) :: info             ! field description
    REAL(wp),              INTENT(in) :: array(:,:,:,:,:) ! output field
    !
    INTEGER :: fileID                       ! File ID
    INTEGER :: varID                        ! Variable ID
    INTEGER :: nmiss = 0
    !
    fileID  = this_list%p%cdiFileID_output
    varID   = info%cdiVarID
    !
    CALL streamWriteVar(fileID, varID, array, nmiss)
    !
  END SUBROUTINE write_var
  !------------------------------------------------------------------------------------------------
  !
  ! deallocate module variables
  !
  SUBROUTINE finish_output

    INTEGER :: i, j

    for_all_var_lists: DO i = 1, nvar_lists    
      IF (var_lists(i)%p%cdiFileID_output >= 0) THEN
        CALL gridDestroy(var_lists(i)%p%cdiCellGridID)
        CALL gridDestroy(var_lists(i)%p%cdiVertGridID)
        CALL gridDestroy(var_lists(i)%p%cdiEdgeGridID)
        DO j = 1, SIZE(var_lists(i)%p%cdiZaxisID)
          IF (var_lists(i)%p%cdiZaxisID(j) /= CDI_UNDEFID) &
            CALL zaxisDestroy(var_lists(i)%p%cdiZaxisID(j))
        ENDDO
        var_lists(i)%p%cdiFileId_output     = CDI_UNDEFID
        var_lists(i)%p%cdiVlistId           = CDI_UNDEFID
        var_lists(i)%p%cdiCellGridID        = CDI_UNDEFID
        var_lists(i)%p%cdiVertGridID        = CDI_UNDEFID
        var_lists(i)%p%cdiEdgeGridID        = CDI_UNDEFID
        var_lists(i)%p%cdiZaxisID(:)        = CDI_UNDEFID
        var_lists(i)%p%cdiTaxisID           = CDI_UNDEFID
        var_lists(i)%p%cdiTimeIndex         = CDI_UNDEFID
      ENDIF
    ENDDO for_all_var_lists
    !
  END SUBROUTINE finish_output
  !------------------------------------------------------------------------------------------------
END MODULE mo_io_output
