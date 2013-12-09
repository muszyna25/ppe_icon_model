! -*- f90 -*-
  ! -------------------------------------------------------------------
  !  CHANNEL INTERFACE FOR MBM OUTPUT
  ! -------------------------------------------------------------------

  !PRIVATE SUBROUTINE main_channel_initialize_gatts
  !PRIVATE SUBROUTINE main_channel_initialize_dims
  !PRIVATE SUBROUTINE main_channel_initialize_reprs

  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_gatts

    USE messy_main_data_bi,       ONLY: basemodver=>modver, basemodstr=>modstr

    USE messy_main_channel,       ONLY: new_attribute, write_attribute &
                                      , AF_RST_CMP, AF_RST_INP, EXP_NAME &
                                      , modstr, modver

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_gatts'
    INTEGER :: status
    CHARACTER (8)   :: ydate
    CHARACTER (10)  :: ytime

    ! example 1: system date and time
    CALL DATE_AND_TIME(ydate, ytime)
    CALL new_attribute(status, 'operating_date_time' &
         , c = ydate(1:8)//' '//ytime(1:6) )
    CALL channel_halt(substr, status)

    ! example 2: experiment name from namelist
    CALL new_attribute(status, 'MESSy_experiment', c=TRIM(EXP_NAME))
    CALL channel_halt(substr, status)

    ! example 3: basemodel version information
    CALL new_attribute(status, 'MESSy_basemodel' &
         , c = basemodstr//' version '//basemodver//&
         &', MPI-Met. Hamburg and Deutscher Wetterdienst Offenbach' )
    CALL channel_halt(substr, status)

    ! example 4: special flag to be checked from restart file
    CALL new_attribute(status, 'restart_test_flag' &
         , i=31415, iflag=AF_RST_CMP)
    CALL channel_halt(substr, status)

    ! diagnostic output 
    CALL write_attribute(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_gatts
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_dims

    USE mo_model_domain,         ONLY: p_patch

    USE mo_grid_config,          ONLY: n_dom

    USE mo_run_config,           ONLY: ntracer

    USE mo_nonhydrostatic_config, ONLY: iadv_rcf, itime_scheme

    USE mo_io_config,            ONLY: inextra_2d

    USE messy_main_timer,        ONLY: delta_time, current_time_step       &
                                  , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                  , YEAR_START, MONTH_START, DAY_START     &
                                  , HOUR_START, MINUTE_START, SECOND_START

    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att

    IMPLICIT NONE
    INTRINSIC :: RANDOM_SEED, ABS, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_channel_initialize_dims'
    INTEGER                             :: status
    INTEGER                             :: nlev, nlevp1, ncells, nverts        &
         &                               , nedges, ntl_prog
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
    INTEGER                             :: i, jg
    CHARACTER(LEN=80)                   :: start_date_string
    REAL(dp)                            :: yyyymmdd
    INTEGER                             :: n

    ! #############
    ! ### TIME ####
    ! #############
    ! ... dimension ...
    CALL new_dimension(status, DIMID_TIME, 'time', 1, .TRUE.)
    CALL channel_halt(substr, status)

    ! ... dimension variable 'time' ...
    CALL add_dimension_variable(status, 'time', 'time', (/ 0.0_dp /))
    CALL channel_halt(substr, status)
    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'long_name', c='time')
    CALL channel_halt(substr, status)
    ! ... with attribute 'units' ...

    WRITE(start_date_string, &
         '("day since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         YEAR_START,MONTH_START,DAY_START,HOUR_START,MINUTE_START,SECOND_START
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'units', c=start_date_string)
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'time', &
         'calendar', c='gregorian')
    CALL channel_halt(substr, status)
    ! - YYYYMMDD
    yyyymmdd = ABS(YEAR)*10000 + MONTH*100 + DAY  &
         + (HOUR*3600 + MINUTE*60 + SECOND)/86400.
    IF (YEAR<0) yyyymmdd = -yyyymmdd
    CALL add_dimension_variable(status, 'time', 'YYYYMMDD', (/ yyyymmdd /))
    CALL channel_halt(substr, status)   
    !
    CALL add_dimension_variable_att(status, 'time', 'YYYYMMDD', &
         'long_name', c='time')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'YYYYMMDD', &
         'units', c='days as %Y%m%d.%f')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'YYYYMMDD', &
         'calendar', c='gregorian')
    CALL channel_halt(substr, status)
    ! - DT
    CALL add_dimension_variable(status, 'time', 'dt', (/ delta_time /))
    CALL channel_halt(substr, status)   
    !
    CALL add_dimension_variable_att(status, 'time', 'dt', &
         'long_name', c='delta_time')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'dt', &
         'units', c='s')
    CALL channel_halt(substr, status)    
    ! - CURRENT TIME STEP
    CALL add_dimension_variable(status, 'time', 'nstep' &
         , (/ REAL(current_time_step, DP) /))
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'time', 'nstep', &
         'long_name', c='current time step')
    CALL channel_halt(substr, status) 

    ! ###############
    ! ### LEVELS ####
    ! ###############

    ALLOCATE(DIMID_LEV_DOM(n_dom)                                             &
         & , DIMID_2LEV_DOM(n_dom)                                            &
         & , DIMID_3LEV_DOM(n_dom)                                            &
         & , DIMID_LEVP1_DOM(n_dom)                                           &
         & , DIMID_NCELLS_DOM(n_dom)                                          &
         & , DIMID_NVERTS_DOM(n_dom)                                          &
         & , DIMID_NEDGES_DOM(n_dom)                                          &
         & , DIMID_NTRACER_DOM(n_dom)                                         &
         & , DIMID_NTLV_DOM(n_dom)                                            &
         & , DIMID_UBCP_DOM(n_dom)                                            &
         & , DIMID_UBCP1_DOM(n_dom)                                           &
         & , DIMID_EXTRA2D_DOM(n_dom)                                         &
         & )
    DIMID_LEV_DOM(:)     = DIMID_UNDEF
    DIMID_2LEV_DOM(:)    = DIMID_UNDEF
    DIMID_3LEV_DOM(:)    = DIMID_UNDEF
    DIMID_LEVP1_DOM(:)   = DIMID_UNDEF
    DIMID_NCELLS_DOM(:)  = DIMID_UNDEF
    DIMID_NVERTS_DOM(:)  = DIMID_UNDEF
    DIMID_NEDGES_DOM(:)  = DIMID_UNDEF
    DIMID_NTRACER_DOM(:) = DIMID_UNDEF
    DIMID_NTLV_DOM(:)    = DIMID_UNDEF
    DIMID_UBCP_DOM(:)    = DIMID_UNDEF
    DIMID_UBCP1_DOM(:)   = DIMID_UNDEF
    DIMID_EXTRA2D_DOM(:) = DIMID_UNDEF

    DO jg=1,n_dom
       nlev = p_patch(jg)%nlev
       nlevp1 = p_patch(jg)%nlevp1
       ncells = p_patch(jg)%n_patch_cells_g
       nverts = p_patch(jg)%n_patch_verts_g
       nedges = p_patch(jg)%n_patch_edges_g
       IF (itime_scheme >= 2) THEN
          ntl_prog = 2
       ELSE
          ntl_prog = 1
       END IF
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_LEV_DOM(jg), 'lev', nlev)
       CALL channel_halt(substr, status)
       !
       ! ... dimension variable ...
       ALLOCATE(array(nlev))
       DO i=1, nlev
          array(i) = REAL(i, DP)
       END DO
       CALL add_dimension_variable(status, 'lev', 'lev', array)
       CALL channel_halt(substr, status)
       DEALLOCATE(array)

       ! ... with attribute 'long_name' ...
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
            'long_name', c='level index')
       CALL channel_halt(substr, status)

       ! ... with attribute 'units' ...
       CALL add_dimension_variable_att(status, 'lev', 'lev', &
            'units', c='level')
       CALL channel_halt(substr, status)
       !
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_3LEV_DOM(jg), 'lev_3', 3)
       CALL channel_halt(substr, status)
       !
       ! ... dimension variable ...
       ALLOCATE(array(3))
       DO i=1, 3
          array(i) = REAL(i, DP)
       END DO
       CALL add_dimension_variable(status, 'lev_3', 'lev_3', array)
       CALL channel_halt(substr, status)
       DEALLOCATE(array)

       ! ... with attribute 'long_name' ...
       CALL add_dimension_variable_att(status, 'lev_3', 'lev_3', &
            'long_name', c='level index')
       CALL channel_halt(substr, status)

       ! ... with attribute 'units' ...
       CALL add_dimension_variable_att(status, 'lev_3', 'lev_3', &
            'units', c='level')
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_2LEV_DOM(jg), 'lev_4', 2)
       CALL channel_halt(substr, status)
       !
       ! ... dimension variable ...
       ALLOCATE(array(2))
       DO i=1, 2
          array(i) = REAL(i, DP)
       END DO
       CALL add_dimension_variable(status, 'lev_4', 'lev_4', array)
       CALL channel_halt(substr, status)
       DEALLOCATE(array)

       ! ... with attribute 'long_name' ...
       CALL add_dimension_variable_att(status, 'lev_4', 'lev_4', &
            'long_name', c='level index')
       CALL channel_halt(substr, status)

       ! ... with attribute 'units' ...
       CALL add_dimension_variable_att(status, 'lev_4', 'lev_4', &
            'units', c='level')
       CALL channel_halt(substr, status)
       !
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_UBCP_DOM(jg), 'ubcp', iadv_rcf+2)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_UBCP1_DOM(jg), 'ubcp1', iadv_rcf+1)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_EXTRA2D_DOM(jg), 'extra_2d', inextra_2d)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_LEVP1_DOM(jg), 'lev_2', nlevp1)
       CALL channel_halt(substr, status)
       !
       ! ... dimension variable ...
       ALLOCATE(array(nlevp1))
       DO i=1, nlevp1
          array(i) = REAL(i, DP)
       END DO
       CALL add_dimension_variable(status, 'lev_2', 'lev_2', array)
       CALL channel_halt(substr, status)
       DEALLOCATE(array)

       ! ... with attribute 'long_name' ...
       CALL add_dimension_variable_att(status, 'lev_2', 'lev_2', &
            'long_name', c='level index')
       CALL channel_halt(substr, status)

       ! ... with attribute 'units' ...
       CALL add_dimension_variable_att(status, 'lev_2', 'lev_2', &
            'units', c='level')
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_NCELLS_DOM(jg), 'ncells', ncells)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_NVERTS_DOM(jg), 'nverts', nverts)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_NEDGES_DOM(jg), 'nedges', nedges)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_NTRACER_DOM(jg), 'ntracer', ntracer)
       CALL channel_halt(substr, status)
       !
       ! ... dimension ...
       CALL new_dimension(status, DIMID_NTLV_DOM(jg), 'ntl', ntl_prog)
       CALL channel_halt(substr, status)
       !
    END DO

    !
    ! ... dimension ...
    CALL new_dimension(status, DIMID_ONE, 'one', 1)
    CALL channel_halt(substr, status)
    !
    !
    ! ... dimension ...
    CALL new_dimension(status, DIMID_TWO, 'two', 2)
    CALL channel_halt(substr, status)
    !
    !
    ! ... dimension ...
    CALL new_dimension(status, DIMID_FIVE, 'five', 5)
    CALL channel_halt(substr, status)
    !

    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_dimension(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_dims
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE main_channel_initialize_reprs

    USE mo_util_string,                  ONLY: int2string

    USE mo_parallel_config,              ONLY: nproma

    USE mo_model_domain,                 ONLY: p_patch

    USE mo_grid_config,                  ONLY: n_dom

    USE mo_run_config,                   ONLY: ntracer

    USE mo_nonhydrostatic_config,        ONLY: iadv_rcf

    USE mo_io_config,                    ONLY: inextra_2d

    USE messy_main_channel_repr,         ONLY: new_representation        &
                                             , write_representation_dc   &
                                             , set_representation_decomp &
                                             , write_representation      &
                                             , AUTO, IRANK               &
                                             , PIOTYPE_SGL, PIOTYPE_IND  &
                                             , PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'main_channel_initialize_reprs'
    INTEGER                        :: status, jg
    INTEGER, DIMENSION(:), POINTER :: scdim => NULL()
    INTEGER                        :: nblks_c, nblks_v, nblks_e, nlev, nlevp1  &
         &                          , ntl_prog

    ALLOCATE(UNSTRUCTURED_CELL_HORIZONTAL_DOM(n_dom)                           &
         & , UNSTRUCTURED_CELL_HOR_2_DOM(n_dom)                                &
         & , UNSTRUCTURED_CELL_HOR_5_DOM(n_dom)                                &
         & , UNSTRUCTURED_CELL_HOR_TRACER_DOM(n_dom)                           &
         & , UNSTRUCTURED_CELL_HOR_UBCP_DOM(n_dom)                             &
         & , UNSTRUCTURED_CELL_HOR_UBCP1_DOM(n_dom)                            &
         & , UNSTRUCTURED_CELL_HOR_EXTRA_DOM(n_dom)                            &
         & , UNSTRUCTURED_CELL_LEV_DOM(n_dom)                                  &
         & , UNSTRUCTURED_CELL_2LEV_DOM(n_dom)                                 &
         & , UNSTRUCTURED_CELL_3LEV_DOM(n_dom)                                 &
         & , UNSTRUCTURED_CELL_LEV_TRACER_DOM(n_dom)                           &
         & , UNSTRUCTURED_CELL_LEV_TLV_DOM(n_dom)                              &
         & , UNSTRUCTURED_CELL_LEVP1_DOM(n_dom)                                &
         & , UNSTRUCTURED_CELL_LEVP1_TRACER_DOM(n_dom)                         &
         & , UNSTRUCTURED_CELL_LEVP1_TLV_DOM(n_dom)                            &
         & )
    UNSTRUCTURED_CELL_HORIZONTAL_DOM(:)   = REPR_UNDEF
    UNSTRUCTURED_CELL_HOR_2_DOM(:)        = REPR_UNDEF
    UNSTRUCTURED_CELL_HOR_5_DOM(:)        = REPR_UNDEF
    UNSTRUCTURED_CELL_HOR_TRACER_DOM(:)   = REPR_UNDEF
    UNSTRUCTURED_CELL_HOR_UBCP_DOM(:)     = REPR_UNDEF
    UNSTRUCTURED_CELL_HOR_UBCP1_DOM(:)    = REPR_UNDEF
    UNSTRUCTURED_CELL_HOR_EXTRA_DOM(:)    = REPR_UNDEF
    UNSTRUCTURED_CELL_LEV_DOM(:)          = REPR_UNDEF
    UNSTRUCTURED_CELL_2LEV_DOM(:)         = REPR_UNDEF
    UNSTRUCTURED_CELL_3LEV_DOM(:)         = REPR_UNDEF
    UNSTRUCTURED_CELL_LEV_TRACER_DOM(:)   = REPR_UNDEF
    UNSTRUCTURED_CELL_LEV_TLV_DOM(:)      = REPR_UNDEF
    UNSTRUCTURED_CELL_LEVP1_DOM(:)        = REPR_UNDEF
    UNSTRUCTURED_CELL_LEVP1_TRACER_DOM(:) = REPR_UNDEF
    UNSTRUCTURED_CELL_LEVP1_TLV_DOM(:)    = REPR_UNDEF

    ALLOCATE(UNSTRUCTURED_VERT_HORIZONTAL_DOM(n_dom)                           &
         & , UNSTRUCTURED_VERT_LEV_DOM(n_dom)                                  &
         & , UNSTRUCTURED_VERT_LEVP1_DOM(n_dom))
    UNSTRUCTURED_VERT_HORIZONTAL_DOM(:) = REPR_UNDEF
    UNSTRUCTURED_VERT_LEV_DOM(:)       = REPR_UNDEF
    UNSTRUCTURED_VERT_LEVP1_DOM(:)     = REPR_UNDEF

    ALLOCATE(UNSTRUCTURED_EDGE_HORIZONTAL_DOM(n_dom)                           &
         & , UNSTRUCTURED_EDGE_LEV_DOM(n_dom)                                  &
         & , UNSTRUCTURED_EDGE_2LEV_DOM(n_dom)                                 &
         & , UNSTRUCTURED_EDGE_3LEV_DOM(n_dom)                                 &
         & , UNSTRUCTURED_EDGE_LEV_TLV_DOM(n_dom)                              &
         & , UNSTRUCTURED_EDGE_LEVP1_DOM(n_dom))
    UNSTRUCTURED_EDGE_HORIZONTAL_DOM(:) = REPR_UNDEF
    UNSTRUCTURED_EDGE_LEV_DOM(:)        = REPR_UNDEF
    UNSTRUCTURED_EDGE_2LEV_DOM(:)       = REPR_UNDEF
    UNSTRUCTURED_EDGE_3LEV_DOM(:)       = REPR_UNDEF
    UNSTRUCTURED_EDGE_LEV_TLV_DOM(:)      = REPR_UNDEF
    UNSTRUCTURED_EDGE_LEVP1_DOM(:)      = REPR_UNDEF

    DO jg=1,n_dom
       nblks_c = p_patch(jg)%nblks_c
       nblks_v = p_patch(jg)%nblks_v
       nblks_e = p_patch(jg)%nblks_e
       nlev    = p_patch(jg)%nlev
       nlevp1  = p_patch(jg)%nlevp1

       CALL new_representation(status, UNSTRUCTURED_CELL_HORIZONTAL_DOM(jg)    &
            & , TRIM('UNSTR_CELL_HOR_D'//int2string(jg,'(I2.2)'))              &
            & , rank = 2, link = 'xx--', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NCELLS_DOM(jg), DIMID_ONE /)          &
            & , ldimlen       = (/ nproma     , nblks_c /)                     &
            & , output_order  = (/ 2, 1 /)                                     &
            & , axis = 'N---'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_HOR_TRACER_DOM(jg)    &
            & , TRIM('UNSTR_CELL_HOR_TR_D'//int2string(jg,'(I2.2)'))           &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            ! op_bk_20131121+
!             & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_ONE            &
!             &                     , DIMID_NTRACER_DOM(jg) /)                   &
!             & , ldimlen       = (/ nproma, nblks_c, ntracer /)                 &
!             & , output_order  = (/ 1, 3, 2/)                                   &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_NTRACER_DOM(jg)&
            &                     , DIMID_ONE /)                               &
            & , ldimlen       = (/ nproma, nblks_c, ntracer /)                 &
            & , output_order  = (/ 3, 2, 1/)                                   &
            ! op_bk_20131121-
            & , axis = '-N--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_HOR_2_DOM(jg)         &
            & , TRIM('UNSTR_CELL_HOR_2_D'//int2string(jg,'(I2.2)'))            &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_TWO            &
            &                     , DIMID_ONE /)                               &
            & , ldimlen       = (/ nproma, nblks_c, 2 /)                       &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = '-N--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_HOR_5_DOM(jg)         &
            & , TRIM('UNSTR_CELL_HOR_5_D'//int2string(jg,'(I2.2)'))            &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_FIVE           &
            &                     , DIMID_ONE /)                               &
            & , ldimlen       = (/ nproma, nblks_c, 5 /)                       &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = '-N--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_HOR_UBCP_DOM(jg)      &
            & , TRIM('UNSTR_CELL_HOR_UBCP_D'//int2string(jg,'(I2.2)'))         &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_UBCP_DOM(jg)   &
            &                     , DIMID_ONE /)                               &
            & , ldimlen       = (/ nproma, nblks_c, iadv_rcf+2 /)              &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = '-N--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_HOR_UBCP1_DOM(jg)     &
            & , TRIM('UNSTR_CELL_HOR_UBCP1_D'//int2string(jg,'(I2.2)'))        &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_UBCP1_DOM(jg)  &
            &                     , DIMID_ONE /)                               &
            & , ldimlen       = (/ nproma, nblks_c, iadv_rcf+1 /)              &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = '-N--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_HOR_EXTRA_DOM(jg)     &
            & , TRIM('UNSTR_CELL_HOR_EXTRA_D'//int2string(jg,'(I2.2)'))        &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_EXTRA2D_DOM(jg)&
            &                     , DIMID_ONE /)                               &
            & , ldimlen       = (/ nproma, nblks_c, inextra_2d /)              &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = '-N--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_LEV_DOM(jg)           &
            & , TRIM('UNSTR_CELL_LEV_D'//int2string(jg,'(I2.2)'))              &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_LEV_DOM(jg)    &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_c /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_2LEV_DOM(jg)          &
            & , TRIM('UNSTR_CELL_2LEV_D'//int2string(jg,'(I2.2)'))             &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_2LEV_DOM(jg)   &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_c /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_3LEV_DOM(jg)          &
            & , TRIM('UNSTR_CELL_3LEV_D'//int2string(jg,'(I2.2)'))             &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_3LEV_DOM(jg)   &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_c /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_LEV_TRACER_DOM(jg)    &
            & , TRIM('UNSTR_CELL_LEV_TR_D'//int2string(jg,'(I2.2)'))           &
            & , rank = 4, link = 'xxxx', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_LEV_DOM(jg)    &
            &                     , DIMID_ONE, DIMID_NTRACER_DOM(jg) /)        &
            & , ldimlen       = (/ nproma, AUTO, nblks_c, AUTO /)              &
            & , output_order  = (/ 3, 4, 2, 1/)                                &
            & , axis = 'Z--N'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_LEV_TLV_DOM(jg)       &
            & , TRIM('UNSTR_CELL_LEV_TL_D'//int2string(jg,'(I2.2)'))           &
            & , rank = 4, link = 'xxxx', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NCELLS_DOM(jg), DIMID_LEV_DOM(jg)    &
            &                     , DIMID_ONE, DIMID_NTLV_DOM(jg) /)           &
            & , ldimlen       = (/ nproma, AUTO, nblks_c, AUTO /)              &
            & , output_order  = (/ 3, 4, 2, 1/)                                &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_LEVP1_DOM(jg)         &
            & , TRIM('UNSTR_CELL_LEVP1_D'//int2string(jg,'(I2.2)'))            &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NCELLS_DOM(jg), DIMID_LEVP1_DOM(jg)   &
            &                    , DIMID_ONE /)                                &
            & , ldimlen       = (/ nproma, AUTO , nblks_c /)                   &
            & , output_order  = (/ 3, 2, 1 /)                                  &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_LEVP1_TRACER_DOM(jg)  &
            & , TRIM('UNSTR_CELL_LEVP1_TR_D'//int2string(jg,'(I2.2)'))         &
            & , rank = 4, link = 'xxxx', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NCELLS_DOM(jg), DIMID_LEVP1_DOM(jg)   &
            &                    , DIMID_ONE, DIMID_NTRACER_DOM(jg) /)         &
            & , ldimlen       = (/ nproma, AUTO , nblks_c, AUTO /)             &
            & , output_order  = (/ 3, 4, 2, 1 /)                               &
            & , axis = 'Z--N'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_CELL_LEVP1_TLV_DOM(jg)     &
            & , TRIM('UNSTR_CELL_LEVP1_TL_D'//int2string(jg,'(I2.2)'))         &
            & , rank = 4, link = 'xxxx', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NCELLS_DOM(jg), DIMID_LEVP1_DOM(jg)   &
            &                    , DIMID_ONE, DIMID_NTLV_DOM(jg) /)            &
            & , ldimlen       = (/ nproma, AUTO , nblks_c, AUTO /)             &
            & , output_order  = (/ 3, 4, 2, 1 /)                               &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_VERT_HORIZONTAL_DOM(jg)    &
            & , TRIM('UNSTR_VERT_HOR_D'//int2string(jg,'(I2.2)'))              &
            & , rank = 2, link = 'xx--', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NVERTS_DOM(jg), DIMID_ONE /)          &
            & , ldimlen       = (/ nproma     , nblks_v /)                     &
            & , output_order  = (/ 2, 1 /)                                     &
            & , axis = 'N---'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_VERT_LEV_DOM(jg)           &
            & , TRIM('UNSTR_VERT_LEV_D'//int2string(jg,'(I2.2)'))              &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NVERTS_DOM(jg), DIMID_LEV_DOM(jg)    &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_v /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_VERT_LEVP1_DOM(jg)         &
            & , TRIM('UNSTR_VERT_LEVP1_D'//int2string(jg,'(I2.2)'))            &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NVERTS_DOM(jg), DIMID_LEVP1_DOM(jg)   &
            &                    , DIMID_ONE /)                                &
            & , ldimlen       = (/ nproma, AUTO , nblks_v /)                   &
            & , output_order  = (/ 3, 2, 1 /)                                  &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_EDGE_HORIZONTAL_DOM(jg)    &
            & , TRIM('UNSTR_EDGE_HOR_D'//int2string(jg,'(I2.2)'))              &
            & , rank = 2, link = 'xx--', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NEDGES_DOM(jg), DIMID_ONE /)          &
            & , ldimlen       = (/ nproma     , nblks_e /)                     &
            & , output_order  = (/ 2, 1 /)                                     &
            & , axis = 'N---'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_EDGE_LEV_DOM(jg)           &
            & , TRIM('UNSTR_EDGE_LEV_D'//int2string(jg,'(I2.2)'))              &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NEDGES_DOM(jg), DIMID_LEV_DOM(jg)    &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_e /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_EDGE_2LEV_DOM(jg)          &
            & , TRIM('UNSTR_EDGE_2LEV_D'//int2string(jg,'(I2.2)'))             &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NEDGES_DOM(jg), DIMID_2LEV_DOM(jg)   &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_e /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_EDGE_3LEV_DOM(jg)          &
            & , TRIM('UNSTR_EDGE_3LEV_D'//int2string(jg,'(I2.2)'))             &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/  DIMID_NEDGES_DOM(jg), DIMID_3LEV_DOM(jg)   &
            &                     , DIMID_ONE/)                                &
            & , ldimlen       = (/ nproma, AUTO, nblks_e /)                    &
            & , output_order  = (/ 3, 2, 1/)                                   &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_EDGE_LEVP1_DOM(jg)         &
            & , TRIM('UNSTR_EDGE_LEVP1_D'//int2string(jg,'(I2.2)'))            &
            & , rank = 3, link = 'xxx-', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NEDGES_DOM(jg), DIMID_LEVP1_DOM(jg)   &
            &                    , DIMID_ONE /)                                &
            & , ldimlen       = (/ nproma, AUTO , nblks_e /)                   &
            & , output_order  = (/ 3, 2, 1 /)                                  &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

       CALL new_representation(status, UNSTRUCTURED_EDGE_LEV_TLV_DOM(jg)       &
            & , TRIM('UNSTR_EDGE_LEV_TL_D'//int2string(jg,'(I2.2)'))           &
            & , rank = 4, link = 'xxxx', dctype = DC_GP                        &
            & , dimension_ids = (/ DIMID_NEDGES_DOM(jg), DIMID_LEV_DOM(jg)     &
            &                    , DIMID_ONE, DIMID_NTLV_DOM(jg) /)            &
            & , ldimlen       = (/ nproma, AUTO , nblks_e, AUTO /)             &
            & , output_order  = (/ 3, 4, 2, 1 /)                               &
            & , axis = 'ZN--'                                                  &
            & , patch_id = jg                                                  &
            & )
       CALL channel_halt(substr, status)

    END DO
    ! #################
    ! DIAGNOSTIC OUTPUT
    ! #################
    CALL write_representation(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_channel_initialize_reprs
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE associate_var_lists_to_channels

    ! ICON
    ! op_bk_20130830+
    USE mo_var_list,      ONLY: var_lists, nvar_lists, get_var
    USE mo_linked_list,   ONLY: t_var_list, t_list_element
    USE mo_var_metadata,  ONLY: t_var_metadata
    USE mo_grid_config,   ONLY: n_dom
    USE mo_run_config,    ONLY: ntracer
    USE mo_nonhydrostatic_config, ONLY: iadv_rcf, itime_scheme

    USE mo_io_config,            ONLY: inextra_2d

    USE mo_cdi_constants, ONLY: GRID_UNSTRUCTURED_CELL                         &
         &                    , GRID_UNSTRUCTURED_VERT                         &
         &                    , GRID_UNSTRUCTURED_EDGE                         &
         &                    , ZA_HYBRID, ZA_HYBRID_HALF, ZA_HYBRID_HALF_HHL  &
         &                    , is_2d_field
    USE mo_model_domain,  ONLY: p_patch
    ! op_bk_20130830-
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,  ONLY: p_parallel_io
    ! MESSy
    USE messy_main_channel, ONLY: new_channel, new_channel_object &
                                , new_attribute                   &
                                , AF_RST_CMP, AF_RST_INP
    USE messy_main_channel_repr, ONLY: t_representation, get_representation    &
         &                           , get_representation_info

    USE messy_main_constants_mem,  ONLY: STRLEN_MEDIUM

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER                      :: i
    ! op_bk_20130830+
    TYPE (t_var_list),     POINTER :: var_list
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: element_1st
    TYPE (t_var_metadata), POINTER :: info
    INTEGER                        :: patch_id
!     TYPE (t_stream),     POINTER :: stream
!     TYPE (list_element), POINTER :: element
!     TYPE (list_element), TARGET  :: element_1st
!     TYPE (memory_info),  POINTER :: info
    ! op_bk_20130830-
    !
    CHARACTER(LEN=*), PARAMETER  :: substr = 'associate_var_lists_to_channels'
    INTEGER                      :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: mem
    INTEGER                      :: reprid
    INTEGER                      :: nblks_c, nblks_v, nblks_e, nlev, nlevp1    &
         &                        , ntl_prog
    LOGICAL                      :: lok
    TYPE(t_representation), POINTER :: repr
    CHARACTER(LEN=STRLEN_MEDIUM)    :: repr_name

    CALL start_message_bi(modstr,'VAR_LISTS TO CHANNELS',substr)

    var_list_loop: DO i=1, nvar_lists

       var_list => var_lists(i)
       patch_id = var_list%p%patch_id

       nblks_c = p_patch(patch_id)%nblks_c
       nblks_v = p_patch(patch_id)%nblks_v
       nblks_e = p_patch(patch_id)%nblks_e
       nlev    = p_patch(patch_id)%nlev
       nlevp1  = p_patch(patch_id)%nlevp1

       IF (itime_scheme >= 2) THEN
          WRITE(*,*) "MESSY: ntl_prog = 2"
          ntl_prog = 2
       ELSE
          WRITE(*,*) "MESSY: ntl_prog = 1"
          ntl_prog = 1
       END IF

       WRITE(*,'(A)') "MESSY: processing channel", TRIM(var_list%p%name)
       IF (TRIM(var_list%p%name) == "nh_state_diag_of_domain_01") THEN
          WRITE(*,*) "SKIPPED!"
          CYCLE
       END IF
       ! ---------------------------
       ! VAR_LIST -> CHANNEL
       ! ---------------------------
       CALL new_channel(status, TRIM(var_list%p%name))
       CALL channel_halt(substr, status)
       ! ---------------------------

       element_1st%next_list_element => var_list%p%first_list_element
       element => element_1st
       element_loop: DO ! loop over elements in linked list
          element => element%next_list_element
          IF (.NOT.ASSOCIATED(element)) EXIT

          info => element%field%info

          ! ---------------------------
          ! ELEMENT -> OBJECT
          ! ---------------------------
          SELECT CASE(info%hgrid)
          CASE(GRID_UNSTRUCTURED_CELL)
             !
             IF (is_2d_field(info%vgrid)) THEN
                SELECT CASE(info%ndims)
                CASE(2)
                   lok = .TRUE.
                   reprid = UNSTRUCTURED_CELL_HORIZONTAL_DOM(patch_id)
                CASE(3)
                   IF (info%used_dimensions(3) == ntracer) THEN
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_HOR_TRACER_DOM(patch_id)
                   ELSE IF ((INDEX(info%name, 'tracer') /= 0) .AND.            &
                        & (info%used_dimensions(3) == 5)) THEN
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_HOR_5_DOM(patch_id)
                   ELSE IF (info%used_dimensions(3) == 2) THEN
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_HOR_2_DOM(patch_id)
                   ELSE IF (info%used_dimensions(3) == iadv_rcf + 2) THEN
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_HOR_UBCP_DOM(patch_id)
                   ELSE IF (info%used_dimensions(3) == iadv_rcf + 1) THEN
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_HOR_UBCP1_DOM(patch_id)
                   ELSE IF (info%used_dimensions(3) == inextra_2d) THEN
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_HOR_EXTRA_DOM(patch_id)
                   ELSE
                      lok = .FALSE.
                   END IF
                CASE DEFAULT
                   lok = .FALSE.
                END SELECT
             ELSE
                SELECT CASE(info%vgrid)
                CASE(ZA_HYBRID)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      IF (info%used_dimensions(2) == nlev) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEV_DOM(patch_id)
                      ELSE IF (info%used_dimensions(2) == 2) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_2LEV_DOM(patch_id)
                      ELSE IF (info%used_dimensions(2) == 3) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_3LEV_DOM(patch_id)
                      ELSE
                         lok = .FALSE.
                      END IF
                   CASE(4)
                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
                           & (info%used_dimensions(4) == ntracer)) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEV_TRACER_DOM(patch_id)
                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEV_TLV_DOM(patch_id)
                      ELSE
                         lok = .FALSE.
                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE(ZA_HYBRID_HALF)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_LEVP1_DOM(patch_id)
                   CASE(4)
                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
                           & (info%used_dimensions(4) == ntracer)) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEVP1_TRACER_DOM(patch_id)
                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEVP1_TLV_DOM(patch_id)
                      ELSE
                         lok = .FALSE.
                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE(ZA_HYBRID_HALF_HHL)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_CELL_LEVP1_DOM(patch_id)
                   CASE(4)
                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
                           & (info%used_dimensions(4) == ntracer)) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEVP1_TRACER_DOM(patch_id)
                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_CELL_LEVP1_TLV_DOM(patch_id)
                      ELSE
                         lok = .FALSE.
                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE DEFAULT
                   lok = .FALSE.
                END SELECT
             END IF
          CASE(GRID_UNSTRUCTURED_VERT)
             !
             IF (is_2d_field(info%vgrid)) THEN
                SELECT CASE(info%ndims)
                CASE(2)
                   lok = .TRUE.
                   reprid = UNSTRUCTURED_VERT_HORIZONTAL_DOM(patch_id)
                CASE(3)
                   lok = .FALSE.
                CASE DEFAULT
                   lok = .FALSE.
                END SELECT
             ELSE
                SELECT CASE(info%vgrid)
                CASE(ZA_HYBRID)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_VERT_LEV_DOM(patch_id)
                   CASE(4)
!                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
!                           & (info%used_dimensions(4) == ntracer)) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_VERT_LEV_TRACER_DOM(patch_id)
!                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_VERT_LEV_TLV_DOM(patch_id)
!                      ELSE
                         lok = .FALSE.
!                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE(ZA_HYBRID_HALF)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_VERT_LEVP1_DOM(patch_id)
                   CASE(4)
!                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
!                           & (info%used_dimensions(4) == ntracer)) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_VERT_LEVP1_TRACER_DOM(patch_id)
!                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_VERT_LEVP1_TLV_DOM(patch_id)
!                      ELSE
                         lok = .FALSE.
!                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE(ZA_HYBRID_HALF_HHL)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_VERT_LEVP1_DOM(patch_id)
                   CASE(4)
!                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
!                           & (info%used_dimensions(4) == ntracer)) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_VERT_LEVP1_TRACER_DOM(patch_id)
!                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_VERT_LEVP1_TLV_DOM(patch_id)
!                      ELSE
                         lok = .FALSE.
!                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE DEFAULT
                   lok = .FALSE.
                END SELECT
             END IF
          CASE(GRID_UNSTRUCTURED_EDGE)
             !
             IF (is_2d_field(info%vgrid)) THEN
                SELECT CASE(info%ndims)
                CASE(2)
                   lok = .TRUE.
                   reprid = UNSTRUCTURED_EDGE_HORIZONTAL_DOM(patch_id)
                CASE(3)
                   lok = .FALSE.
                CASE DEFAULT
                   lok = .FALSE.
                END SELECT
             ELSE
                SELECT CASE(info%vgrid)
                CASE(ZA_HYBRID)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      IF (info%used_dimensions(2) == nlev) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_EDGE_LEV_DOM(patch_id)
                      ELSE IF (info%used_dimensions(2) == 2) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_EDGE_2LEV_DOM(patch_id)
                      ELSE IF (info%used_dimensions(2) == 3) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_EDGE_3LEV_DOM(patch_id)
                      ELSE
                         lok = .FALSE.
                      END IF
                   CASE(4)
                      IF (info%used_dimensions(4) == ntl_prog) THEN
                         lok = .TRUE.
                         reprid = UNSTRUCTURED_EDGE_LEV_TLV_DOM(patch_id)
                      ELSE
                         lok = .FALSE.
                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE(ZA_HYBRID_HALF)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_EDGE_LEVP1_DOM(patch_id)
                   CASE(4)
!                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
!                           & (info%used_dimensions(4) == ntracer)) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_EDGE_LEVP1_TRACER_DOM(patch_id)
!                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_EDGE_LEVP1_TLV_DOM(patch_id)
!                      ELSE
                         lok = .FALSE.
!                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE(ZA_HYBRID_HALF_HHL)
                   SELECT CASE (info%ndims)
                   CASE(3)
                      lok = .TRUE.
                      reprid = UNSTRUCTURED_EDGE_LEVP1_DOM(patch_id)
                   CASE(4)
!                      IF ((INDEX(info%name, 'tracer') /= 0).AND.               &
!                           & (info%used_dimensions(4) == ntracer)) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_EDGE_LEVP1_TRACER_DOM(patch_id)
!                      ELSE IF (info%used_dimensions(4) == ntl_prog) THEN
!                         lok = .TRUE.
!                         reprid = UNSTRUCTURED_EDGE_LEVP1_TLV_DOM(patch_id)
!                      ELSE
                         lok = .FALSE.
!                      END IF
                   CASE DEFAULT
                      lok = .FALSE.
                   END SELECT
                CASE DEFAULT
                   lok = .FALSE.
                END SELECT
             END IF
          CASE DEFAULT
             lok = .FALSE.
          END SELECT

!!!!!!
             IF (lok) THEN

                IF ((info%allocated) .AND. ASSOCIATED(element%field%r_ptr)) THEN
!                   CALL get_var(var_list, TRIM(info%name), mem)
                   mem => element%field%r_ptr(:,:,:,:,1)
                ELSE
                   CYCLE
                END IF

                CALL new_channel_object(status, TRIM(var_list%p%name) &
                     , TRIM(info%name)                            &
                     , mem = mem                                  &
                     , reprid = reprid                            &
                     , lrestreq = .NOT. info%lrestart_cont .OR. info%lrestart )
                ! op_bk_20130909+
!                 IF ( status == 0 ) THEN
!                    CALL get_representation_info(status,'',id=reprid,name=repr_name)
!                    WRITE(*,'(A24, A32, A20, A25)') "MESSY: CHANNEL OBJECT "    &
!                         & , TRIM(var_list%p%name), TRIM(info%name), TRIM(repr_name)
!                    WRITE(*,'(5I4, A4)') info%used_dimensions, " OK."
!                 ELSE
!                    CALL get_representation_info(status,'',id=reprid,name=repr_name)
!                    WRITE(*,'(A24, A32, A20, A25)') "MESSY: CHANNEL OBJECT "    &
!                         & , TRIM(var_list%p%name), TRIM(info%name), TRIM(repr_name)
!                    WRITE(*,'(5I4, A7)') info%used_dimensions, " ERROR."
!                 END IF
                ! op_bk_20130909-
                CALL channel_halt(substr, status)
             
             IF (TRIM(info%cf%long_name) /= '') THEN
                CALL new_attribute(status &
                     , TRIM(var_list%p%name), TRIM(info%name), 'long_name' &
                     , c=TRIM(info%cf%long_name), iflag=AF_RST_CMP)
                CALL channel_halt(substr, status)
             END IF
             IF (TRIM(info%cf%units) /= '') THEN
                CALL new_attribute(status &
                     , TRIM(var_list%p%name), TRIM(info%name), 'units' &
                     , c=TRIM(info%cf%units), iflag=AF_RST_CMP)
                CALL channel_halt(substr, status)
             END IF

             CALL get_representation(status, reprid, repr)
             CALL channel_halt(substr, status)
             CALL new_attribute(status &
                  , TRIM(var_list%p%name), TRIM(info%name), 'representation' &
                  , c=TRIM(repr%name), iflag=AF_RST_CMP)
             CALL channel_halt(substr, status)             

          ELSE

             IF (p_parallel_io) THEN
                WRITE(*,*) 'WARNING: stream/element ''' &
                     ,TRIM(var_list%p%name),'/', TRIM(info%name) &
                     ,''' cannot be associated to channel!'
                WRITE(*,'(A3, 5I4, A3)') "(/ ", info%used_dimensions, " /)."
             END IF

             END IF

          END DO element_loop
       END DO var_list_loop

    CALL end_message_bi(modstr,'STREAMS TO CHANNELS',substr)

  END SUBROUTINE associate_var_lists_to_channels
  ! -------------------------------------------------------------------
