!>
!! Subroutines needed for the time dependent SST and Sea Ice fraction
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!!
!! @author Pilar Ripodas, DWD, Offenbach (2012-12)
!!
!!
!! @par Revision History
!! Initial release by  Pilar Ripodas, DWD, Offenbach (2012-12)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_td_ext_data

  USE mo_kind,               ONLY: wp
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_io_units,           ONLY: filename_max
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_mpi,                ONLY: p_comm_work_test, p_comm_work
#ifdef NOMPI
  USE mo_mpi,                 ONLY: my_process_is_mpi_all_seq
#endif
  USE mo_io_config,           ONLY: default_read_method
  USE mo_read_interface,      ONLY: openInputFile, closeFile, on_cells, &
    &                               t_stream_id, read_2D_1time
  USE mo_datetime,            ONLY: t_datetime, month2hour
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell_int
  USE mo_grid_config,         ONLY: n_dom
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_seaice_nwp,          ONLY: frsi_min

  USE mo_extpar_config,       ONLY: generate_td_filename
  USE mo_lnd_nwp_config,      ONLY: sst_td_filename, ci_td_filename

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: set_actual_td_ext_data
  PUBLIC  :: read_td_ext_data_file

CONTAINS

  !>
  !! Set the SST and sea ice fraction to the values corresponding to datetime
  !!
  !! The SST and sea ice fraction are set to the values of the corresponding
  !!  day and month. Values are interpolated from  climatological monthly
  !!  means (ext_data_mode = 2) or from the actual monthly means (ext_data_mode = 3).
  !!  Another option not yet implemented (ext_data_mode = 4) will set SST and
  !!  sea ice cover to the actual (day, month, year) daily mean
  !!
  !! @par Revision History
  !! Developed  by P. Ripodas (2012-12)
  !!
  SUBROUTINE set_actual_td_ext_data (lread, datetime,datetime_old,ext_data_mode,  &
                                  &  p_patch, ext_data, p_lnd_state)

    LOGICAL , INTENT(IN)          :: lread !force the read of the ext dat files for sstice_mode=3
    TYPE(t_datetime), INTENT(IN)  :: datetime, datetime_old
    INTEGER,INTENT(IN)            :: ext_data_mode
    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    TYPE(t_lnd_state), INTENT(INOUT)     :: p_lnd_state(:)

    INTEGER                       :: month1, month2, year1, year2
    INTEGER                       :: m1, m2
    REAL (wp)                     :: pw1, pw2, pw2_old
    INTEGER                       :: jg, jb, jc, i_startidx, i_endidx
    INTEGER                       :: i_nchdom, i_rlstart, i_rlend
    INTEGER                       :: i_startblk, i_endblk

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_td_ext_data:set_actual_td_ext_data  '

    !---------------------------------------------------------------

      SELECT CASE (ext_data_mode)

       CASE (2) !SST and sea ice fraction updated based
                !  on the climatological monthly values
        CALL month2hour (datetime, month1, month2, pw2 )
        pw1 = 1._wp - pw2

        DO jg = 1, n_dom

           i_nchdom  = MAX(1,p_patch(jg)%n_childdom)
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk, i_rlstart, i_rlend)
           i_rlstart = grf_bdywidth_c+1
           i_rlend   = min_rlcell_int

           i_startblk = p_patch(jg)%cells%start_blk(i_rlstart,1)
           i_endblk   = p_patch(jg)%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
          DO jb=i_startblk, i_endblk

            CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
               & i_startidx, i_endidx, i_rlstart, i_rlend)

            DO jc = i_startidx, i_endidx
             IF (ext_data(jg)%atm_td%sst_m(jc,jb,1) >= 0._wp .AND.   &
                & ext_data(jg)%atm_td%sst_m(jc,jb,2) >= 0._wp ) THEN
               p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) =                   &
                            pw1*ext_data(jg)%atm_td%sst_m(jc,jb,month1)    &
                &         + pw2*ext_data(jg)%atm_td%sst_m(jc,jb,month2)
               p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) =                  &
                            pw1*ext_data(jg)%atm_td%fr_ice_m(jc,jb,month1) &
                &         + pw2*ext_data(jg)%atm_td%fr_ice_m(jc,jb,month2)

               IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
                 p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb)= 0._wp
               ELSE IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) > 1._wp-frsi_min ) THEN
                 p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb)= 1._wp
               END IF
              END IF
            ENDDO
          ENDDO
!$OMP END DO
!$OMP END PARALLEL

        END DO ! jg


       CASE (3) !SST and sea ice fraction updated based
                !  on the actual monthly values
        CALL month2hour (datetime_old, m1, m2, pw2_old )
        CALL month2hour (datetime, month1, month2, year1, year2, pw2 )

        WRITE( message_text,'(a,5i6,f10.5)') 'sst ci interp,',datetime%day,month1,month2,year1,year2,pw2
        CALL message  (routine, TRIM(message_text))
        IF (m1 /= month1 .OR. lread ) THEN
          CALL read_td_ext_data_file (month1,month2,year1,year2,p_patch(1:),ext_data)
        END IF
        pw1 = 1._wp - pw2

        DO jg = 1, n_dom

           i_nchdom  = MAX(1,p_patch(jg)%n_childdom)
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk, i_rlstart, i_rlend)
           i_rlstart = grf_bdywidth_c+1
           i_rlend   = min_rlcell_int

           i_startblk = p_patch(jg)%cells%start_blk(i_rlstart,1)
           i_endblk   = p_patch(jg)%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
          DO jb=i_startblk, i_endblk

            CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
               & i_startidx, i_endidx, i_rlstart, i_rlend)

            DO jc = i_startidx, i_endidx
              IF (ext_data(jg)%atm_td%sst_m(jc,jb,1) < 0._wp .AND.   &
                & ext_data(jg)%atm_td%sst_m(jc,jb,2) >= 0._wp )  THEN
                 WRITE( message_text,'(a,2g10.5)') 'something is wrong,', &
                 & ext_data(jg)%atm_td%sst_m(jc,jb,1), ext_data(jg)%atm_td%sst_m(jc,jb,2)
                 CALL message  (routine, TRIM(message_text))
              ELSE IF (ext_data(jg)%atm_td%sst_m(jc,jb,1) >= 0._wp .AND.   &
                & ext_data(jg)%atm_td%sst_m(jc,jb,2) >= 0._wp ) THEN
               p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) =              &
                            pw1*ext_data(jg)%atm_td%sst_m(jc,jb,1)    &
                &         + pw2*ext_data(jg)%atm_td%sst_m(jc,jb,2)
               p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) =            &
                            pw1*ext_data(jg)%atm_td%fr_ice_m(jc,jb,1) &
                &         + pw2*ext_data(jg)%atm_td%fr_ice_m(jc,jb,2)

               IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
                 p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb)= 0._wp
               ELSE IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) > 1._wp-frsi_min ) THEN
                 p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb)= 1._wp
               END IF
              END IF
            ENDDO
          ENDDO
!$OMP END DO
!$OMP END PARALLEL

        END DO ! jg

       CASE (4) !SST and sea ice fraction updated based
                !  on the actual daily values
        !Not implemented
        WRITE( message_text,'(a)') 'ext_data_mode == 4 not yet implemented '
        CALL finish  (routine, TRIM(message_text))

      END SELECT



  END SUBROUTINE set_actual_td_ext_data
!-----------------------------------------------------------------------

  !>
  !! <Short description of the subroutine for listings and indices>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE read_td_ext_data_file (m1,m2,y1,y2,p_patch,ext_data)
  !>
  !! <Short description of the function for listings and indices>
  !!
  !! <Describe the purpose of the function and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!

    INTEGER, INTENT(IN)                  :: m1,m2,y1,y2 ! month and year of the
                                               ! nearest months to datetime
    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(LEN=filename_max) :: extpar_file
    INTEGER                     :: jg
    TYPE(t_stream_id)           :: stream_id
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = 'mo_td_ext_data:read_td_ext_data_file:'

!-----------------------------------------------------------------------
   ! extpar_td_filename = "<path>extpar_<year>_<month>_<gridfile>"

   ! Set the months needed to interpolate the ext_para to the actual day

      DO jg = 1,n_dom
        IF (p_patch(jg)%geometry_info%cell_type == 6) THEN ! hexagonal grid

          CALL finish(TRIM(ROUTINE),&
            & 'Hexagonal grid is not supported, yet.')

        ENDIF
        WRITE( message_text,'(a,2i6,a)') 'n_dom, jg, grid_filename ',n_dom, &
          &                              jg, TRIM(p_patch(jg)%grid_filename)
        CALL message  (routine, TRIM(message_text))

        !! READ SST files
        extpar_file = generate_td_filename(sst_td_filename,                   &
          &                                getModelBaseDir(),                 &
          &                                TRIM(p_patch(jg)%grid_filename),   &
          &                                m1,y1                   )

        CALL message  (routine, TRIM(extpar_file))
        stream_id = openInputFile(extpar_file, p_patch(jg), default_read_method)
        CALL read_2D_1time(stream_id, on_cells, 'SST', &
          &          ext_data(jg)%atm_td%sst_m(:,:,1))
        CALL closeFile(stream_id)

        extpar_file = generate_td_filename(sst_td_filename,                   &
          &                                getModelBaseDir(),                 &
          &                                TRIM(p_patch(jg)%grid_filename),   &
          &                                m2, y2                   )
        CALL message  (routine, TRIM(extpar_file))
        stream_id = openInputFile(extpar_file, p_patch(jg), default_read_method)
        CALL read_2D_1time(stream_id, on_cells, 'SST', &
          &          ext_data(jg)%atm_td%sst_m(:,:,2))
        CALL closeFile(stream_id)

        !! READ CI files

        extpar_file = generate_td_filename(ci_td_filename,                    &
          &                                getModelBaseDir(),                 &
          &                                TRIM(p_patch(jg)%grid_filename),   &
          &                                m1,y1                   )
        CALL message  (routine, TRIM(extpar_file))
        stream_id = openInputFile(extpar_file, p_patch(jg), default_read_method)
        CALL read_2D_1time(stream_id, on_cells, 'CI', &
          &          ext_data(jg)%atm_td%fr_ice_m(:,:,1))
        CALL closeFile(stream_id)

        extpar_file = generate_td_filename(ci_td_filename,                    &
          &                                getModelBaseDir(),                 &
          &                             TRIM(p_patch(jg)%grid_filename),      &
          &                             m2,y2                   )
        CALL message  (routine, TRIM(extpar_file))
        stream_id = openInputFile(extpar_file, p_patch(jg), default_read_method)
        CALL read_2D_1time(stream_id, on_cells, 'CI', &
          &          ext_data(jg)%atm_td%fr_ice_m(:,:,2))
        CALL closeFile(stream_id)
      ENDDO  ! jg

  END SUBROUTINE read_td_ext_data_file

!-----------------------------------------------------------------------

END MODULE mo_td_ext_data

