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
MODULE mo_td_ext_data

  USE mo_kind,               ONLY: wp
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_io_units,           ONLY: filename_max
  USE mo_master_nml,         ONLY: model_base_dir
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
#ifdef NOMPI
  USE mo_mpi,                 ONLY: my_process_is_mpi_all_seq
#endif
  USE mo_util_netcdf,         ONLY: read_netcdf_data, nf
  USE mo_datetime,            ONLY: t_datetime, month2hour
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rlcell, min_rlcell_int
  USE mo_grid_config,         ONLY: n_dom
  USE mo_util_string,         ONLY: MAX_STRING_LEN, t_keyword_list,   &
                                 &  associate_keyword, with_keywords
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_parallel_config,     ONLY: p_test_run
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_seaice_nwp,          ONLY: frsi_min

  IMPLICIT NONE

  ! required for reading external data
  INCLUDE 'netcdf.inc'

  PRIVATE

  PUBLIC  :: set_actual_td_ext_data,  read_td_ext_data_file

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


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
  !! Developed  by P. R´ipodas (2012-12)
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
    INTEGER                       :: mpi_comm
    REAL (wp)                     :: pw1, pw2, pw2_old
    INTEGER                       :: jg, jb, jc, i_startidx, i_endidx
    INTEGER                       :: i_nchdom, i_rlstart, i_rlend
    INTEGER                       :: i_startblk, i_endblk

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_td_ext_data:set_actual_td_ext_data  '

    !---------------------------------------------------------------
 
      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      SELECT CASE (ext_data_mode)

       CASE (2) !SST and sea ice fraction updated based 
                !  on the climatological monthly values
        CALL month2hour (datetime, month1, month2, pw2 )
        pw1 = 1._wp - pw2

        DO jg = 1, n_dom
     
           i_nchdom  = MAX(1,p_patch(jg)%n_childdom)
!$OMP PARALLEL PRIVATE(i_startblk,i_endblk, i_rlstart, i_rlend)     
           i_rlstart = 1
           i_rlend   = min_rlcell

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
           i_rlstart = 1
           i_rlend   = min_rlcell

           i_startblk = p_patch(jg)%cells%start_blk(i_rlstart,1)
           i_endblk   = p_patch(jg)%cells%end_blk(i_rlend,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
          DO jb=i_startblk, i_endblk

            CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
               & i_startidx, i_endidx, i_rlstart, i_rlend)

            DO jc = i_startidx, i_endidx
              IF (ext_data(jg)%atm_td%sst_m(jc,jb,1) < 0._wp .AND.   &
                & ext_data(jg)%atm_td%sst_m(jc,jb,2) >= 0._wp )  THEN
                 WRITE( message_text,'(a,2g10.5)') 'something ist wrong,', &
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
    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(LEN=filename_max) :: extpar_td_filename,extpar_file 
    INTEGER                     :: ncid, jg, mpi_comm
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = 'mo_td_ext_data:read_td_ext_data_file:'

!-----------------------------------------------------------------------
   ! extpar_td_filename = "<path>extpar_<year>_<month>_<gridfile>"

   ! Set the months needed to interpolate the ext_para to the actual day

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

      DO jg = 1,n_dom
        IF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid

          CALL finish(TRIM(ROUTINE),&
            & 'Hexagonal grid is not supported, yet.')

        ENDIF
        WRITE( message_text,'(a,2i6,a)') 'n_dom, jg, grid_filename ',n_dom, jg, TRIM(p_patch(jg)%grid_filename)
        CALL message  (routine, TRIM(message_text))
        !! READ SST files

        extpar_td_filename = "<path>SST_<year>_<month>_<gridfile>"
        IF(my_process_is_stdio()) THEN

          extpar_file = generate_td_filename(extpar_td_filename,                &
            &                             model_base_dir,                    &
            &                             TRIM(p_patch(jg)%grid_filename),   &
            &                             y1,m1                   )

          CALL message  (routine, TRIM(extpar_file))
          CALL nf( nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine )

        ENDIF
      
        !
        CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm_td%sst_m(:,:,1))

        IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

        IF(my_process_is_stdio()) THEN

          extpar_file = generate_td_filename(extpar_td_filename,                &
            &                             model_base_dir,                    &
            &                             TRIM(p_patch(jg)%grid_filename),   &
            &                             y2,m2                   )
          CALL message  (routine, TRIM(extpar_file))
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine)

        ENDIF
        CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm_td%sst_m(:,:,2))

        IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)



        !! READ CI files

        extpar_td_filename = "<path>CI_<year>_<month>_<gridfile>"
        IF(my_process_is_stdio()) THEN

          extpar_file = generate_td_filename(extpar_td_filename,                &
            &                             model_base_dir,                    &
            &                             TRIM(p_patch(jg)%grid_filename),   &
            &                             y1,m1                   )
          CALL message  (routine, TRIM(extpar_file))
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine)

        ENDIF
      
        !
        CALL read_netcdf_data (ncid, 'CI', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm_td%fr_ice_m(:,:,1))
        IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)

        IF(my_process_is_stdio()) THEN

          extpar_file = generate_td_filename(extpar_td_filename,                &
            &                             model_base_dir,                    &
            &                             TRIM(p_patch(jg)%grid_filename),   &
            &                             y2,m2                   )
          CALL message  (routine, TRIM(extpar_file))
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid), routine)

        ENDIF
        CALL read_netcdf_data (ncid, 'CI', p_patch(jg)%n_patch_cells_g, &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     ext_data(jg)%atm_td%fr_ice_m(:,:,2))

        IF( my_process_is_stdio()) CALL nf(nf_close(ncid), routine)
      ENDDO  ! jg

  END SUBROUTINE read_td_ext_data_file

!-----------------------------------------------------------------------
  FUNCTION generate_td_filename(extpar_filename, model_base_dir, grid_filename,year,month) &
    &  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)   :: extpar_filename, &
      &                               model_base_dir,  &
      &                               grid_filename
    INTEGER, INTENT(IN)            :: year,month
    CHARACTER(len=MAX_STRING_LEN)  :: syear,smonth
    CHARACTER(len=MAX_STRING_LEN)  :: result_str
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    WRITE(syear, '(i4.4)') year
    WRITE(smonth,'(i2.2)') month

    CALL associate_keyword("<path>",     TRIM(model_base_dir), keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),  keywords)
    CALL associate_keyword("<year>", TRIM(syear),  keywords)
    CALL associate_keyword("<month>", TRIM(smonth),  keywords)
    ! replace keywords in "extpar_filename", which is by default
    ! extpar_filename = "<path>extpar_<year>_<month>_<gridfile>"
    result_str = TRIM(with_keywords(keywords, TRIM(extpar_filename)))

  END FUNCTION generate_td_filename
END MODULE mo_td_ext_data

