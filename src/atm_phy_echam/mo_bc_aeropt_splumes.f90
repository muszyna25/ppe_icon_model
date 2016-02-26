!>
!! @brief Read and apply monthly aerosol optical properties of S. Kinne
!! from yearly files.
!!
!! @author B. Stevens, K. Peters, J.S. Rast (MPI-M)
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_bc_aeropt_splumes

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mo_read_interface,       ONLY: openInputFile, read_1D, &
                                   & read_bcast_real_2D, read_bcast_real_3D, &
                                   & closeFile
!!$, on_cells, &
!!$    &                                t_stream_id, read_0D_real, read_3D_time

  IMPLICIT NONE

  PRIVATE
  PUBLIC                  :: sp_setup

  INTEGER, PARAMETER      ::     &
       nplumes   = 9            ,& !< Number of plumes
       nfeatures = 2            ,& !< Number of features per plume
       ntimes    = 52           ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251             !< Number of years of available forcing
  REAL(wp), POINTER ::                    &
       plume_lat   (:)  ,& !< (nplumes) latitude where plume maximizes
       plume_lon   (:)  ,& !< (nplumes) longitude where plume maximizes
       beta_a      (:)  ,& !< (nplumes) parameter a for beta function 
                           !< vertical profile
       beta_b      (:)  ,& !< (nplumes) parameter b for beta function 
                           !< vertical profile
       aod_spmx    (:)  ,& !< (nplumes) aod at 550 for simple plume (maximum)
       aod_fmbg    (:)  ,& !< (nplumes) aod at 550 for fine mode 
                           !< background (for twomey effect)
       asy550      (:)  ,& !< (nplumes) asymmetry parameter for plume at 550nm
       ssa550      (:)  ,& !< (nplumes) single scattering albedo for 
                           !< plume at 550nm
       angstrom    (:)  ,& !< (nplumes) angstrom parameter for plume 
       sig_lon_E   (:,:),& !< (nfeatures,nplumes) Eastward extent of 
                           !< plume feature
       sig_lon_W   (:,:),& !< (nfeatures,nplumes) Westward extent of 
                           !< plume feature
       sig_lat_E   (:,:),& !< (nfeatures,nplumes) Southward extent of 
                           !< plume feature
       sig_lat_W   (:,:),& !< (nfeatures,nplumes) Northward extent of 
                           !< plume feature
       theta       (:,:),& !< (nfeatures,nplumes) Rotation angle of feature
       ftr_weight  (:,:),& !< (nfeatures,nplumes) Feature weights = 
                           !< (nfeatures + 1) to account for BB background
       time_weight (:,:),& !< (nfeatures,nplumes) Time-weights = 
                           !< (nfeatures +1) to account for BB background
       year_weight (:,:)    ,& !< (nyear,nplumes) Yearly weight for plume
       ann_cycle   (:,:,:)     !< (nfeatures,ntimes,nplumes) annual cycle for feature
  CHARACTER(LEN=256)       :: cfname
  LOGICAL                  :: sp_initialized

  CONTAINS

  ! -----------------------------------------------------------------
  ! SP_SETUP:  This subroutine should be called at initialization to 
  !            read the netcdf data that describes the simple plume
  !            climatology.  The information needs to be either read 
  !            by each processor or distributed to processors.
  !
  SUBROUTINE sp_setup
    !
    ! ---------- 
    !
    INTEGER           :: ifile_id
    CHARACTER(len=32) :: ci_length, cj_length

    cfname='MACv2.0-SP_v1-beta.nc'
    ifile_id=openInputFile(cfname)

    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lat',&
                       & return_pointer=plume_lat, file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='plume_lon',&
                       & return_pointer=plume_lon, file_name=cfname,         &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_a',   &
                       & return_pointer=beta_a, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='beta_b',   &
                       & return_pointer=beta_b, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_spmx', &
                       & return_pointer=aod_spmx, file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='aod_fmbg', &
                       & return_pointer=aod_fmbg, file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='ssa550',   &
                       & return_pointer=ssa550, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='asy550',   &
                       & return_pointer=asy550, file_name=cfname,            &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_1d_wrapper(ifile_id=ifile_id,        variable_name='angstrom', &
                       & return_pointer=angstrom, file_name=cfname,          &
                       & variable_dimls=(/nplumes/),                         &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_W',&
                       & return_pointer=sig_lat_W, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lat_E',&
                       & return_pointer=sig_lat_E, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_W',&
                       & return_pointer=sig_lon_W, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='sig_lon_E',&
                       & return_pointer=sig_lon_E, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='theta',    &
                       & return_pointer=theta, file_name=cfname,             &
                       & variable_dimls=(/nfeatures,nplumes/),               &
                       & module_name='mo_bc_aeropt_splumes',                 &
                       & sub_prog_name='sp_setup'                            )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='ftr_weight',&
                       & return_pointer=ftr_weight, file_name=cfname,         &
                       & variable_dimls=(/nfeatures,nplumes/),                &
                       & module_name='mo_bc_aeropt_splumes',                  &
                       & sub_prog_name='sp_setup'                             )
    CALL read_2d_wrapper(ifile_id=ifile_id,        variable_name='year_weight',&
                       & return_pointer=year_weight, file_name=cfname,         &
                       & variable_dimls=(/nyears,nplumes/),                    &
                       & module_name='mo_bc_aeropt_splumes',                   &
                       & sub_prog_name='sp_setup'                              )
    CALL read_3d_wrapper(ifile_id=ifile_id,        variable_name='ann_cycle', &
                       & return_pointer=ann_cycle, file_name=cfname,          &
                       & variable_dimls=(/nfeatures,ntimes,nplumes/),         &
                       & module_name='mo_bc_aeropt_splumes',                  &
                       & sub_prog_name='sp_setup'                             )
    CALL closeFile(ifile_id)
    sp_initialized = .TRUE.
    RETURN
  END SUBROUTINE SP_SETUP

  SUBROUTINE read_1d_wrapper(ifile_id,                 variable_name,        &
                           & return_pointer,           file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which 
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable 
                                                    !< to be read
    REAL(wp), POINTER,INTENT(out)  :: return_pointer(:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file 
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(1)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module 
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling 
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length, cj_length
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_1D(file_id=ifile_id,         variable_name=variable_name,      &
                 return_pointer=return_pointer                               )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(return_pointer,1)/=variable_dimls(1)) THEN
         WRITE(ci_length,*) SIZE(return_pointer,1)
         WRITE(cj_length,*) variable_dimls(1)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('// &
                & TRIM(ADJUSTL(cj_length))//') has wrong dimension length '// &
                & TRIM(ADJUSTL(ci_length))//' in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_1d_wrapper
  SUBROUTINE read_2d_wrapper(ifile_id,                 variable_name,        &
                           & return_pointer,           file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which 
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable 
                                                    !< to be read
    REAL(wp), POINTER,INTENT(out)  :: return_pointer(:,:) !< values of variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file 
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(2)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module 
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling 
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(2), cj_length(2)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_2D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  return_pointer=return_pointer                    )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(return_pointer,1)/=variable_dimls(1) .OR. &
         & SIZE(return_pointer,2)/=variable_dimls(2)) THEN
         WRITE(ci_length(1),*) SIZE(return_pointer,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(return_pointer,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_2d_wrapper
  SUBROUTINE read_3d_wrapper(ifile_id,                 variable_name,        &
                           & return_pointer,           file_name,            &
                           & variable_dimls,           module_name,          &
                           & sub_prog_name                                   )
    INTEGER, INTENT(in)            :: ifile_id      !< file id from which 
                                                    !< variable is read
    CHARACTER(LEN=*),INTENT(in)    :: variable_name !< name of variable 
                                                    !< to be read
    REAL(wp), POINTER,INTENT(out)  :: return_pointer(:,:,:) !< values of 
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: file_name     !< file name of file 
                                                    !< contain respective var.
    INTEGER, INTENT(in),OPTIONAL   :: variable_dimls(3)!< dimension length of
                                                    !< variable
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: module_name!< name of module 
                                           !< containing calling subprogr.
    CHARACTER(LEN=*),INTENT(in),OPTIONAL:: sub_prog_name!< name of calling 
                                                        !< subprogram
    CHARACTER(LEN=32)                   :: ci_length(3), cj_length(3)
    CHARACTER(LEN=1024)                 :: message1, message2

    CALL read_bcast_REAL_3D(file_id=ifile_id,                                &
                         &  variable_name=variable_name,                     &
                         &  return_pointer=return_pointer                    )
    IF (PRESENT(variable_dimls)) THEN
       IF (SIZE(return_pointer,1)/=variable_dimls(1) .OR. &
         & SIZE(return_pointer,2)/=variable_dimls(2) .OR. &
         & SIZE(return_pointer,3)/=variable_dimls(3)) THEN
         WRITE(ci_length(1),*) SIZE(return_pointer,1)
         WRITE(cj_length(1),*) variable_dimls(1)
         WRITE(ci_length(2),*) SIZE(return_pointer,2)
         WRITE(cj_length(2),*) variable_dimls(2)
         WRITE(ci_length(3),*) SIZE(return_pointer,3)
         WRITE(cj_length(3),*) variable_dimls(3)
         IF (PRESENT(sub_prog_name)) THEN
           message1=TRIM(ADJUSTL(sub_prog_name))//' of'
         ELSE
           message1='Unknown subprogram of '
         END IF
         IF (PRESENT(module_name)) THEN
           message1=TRIM(ADJUSTL(message1))//' '//TRIM(ADJUSTL(module_name))
         ELSE
           message1=TRIM(ADJUSTL(message1))//' unknown module'
         END IF
         message2=TRIM(ADJUSTL(variable_name))//'('//&
                & TRIM(ADJUSTL(cj_length(1)))//','//&
                & TRIM(ADJUSTL(cj_length(2)))//','//&
                & TRIM(ADJUSTL(cj_length(3)))//&
                & ') has wrong dimension length ('//&
                & TRIM(ADJUSTL(ci_length(1)))//','//&
                & TRIM(ADJUSTL(ci_length(2)))//','//&
                & TRIM(ADJUSTL(ci_length(3)))//&
                & ') in'
         IF (PRESENT(file_name)) THEN
           message2=TRIM(ADJUSTL(message2))//' '//TRIM(ADJUSTL(file_name))
         ELSE
           message2=TRIM(ADJUSTL(message2))//' unknown file'
         END IF
         WRITE(0,*) TRIM(ADJUSTL(message1))
         WRITE(0,*) TRIM(ADJUSTL(message2))
         CALL finish(TRIM(ADJUSTL(message1)),TRIM(ADJUSTL(message2)))
       END IF
    END IF
  END SUBROUTINE read_3d_wrapper
  
END MODULE mo_bc_aeropt_splumes
