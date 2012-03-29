!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Guenther Zaengl, DWD
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial Revision by Daniel Reinert, DWD, 2012-03-22
!! - some IO-routines, which might be of future use, moved here from 
!!   the outdated output module mo_io_vlist. 
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
MODULE mo_io_util

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message
  USE mo_mpi,                   ONLY: my_process_is_mpi_workroot, my_process_is_stdio, &
    &                                 my_process_is_mpi_test, process_mpi_all_test_id, &
    &                                 process_mpi_all_workroot_id, p_send, p_recv,     &
    &                                 my_process_is_mpi_seq
  USE mo_impl_constants,        ONLY: SUCCESS
  USE mo_model_domain,          ONLY: t_patch
  USE mo_communication,         ONLY: t_comm_pattern, exchange_data
  USE mo_parallel_config,       ONLY: nproma, p_test_run


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gather_array1 
  PUBLIC :: gather_array2

  PUBLIC :: t_outvar_desc
  PUBLIC :: t_collected_var_ptr
  PUBLIC :: outvar_desc

  PUBLIC :: GATHER_C
  PUBLIC :: GATHER_E
  PUBLIC :: GATHER_V

  PUBLIC :: max_outvars
  PUBLIC :: max_gridlevs
  PUBLIC :: num_output_vars


  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  INTEGER, PARAMETER :: GATHER_C      = 1
  INTEGER, PARAMETER :: GATHER_E      = 2
  INTEGER, PARAMETER :: GATHER_V      = 3

  INTEGER, PARAMETER :: max_outvars  = 250 ! max. number of output variables
  INTEGER, PARAMETER :: max_gridlevs = 12  ! max. number of grid levels

  ! Descriptions of output variables
  INTEGER :: num_output_vars(max_gridlevs)

  TYPE t_outvar_desc

    INTEGER ::           TYPE ! GATHER_C, GATHER_E, GATHER_V
    INTEGER ::           nlev
    CHARACTER(LEN=80) :: name

  END TYPE t_outvar_desc

  !> pointer type (avoids clash of INTENT and POINTER attributes)
  TYPE t_collected_var_ptr
    REAL(wp), POINTER :: ptr(:,:,:)
  END TYPE t_collected_var_ptr


  TYPE(t_outvar_desc) :: outvar_desc(max_outvars, max_gridlevs)


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! For output, the fields have to be provided in the form:.
  !!
  !! For output, the fields have to be provided in the form:
  !! field(vector_length,nlevs).
  !! Therefore, the reshaping performed for optimization should be reversed.
  !! This is the version for a horizontal field.
  !!
  !! The optional parameter "out_field_2d" allows the calling procedure
  !! to work with the collected 2d variable (not only with the RESHAPEd
  !! variable "out_field"). However, in this case the caller is in
  !! charge of deallocating this field.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2008-11-19)
  !!
  SUBROUTINE gather_array1(typ,p_patch,in_field,out_field,opt_name,opt_out_field_2d)
    !
    INTEGER,                   INTENT(in)    :: typ
    TYPE(t_patch), TARGET,     INTENT(in)    :: p_patch
    ! 2dimensional input field (nblks,nproma)
    REAL(wp),                  INTENT(in)    :: in_field(:,:)
    ! one vector line version of the input
    REAL(wp),                  INTENT(inout) :: out_field(:)
    CHARACTER(LEN=*), OPTIONAL,INTENT(IN)    :: opt_name
    ! quasi-2d output (nproma,1,nblks)
    TYPE(t_collected_var_ptr), INTENT(OUT), OPTIONAL :: opt_out_field_2d

    REAL(wp), ALLOCATABLE :: out_field2(:,:)

    !-----------------------------------------------------------------------

    IF (my_process_is_stdio()) THEN
      ALLOCATE(out_field2(UBOUND(out_field,1),1))
    ELSE
      ALLOCATE(out_field2(0,0))
    ENDIF
    IF (PRESENT(opt_name)) THEN
      CALL gather_array2(typ,p_patch,&
                         RESHAPE(in_field,(/UBOUND(in_field,1),1,UBOUND(in_field,2)/)), &
                         out_field2, opt_name, opt_out_field_2d)
    ELSE
      CALL gather_array2(typ,p_patch,&
                         RESHAPE(in_field,(/UBOUND(in_field,1),1,UBOUND(in_field,2)/)), &
                         out_field2, opt_out_field_3d=opt_out_field_2d)
    ENDIF

    IF(my_process_is_stdio()) THEN
      out_field(:) = out_field2(:,1)
    ENDIF

    DEALLOCATE(out_field2)

  END SUBROUTINE gather_array1


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! For output, the fields have to be provided in the form:
  !!   field(vector_length,nlevs).
  !! Therefore, the reshaping performed for optimization should be reversed.
  !! This is the version for a 3-d field
  !!
  !! The optional parameter "out_field_3d" allows the calling procedure
  !! to work with the collected 3d variable (not only with the RESHAPEd
  !! variable "out_field"). However, in this case the caller is in
  !! charge of deallocating this field.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M, (2008-11-19)
  !!
  SUBROUTINE gather_array2(typ,p_patch,in_field,out_field, opt_name, opt_out_field_3d)

    INTEGER,     INTENT(in)                  :: typ
    TYPE(t_patch), INTENT(in), TARGET        :: p_patch
    ! 3d input (nproma,nlev,nblks)
    REAL(wp),    INTENT(in)                  :: in_field(:,:,:)
    ! 2d output (length,nlev)
    REAL(wp),    INTENT(inout)               :: out_field(:,:)
    CHARACTER(LEN=*), OPTIONAL,INTENT(IN)    :: opt_name
    ! 3d output (nproma,nlev,nblks)
    TYPE(t_collected_var_ptr), INTENT(OUT), OPTIONAL :: opt_out_field_3d

    ! local variables
    REAL(wp), POINTER             :: tmp_field(:,:,:)
    INTEGER                       :: nblks, npromz, jb, jl, jk, jend, &
      &                              dim1, dim3, ierrstat
    INTEGER                       :: isize_out, isize_lev               ! array size of output
    TYPE(t_comm_pattern), POINTER :: p_comm_pat

    !-----------------------------------------------------------------------

    IF(UBOUND(in_field,1) /= nproma) THEN
      CALL finish('mo_io_vlist/gather_array2','Illegal 1st array dimension')
    ENDIF
!     IF(p_io/=p_test_pe .AND. p_io/=p_work_pe0) THEN ! Safety check only
!       CALL finish('mo_io_vlist/gather_array2','Illegal I/O PE number for this routine')
!     ENDIF

    IF(typ == GATHER_C) THEN

      IF(UBOUND(in_field,3) /= p_patch%nblks_c) &
        CALL finish('mo_io_vlist/gather_array2','Illegal 3rd array dimension')
      dim1 = nproma
      dim3 = (p_patch%n_patch_cells_g-1)/nproma+1

      p_comm_pat => p_patch%comm_pat_gather_c
      nblks      =  p_patch%nblks_c
      npromz     =  p_patch%npromz_c

    ELSE IF(typ == GATHER_E) THEN

      IF(UBOUND(in_field,3) /= p_patch%nblks_e) &
        CALL finish('mo_io_vlist/gather_array2','Illegal 3rd array dimension')
      dim1 = nproma
      dim3 = (p_patch%n_patch_edges_g-1)/nproma+1

      p_comm_pat => p_patch%comm_pat_gather_e
      nblks      =  p_patch%nblks_e
      npromz     =  p_patch%npromz_e

    ELSE IF(typ == GATHER_V) THEN

      IF(UBOUND(in_field,3) /= p_patch%nblks_v) &
        CALL finish('mo_io_vlist/gather_array2','Illegal 3rd array dimension')
      dim1 = nproma
      dim3 = (p_patch%n_patch_verts_g-1)/nproma+1

      p_comm_pat => p_patch%comm_pat_gather_v
      nblks      =  p_patch%nblks_v
      npromz     =  p_patch%npromz_v

    ELSE

      CALL finish('mo_io_vlist/gather_array2','Illegal type parameter')

      ! To get rid of compiler warnings (by gcc) about variables which may be used uninitialized,
      ! define these varaibles also here. They are not used since the "finish" above stops
      ! the model integration.
      p_comm_pat => p_patch%comm_pat_gather_c
      nblks      =  p_patch%nblks_c
      npromz     =  p_patch%npromz_c

    ENDIF

    IF (PRESENT(opt_out_field_3d)) THEN
      ALLOCATE(opt_out_field_3d%ptr(dim1, UBOUND(in_field,2), dim3), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) THEN
        CALL finish ('mo_io_vlist/gather_array2', 'allocation failed')
      ENDIF
      tmp_field => opt_out_field_3d%ptr
    ELSE
      ALLOCATE(tmp_field(dim1, UBOUND(in_field,2), dim3))
    END IF
    tmp_field(:,:,:)=0.0_wp

    IF(p_test_run) THEN
      IF(.NOT. my_process_is_mpi_test()) THEN
        ! Gather all data on process_mpi_all_workroot_id and send it to process_mpi_test_id for verification
        CALL exchange_data(p_comm_pat, RECV=tmp_field, SEND=in_field)
        IF(my_process_is_mpi_workroot()) CALL p_send(tmp_field, process_mpi_all_test_id, 1)
      ELSE
        ! Receive result from parallel worker PEs and check for correctness
        CALL p_recv(tmp_field, process_mpi_all_workroot_id, 1)
        DO jb = 1, nblks
          jend = nproma
          IF(jb==nblks) jend = npromz
          DO jl = 1, jend
            IF(ANY(tmp_field(jl,:,jb) /= in_field(jl,:,jb))) THEN
                IF (PRESENT(opt_name)) THEN
                  WRITE(0,*)'Error ',TRIM(opt_name),jl,jb ,tmp_field(jl,:,jb),in_field(jl,:,jb)
                ELSE
                  WRITE(0,*)'Error ',jl,jb !,tmp_field(jl,:,jb),in_field(jl,:,jb)
               ENDIF
              CALL message('mo_io_vlist/gather_array2','Sync error test PE/worker PEs')
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ELSE
      IF(my_process_is_mpi_seq()) THEN
        ! We are running on 1 PE. Thus, just copy in_field.
        DO jb= 1, nblks
          jend = nproma
          IF(jb==nblks) jend = npromz
          DO jl = 1, jend
            tmp_field(jl,:,jb) = in_field(jl,:,jb)
          ENDDO
        ENDDO
      ELSE
        ! Gather all data on process_mpi_all_workroot_id
        CALL exchange_data(p_comm_pat, RECV=tmp_field, SEND=in_field)
      ENDIF
    ENDIF

    IF(my_process_is_stdio()) THEN
      isize_out = SIZE(out_field,1)
      isize_lev = SIZE(in_field,2)

      DO jk = 1, isize_lev
        out_field(:,jk) = RESHAPE(tmp_field(:,jk,:),(/isize_out/))
      ENDDO
    ENDIF

    IF (.NOT. PRESENT(opt_out_field_3d)) &
      DEALLOCATE(tmp_field)

  END SUBROUTINE gather_array2


END MODULE mo_io_util

