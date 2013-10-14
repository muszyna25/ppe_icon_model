!>
!!               The module <i>mo_communication</i>.
!!
!!               The module <i>mo_communication</i>
!! provides functionality for boundary exchange and global sums.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
MODULE mo_sync
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------


USE mo_kind,               ONLY: wp, dp, i8
USE mo_exception,          ONLY: finish, message, message_text
USE mo_model_domain,       ONLY: t_patch
USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info
USE mo_parallel_config,    ONLY: nproma
USE mo_run_config,         ONLY: msg_level
USE mo_math_constants,     ONLY: pi
USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int, min_rlvert_int, max_dom
USE mo_impl_constants_grf, ONLY: grf_bdywidth_c, grf_bdywidth_e
USE mo_io_units,           ONLY: find_next_free_unit, filename_max
USE mo_mpi,                ONLY: p_pe, p_bcast, p_sum, p_max, p_min, p_send, p_recv,               &
  &                              p_comm_work_test,  p_comm_work, p_n_work, my_process_is_mpi_test, &
  &                              get_my_mpi_all_id, process_mpi_all_test_id,                       &
  &                              my_process_is_mpi_parallel, p_work_pe0,p_pe_work, push_glob_comm, &
  &                              pop_glob_comm, get_glob_proc0, comm_lev, glob_comm, comm_proc0,   &
  &                              p_gather, p_gatherv
USE mo_parallel_config, ONLY:p_test_run,   &
  & n_ghost_rows, l_log_checks, l_fast_sum
USE mo_communication,      ONLY: exchange_data, exchange_data_4de3,            &
                                 exchange_data_mult, t_comm_pattern,           &
                                 blk_no, idx_no, idx_1d

USE mo_timer,           ONLY: timer_start, timer_stop, activate_sync_timers, &
  & timer_global_sum, timer_omp_global_sum, timer_ordglb_sum, timer_omp_ordglb_sum
 

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

!modules interface-------------------------------------------
!subroutines
PUBLIC :: sync_patch_array, check_patch_array, sync_idx,              &
          global_sum_array, omp_global_sum_array,                     &
          global_sum_array2, global_sum_array3,                       &
          sync_patch_array_mult, global_min, global_max,              &
          sync_patch_array_4de3, decomposition_statistics,            &
          enable_sync_checks, disable_sync_checks,                    &
          cumulative_sync_patch_array, complete_cumulative_sync

!
!variables

INTEGER, PARAMETER, PUBLIC :: SYNC_C = 1
INTEGER, PARAMETER, PUBLIC :: SYNC_E = 2
INTEGER, PARAMETER, PUBLIC :: SYNC_V = 3
INTEGER, PARAMETER, PUBLIC :: SYNC_C1 = 4

INTERFACE sync_patch_array
  MODULE PROCEDURE sync_patch_array_r2
  MODULE PROCEDURE sync_patch_array_r3
  MODULE PROCEDURE sync_patch_array_i2
  MODULE PROCEDURE sync_patch_array_i3
END INTERFACE

INTERFACE check_patch_array
  MODULE PROCEDURE check_patch_array_2
  MODULE PROCEDURE check_patch_array_3
  MODULE PROCEDURE check_patch_array_4
END INTERFACE

INTERFACE global_min
  MODULE PROCEDURE global_min_0d
  MODULE PROCEDURE global_min_1d
END INTERFACE

INTERFACE global_max
  MODULE PROCEDURE global_max_0d
  MODULE PROCEDURE global_max_1d
END INTERFACE

INTERFACE global_sum_array
  MODULE PROCEDURE global_sum_array_0d
  MODULE PROCEDURE global_sum_array_0di
  MODULE PROCEDURE global_sum_array_1d
  MODULE PROCEDURE global_sum_array_2d
  MODULE PROCEDURE global_sum_array_3d
END INTERFACE

INTERFACE omp_global_sum_array
  MODULE PROCEDURE omp_global_sum_array_1d
  MODULE PROCEDURE omp_global_sum_array_2d
END INTERFACE


! Unit for logging sync errors
INTEGER, SAVE :: log_unit = -1

! Flag if sync checks are enabled when sync_patch_array et al is called
LOGICAL, SAVE :: do_sync_checks = .TRUE.

!-------------------------------------------------------------------------

!> Type definition for "cumulative syncs": These are boundary
!  exchanges for 3D cell-based fields, collecting as many fields as
!  possible before actually performing the sync.
TYPE t_cumulative_sync
  REAL(wp),      POINTER :: f3d(:,:,:)
  TYPE(t_patch), POINTER :: p_patch
END TYPE t_cumulative_sync

!> max. no. fields that can be handled by "sync_patch_array_mult"
INTEGER, PARAMETER :: MAX_CUMULATIVE_SYNC = 5
!> No. of pending cumulative sync; dim (typ,patch_id,fieldno)
INTEGER :: ncumul_sync(4,max_dom) = 0
!> List of cumulative sync fields:
TYPE(t_cumulative_sync) :: cumul_sync(4,max_dom,MAX_CUMULATIVE_SYNC)


CONTAINS


!-------------------------------------------------------------------------
!
!> enable_sync_checks: Enables sync checks when sync_patch_array et al is called

SUBROUTINE enable_sync_checks()
  do_sync_checks = .TRUE.
END SUBROUTINE enable_sync_checks
!-------------------------------------------------------------------------
!
!> disable_sync_checks: Disables sync checks when sync_patch_array et al is called

SUBROUTINE disable_sync_checks()
  do_sync_checks = .FALSE.
END SUBROUTINE disable_sync_checks


!-------------------------------------------------------------------------
!> Does boundary exchange for a 3-D REAL array.
!
!  @par Revision History
!  Initial version by Rainer Johanni, Nov 2009
!
SUBROUTINE sync_patch_array_r3(typ, p_patch, arr)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), INTENT(IN)    :: p_patch
   REAL(wp),      INTENT(INOUT) :: arr(:,:,:)

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) CALL check_patch_array_3(typ, p_patch, arr, 'sync')

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
      IF(typ == SYNC_C) THEN
         CALL exchange_data(p_patch%comm_pat_c, arr)
      ELSE IF(typ == SYNC_E) THEN
         CALL exchange_data(p_patch%comm_pat_e, arr)
      ELSE IF(typ == SYNC_V) THEN
         CALL exchange_data(p_patch%comm_pat_v, arr)
      ELSE IF(typ == SYNC_C1) THEN
         CALL exchange_data(p_patch%comm_pat_c1, arr)
      ELSE
         CALL finish('sync_patch_array','Illegal type parameter')
      ENDIF
   ENDIF
END SUBROUTINE sync_patch_array_r3


!-------------------------------------------------------------------------
!> Does boundary exchange for a 3-D INTEGER array.
!
!  @par Revision History
!  Initial version by Rainer Johanni, Nov 2009
!
!  @note This implementation does not perform a consistency check
!        (p_test_run)!
!
SUBROUTINE sync_patch_array_i3(typ, p_patch, arr)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), INTENT(IN)    :: p_patch
   INTEGER,       INTENT(INOUT) :: arr(:,:,:)

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
      IF(typ == SYNC_C) THEN
         CALL exchange_data(p_patch%comm_pat_c, arr)
      ELSE IF(typ == SYNC_E) THEN
         CALL exchange_data(p_patch%comm_pat_e, arr)
      ELSE IF(typ == SYNC_V) THEN
         CALL exchange_data(p_patch%comm_pat_v, arr)
      ELSE IF(typ == SYNC_C1) THEN
         CALL exchange_data(p_patch%comm_pat_c1, arr)
      ELSE
         CALL finish('sync_patch_array','Illegal type parameter')
      ENDIF
   ENDIF
END SUBROUTINE sync_patch_array_i3


!-------------------------------------------------------------------------
!> Does boundary exchange for a 2-D REAL array.
!
!  @par Revision History
!  Initial version by Rainer Johanni, Nov 2009
!
SUBROUTINE sync_patch_array_r2(typ, p_patch, arr)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), INTENT(IN)    :: p_patch
   REAL(wp),      INTENT(INOUT) :: arr(:,:)
   ! local variable
   REAL(wp), ALLOCATABLE :: arr3(:,:,:)

   ALLOCATE(arr3(UBOUND(arr,1), 1, UBOUND(arr,2)))
   arr3(:,1,:) = arr(:,:)
   CALL sync_patch_array_r3(typ, p_patch, arr3)
   arr(:,:) = arr3(:,1,:)
   DEALLOCATE(arr3)
END SUBROUTINE sync_patch_array_r2


!-------------------------------------------------------------------------
!> Does boundary exchange for a 2-D INTEGER array.
!
!  @par Revision History
!  Initial version by Rainer Johanni, Nov 2009
!
!  @note This implementation does not perform a consistency check
!        (p_test_run)!
!
SUBROUTINE sync_patch_array_i2(typ, p_patch, arr)
   INTEGER,       INTENT(IN)    :: typ
   TYPE(t_patch), INTENT(IN)    :: p_patch
   INTEGER,       INTENT(INOUT) :: arr(:,:)
   ! local variable
   INTEGER, ALLOCATABLE :: arr3(:,:,:)

   ALLOCATE(arr3(UBOUND(arr,1), 1, UBOUND(arr,2)))
   arr3(:,1,:) = arr(:,:)
   CALL sync_patch_array_i3(typ, p_patch, arr3)
   arr(:,:) = arr3(:,1,:)
   DEALLOCATE(arr3)
END SUBROUTINE sync_patch_array_i2


!-------------------------------------------------------------------------
!>
!! Does boundary exchange for up to 5 3D cell-based fields and/or a 4D field.
!!
!! @par Revision History
!! Optimized version by Guenther Zaengl, Apr 2010, based on routines
!! developed by Rainer Johanni
!!
SUBROUTINE sync_patch_array_mult(typ, p_patch, nfields, f3din1, f3din2, f3din3, &
                                 f3din4, f3din5, f4din)

   INTEGER, INTENT(IN)             :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   INTEGER,     INTENT(IN)         :: nfields

   REAL(wp), OPTIONAL, INTENT(INOUT) ::  f3din1(:,:,:), f3din2(:,:,:), f3din3(:,:,:), &
                                         f3din4(:,:,:), f3din5(:,:,:), f4din(:,:,:,:)

   REAL(wp), ALLOCATABLE :: arr3(:,:,:)
   TYPE(t_comm_pattern), POINTER :: p_pat
   INTEGER :: i
   INTEGER :: ndim2tot ! Sum of second dimensions over all input fields

!-----------------------------------------------------------------------

    IF(typ == SYNC_C) THEN
      p_pat => p_patch%comm_pat_c
    ELSE IF(typ == SYNC_E) THEN
      p_pat => p_patch%comm_pat_e
    ELSE IF(typ == SYNC_V) THEN
      p_pat => p_patch%comm_pat_v
    ELSE IF(typ == SYNC_C1) THEN
      p_pat => p_patch%comm_pat_c1
    ENDIF

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
     IF (PRESENT(f4din)) THEN
       ALLOCATE(arr3(UBOUND(f4din,1), UBOUND(f4din,2), UBOUND(f4din,3)))
       DO i = 1, SIZE(f4din,4)
         arr3(:,:,:) = f4din(:,:,:,i)
         CALL check_patch_array_3(typ, p_patch, arr3, 'sync')
       ENDDO
       DEALLOCATE(arr3)
     ENDIF
     IF (PRESENT(f3din1)) CALL check_patch_array_3(typ, p_patch, f3din1, 'sync')
     IF (PRESENT(f3din2)) CALL check_patch_array_3(typ, p_patch, f3din2, 'sync')
     IF (PRESENT(f3din3)) CALL check_patch_array_3(typ, p_patch, f3din3, 'sync')
     IF (PRESENT(f3din4)) CALL check_patch_array_3(typ, p_patch, f3din4, 'sync')
     IF (PRESENT(f3din5)) CALL check_patch_array_3(typ, p_patch, f3din5, 'sync')
   ENDIF

   ! Boundary exchange for work PEs
   IF(my_process_is_mpi_parallel()) THEN
     IF (PRESENT(f4din)) THEN
       ndim2tot = SIZE(f4din,4)*SIZE(f4din,2)
     ELSE
       ndim2tot = 0
     ENDIF
     IF (PRESENT(f3din1)) ndim2tot = ndim2tot+SIZE(f3din1,2)
     IF (PRESENT(f3din2)) ndim2tot = ndim2tot+SIZE(f3din2,2)
     IF (PRESENT(f3din3)) ndim2tot = ndim2tot+SIZE(f3din3,2)
     IF (PRESENT(f3din4)) ndim2tot = ndim2tot+SIZE(f3din4,2)
     IF (PRESENT(f3din5)) ndim2tot = ndim2tot+SIZE(f3din5,2)

     CALL exchange_data_mult(p_pat, nfields, ndim2tot, recv1=f3din1, recv2=f3din2, &
                             recv3=f3din3, recv4=f3din4, recv5=f3din5, recv4d=f4din)
   ENDIF

END SUBROUTINE sync_patch_array_mult

!>
!! Does boundary exchange for a 4D field for which the extra dimension
!! is on the third index.
!!
!! @par Revision History
!! Optimized version by Guenther Zaengl, Apr 2010, based on routines
!! developed by Rainer Johanni
!!
SUBROUTINE sync_patch_array_4de3(typ, p_patch, nfields, f4din)

   INTEGER, INTENT(IN)             :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   INTEGER,     INTENT(IN)         :: nfields

   REAL(wp), INTENT(INOUT) :: f4din(:,:,:,:)

   REAL(wp), ALLOCATABLE :: arr3(:,:,:)
   TYPE(t_comm_pattern), POINTER :: p_pat
   INTEGER :: i, ndim2tot

!-----------------------------------------------------------------------

    IF(typ == SYNC_C) THEN
      p_pat => p_patch%comm_pat_c
    ELSE IF(typ == SYNC_E) THEN
      p_pat => p_patch%comm_pat_e
    ELSE IF(typ == SYNC_V) THEN
      p_pat => p_patch%comm_pat_v
    ELSE IF(typ == SYNC_C1) THEN
      p_pat => p_patch%comm_pat_c1
    ENDIF

   ! If this is a verification run, check consistency before doing boundary exchange
   IF (p_test_run .AND. do_sync_checks) THEN
     ALLOCATE(arr3(UBOUND(f4din,1), UBOUND(f4din,2), UBOUND(f4din,4)))
     DO i = 1, nfields
       arr3(:,:,:) = f4din(:,:,i,:)
       CALL check_patch_array_3(typ, p_patch, arr3, 'sync')
     ENDDO
     DEALLOCATE(arr3)
   ENDIF

   ! Boundary exchange for work PEs
!   IF(p_nprocs /= 1 .AND. p_pe /= p_test_pe) THEN
     IF(my_process_is_mpi_parallel()) THEN
     IF (nfields/=UBOUND(f4din,3)) &
       CALL finish('sync_patch_array_4de3','inconsistent arguments')
     ndim2tot = nfields*SIZE(f4din,2)
     CALL exchange_data_4de3(p_pat, nfields, ndim2tot, recv=f4din)
   ENDIF

END SUBROUTINE sync_patch_array_4de3


!-------------------------------------------------------------------------
!
!

!>
!! In a verification run, this routine checks the consistency of an array, i.
!!
!! e. if the parts on the worker PEs are identical with the data on
!! the verification PE.
!! For a non-verification run it just does nothing.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE check_patch_array_3(typ, p_patch, arr, opt_varname)

!

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   CHARACTER*(*), INTENT(IN), OPTIONAL :: opt_varname

!  Please note: this is also an input parameter in reality, but for using it
!  as an argument of p_bcast it has to be declared INTENT(INOUT)
   REAL(wp), INTENT(INOUT) :: arr(:,:,:)

   REAL(wp), ALLOCATABLE:: arr_g(:,:,:)
   INTEGER :: j, jb, jl, jb_g, jl_g, n, ndim2, ndim3, nblks_g, flag
   INTEGER :: ityp, ndim, ndim_g
   INTEGER :: nerr(0:n_ghost_rows)
   INTEGER, POINTER :: p_glb_index(:)=>NULL(), p_decomp_domain(:,:)=>NULL()

   CHARACTER(len=256) :: varname, cfmt

   CHARACTER(filename_max) :: log_file
   REAL(wp) :: absmax
   LOGICAL :: sync_error

   ityp   = -1
   ndim   = -1
   ndim_g = -1

!-----------------------------------------------------------------------

   IF(.NOT. p_test_run) RETURN ! This routine is only effective in a verification run

   IF(PRESENT(opt_varname)) THEN
      varname = opt_varname
   ELSE
      varname = ''
   ENDIF

   ! Check dimensions of arr, determine if this is an cell/edge/vert array

   IF(UBOUND(arr,1) /= nproma) THEN
      CALL finish('sync_patch_array','first dimension /= nproma')
   ENDIF

   ndim2 = UBOUND(arr,2)
   ndim3 = UBOUND(arr,3)

   IF(typ == SYNC_C .OR. typ == SYNC_C1) THEN
      ndim   = p_patch%n_patch_cells
      ndim_g = p_patch%n_patch_cells_g
      p_glb_index => p_patch%cells%decomp_info%glb_index
      p_decomp_domain => p_patch%cells%decomp_info%decomp_domain
      ityp = typ
   ELSE IF(typ == SYNC_E) THEN
      ndim   = p_patch%n_patch_edges
      ndim_g = p_patch%n_patch_edges_g
      p_glb_index => p_patch%edges%decomp_info%glb_index
      p_decomp_domain => p_patch%edges%decomp_info%decomp_domain
      ityp = typ
   ELSE IF(typ == SYNC_V) THEN
      ndim   = p_patch%n_patch_verts
      ndim_g = p_patch%n_patch_verts_g
      p_glb_index => p_patch%verts%decomp_info%glb_index
      p_decomp_domain => p_patch%verts%decomp_info%decomp_domain
      ityp = typ
   ELSE IF(typ == 0) THEN
      ! typ == 0 may be set for quick checks without knowing the type of the array.
      ! It may only be used if the array is correctly dimensioned.

      IF(ndim3 == p_patch%nblks_c) THEN
         ndim   = p_patch%n_patch_cells
         ndim_g = p_patch%n_patch_cells_g
         p_glb_index => p_patch%cells%decomp_info%glb_index
         p_decomp_domain => p_patch%cells%decomp_info%decomp_domain
         ityp = SYNC_C
      ELSE IF(ndim3 == p_patch%nblks_e) THEN
         ndim   = p_patch%n_patch_edges
         ndim_g = p_patch%n_patch_edges_g
         p_glb_index => p_patch%edges%decomp_info%glb_index
         p_decomp_domain => p_patch%edges%decomp_info%decomp_domain
         ityp = SYNC_E
      ELSE IF(ndim3 == p_patch%nblks_v) THEN
         ndim   = p_patch%n_patch_verts
         ndim_g = p_patch%n_patch_verts_g
         p_glb_index => p_patch%verts%decomp_info%glb_index
         p_decomp_domain => p_patch%verts%decomp_info%decomp_domain
         ityp = SYNC_V
      ELSE
         CALL finish('check_patch_array','typ==0 but unknown blocksize of array')
      ENDIF
   ELSE
      CALL finish('sync_patch_array','Illegal type parameter')
   ENDIF

   ! Actually do the check.
   ! The test PE broadcasts its full array, the other check if their section matches.

   nblks_g = (ndim_g-1)/nproma+1

   IF(get_my_mpi_all_id() == process_mpi_all_test_id) THEN

     ! the test PE may also have reordered global indices. create a
     ! temporary array in the correct order:
     ALLOCATE(arr_g(nproma,ndim2,nblks_g))
     DO j=1,ndim
         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch

         jb_g = blk_no(p_glb_index(j)) ! Block index in global patch
         jl_g = idx_no(p_glb_index(j)) ! Line  index in global patch
         arr_g(jl_g,1:ndim2,jb_g) = arr(jl,1:ndim2,jb)
     END DO

     IF(comm_lev==0) THEN
       CALL p_bcast(arr_g(:,:,1:nblks_g), process_mpi_all_test_id, comm=p_comm_work_test)
     ELSE
       CALL p_send(arr_g(:,:,1:nblks_g),comm_proc0(comm_lev)+p_work_pe0,1)
     ENDIF
     DEALLOCATE(arr_g)

   ELSE

      ALLOCATE(arr_g(nproma,ndim2,nblks_g))
      IF(comm_lev==0) THEN
         CALL p_bcast(arr_g(:,:,1:nblks_g), process_mpi_all_test_id, comm=p_comm_work_test)
      ELSE
         IF(p_pe_work==comm_proc0(comm_lev)) &
           & CALL p_recv(arr_g(:,:,1:nblks_g), process_mpi_all_test_id, 1)
         CALL p_bcast(arr_g(:,:,1:nblks_g),0,comm=glob_comm(comm_lev))
      ENDIF

      ! Count errors in the inner domain and the different ghost rows

      nerr(:) = 0
      absmax = 0.0_wp
      sync_error = .FALSE.

      DO j = 1, ndim

         jb = blk_no(j) ! Block index in distributed patch
         jl = idx_no(j) ! Line  index in distributed patch

         jb_g = blk_no(p_glb_index(j)) ! Block index in global patch
         jl_g = idx_no(p_glb_index(j)) ! Line  index in global patch

         flag = p_decomp_domain(jl,jb)

         ! Safety measures only:
         flag = MAX(flag,0)
         flag = MIN(flag,UBOUND(nerr,1))

         DO n=1,ndim2
            IF(arr(jl,n,jb) /= arr_g(jl_g,n,jb_g)) THEN
               nerr(flag) = nerr(flag)+1
               IF(flag==0) THEN
                  ! Real sync error detected
                  sync_error = .TRUE.
                  absmax = MAX(absmax,ABS(arr(jl,n,jb) - arr_g(jl_g,n,jb_g)))
                  IF (l_log_checks) &
                     WRITE(log_unit,'(a,5i7,3e18.10)') 'sync error location:',&
                       jb,jl,jb_g,jl_g,n,arr(jl,n,jb),arr_g(jl_g,n,jb_g),    &
                       ABS(arr(jl,n,jb)-arr_g(jl_g,n,jb_g))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF(l_log_checks) THEN

         IF(log_unit<0) THEN
            WRITE(log_file,'(''log'',i4.4,''.txt'')') p_pe
            log_unit = find_next_free_unit(10,99)
            OPEN(log_unit,FILE=log_file)
         ENDIF

         n = n_ghost_rows
         WRITE(cfmt,'(a,i3,a)') '(',n+1,'i8,'' '',a)'

         IF(ALL(arr == 0.0_wp)) THEN
            WRITE(log_unit,cfmt) nerr(0:n),TRIM(varname)//': ALL 0 !!!'
         ELSE
            WRITE(log_unit,cfmt) nerr(0:n),TRIM(varname)
         ENDIF
         IF(absmax > 0.0_wp) WRITE(log_unit,*) 'Max abs inner err:',absmax
      ENDIF

      ! Terminate the programm if the array is out of sync

      IF(sync_error) THEN
        IF(l_log_checks) THEN
          CLOSE (log_unit)
        ENDIF
        CALL finish('sync_patch_array','Out of sync detected!')
      ENDIF

      DEALLOCATE(arr_g)

   ENDIF

END SUBROUTINE check_patch_array_3
!-------------------------------------------------------------------------
!
!

!>
!! 2-D Interface to check_patch_array.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE check_patch_array_2(typ, p_patch, arr, opt_varname)

!

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   REAL(wp), INTENT(IN)    :: arr(:,:)
   CHARACTER*(*), INTENT(IN), OPTIONAL :: opt_varname

   REAL(wp), ALLOCATABLE :: arr3(:,:,:)
!-----------------------------------------------------------------------

   IF(.NOT. p_test_run) RETURN ! This routine is only effective in a verification run

   ALLOCATE(arr3(UBOUND(arr,1), 1, UBOUND(arr,2)))
   arr3(:,1,:) = arr(:,:)

   IF(PRESENT(opt_varname)) THEN
      CALL check_patch_array_3(typ, p_patch, arr3, opt_varname)
   ELSE
      CALL check_patch_array_3(typ, p_patch, arr3)
   ENDIF

   DEALLOCATE(arr3) ! NB: No back-copy here!

END SUBROUTINE check_patch_array_2
!-------------------------------------------------------------------------
!
!

!>
!! 4-D Interface to check_patch_array.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
SUBROUTINE check_patch_array_4(typ, p_patch, arr, opt_varname)

!

   INTEGER, INTENT(IN)     :: typ
   TYPE(t_patch), INTENT(IN), TARGET :: p_patch
   REAL(wp), INTENT(INOUT) :: arr(:,:,:,:)
   CHARACTER(len=*), INTENT(IN), OPTIONAL :: opt_varname

   CHARACTER(len=256) :: new_var
   INTEGER :: jt
!-----------------------------------------------------------------------

   IF(.NOT. p_test_run) RETURN ! This routine is only effective in a verification run

   DO jt=1,UBOUND(arr,4)
      IF(PRESENT(opt_varname)) THEN
         WRITE(new_var,'(a,''['',i2,'']'')') opt_varname, jt
         CALL check_patch_array_3(typ, p_patch, arr(:,:,:,jt),TRIM(new_var))
      ELSE
         CALL check_patch_array_3(typ, p_patch, arr(:,:,:,jt))
      ENDIF
   ENDDO

END SUBROUTINE check_patch_array_4
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!>
!! Syncs an idx/blk pair of arrays
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Oct 2011

SUBROUTINE sync_idx(type_arr, type_idx, p_patch, idx, blk, opt_remap)

  INTEGER, INTENT(IN) :: type_arr, type_idx
  TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
  INTEGER, INTENT(INOUT) :: idx(:,:), blk(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: opt_remap

  INTEGER :: nblks, n_idx, n_idx_g, jb, jl, i_l, i_g
  LOGICAL :: remap
  REAL(wp), ALLOCATABLE :: z_idx(:,:)
  TYPE(t_grid_domain_decomp_info), POINTER :: decomp_info

  ! opt_remap: Flag if index values pointing outside local domain
  ! should be remapped to values within local domain

  IF(PRESENT(opt_remap)) THEN
    remap = opt_remap
  ELSE
    remap = .TRUE.
  ENDIF

  IF(type_arr == SYNC_C) THEN
    nblks = p_patch%nblks_c
  ELSEIF(type_arr == SYNC_E) THEN
    nblks = p_patch%nblks_e
  ELSEIF(type_arr == SYNC_V) THEN
    nblks = p_patch%nblks_v
  ELSE
    CALL finish('sync_idx','Unsupported type_arr')
  ENDIF

  IF(type_idx == SYNC_C) THEN
    decomp_info => p_patch%cells%decomp_info
    n_idx = p_patch%n_patch_cells
    n_idx_g = p_patch%n_patch_cells_g
  ELSEIF(type_idx == SYNC_E) THEN
    decomp_info => p_patch%edges%decomp_info
    n_idx = p_patch%n_patch_edges
    n_idx_g = p_patch%n_patch_edges_g
  ELSEIF(type_idx == SYNC_V) THEN
    decomp_info => p_patch%verts%decomp_info
    n_idx = p_patch%n_patch_verts
    n_idx_g = p_patch%n_patch_verts_g
  ELSE
    CALL finish('sync_idx','Unsupported type_idx')
  ENDIF


  if(ubound(idx,1) /= nproma) CALL finish('sync_idx','ubound(idx,1) /= nproma')
  if(ubound(idx,2) /= nblks)  CALL finish('sync_idx','ubound(idx,2) /= nblks')

  ALLOCATE(z_idx(nproma,nblks))
  z_idx = 0._wp

  ! Set z_idx with the global 1D-index of all points

  DO jb = 1, nblks
    DO jl = 1, nproma

      i_l = idx_1d(idx(jl,jb),blk(jl,jb))

      IF(i_l <= 0 .or. i_l > n_idx) THEN
        z_idx(jl,jb) = 0._wp
      ELSE
        z_idx(jl,jb) = REAL(decomp_info%glb_index(i_l),wp)
      ENDIF

    END DO
  END DO

  ! Sync z_idx
  CALL sync_patch_array(type_arr, p_patch, z_idx)

  ! Set all points with local index corresponding to z_idx
  DO jb = 1, nblks
    DO jl = 1, nproma

      i_g = INT(z_idx(jl,jb))

      IF(i_g <= 0 .or. i_g > n_idx_g) THEN
        idx(jl,jb) = idx_no(0)
        blk(jl,jb) = blk_no(0)
      ELSE
        i_l = decomp_info%loc_index(i_g)
        ! Determine what to do with nonlocal values
        IF(i_l<0) THEN
          IF(remap) THEN
            ! Remap index to positive value like in mo_subdivision/remap_index
            i_l = MAX(ABS(i_l)-1,1)
          ELSE
            ! Set it to negative global index like in mo_subdivision/get_local_idx_blk
            i_l = -i_g
          ENDIF
        ENDIF
        idx(jl,jb) = idx_no(i_l)
        blk(jl,jb) = blk_no(i_l)
      ENDIF

    END DO
  END DO

END SUBROUTINE sync_idx

!-------------------------------------------------------------------------
!>
!! Calculates the global sum of an integer scalar.
!! This routine should be called outside an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION global_sum_array_0di (zfield) RESULT (global_sum)

  INTEGER,           INTENT(in) :: zfield
  INTEGER                       :: global_sum
  REAL(wp)                      :: z_aux, z_auxs

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------
  IF (activate_sync_timers) CALL timer_start(timer_global_sum)

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  z_aux =  REAL(zfield,wp)
  z_auxs = p_sum(z_aux, comm=p_comm_glob)

  global_sum = NINT(z_auxs)
  
  IF (activate_sync_timers) CALL timer_stop(timer_global_sum)

END FUNCTION global_sum_array_0di

!-------------------------------------------------------------------------
!>
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine shuold be called outside an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION global_sum_array_0d (zfield) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1), z_aux(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  z_aux(1) = zfield

  IF(l_fast_sum) THEN
    global_sum = simple_sum(z_aux, SIZE(z_aux), p_comm_glob)
  ELSE
    global_sum = order_insensit_ieee64_sum(z_aux, SIZE(z_aux), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF


END FUNCTION global_sum_array_0d

!-------------------------------------------------------------------------
!>
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine shuold be called outside an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION global_sum_array_1d (zfield) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF
  

END FUNCTION global_sum_array_1d

!-------------------------------------------------------------------------
!>
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine shuold be called outside an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION global_sum_array_2d (zfield) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:, :)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF
  

END FUNCTION global_sum_array_2d

!-------------------------------------------------------------------------
!>
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine shuold be called outside an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION global_sum_array_3d (zfield) RESULT (global_sum)

  REAL(wp),          INTENT(in) :: zfield(:,:,:)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
  ENDIF

END FUNCTION global_sum_array_3d
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!>
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine shuold be called from within an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION omp_global_sum_array_1d (zfield) RESULT (global_sum)

!
  REAL(wp),          INTENT(in) :: zfield(:)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = omp_simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = omp_order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
!$OMP BARRIER
!$OMP MASTER
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  ENDIF

END FUNCTION omp_global_sum_array_1d

!-------------------------------------------------------------------------
!>
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine shuold be called from within an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION omp_global_sum_array_2d (zfield) RESULT (global_sum)

!
  REAL(wp),          INTENT(in) :: zfield(:, :)
  REAL(wp)                      :: global_sum
  REAL(wp)                      :: sum_on_testpe(1)

  INTEGER :: p_comm_glob
!-----------------------------------------------------------------------

  IF(comm_lev==0) THEN
    p_comm_glob = p_comm_work
  ELSE
    p_comm_glob = glob_comm(comm_lev)
  ENDIF

  IF(l_fast_sum) THEN
    global_sum = omp_simple_sum(zfield, SIZE(zfield), p_comm_glob)
  ELSE
    global_sum = omp_order_insensit_ieee64_sum(zfield, SIZE(zfield), p_comm_glob)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) THEN
!$OMP BARRIER
!$OMP MASTER
    IF(l_fast_sum) THEN
      CALL check_result( (/ global_sum /), 'global_sum_array', sum_on_testpe)
      global_sum = sum_on_testpe(1)
    ELSE
      CALL check_result( (/ global_sum /), 'global_sum_array')
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
  ENDIF

END FUNCTION omp_global_sum_array_2d


!-------------------------------------------------------------------------
!
!

!>
!! Calculates the global sum of zfield and checks for consistency.
!!
!! Calculates the global sum of zfield and checks for consistency
!! when doing a verification run.
!! This routine has to be called from outside an OMP parallel Region!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! NEC-optimized version by Guenther Zaengl, Dec 2009
!!
FUNCTION global_sum_array2 (zfield) RESULT (global_sum)

!
   REAL(wp), INTENT(in) :: zfield(:,:)
   REAL(wp)             :: global_sum

   INTEGER(i8)       :: itmp(2), isum(2), ival1, ival2
   INTEGER(i8)       :: ilsum(SIZE(zfield,2),2)
   INTEGER           :: i, j, iexp, nblks, nlen
   REAL(dp)          :: abs_max
   REAL(dp)          :: fact, r_fact, rval
   REAL(dp)          :: max_aux(SIZE(zfield,2))

   REAL(dp), PARAMETER :: two_40 = 1099511627776._dp ! 2.**40
   REAL(dp), PARAMETER :: r_two_40 = 1._dp/two_40

#if defined (__SX__) || defined (__PGI)
   INTEGER(i8) :: mask40
   DATA mask40 / z'000000ffffffffff' / ! last 40 bits set
#else
   INTEGER(i8), PARAMETER :: mask40 = INT(z'000000ffffffffff',i8)
#endif
   INTEGER :: p_comm_glob

!-----------------------------------------------------------------------

   ilsum = 0_i8
   nblks = SIZE(zfield,2)
   nlen  = SIZE(zfield,1)

   IF(comm_lev==0) THEN
     p_comm_glob = p_comm_work
   ELSE
     p_comm_glob = glob_comm(comm_lev)
   ENDIF

   ! Get the maximum absolute value of all numbers.
!$OMP PARALLEL PRIVATE(rval)
!$OMP DO PRIVATE(j)
   DO j=1,nblks
      max_aux(j) = MAXVAL(ABS(zfield(:,j)))
   ENDDO
!$OMP END DO

!$OMP MASTER
   rval = MAXVAL(max_aux)
   abs_max = p_max(rval, comm=p_comm_glob)

   ! Get the exponent of abs_max for scaling
   iexp = EXPONENT(abs_max)

   ! Calculate a factor for scaling the input numbers
   ! so that the maximum absolute value of a scaled number
   ! is below 2**40

   fact = SCALE(1._dp,40-iexp) ! same as 2**(40-iexp)
   r_fact = SCALE(1._dp,iexp-40) ! 1./fact

!$OMP END MASTER
!$OMP BARRIER

   ! Sum up all numbers as scaled integers

!$OMP DO PRIVATE(i, j, rval, ival1, ival2)
   DO j=1,nblks
     DO i=1,nlen
       ! Scale number into range -2**40 < rval < 2**40
       ! and store integer part in ival1

       rval = zfield(i,j)*fact
       ival1 = INT(rval,i8)

       ! Scale fraction by 2**40 and store integer part in ival2

       ival2 = INT((rval - REAL(ival1,dp))*two_40,i8)

       ! Sum up ival1 and ival2; since we are using 8-byte integers
       ! for the sum there are no problems with overflow

       ilsum(j,1) = ilsum(j,1) + ival1
       ilsum(j,2) = ilsum(j,2) + ival2
     ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   ! If the exponent is too small, return 0 in order to avoid
   ! problems with overflow below

   IF(iexp < -980) THEN
      global_sum = 0._dp
      RETURN
   ENDIF

   itmp(1) = SUM(ilsum(:,1))
   itmp(2) = SUM(ilsum(:,2))
   isum = p_sum(itmp, comm=p_comm_glob)

   ! Scale integer numbers back to real numbers and add them.
   ! For safety, we use only positive INTEGERS < 2**40 when converting to REAL

    IF (isum(1) >= 0_i8) THEN
       ival1 = ISHFT(isum(1),-40)
       ival2 = IAND (isum(1),mask40)
       global_sum = (REAL(ival1,dp)*r_fact)*two_40 + REAL(ival2,dp)*r_fact
    ELSE
       ival1 = ISHFT(ABS(isum(1)),-40)
       ival2 = IAND (ABS(isum(1)),mask40)
       global_sum = (REAL(ival1,dp)*r_fact)*two_40 + REAL(ival2,dp)*r_fact
       global_sum = -global_sum
    ENDIF

    IF (isum(2) >= 0_i8) THEN
       ival1 = ISHFT(isum(2),-40)
       ival2 = IAND (isum(2),mask40)
       global_sum = global_sum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_40
    ELSE
       ival1 = ISHFT(ABS(isum(2)),-40)
       ival2 = IAND (ABS(isum(2)),mask40)
       global_sum = global_sum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_40
    ENDIF

   IF(p_test_run .AND. do_sync_checks) CALL check_result( (/ global_sum /), 'global_sum_array2')

END FUNCTION global_sum_array2


!>
!! Calculates the global sum of 3D and/or 4D input fields and checks for
!! consistency when doing a verification run.
!! This routine has to be called from outside an OMP parallel Region!
!! If ldiff=.true., a difference or ratio between f??in and f??d is computed
!! nfields specifies the number of 3D input fields, or pairs of input fields if
!! ldiff=.true. (when providing a 4D field as input, it counts as "n" 3D fields)
!! diffmask specifies if differences or ratios are to be computed.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Optimized version by Guenther Zaengl, Apr 2010
!!
FUNCTION global_sum_array3 (nfields,ldiff,f3din,f3dd,f3din2,f3dd2,f4din,f4dd,diffmask) &
         RESULT (global_sum)

   REAL(wp), INTENT(in), TARGET           :: f3din(:,:,:)
   REAL(wp), INTENT(in), TARGET, OPTIONAL :: f3dd(:,:,:), f3din2(:,:,:), f3dd2(:,:,:), &
                                             f4din(:,:,:,:), f4dd(:,:,:,:)

   INTEGER,  INTENT(in)           :: nfields  ! Total number of 3D input fields
   LOGICAL,  INTENT(in)           :: ldiff    ! .TRUE.: compute differences or ratios (depending
                                              ! on diffmask) between f??in and f??d
   INTEGER,  INTENT(in), OPTIONAL :: diffmask(nfields) ! 1: compute differences, 2: compute ratios

   REAL(wp)                       :: global_sum(nfields*SIZE(f3din,2))

   INTEGER(i8)       :: ival1, ival2
   INTEGER           :: i, j, k, kk, n, nblks, nlen, nlev, nblks2, nfldtot, &
                        nlevtot, icnt, icnt2, nnblks(2*nfields)
   REAL(dp)          :: rlog2, rval1

   INTEGER(i8), ALLOCATABLE :: itmp(:,:),isum(:,:)
   INTEGER    , ALLOCATABLE :: iexp(:)
   REAL(dp),    ALLOCATABLE :: fact(:),r_fact(:),rval(:),abs_max(:)
   REAL(wp),    ALLOCATABLE :: aux_sum(:)

   TYPE t_fieldptr
     REAL(wp), POINTER :: fld(:,:,:)
   END TYPE t_fieldptr

   TYPE(t_fieldptr), ALLOCATABLE :: ff(:)

   REAL(dp), PARAMETER :: two_40 = 1099511627776._dp ! 2.**40
   REAL(dp), PARAMETER :: r_two_40 = 1._dp/two_40

#if defined (__SX__) || defined (__PGI)
   INTEGER(i8) :: mask40
   DATA mask40 / z'000000ffffffffff' / ! last 40 bits set
#else
   INTEGER(i8), PARAMETER :: mask40 = INT(z'000000ffffffffff',i8)
#endif
   INTEGER :: p_comm_glob

!-----------------------------------------------------------------------

   IF (ldiff .AND. .NOT.(PRESENT(diffmask))) THEN
     CALL finish('global_sum_array3','ldiff=.TRUE. requires the presence of diffmask')
   ENDIF

   IF(comm_lev==0) THEN
     p_comm_glob = p_comm_work
   ELSE
     p_comm_glob = glob_comm(comm_lev)
   ENDIF

   rlog2 = LOG(2._dp)
   nlen  = SIZE(f3din,1)
   nlev  = SIZE(f3din,2)
   nblks = SIZE(f3din,3)
   IF (ldiff) THEN
     nblks2 = SIZE(f3dd,3)
     nfldtot = 2*nfields
   ELSE
     nblks2 = 0
     nfldtot = nfields
   ENDIF
   nlevtot           = nlev*nfldtot
   nnblks(1:nfields) = nblks

   ALLOCATE (ff(nfldtot),itmp(nlevtot,2),isum(nlevtot,2),iexp(nlevtot),fact(nlevtot),&
             r_fact(nlevtot),rval(nlevtot),abs_max(nlevtot),aux_sum(nlevtot))

   itmp = 0_i8

   ! Set pointers to input fields
   icnt = 1
   icnt2 = 1
   ff(icnt)%fld => f3din
   IF (PRESENT(f3din2)) THEN
     icnt = icnt + 1
     icnt2 = icnt2 + 1
     ff(icnt)%fld => f3din2
   ENDIF
   IF (PRESENT(f4din)) THEN
     DO i = 1, nfields-icnt2
       icnt = icnt + 1
       ff(icnt)%fld => f4din(:,:,:,i)
     ENDDO
   ENDIF

   IF (ldiff) THEN
     icnt = icnt + 1
     icnt2 = 1
     ff(icnt)%fld => f3dd
     nnblks(icnt) = nblks2
     IF (PRESENT(f3dd2)) THEN
       icnt = icnt + 1
       icnt2 = icnt2 + 1
       ff(icnt)%fld => f3dd2(:,:,:)
       nnblks(icnt) = nblks2
     ENDIF
     IF (PRESENT(f4dd)) THEN
       DO i = 1, nfields-icnt2
         icnt = icnt + 1
         ff(icnt)%fld => f4dd(:,:,:,i)
         nnblks(icnt) = nblks2
       ENDDO
     ENDIF
   ENDIF

! OpenMP parallelization is done over the vertical levels in this routine
!$OMP PARALLEL
!$OMP DO PRIVATE(k,kk,n)
   ! Get the maximum absolute value of all field elements of a model level.
   DO k=1,nlevtot
     n = INT((k-1)/nlev)+1
     kk = k-(n-1)*nlev
     rval(k) = MAXVAL(ABS(ff(n)%fld(:,kk,:)))
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   abs_max = p_max(rval, comm=p_comm_glob)

   DO k=1,nlevtot
     abs_max(k) = MAX(abs_max(k),1.d-250)
     ! Get the exponent of abs_max for scaling
     iexp(k) = INT(LOG(abs_max(k))/rlog2)
     ! Calculate a factor for scaling the input numbers
     ! so that the maximum absolute value of a scaled number
     ! is below 2**40
     fact(k) = 2._wp**(40._wp-REAL(iexp(k),wp))
     r_fact(k) = 1._dp/fact(k)
   ENDDO

!$OMP PARALLEL
!$OMP DO PRIVATE(i, j, k, kk, n, ival1, rval1)
   DO k=1,nlevtot
     n = INT((k-1)/nlev)+1
     kk = k-(n-1)*nlev
     ! Sum up all numbers as scaled integers
!CDIR UNROLL=4
     DO j=1,nnblks(n)
       DO i=1,nlen
       ! Scale number into range -2**40 < rval < 2**40
       ! and store integer part in ival1

         rval1 = ff(n)%fld(i,kk,j)*fact(k)
         ival1 = INT(rval1,i8)

         ! Sum up ival1 and ival2; since we are using 8-byte integers
         ! for the sum there are no problems with overflow

         itmp(k,1) = itmp(k,1) + ival1
         itmp(k,2) = itmp(k,2) + INT((rval1 - REAL(ival1,dp))*two_40,i8)
       ENDDO
     ENDDO

   ENDDO
!$OMP END DO
!$OMP END PARALLEL

   isum(:,1) = p_sum(itmp(:,1), comm=p_comm_glob)
   isum(:,2) = p_sum(itmp(:,2), comm=p_comm_glob)

   DO k=1,nlevtot

     ! If the exponent is too small, return 0 in order to avoid
     ! problems with overflow below
     IF(iexp(k) < -980) THEN
       aux_sum(k) = 0._dp
       CYCLE
     ENDIF

     ! Scale integer numbers back to real numbers and add them.
     ! For safety, we use only positive INTEGERS < 2**40 when converting to REAL

     IF (isum(k,1) >= 0_i8) THEN
       ival1 = ISHFT(isum(k,1),-40)
       ival2 = IAND (isum(k,1),mask40)
       aux_sum(k) = (REAL(ival1,dp)*r_fact(k))*two_40 + REAL(ival2,dp)*r_fact(k)
     ELSE
       ival1 = ISHFT(ABS(isum(k,1)),-40)
       ival2 = IAND (ABS(isum(k,1)),mask40)
       aux_sum(k) = (REAL(ival1,dp)*r_fact(k))*two_40 + REAL(ival2,dp)*r_fact(k)
       aux_sum(k) = -aux_sum(k)
     ENDIF

     IF (isum(k,2) >= 0_i8) THEN
       ival1 = ISHFT(isum(k,2),-40)
       ival2 = IAND (isum(k,2),mask40)
       aux_sum(k) = aux_sum(k) + (REAL(ival1,dp)*r_fact(k)) + &
                    (REAL(ival2,dp)*r_fact(k))*r_two_40
     ELSE
       ival1 = ISHFT(ABS(isum(k,2)),-40)
       ival2 = IAND (ABS(isum(k,2)),mask40)
       aux_sum(k) = aux_sum(k) - (REAL(ival1,dp)*r_fact(k)) - &
                    (REAL(ival2,dp)*r_fact(k))*r_two_40
     ENDIF

   ENDDO

   IF (ldiff) THEN ! Compute difference or ratio of global sums
     n = nlevtot/2
     DO k=1,n
       i = INT((k-1)/nlev)+1
       IF (diffmask(i) == 1) THEN
         global_sum(k) = aux_sum(k) - aux_sum(k+n)
       ELSE
         global_sum(k) = MIN(1.001_wp, aux_sum(k) / MAX(1.e-50_wp,aux_sum(k+n)))
         global_sum(k) = MAX(0.999_wp, global_sum(k))
       ENDIF
     ENDDO
   ELSE
     n = nlevtot
     global_sum(1:n) = aux_sum(1:n)
   ENDIF

   IF(p_test_run .AND. do_sync_checks) CALL check_result(global_sum, 'global_sum_array3')

   DEALLOCATE (ff,itmp,isum,iexp,fact,r_fact,rval,abs_max,aux_sum)

END FUNCTION global_sum_array3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> global_min/global_max family:
!! Calculate global min/max (using current communicator) and compare
!! result to result on test PE (in case of a test run)
FUNCTION global_min_0d(zfield) RESULT(global_min)

  REAL(wp), INTENT(IN) :: zfield
  REAL(wp) :: global_min

  IF(comm_lev==0) THEN
    global_min = p_min(zfield, comm=p_comm_work)
  ELSE
    global_min = p_min(zfield, comm=glob_comm(comm_lev))
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( (/ global_min /), 'global_min' )

END FUNCTION global_min_0d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
FUNCTION global_min_1d(zfield) RESULT(global_min)

  REAL(wp), INTENT(IN) :: zfield(:)
  REAL(wp) :: global_min(SIZE(zfield))

  IF(comm_lev==0) THEN
    global_min = p_min(zfield, comm=p_comm_work)
  ELSE
    global_min = p_min(zfield, comm=glob_comm(comm_lev))
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( global_min, 'global_min' )

END FUNCTION global_min_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! @param[out]   proc_id  (Optional:) PE number of maximum value
! @param[inout] keyval   (Optional:) additional meta information
! @param[in]    iroot    (Optional:) root PE, otherwise we perform an
!                                    ALL-TO-ALL operation
!
! The parameter @p keyval can be used to communicate
! additional data on the maximum value, e.g., the level
! index where the maximum occurred.
!
FUNCTION global_max_0d(zfield, proc_id, keyval, iroot, icomm) RESULT(global_max)

  REAL(wp), INTENT(IN) :: zfield
  INTEGER, OPTIONAL, INTENT(inout) :: proc_id
  INTEGER, OPTIONAL, INTENT(inout) :: keyval
  INTEGER, OPTIONAL, INTENT(in)    :: iroot                !< rank of collecting PE
  INTEGER, OPTIONAL, INTENT(in)    :: icomm                !< MPI communicator
  REAL(wp) :: global_max
  INTEGER  :: pcomm

  IF(comm_lev==0) THEN
    pcomm=p_comm_work
  ELSE
    pcomm=glob_comm(comm_lev)
  END IF
  IF (PRESENT(icomm)) pcomm = icomm

  IF (p_test_run) THEN ! all-to-all communication required
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm)
  ELSE
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm, root=iroot)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( (/ global_max /), 'global_max' )

END FUNCTION global_max_0d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
! @param[out]   proc_id  (Optional:) PE number of maximum value
! @param[inout] keyval   (Optional:) additional meta information
! @param[in]    iroot    (Optional:) root PE, otherwise we perform an
!                                    ALL-TO-ALL operation
!
! The parameter @p keyval can be used to communicate
! additional data on the maximum value, e.g., the level
! index where the maximum occurred.
!
FUNCTION global_max_1d(zfield, proc_id, keyval, iroot, icomm) RESULT(global_max)

  REAL(wp), INTENT(IN) :: zfield(:)
  INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
  INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
  INTEGER, OPTIONAL, INTENT(in)    :: iroot                !< rank of collecting PE
  INTEGER, OPTIONAL, INTENT(in)    :: icomm                !< MPI communicator
  REAL(wp) :: global_max(SIZE(zfield))
  INTEGER  :: pcomm

  IF(comm_lev==0) THEN
    pcomm=p_comm_work
  ELSE
    pcomm=glob_comm(comm_lev)
  END IF
  IF (PRESENT(icomm)) pcomm = icomm

  IF (p_test_run) THEN ! all-to-all communication required
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm)
  ELSE
    global_max = p_max(zfield, proc_id=proc_id, keyval=keyval, comm=pcomm, root=iroot)
  ENDIF

  IF(p_test_run .AND. do_sync_checks) CALL check_result( global_max, 'global_max' )

END FUNCTION global_max_1d

!-------------------------------------------------------------------------------

!> Checks if res is identical on Test PE and working PEs.
!! If res_on_testpe is present, it contains the result on the test PE on exit,
!! otherwise the routine finishes if the results are not identical.

SUBROUTINE check_result(res, routine, res_on_testpe)

  REAL(wp), INTENT(IN) :: res(:)
  CHARACTER(len=*), INTENT(IN) :: routine
  REAL(wp), INTENT(out), OPTIONAL :: res_on_testpe(:)

  REAL(wp) :: aux(SIZE(res))
  INTEGER :: k
  LOGICAL :: out_of_sync


  aux(:) = 0.0_wp ! Safety only
  IF(my_process_is_mpi_test()) aux(:) = res(:)

  IF(comm_lev==0) THEN
    CALL p_bcast(aux, process_mpi_all_test_id, comm=p_comm_work_test)
  ELSE
    IF(get_my_mpi_all_id() == process_mpi_all_test_id) THEN
      CALL p_send(aux, comm_proc0(comm_lev)+p_work_pe0, 1)
    ELSE
      IF(p_pe_work==comm_proc0(comm_lev)) CALL p_recv(aux, process_mpi_all_test_id, 1)
      CALL p_bcast(aux, 0, comm=glob_comm(comm_lev))
    ENDIF
  ENDIF

  out_of_sync = .FALSE.
  DO k = 1, SIZE(res)
    IF( .NOT. my_process_is_mpi_test() .AND. l_log_checks .AND. log_unit>0) &
      & WRITE(log_unit,'(a,2g25.18,a,g25.18)') routine,aux(k),res(k),' Error: ',ABS(aux(k)-res(k))
    IF(PRESENT(res_on_testpe)) THEN
      res_on_testpe(k) = aux(k)
    ELSE
      ! Check if result is identical
      IF(aux(k)/=res(k)) out_of_sync = .TRUE.
    ENDIF
  ENDDO

  IF(out_of_sync) CALL finish(routine, 'Result out of sync')

END SUBROUTINE check_result


!-------------------------------------------------------------------------------
!>
!! This routine calculates the sum of an array of IEEE 64 bit
!! floating point values in an order insensitve way.
!! This is done by calculating the sum in INTEGER arithmetic,
!! so it is always exactly the same for any permutation of the numbers
!! (even among several processors).
!! Since 60 bits are used for the mantissa of each operand (and even more
!! for the accumulator), this sum should be almost always more precise
!! than naivly summing up the operands.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called from a parallel region!!!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION omp_order_insensit_ieee64_sum(vals, num_vals, mpi_comm) RESULT(global_sum)

!
   INTEGER  :: num_vals, mpi_comm
   REAL(dp) :: vals(num_vals)

   REAL(dp) :: global_sum

   INTEGER(i8)       :: itmp(2),  ival1, ival2
   INTEGER(i8), SAVE :: isum(2) ! This must be a shared variable
   INTEGER           :: i, iexp
   REAL(dp), SAVE    :: abs_max ! This must be a shared variable
   REAL(dp)          :: fact, r_fact, rval

   REAL(dp), PARAMETER :: two_30 = 1073741824._dp ! 2.**30
   REAL(dp), PARAMETER :: r_two_30 = 1._dp/two_30

#if defined (__SX__) || defined (__PGI)
   INTEGER(i8) :: mask30
   DATA mask30 / z'000000003fffffff' / ! last 30 bits set
#else
   INTEGER(i8), PARAMETER :: mask30 = INT(z'000000003fffffff',i8)
#endif

!-----------------------------------------------------------------------
!    IF (activate_sync_timers) CALL timer_start(timer_omp_ordglb_sum)

   ! Set shared variables in a MASTER region

!$OMP MASTER
   abs_max = 0._dp
   isum(:) = 0_i8
!$OMP END MASTER
!$OMP BARRIER

   ! Get the maximum absolute value of all numbers.

!$OMP DO PRIVATE(i), REDUCTION(MAX:abs_max)
   DO i=1,num_vals
      abs_max = MAX(abs_max, ABS(vals(i)))
   ENDDO
!$OMP END DO

!$OMP MASTER
   rval = abs_max
   abs_max = p_max(rval, comm=mpi_comm)
!$OMP END MASTER
!$OMP BARRIER

   ! If abs_max is 0, all input values are 0
   ! and we are done

   IF(abs_max == 0.0_dp) THEN
      global_sum = 0._dp
      RETURN
   ENDIF

   ! Get the exponent of abs_max for scaling

   iexp = EXPONENT(abs_max)

   ! If the exponent is too small, return 0 in order to avoid
   ! problems with overflow below

   IF(iexp < -980) THEN
      global_sum = 0._dp
      RETURN
   ENDIF

   ! Calculate a factor for scaling the input numbers
   ! so that the maximum absolute value of a scaled number
   ! is below 2**30

   fact = SCALE(1._dp,30-iexp) ! same as 2**(30-iexp)
   r_fact = SCALE(1._dp,iexp-30) ! 1./fact

   ! Sum up all numbers as scaled integers

!$OMP DO PRIVATE(i, rval, ival1, ival2), REDUCTION(+:isum)
   DO i=1,num_vals

      ! Scale number into range -2**30 < rval < 2**30
      ! and store integer part in ival1

      rval = vals(i)*fact
      ival1 = INT(rval,i8)

      ! Scale fraction by 2**30 and store integer part in ival2

      ival2 = INT((rval - REAL(ival1,dp))*two_30,i8)

      ! Sum up ival1 and ival2; since we are using 8-byte integers
      ! for the sum there are no problems with overflow

      isum(1) = isum(1) + ival1
      isum(2) = isum(2) + ival2

   ENDDO
!$OMP END DO

!$OMP MASTER
   itmp = isum
   isum = p_sum(itmp, comm=mpi_comm)
!$OMP END MASTER
!$OMP BARRIER

   ! Scale integer numbers back to real numbers and add them.
   ! For safety, we use only positive INTEGERS < 2**30 when converting to REAL

    IF(isum(1) >= 0_i8)THEN
       ival1 = ISHFT(isum(1),-30)
       ival2 = IAND (isum(1),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
    ELSE
       ival1 = ISHFT(ABS(isum(1)),-30)
       ival2 = IAND (ABS(isum(1)),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
       global_sum = -global_sum
    ENDIF

    IF(isum(2) >= 0_i8)THEN
       ival1 = ISHFT(isum(2),-30)
       ival2 = IAND (isum(2),mask30)
       global_sum = global_sum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_30
    ELSE
       ival1 = ISHFT(ABS(isum(2)),-30)
       ival2 = IAND (ABS(isum(2)),mask30)
       global_sum = global_sum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_30
    ENDIF
   
!    IF (activate_sync_timers) CALL timer_stop(timer_omp_ordglb_sum)

END FUNCTION omp_order_insensit_ieee64_sum
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!>
!! This routine calculates the sum of an array of IEEE 64 bit
!! floating point values in an order insensitve way.
!! This is done by calculating the sum in INTEGER arithmetic,
!! so it is always exactly the same for any permutation of the numbers
!! (even among several processors).
!! Since 60 bits are used for the mantissa of each operand (and even more
!! for the accumulator), this sum should be almost always more precise
!! than naivly summing up the operands.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called outside an omp parallel region!!!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION order_insensit_ieee64_sum(vals, num_vals, mpi_comm) RESULT(global_sum)

!
   INTEGER  :: num_vals, mpi_comm
   REAL(dp) :: vals(num_vals)

   REAL(dp) :: global_sum

   INTEGER(i8)       :: itmp(2),  ival1, ival2
   INTEGER(i8), SAVE :: isum(2) ! This must be a shared variable
   INTEGER           :: i, iexp
   REAL(dp), SAVE    :: abs_max ! This must be a shared variable
   REAL(dp)          :: fact, r_fact, rval

   REAL(dp), PARAMETER :: two_30 = 1073741824._dp ! 2.**30
   REAL(dp), PARAMETER :: r_two_30 = 1._dp/two_30

#if defined (__SX__) || defined (__PGI)
   INTEGER(i8) :: mask30
   DATA mask30 / z'000000003fffffff' / ! last 30 bits set
#else
   INTEGER(i8), PARAMETER :: mask30 = INT(z'000000003fffffff',i8)
#endif

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_ordglb_sum)

   ! Set shared variables in a MASTER region
   abs_max = 0._dp
   isum(:) = 0_i8
   ! Get the maximum absolute value of all numbers.
   DO i=1,num_vals
      abs_max = MAX(abs_max, ABS(vals(i)))
   ENDDO
   rval = abs_max
   abs_max = p_max(rval, comm=mpi_comm)
   ! If abs_max is 0, all input values are 0
   ! and we are done
   IF(abs_max == 0.0_dp) THEN
      global_sum = 0._dp
      IF (activate_sync_timers) CALL timer_stop(timer_ordglb_sum)
      RETURN
   ENDIF

   ! Get the exponent of abs_max for scaling

   iexp = EXPONENT(abs_max)

   ! If the exponent is too small, return 0 in order to avoid
   ! problems with overflow below

   IF(iexp < -980) THEN
      global_sum = 0._dp
      IF (activate_sync_timers) CALL timer_stop(timer_ordglb_sum)
      RETURN
   ENDIF

   ! Calculate a factor for scaling the input numbers
   ! so that the maximum absolute value of a scaled number
   ! is below 2**30

   fact = SCALE(1._dp,30-iexp) ! same as 2**(30-iexp)
   r_fact = SCALE(1._dp,iexp-30) ! 1./fact

   ! Sum up all numbers as scaled integers
   DO i=1,num_vals

      ! Scale number into range -2**30 < rval < 2**30
      ! and store integer part in ival1

      rval = vals(i)*fact
      ival1 = INT(rval,i8)

      ! Scale fraction by 2**30 and store integer part in ival2

      ival2 = INT((rval - REAL(ival1,dp))*two_30,i8)

      ! Sum up ival1 and ival2; since we are using 8-byte integers
      ! for the sum there are no problems with overflow

      isum(1) = isum(1) + ival1
      isum(2) = isum(2) + ival2

   ENDDO
   itmp = isum
   isum = p_sum(itmp, comm=mpi_comm)
   
   ! Scale integer numbers back to real numbers and add them.
   ! For safety, we use only positive INTEGERS < 2**30 when converting to REAL

    IF(isum(1) >= 0_i8)THEN
       ival1 = ISHFT(isum(1),-30)
       ival2 = IAND (isum(1),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
    ELSE
       ival1 = ISHFT(ABS(isum(1)),-30)
       ival2 = IAND (ABS(isum(1)),mask30)
       global_sum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
       global_sum = -global_sum
    ENDIF

    IF(isum(2) >= 0_i8)THEN
       ival1 = ISHFT(isum(2),-30)
       ival2 = IAND (isum(2),mask30)
       global_sum = global_sum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_30
    ELSE
       ival1 = ISHFT(ABS(isum(2)),-30)
       ival2 = IAND (ABS(isum(2)),mask30)
       global_sum = global_sum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_30
    ENDIF
   
    IF (activate_sync_timers) CALL timer_stop(timer_ordglb_sum)

END FUNCTION order_insensit_ieee64_sum
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!>
!! This routine calculates the sum of an array in the
!! straightforward way in parallel.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called outside a parallel omp region!!!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION simple_sum(vals, num_vals, mpi_comm) RESULT(global_sum)

!
   INTEGER  :: num_vals, mpi_comm
   REAL(dp) :: vals(num_vals)

   REAL(dp) :: global_sum

   INTEGER :: i
   REAL(dp), SAVE :: s, res

!-----------------------------------------------------------------------
   IF (activate_sync_timers) CALL timer_start(timer_global_sum)

   s = 0._dp
   ! Sum up all numbers
   DO i=1,num_vals
      s = s + vals(i)
   ENDDO
   res = p_sum(s, comm=mpi_comm)

   global_sum = res
  
   IF (activate_sync_timers) CALL timer_stop(timer_global_sum)

END FUNCTION simple_sum
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!>
!! This routine calculates the sum of an array in the
!! straightforward way in parallel.
!! ATTENTION: When compiled with OpenMP in effect, this routine
!! should be called from a parallel region!!!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
FUNCTION omp_simple_sum(vals, num_vals, mpi_comm) RESULT(global_sum)

!
   INTEGER  :: num_vals, mpi_comm
   REAL(dp) :: vals(num_vals)

   REAL(dp) :: global_sum

   INTEGER :: i
   REAL(dp), SAVE :: s, res

!-----------------------------------------------------------------------
  IF (activate_sync_timers) CALL timer_start(timer_omp_global_sum)

   ! Set shared variables in a MASTER region

!$OMP MASTER
   s = 0._dp
!$OMP END MASTER
!$OMP BARRIER

   ! Sum up all numbers

!$OMP DO PRIVATE(i), REDUCTION(+:s)
   DO i=1,num_vals
      s = s + vals(i)
   ENDDO
!$OMP END DO

!$OMP MASTER
   res = p_sum(s, comm=mpi_comm)
!$OMP END MASTER
!$OMP BARRIER

   global_sum = res
   
   IF (activate_sync_timers) CALL timer_stop(timer_omp_global_sum)

END FUNCTION omp_simple_sum
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
! exact_ieee64_sum is currently unused !!!
! Please note: If it should ever be reactivated, it needs an
! MPI communicator as additinal input argument.
! It currently works only on MPI_COMM_WORLD.
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
!

!>
!! This routine calculates the exact sum of an array of IEEE 64 bit.
!!
!! This routine calculates the exact sum of an array of IEEE 64 bit
!! floating point values rounded to the nearest representable 64 bit number.
!! Since the sum is calculated exactly in INTEGER arithmetic, the result
!! for any permutation of the numbers (even among several processors)
!! is always exactly the same.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!!
!SUBROUTINE exact_ieee64_sum(vals, num_vals, res, mpi_comm)
SUBROUTINE exact_ieee64_sum(vals, num_vals, res)

!

   IMPLICIT NONE

   INTEGER  :: num_vals
!   INTEGER  :: mpi_comm
   REAL(dp) :: vals(num_vals), res

   INTEGER, PARAMETER :: nacc = 66 ! Length of accumulator in 32 bit units

   INTEGER(i8) :: iacc(nacc+1), iacc_send(nacc+1), ix, i1, i2, i3, ival
   INTEGER     :: iexp, ipos
   INTEGER     :: i, isig, ish
!   INTEGER     :: mpi_err

#if defined (__SX__) || defined (__PGI)
   INTEGER(i8) :: ibit32
   INTEGER(i8) :: ione
   INTEGER(i8) :: mask32
   INTEGER(i8) :: mask52
   INTEGER(i8) :: signan

   DATA ibit32 / z'0000000100000000' /
   DATA ione   / z'0000000000000001' / ! 1 as an 8 byte integer
   DATA mask32 / z'00000000ffffffff' / ! last 32 bits set
   DATA mask52 / z'000fffffffffffff' / ! last 52 bits set
   DATA signan / z'7ff7ffffffffffff' / ! signalling nan
#else
   INTEGER(i8), PARAMETER :: ibit32 = INT(z'0000000100000000',i8)
   INTEGER(i8), PARAMETER :: ione   = INT(z'0000000000000001',i8)
   INTEGER(i8), PARAMETER :: mask32 = INT(z'00000000ffffffff',i8)
   INTEGER(i8), PARAMETER :: mask52 = INT(z'000fffffffffffff',i8)
   INTEGER(i8), PARAMETER :: signan = INT(z'7ff7ffffffffffff',i8)
#endif

!-----------------------------------------------------------------------

   iacc(:) = 0_i8

   ! Add all input numbers to the accumulator

   DO i=1,num_vals

      ix = TRANSFER(vals(i),ix)

      ! Get exponent

      iexp = INT(IBITS(ix,52,11))

      ! If iexp==0 this is +/-0 or a denormalized number (we treat that as 0)

      IF (iexp == 0) CYCLE

      ! If iexp is 2047, there is Inf or NaN among the input,
      ! set iacc(nacc+1) (as a flag) and terminate summing.

      IF(iexp == 2047) THEN
         iacc(nacc+1) = 1_i8
         EXIT
      ENDIF

      ! Get mantissa, we have to add the implicit leading 1 (by setting bit 52)

      ival = IBSET(IAND(ix,mask52),52)

      ! Shift mantissa to the right position and add it to the accumulator

      ipos = iexp/32
      ish  = MOD(iexp,32)

      i1 = IAND(ISHFT(ival,ish   ),mask32)
      i2 = IAND(ISHFT(ival,ish-32),mask32)
      i3 = IAND(ISHFT(ival,ish-64),mask32)

      IF(BTEST(ix,63)) THEN
         iacc(ipos+1) = iacc(ipos+1) - i1
         iacc(ipos+2) = iacc(ipos+2) - i2
         iacc(ipos+3) = iacc(ipos+3) - i3
      ELSE
         iacc(ipos+1) = iacc(ipos+1) + i1
         iacc(ipos+2) = iacc(ipos+2) + i2
         iacc(ipos+3) = iacc(ipos+3) + i3
      ENDIF

   ENDDO

   ! Get global sum over iacc

#ifndef NOMPI
   iacc_send(:) = iacc(:)
!   CALL mpi_allreduce(iacc_send, iacc, nacc+1, MPI_INTEGER8, MPI_SUM, mpi_comm, mpi_err)
   iacc = p_sum(iacc_send)
#endif

   ! If Inf or NaN is among the input, return signalling Nan

   IF (iacc(nacc+1) > 0_i8) THEN
      res = TRANSFER(signan, res)
      RETURN
   ENDIF

   ! Normalize accumulator transferring carries to the next higher word

   DO i=1,nacc-1
      iacc(i+1) = iacc(i+1) + IAND(ISHFT(iacc(i),-32),mask32)
      IF(BTEST(iacc(i),63)) iacc(i+1) = iacc(i+1) - ibit32
      iacc(i) = IAND(iacc(i),mask32)
   ENDDO

   ! If the result is negative invert it and store sign

   isig = 0
   IF (iacc(nacc) < 0_i8) THEN
      isig = 1
      iacc(1:nacc) = IAND(NOT(iacc(1:nacc)),mask32)
      iacc(1) = iacc(1) + 1_i8
      DO i=1,nacc-1
         IF(BTEST(iacc(i),32)) THEN
            iacc(i+1) = iacc(i+1) + 1_i8
            iacc(i) = IAND(iacc(i),mask32)
         ELSE
            EXIT
         ENDIF
      ENDDO
   ENDIF

   ! Search highest word in accumulator not equal to 0

   DO i=nacc,1,-1
      IF (iacc(i) /= 0_i8) EXIT
   ENDDO

   IF (i == 0) THEN
      ! Result is exactly 0
      res = 0.0_dp
      RETURN
   ENDIF

   ! Get high order 64 bits of accumulator ...

   IF(i==1) THEN
      ival = ISHFT(iacc(i),32)
   ELSE
      ival = IOR(ISHFT(iacc(i),32),iacc(i-1))
   ENDIF
   iexp = (i-2)*32 ! This will be changed by shifting below

   ! ... and shift until first bit is 1

   DO WHILE(.NOT. BTEST(ival,63))
      ival = ISHFT(ival,1)
      IF(i>2) THEN
         IF(BTEST(iacc(i-2),31)) ival = IOR(ival,ione)
         iacc(i-2) = ISHFT(iacc(i-2),1)
      ENDIF
      iexp = iexp-1
   ENDDO
   iexp = iexp + 11 ! We have shifted 11 bits too far to the left

   ! Round result

   ival = ival + 1024_i8 ! 11 bits at the right will be discarded
   ! If the above operation caused an overflow, all relevant 53 bits of the
   ! result were set and are now zero. We don't need to shift again,
   ! but we have to increase the exponent:
   IF(.NOT.BTEST(ival,63)) iexp = iexp+1

   ! Check for underflow

   IF(iexp<=0) THEN
      res = 0.0_dp
      RETURN
   ENDIF

   ! Check for overflow

   IF(iexp>=2047) THEN
      iexp = 2047 ! This will result in +/-Inf in the result
      ival = 0_i8
   ENDIF

   ! Assemble result

   ix = IBITS(ival,11,52)
   ix = IOR(ix,ISHFT(INT(iexp,i8),52))
   IF(isig==1) ix = IBSET(ix,63)

   res = TRANSFER(ix,res)

END SUBROUTINE exact_ieee64_sum

!>
!! Computes and prints summary information on the domain decomposition
!!
!!
!! @par Revision History
!! Initial version by Guenther Zaengl, Nov 2011
!!
SUBROUTINE decomposition_statistics(p_patch)

   TYPE(t_patch), INTENT(INOUT) :: p_patch

   REAL(wp) :: cellstat(0:6),edgestat(0:6),vertstat(0:5), csmax(0:5),csmin(0:5),csavg(0:6), &
               esmax(0:6),esmin(0:6),esavg(0:6),vsmax(0:5),vsmin(0:5),vsavg(0:5),avglat,avglon
   INTEGER  :: i_nchdom, i, i_pe, max_nprecv, i1, i2, i1m, i2m
   INTEGER,  ALLOCATABLE :: nprecv_buf(:),displs(:),recvlist_buf(:)
   REAL(wp), ALLOCATABLE :: avglat_buf(:),avglon_buf(:)

!-----------------------------------------------------------------------

   i_nchdom = MAX(1,p_patch%n_childdom)

   cellstat(0) = REAL(nproma*(p_patch%cells%end_blk(grf_bdywidth_c,1)-1) + &
                      p_patch%cells%end_idx(grf_bdywidth_c,1),wp)
   cellstat(1) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int,i_nchdom),wp) - cellstat(0)
   cellstat(2) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int-1,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int-1,i_nchdom),wp) - cellstat(0)
   cellstat(3) = REAL(nproma*(p_patch%cells%end_blk(min_rlcell_int-2,i_nchdom)-1) + &
                      p_patch%cells%end_idx(min_rlcell_int-2,i_nchdom),wp) - cellstat(0)
   cellstat(4) = REAL(p_patch%comm_pat_c%np_send,wp)
   cellstat(5) = REAL(p_patch%comm_pat_c%np_recv,wp)

   ! The purpose of this is to compute average quantities only over those PEs
   ! that actually contain grid points of a given model domain (relevant in the case
   ! of processor splitting)
   IF (cellstat(1) > 0._wp) THEN
     cellstat(6) = 1._wp
   ELSE
     cellstat(6) = 0._wp
   ENDIF

   edgestat(0) = REAL(nproma*(p_patch%edges%end_blk(grf_bdywidth_e,1)-1) + &
                      p_patch%edges%end_idx(grf_bdywidth_e,1),wp)
   edgestat(1) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int,i_nchdom),wp) - edgestat(0)
   edgestat(2) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int-1,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int-1,i_nchdom),wp) - edgestat(0)
   edgestat(3) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int-2,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int-2,i_nchdom),wp) - edgestat(0)
   edgestat(4) = REAL(nproma*(p_patch%edges%end_blk(min_rledge_int-3,i_nchdom)-1) + &
                      p_patch%edges%end_idx(min_rledge_int-3,i_nchdom),wp) - edgestat(0)
   edgestat(5) = REAL(p_patch%comm_pat_e%np_send,wp)
   edgestat(6) = REAL(p_patch%comm_pat_e%np_recv,wp)

   vertstat(0) = REAL(nproma*(p_patch%verts%end_blk(grf_bdywidth_c,1)-1) + &
                      p_patch%verts%end_idx(grf_bdywidth_c,1),wp)
   vertstat(1) = REAL(nproma*(p_patch%verts%end_blk(min_rlvert_int,i_nchdom)-1) + &
                      p_patch%verts%end_idx(min_rlvert_int,i_nchdom),wp) - vertstat(0)
   vertstat(2) = REAL(nproma*(p_patch%verts%end_blk(min_rlvert_int-1,i_nchdom)-1) + &
                      p_patch%verts%end_idx(min_rlvert_int-1,i_nchdom),wp) - vertstat(0)
   vertstat(3) = REAL(nproma*(p_patch%verts%end_blk(min_rlvert_int-2,i_nchdom)-1) + &
                      p_patch%verts%end_idx(min_rlvert_int-2,i_nchdom),wp) - vertstat(0)
   vertstat(4) = REAL(p_patch%comm_pat_v%np_send,wp)
   vertstat(5) = REAL(p_patch%comm_pat_v%np_recv,wp)

   ! Question: how can I exclude PEs containing zero grid points of a model domain
   ! from global minimum computation?
   csmax = global_max(cellstat(0:5))
   csmin = global_min(cellstat(0:5))
   DO i = 0, 6
     csavg(i) = global_sum_array(cellstat(i))
   ENDDO
   csavg(0:5) = csavg(0:5)/MAX(1._wp,csavg(6))
   esmax = global_max(edgestat)
   esmin = global_min(edgestat)
   DO i = 0, 6
     esavg(i) = global_sum_array(edgestat(i))/MAX(1._wp,csavg(6))
   ENDDO
   vsmax = global_max(vertstat)
   vsmin = global_min(vertstat)
   DO i = 0, 5
     vsavg(i) = global_sum_array(vertstat(i))/MAX(1._wp,csavg(6))
   ENDDO

   WRITE(message_text,'(a,i4)') 'grid ',p_patch%id
     CALL message('Information on domain decomposition',TRIM(message_text))
   WRITE(message_text,'(i6)') NINT(csavg(6))
     CALL message('Number of compute PEs used for this grid',TRIM(message_text))
   IF (p_patch%id > 1) THEN
     WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(0)),NINT(csmin(0)),csavg(0)
      CALL message('#   lateral boundary cells', TRIM(message_text))
   ENDIF  
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(1)),NINT(csmin(1)),csavg(1)
     CALL message('#         prognostic cells', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(2)),NINT(csmin(2)),csavg(2)
     CALL message('# cells up to halo level 1', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(3)),NINT(csmin(3)),csavg(3)
     CALL message('# cells up to halo level 2', TRIM(message_text))

   IF (p_patch%id > 1) THEN
     WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(0)),NINT(esmin(0)),esavg(0)
       CALL message('#   lateral boundary edges', TRIM(message_text))
   ENDIF
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(1)),NINT(esmin(1)),esavg(1)
     CALL message('#         prognostic edges', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(2)),NINT(esmin(2)),esavg(2)
     CALL message('# edges up to halo level 1', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(3)),NINT(esmin(3)),esavg(3)
     CALL message('# edges up to halo level 2', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(4)),NINT(esmin(4)),esavg(4)
     CALL message('# edges up to halo level 3', TRIM(message_text))

   IF (p_patch%id > 1) THEN
     WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(0)),NINT(vsmin(0)),vsavg(0)
       CALL message('#   lateral boundary verts', TRIM(message_text))
   ENDIF
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(1)),NINT(vsmin(1)),vsavg(1)
     CALL message('#         prognostic verts', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(2)),NINT(vsmin(2)),vsavg(2)
     CALL message('# verts up to halo level 1', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(3)),NINT(vsmin(3)),vsavg(3)
     CALL message('# verts up to halo level 2', TRIM(message_text))
     CALL message('','')
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(4)),NINT(csmin(4)),csavg(4)
     CALL message('# send PEs (cells)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(csmax(5)),NINT(csmin(5)),csavg(5)
     CALL message('# recv PEs (cells)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(5)),NINT(esmin(5)),esavg(5)
     CALL message('# send PEs (edges)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(esmax(6)),NINT(esmin(6)),esavg(6)
     CALL message('# recv PEs (edges)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(4)),NINT(vsmin(4)),vsavg(4)
     CALL message('# send PEs (verts)', TRIM(message_text))
   WRITE(message_text,'(a,2i7,f10.2)') 'max/min/avg ',NINT(vsmax(5)),NINT(vsmin(5)),vsavg(5)
     CALL message('# recv PEs (verts)', TRIM(message_text))
     CALL message('','')

   ! Stop if the current model domain lives on only one processor. MPI communication
   ! will not work in this case
   IF (.NOT. p_test_run .AND. NINT(csavg(6)) <= 1) &
     CALL finish('Bad use of processor splitting','This grid is processed by only one PE')


   IF (msg_level >= 14) THEN ! add list of cell recv PEs for each PE

     WRITE(message_text,'(a,i4)') 'grid ',p_patch%id
     CALL message('List of receive PEs (cells)', TRIM(message_text))
     max_nprecv = NINT(csmax(5))

     ! Compute average latitude and longitude over prognostic grid points (including nest boundary)
     avglat = 0._wp
     avglon = 0._wp
     i2m = p_patch%cells%end_blk(min_rlcell_int,i_nchdom)
     DO i2 = 1, i2m
       i1m = nproma
       if (i2 == i2m) i1m = p_patch%cells%end_idx(min_rlcell_int,i_nchdom)
       DO i1 = 1, i1m
         avglat = avglat + p_patch%cells%center(i1,i2)%lat
         avglon = avglon + p_patch%cells%center(i1,i2)%lon
       ENDDO
     ENDDO

     IF (cellstat(0)+cellstat(1) > 0._wp) THEN
       avglat = avglat/(cellstat(0)+cellstat(1))*180._wp/pi
       avglon = avglon/(cellstat(0)+cellstat(1))*180._wp/pi
     ELSE
       avglat = 0._wp
       avglon = 0._wp
     ENDIF

     ALLOCATE(nprecv_buf(p_n_work),displs(p_n_work),recvlist_buf(max_nprecv*p_n_work), &
              avglat_buf(p_n_work),avglon_buf(p_n_work))

     CALL p_gather(p_patch%comm_pat_c%np_recv, nprecv_buf, 0, p_comm_work)

     displs(:) = (/ ( (i_pe-1)*max_nprecv, i_pe=1, p_n_work) /)

     CALL p_gatherv(p_patch%comm_pat_c%pelist_recv, p_patch%comm_pat_c%np_recv,  &
                    recvlist_buf, nprecv_buf, displs, 0, p_comm_work)

     CALL p_gather(avglat, avglat_buf, 0, p_comm_work)
     CALL p_gather(avglon, avglon_buf, 0, p_comm_work)

     IF (p_pe == 0) THEN
       WRITE(0,'(a)') 'PE, avg lat/lon, list of receive PEs'
       DO i_pe = 1, p_n_work
         i1 = (i_pe-1)*max_nprecv+1
         i2 = i1-1+MIN(18,nprecv_buf(i_pe))
         IF (nprecv_buf(i_pe) > 0) WRITE(0,'(i6,a,2f8.2,18i6)') i_pe-1,':', &
           avglat_buf(i_pe),avglon_buf(i_pe),recvlist_buf(i1:i2)
       ENDDO
     ENDIF

     DEALLOCATE(nprecv_buf,displs,recvlist_buf,avglat_buf,avglon_buf)

   ENDIF

END SUBROUTINE decomposition_statistics

!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!> Does boundary exchange for 3D cell-based fields, collecting as many
!  fields as possible before actually performing the sync.
!
!  @par Revision History
!  Initial revision : F. Prill, DWD (2013-02-08)
! 
SUBROUTINE cumulative_sync_patch_array(typ, p_patch, f3d)
  INTEGER,       INTENT(IN)            :: typ
  TYPE(t_patch), INTENT(IN),    TARGET :: p_patch
  REAL(wp),      INTENT(INOUT), TARGET :: f3d(:,:,:)
  ! local variables
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_sync:cumulative_sync_patch_array")
  INTEGER :: idx

  ! add pointer to list of cumulative sync fields:
  idx = ncumul_sync(typ,p_patch%id) + 1
  IF (idx > MAX_CUMULATIVE_SYNC) CALL finish(routine, "Internal error!")
  ncumul_sync(typ, p_patch%id)            =  idx
  cumul_sync(typ, p_patch%id,idx)%f3d     => f3d
  cumul_sync(typ, p_patch%id,idx)%p_patch => p_patch
  IF (idx == MAX_CUMULATIVE_SYNC) THEN
!CDIR NOIEXPAND
    CALL complete_cumulative_sync(typ, p_patch%id)
  END IF
END SUBROUTINE cumulative_sync_patch_array


!-------------------------------------------------------------------------
!> If there are any pending "cumulative sync" operations: Complete them!
!
!  @par Revision History
!  Initial revision : F. Prill, DWD (2013-02-08)
! 
RECURSIVE SUBROUTINE complete_cumulative_sync(opt_typ, opt_patch_id)
  INTEGER, OPTIONAL, INTENT(IN) :: opt_typ
  INTEGER, OPTIONAL, INTENT(IN) :: opt_patch_id
  ! local variables
  CHARACTER(*), PARAMETER :: routine = TRIM("mo_sync:complete_cumulative_sync")
  TYPE(t_patch), POINTER :: p_patch
  INTEGER :: i

  ! first, handle the case that this routine has been called for all
  ! sync types or/and all patches:
  IF (.NOT. PRESENT(opt_typ)) THEN
    DO i=1,4
!CDIR NOIEXPAND
      CALL complete_cumulative_sync(i, opt_patch_id)
    END DO
    RETURN
  END IF
  IF (.NOT. PRESENT(opt_patch_id)) THEN
    DO i=1,max_dom
!CDIR NOIEXPAND
      CALL complete_cumulative_sync(opt_typ,i)
    END DO
    RETURN
  END IF
  ! now, call "sync_patch_array_mult"
  IF (PRESENT(opt_typ) .AND. (PRESENT(opt_patch_id))) THEN
    p_patch => cumul_sync(opt_typ,opt_patch_id,1)%p_patch
    SELECT CASE(ncumul_sync(opt_typ,opt_patch_id))
    CASE (0)
      RETURN
    CASE (1)
      CALL sync_patch_array_mult(opt_typ, p_patch, 1, cumul_sync(opt_typ,opt_patch_id,1)%f3d)
    CASE (2)
      CALL sync_patch_array_mult(opt_typ, p_patch, 2, cumul_sync(opt_typ,opt_patch_id, 1)%f3d, &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d)
    CASE (3)
      CALL sync_patch_array_mult(opt_typ, p_patch, 3, cumul_sync(opt_typ,opt_patch_id,1)%f3d,  &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d, cumul_sync(opt_typ,opt_patch_id,3)%f3d)
    CASE (4)
      CALL sync_patch_array_mult(opt_typ, p_patch, 4, cumul_sync(opt_typ,opt_patch_id,1)%f3d,  &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d, cumul_sync(opt_typ,opt_patch_id,3)%f3d, &
        &                        cumul_sync(opt_typ,opt_patch_id,4)%f3d)
    CASE (5)
      CALL sync_patch_array_mult(opt_typ, p_patch, 5, cumul_sync(opt_typ,opt_patch_id,1)%f3d,  &
        &                        cumul_sync(opt_typ,opt_patch_id,2)%f3d, cumul_sync(opt_typ,opt_patch_id,3)%f3d, &
        &                        cumul_sync(opt_typ,opt_patch_id,4)%f3d, cumul_sync(opt_typ,opt_patch_id,5)%f3d)
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT
    ncumul_sync(opt_typ,opt_patch_id) = 0 ! reset sync counter
  END IF
END SUBROUTINE complete_cumulative_sync


END MODULE mo_sync

