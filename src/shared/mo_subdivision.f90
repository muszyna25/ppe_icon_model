#ifdef __xlC__
@PROCESS nosmp
@PROCESS NOOPTimize
#endif

! #ifdef __xlC__
! @PROCESS HOT
! #endif
#ifdef __PGI
!pgi$g opt=1
#endif
!>
!!               This module provides all routines for dividing patches.
!!
!!               This module provides all routines for dividing patches
!! (including interpolation state) and setting up communication.
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
!! $Id: n/a$
!!
MODULE mo_subdivision
  ! If METIS is installed, uncomment the following line
  ! (or better adjust configure to recognize that)
  !#define HAVE_METIS
  !
  !-------------------------------------------------------------------------
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !-------------------------------------------------------------------------
  !
  USE mo_kind,               ONLY: wp, i8
  USE mo_impl_constants,     ONLY: success, min_rlcell, max_rlcell,  &
    & min_rledge, max_rledge, min_rlvert, max_rlvert,                &
    & min_rlcell_int, min_rledge_int, min_rlvert_int, max_hw, max_dom, max_phys_dom
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish, message, message_text,    &
    &                              get_filename_noext

  USE mo_run_config,         ONLY: ltransport, msg_level
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_io_units,           ONLY: find_next_free_unit, filename_max
  USE mo_model_domain,       ONLY: t_patch, t_grid_cells,       &
    &                              p_patch,                     &
    &                              p_patch_local_parent,        &
    &                              t_phys_patch, p_phys_patch
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_mpi,                ONLY: p_bcast, p_send, p_recv, p_sum, p_max, p_min, &
    &                              proc_split
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_UNDEFINED, MPI_COMM_NULL
#endif
  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_test, &
    & my_process_is_mpi_seq, process_mpi_all_test_id, process_mpi_all_workroot_id, &
    & my_process_is_mpi_workroot, p_pe_work, p_n_work,                  &
    & get_my_mpi_all_id, my_process_is_mpi_parallel

  USE mo_parallel_config,       ONLY:  nproma, p_test_run, &
    & division_method, division_file_name, n_ghost_rows, div_from_file, div_geometric
    
#ifdef HAVE_METIS
  USE mo_parallel_config,    ONLY: div_metis
#endif
  USE mo_communication,      ONLY: setup_comm_pattern, blk_no, idx_no, idx_1d
  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
    & grf_bdyintp_end_c, grf_bdyintp_end_e, grf_fbk_start_c, grf_fbk_start_e, &
    & grf_bdywidth_c, grf_bdywidth_e, grf_nudgintp_start_c, grf_nudgintp_start_e
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, patch_weight, n_phys_dom
  USE mo_model_domimp_patches,ONLY: allocate_patch, deallocate_basic_patch, deallocate_patch
  USE mo_dump_restore,        ONLY: dump_all_domain_decompositions

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: decompose_domain, copy_processor_splitting
  PUBLIC :: set_patch_communicators
  PUBLIC :: finalize_decomposition
  PUBLIC :: setup_phys_patches
  PUBLIC :: set_comm_pat_gather
  PUBLIC :: complete_parallel_setup

  ! pointers to the work patches
  TYPE(t_patch), POINTER :: wrk_p_patch
  TYPE(t_patch), POINTER :: wrk_p_patch_g, wrk_p_parent_patch_g
  TYPE(t_patch), POINTER :: wrk_divide_patch

  !-------------------------------------------------------------------------
  ! Definition of local parent patches
  ! For any given patch p_patch(jg) and jgp = p_patch(jg)%parent_id,
  ! p_patch_local_parent(jg) has the same resolution as p_patch(jgp)
  ! but it covers only the area of p_patch(jgp) which is covered by its child p_patch(jg)
  ! and it is divided in the same manner as p_patch(jg).
  ! Please note that p_patch_local_parent(1) is undefined if n_dom_start = 1

  ! Please note: The definitions of the local parents are now at the same locations
  ! as the definitions of the respective patch or state
  !-------------------------------------------------------------------------

  ! Private flag if patch should be divided for radiation calculation
  LOGICAL :: divide_for_radiation = .FALSE.

  ! number of grid points of different categories:
  ! lateral points, interior points, nested points, halo points
  INTEGER(i8), PUBLIC :: npts_local(0:max_dom, 4)

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!  In case of a test run: Copies processor splitting to test PE
  SUBROUTINE copy_processor_splitting(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER :: ibuf(n_dom_start:n_dom,2)

    IF(.NOT. p_test_run) RETURN ! Nothing to do

    IF(my_process_is_mpi_workroot()) THEN
      CALL p_send(proc_split, process_mpi_all_test_id, 1)
      ibuf(:,1) = p_patch(:)%n_proc
      ibuf(:,2) = p_patch(:)%proc0
      CALL p_send(ibuf, process_mpi_all_test_id, 2)
    ENDIF

    IF(my_process_is_mpi_test()) THEN
      CALL p_recv(proc_split, process_mpi_all_workroot_id, 1)
      CALL p_recv(ibuf, process_mpi_all_workroot_id, 2)
      p_patch(:)%n_proc = ibuf(:,1)
      p_patch(:)%proc0  = ibuf(:,2)
    ENDIF

  END SUBROUTINE copy_processor_splitting
  !-------------------------------------------------------------------------
  !>
  !!  Sets the communicators in the patches if these have been read.
  SUBROUTINE set_patch_communicators(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER jc, jgc, jg, jgp, n_proc_total, comm, mpierr
    INTEGER, ALLOCATABLE :: patch_no(:)

    IF (my_process_is_mpi_seq()) THEN
      p_patch(:)%comm = p_comm_work
      RETURN
    ENDIF

#ifndef NOMPI
    ! Default if processor set is not split

    p_patch(:)%comm   = p_comm_work
    p_patch(:)%rank   = p_pe_work

    proc_split = .FALSE.

    IF(p_patch(1)%n_childdom <= 1) RETURN ! No splitting for 0 or 1 childs

    ! Check if the processor set is split for childs of root

    ALLOCATE(patch_no(0:p_n_work-1))

    n_proc_total = 0
    patch_no(:) = 0
    DO jc = 1, p_patch(1)%n_childdom
      jgc = p_patch(1)%child_id(jc)
      n_proc_total = n_proc_total + p_patch(jgc)%n_proc
      patch_no(p_patch(jgc)%proc0 : p_patch(jgc)%proc0+p_patch(jgc)%n_proc-1) = jc
    ENDDO

    ! if any processor has no patch assigned, this is an error

    IF(ANY(patch_no == 0)) &
      CALL finish('set_patch_communicators','Unknown patch split mode (1)')

    IF(n_proc_total == p_n_work) THEN

      proc_split = .TRUE.

      ! Split communicator among childs of root patch

      CALL MPI_Comm_split(p_comm_work, patch_no(p_pe_work), p_pe_work, comm, mpierr)

      ! Set comm and rank for childs of root patch

      DO jc = 1, p_patch(1)%n_childdom
        jgc = p_patch(1)%child_id(jc)
        IF(patch_no(p_pe_work) == jc) THEN
          p_patch(jgc)%comm = comm
          CALL MPI_Comm_rank(comm, p_patch(jgc)%rank, mpierr)
        ELSE
          p_patch(jgc)%comm = MPI_COMM_NULL ! We should never use this comm
          p_patch(jgc)%rank = -1
        ENDIF
      ENDDO

      ! and for deeper level descandants

      DO jg = 2, n_dom

        jgp = p_patch(jg)%parent_id

        IF(jgp /= 1) THEN
          p_patch(jg)%comm   = p_patch(jgp)%comm
          p_patch(jg)%rank   = p_patch(jgp)%rank
        ENDIF

      ENDDO

    ELSEIF(n_proc_total /= p_patch(1)%n_childdom * p_n_work) THEN

      ! If the processor is not split and n_proc_totalindicates that not
      ! every processor is working on every patch, this is an error
      CALL finish('set_patch_communicators','Unknown patch split mode (2)')

    ENDIF

    DEALLOCATE(patch_no)
#endif

  END SUBROUTINE set_patch_communicators

  !------------------------------------------------------------------
  !>
  !!  Divide patches and interpolation states for mpi parallel runs.
  SUBROUTINE decompose_domain( p_patch_global, opt_n_procs )

    TYPE(t_patch), INTENT(INOUT), TARGET :: p_patch_global(n_dom_start:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_n_procs

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_subdivision:decompose_domain")

    CHARACTER(len=20), PARAMETER, DIMENSION(4) :: summary = &
      & (/ "lateral grid points " , &
      &    "interior grid points", &
      &    "nested grid points  ",   &
      &    "halo grid points    " /)

    ! Local variables:
    INTEGER :: mpierr
    INTEGER :: ist, jg, jgp, jc, jgc, comm, n, &
      &        l1, i_nchdom, n_procs_decomp
    INTEGER :: nprocs(p_patch_global(1)%n_childdom)
    INTEGER, ALLOCATABLE :: cell_owner(:), cell_owner_p(:)
    REAL(wp) :: weight(p_patch_global(1)%n_childdom)
    INTEGER(i8) :: npts_global(4)
    ! (Optional:) Print a detailed summary on model grid points
    LOGICAL :: l_detailed_summary
    TYPE(t_patch), ALLOCATABLE, TARGET :: p_patch_out(:), p_patch_lp_out(:)

    CALL message(routine, 'start of domain decomposition')

    ! If opt_n_procs is set, this routine should be called from a single processor run
    ! and it splits the domain into opt_n_procs parts.
    ! Otherwise, the call should be from a parallel run and the domain is split
    ! according to the number of processors

    IF(present(opt_n_procs)) THEN
      n_procs_decomp = opt_n_procs
      IF(my_process_is_mpi_parallel()) &
        CALL finish(routine, 'Call with opt_n_procs in parallel mode')
    ELSE
      n_procs_decomp = p_n_work
    ENDIF


    ! Check if the processor set should be split for patches of the 1st generation.
    ! This is the case if patch_weight > 0 for at least one of the root's childs.
    ! For clearness, we require that patch_weight=0 for any other patch

    ! Get weights for 1st generation patches
    proc_split = .FALSE.
    DO jc = 1, p_patch_global(1)%n_childdom
      jgc = p_patch_global(1)%child_id(jc)
      weight(jc) = patch_weight(jgc)
      IF(weight(jc) > 0._wp) proc_split = .TRUE.
    ENDDO

    ! Check if weights for other patches are 0 (for clearness only)
    IF(patch_weight(1) /= 0._wp) &
      CALL finish(routine,'Weight for root patch must be 0')
    DO jg = 2, n_dom
      jgp = p_patch_global(jg)%parent_id
      IF(jgp /= 1 .AND. patch_weight(jg) > 0._wp) &
        CALL finish(routine,'Weight for higher level patch must be 0')
    ENDDO

    IF(proc_split) THEN

      IF(p_pe_work==0) WRITE(0,*) 'Splitting processor grid for first level patches'
      IF(p_pe_work==0) WRITE(0,'(a,10f12.3)') 'Weights for first level patches:',weight(:)

      ! In this case, the working processor set must be at least as big
      ! as the number of childs of the root patch
      IF(p_patch_global(1)%n_childdom > n_procs_decomp) &
        CALL finish(routine,'Too few procs for processor grid splitting')

      ! Get the number of procs per patch according to weight(:).
      ! Every patch gets at least 1 proc (of course!):
      nprocs(:) = 1

      ! The remaining procs are divided among patches similar to
      ! the d'Hondt method for elections

      DO n = p_patch_global(1)%n_childdom+1, n_procs_decomp
        jg = MAXLOC(weight(:)/REAL(nprocs(:)+1,wp),1)
        nprocs(jg) = nprocs(jg)+1
      ENDDO

      IF(p_pe_work==0) THEN
        WRITE(0,*) 'Processor splitting:'
        DO jc = 1, p_patch_global(1)%n_childdom
          jgc =  p_patch_global(1)%child_id(jc)
          WRITE(0,'(a,i0,a,f10.3,a,i0,a,i0)')                                              &
            &   'Patch ',jgc,' weight ',weight(jc),' gets ',nprocs(jc),' of ',n_procs_decomp
          IF (nprocs(jc) <= 1) &
            CALL finish(routine,'Processor splitting requires at least 2 PEs per patch')
        ENDDO
      ENDIF

      ! Set proc0, n_proc for all patches ...

      ! ... for the root patch and patch 0 if it exists

      p_patch(n_dom_start:1)%n_proc = n_procs_decomp
      p_patch(n_dom_start:1)%proc0  = 0

      ! ... for 1st generation childs

      n = 0
      DO jc = 1, p_patch_global(1)%n_childdom
        jgc = p_patch_global(1)%child_id(jc)
        p_patch(jgc)%proc0  = n
        p_patch(jgc)%n_proc = nprocs(jc)
        n = n + nprocs(jc)
      ENDDO

      ! ... for deeper level descandants

      DO jg = 2, n_dom

        jgp = p_patch_global(jg)%parent_id

        IF(jgp /= 1) THEN
          p_patch(jg)%n_proc = p_patch(jgp)%n_proc
          p_patch(jg)%proc0  = p_patch(jgp)%proc0
        ENDIF

      ENDDO

    ELSE

      ! No splitting, proc0, n_proc are identical for all patches

      IF(p_pe_work==0) PRINT *,'No splitting of processor grid'
      p_patch(:)%n_proc = n_procs_decomp
      p_patch(:)%proc0  = 0

    ENDIF

#ifdef NOMPI
    p_patch(:)%comm = 0
#else
    p_patch(:)%comm = MPI_COMM_NULL
#endif
    p_patch(:)%rank = -1


    IF(PRESENT(opt_n_procs)) THEN
      ! Allocate p_patch_out, p_patch_lp_out so that they can keep
      ! all domain decompositions for one patch
      ALLOCATE(p_patch_out(opt_n_procs))
      ALLOCATE(p_patch_lp_out(opt_n_procs))
    ELSE
      ! p_patch is already allocated, p_patch_local_parent needs allocation
      ALLOCATE(p_patch_local_parent(n_dom_start+1:n_dom))
    ENDIF

    ! Divide patches

    DO jg = n_dom_start, n_dom

      jgp = p_patch_global(jg)%parent_id

      IF(jg == n_dom_start) THEN
        NULLIFY(wrk_p_parent_patch_g) ! Must be NULL for global patch
      ELSE
        wrk_p_parent_patch_g => p_patch_global(jgp)
      ENDIF

      wrk_p_patch_g => p_patch_global(jg)

      ! Set division method, divide_for_radiation is only used for patch 0

      divide_for_radiation = (jg == 0)

      ! Here comes the actual domain decomposition:
      ! Every cells gets assigned an owner.

      ALLOCATE(cell_owner(p_patch_global(jg)%n_patch_cells))
      CALL divide_patch_cells(p_patch(jg)%n_proc, p_patch(jg)%proc0, cell_owner)

      IF(jg > n_dom_start) THEN
        ! Assign the cell owners of the current patch to the parent cells
        ! for the construction of the local parent:
        ALLOCATE(cell_owner_p(p_patch_global(jgp)%n_patch_cells))
        CALL divide_parent_cells(p_patch_global(jg),cell_owner,cell_owner_p)
      ENDIF

      ! Please note: Previously, for jg==0 no ghost rows were set.
      ! Currently, we need ghost rows for jg==0 also for dividing the int state and grf state
      ! Have still to check if int state/grf state is needed at all for jg==0,
      ! if this is not the case, the ghost rows can be dropped again.

      IF(PRESENT(opt_n_procs)) THEN

        ! Calculate all basic patches and output them at once.
        ! This is done in this way (calculating all before output) since
        ! we need the biggest dimension of all patches before real NetCDF output starts.
        ! All arrays in the patch which are not scaling (i.e. having global dimensions)
        ! are discarded so that the stored patches shouldn't need much more space
        ! than the global patch descriptions.

        ! If the storage for the patches should get a problem nontheless, there has to be found
        ! a way to output them before all are calculated.

        ! Do all domain decompositions
        DO n = 1, opt_n_procs

          WRITE(0,'(2(a,i0))') 'Dividing patch ',jg,' for proc ',n-1
          p_patch_out(n)%n_proc = p_patch(jg)%n_proc
          p_patch_out(n)%proc0  = p_patch(jg)%proc0
          wrk_p_patch_g => p_patch_global(jg)
          wrk_p_patch   => p_patch_out(n)
          CALL divide_patch(cell_owner, n_ghost_rows, .TRUE., n-1)
          CALL discard_large_arrays(p_patch_out(n), n)

          IF(jg > n_dom_start) THEN
            ! Divide local parent
            wrk_p_patch_g => p_patch_global(jgp)
            wrk_p_patch   => p_patch_lp_out(n)
            CALL divide_patch(cell_owner_p, 1, .FALSE., n-1)
            CALL discard_large_arrays(p_patch_lp_out(n), n)
          ENDIF

        ENDDO

        ! Dump domain decompositions to NetCDF
        IF(jg > n_dom_start) THEN
          CALL dump_all_domain_decompositions(p_patch_out, p_patch_lp_out)
        ELSE
          CALL dump_all_domain_decompositions(p_patch_out)
        ENDIF

        ! Deallocate patch arrays
        DO n = 1, opt_n_procs
          CALL deallocate_patch(p_patch_out(n))
          IF(jg > n_dom_start) CALL deallocate_patch(p_patch_lp_out(n))
        ENDDO

      ELSE

        wrk_p_patch_g => p_patch_global(jg)
        wrk_p_patch   => p_patch(jg)
        CALL divide_patch(cell_owner, n_ghost_rows, .TRUE., p_pe_work)

        IF(jg > n_dom_start) THEN
          wrk_p_patch_g => p_patch_global(jgp)
          wrk_p_patch   => p_patch_local_parent(jg)
          CALL divide_patch(cell_owner_p, 1, .FALSE., p_pe_work)
        ENDIF

      ENDIF

      DEALLOCATE(cell_owner)
      IF(jg > n_dom_start) DEALLOCATE(cell_owner_p)

    ENDDO


    ! The global patches may be discarded now

    DO jg = n_dom_start, n_dom
      CALL deallocate_basic_patch( p_patch_global(jg) )
    ENDDO

    ! (Optional:)
    ! Print a detailed summary on model grid points
    l_detailed_summary = (msg_level >= 16) .AND. .NOT.PRESENT(opt_n_procs)

    IF (l_detailed_summary) THEN
      WRITE (message_text,*) "PE ", get_my_mpi_all_id(), &
        &                    "Detailed grid summary (cells)"
      CALL message(routine, TRIM(message_text))

      DO jg = n_dom_start, n_dom
        ! count grid points for this PE:
        i_nchdom     = MAX(1,p_patch(jg)%n_childdom)
        ! local, lateral grid points
        npts_local(jg,1)    = count_entries(  &
          &                   p_patch(jg)%cells%start_blk(1,1),         &
          &                   p_patch(jg)%cells%start_idx(1,1),         &
          &                   p_patch(jg)%cells%end_blk(max_rlcell,1),  &
          &                   p_patch(jg)%cells%end_idx(max_rlcell,1) )
        ! local, interior grid points
        npts_local(jg,2)    = count_entries(  &
          &                   p_patch(jg)%cells%start_blk(0,1), &
          &                   p_patch(jg)%cells%start_idx(0,1), &
          &                   p_patch(jg)%cells%end_blk(0,1),   &
          &                   p_patch(jg)%cells%end_idx(0,1) )
        ! local, nested grid points:
        npts_local(jg,3)    = count_entries(  &
          &                   p_patch(jg)%cells%start_blk(-1,1), &
          &                   p_patch(jg)%cells%start_idx(-1,1), &
          &                   p_patch(jg)%cells%end_blk(min_rlcell_int,i_nchdom),        &
          &                   p_patch(jg)%cells%end_idx(min_rlcell_int,i_nchdom) )
        ! local, halo grid points:
        npts_local(jg,4)    = count_entries(  &
          &                   p_patch(jg)%cells%start_blk(min_rlcell_int-1,1), &
          &                   p_patch(jg)%cells%start_idx(min_rlcell_int-1,1), &
          &                   p_patch(jg)%cells%end_blk(min_rlcell,i_nchdom),  &
          &                   p_patch(jg)%cells%end_idx(min_rlcell,i_nchdom) )
        ! sum up over all PEs (collective operation):
        npts_global(:) = p_sum(npts_local(jg,:), p_comm_work)

        WRITE (message_text,'(A8,i4)') "patch # ", jg
        CALL message(routine, TRIM(message_text))
        DO l1=1,4
          WRITE (message_text, '(A25,i6)') ">   "//summary(l1)//":", npts_local(jg,l1)
          CALL message(routine, TRIM(message_text))
        END DO
        WRITE (message_text,'(A20,i4)') "global values, patch # ", jg
        CALL message(routine, TRIM(message_text))
        DO l1=1,4
          WRITE (message_text, '(A25,i6)') ">   "//summary(l1)//":", npts_global(l1)
          CALL message(routine, TRIM(message_text))
        END DO

      END DO
    END IF

  END SUBROUTINE decompose_domain

  !-----------------------------------------------------------------------------
  !> discard_large_arrays:
  !! Deallocate large arrays in patch which are not output to NetCDF anyways

  SUBROUTINE discard_large_arrays(p, n)

    TYPE(t_patch), INTENT(INOUT) :: p
    INTEGER, INTENT(IN) :: n

    ! loc_index is never output to NetCDF

    DEALLOCATE(p%cells%loc_index)
    DEALLOCATE(p%edges%loc_index)
    DEALLOCATE(p%verts%loc_index)
    ! Allocate it again so that we don't get problems during patch destruction
    ALLOCATE(p%cells%loc_index(1))
    ALLOCATE(p%edges%loc_index(1))
    ALLOCATE(p%verts%loc_index(1))

    ! owner_g is identical everywhere and output only for the first patch

    IF(n>1) THEN
      DEALLOCATE(p%cells%owner_g)
      DEALLOCATE(p%edges%owner_g)
      DEALLOCATE(p%verts%owner_g)
      ! Allocate it again
      ALLOCATE(p%cells%owner_g(1))
      ALLOCATE(p%edges%owner_g(1))
      ALLOCATE(p%verts%owner_g(1))
    ENDIF

  END SUBROUTINE discard_large_arrays

  !-----------------------------------------------------------------------------

  SUBROUTINE complete_parallel_setup

    INTEGER :: jg, jgp

    DO jg = n_dom_start, n_dom

      jgp = p_patch(jg)%parent_id

      ! Set communication patterns for boundary exchange
      CALL set_comm_pat_bound_exch(p_patch(jg))

      ! Set communication patterns for gathering on proc 0
      CALL set_comm_pat_gather(p_patch(jg))

      CALL set_owner_mask(p_patch(jg))
     
     ! Fill the owner_local value
      ! this is done in the set_owner_mask
      ! CALL fill_owner_local(p_patch(jg))

      IF(jg == n_dom_start) THEN

        ! parent_idx/blk is set to 0 since it just doesn't exist,
        ! child_idx/blk is set to 0 since it makes sense only on the local parent
        p_patch(jg)%cells%parent_idx = 0
        p_patch(jg)%cells%parent_blk = 0
        p_patch(jg)%cells%child_idx  = 0
        p_patch(jg)%cells%child_blk  = 0
        p_patch(jg)%edges%parent_idx = 0
        p_patch(jg)%edges%parent_blk = 0
        p_patch(jg)%edges%child_idx  = 0
        p_patch(jg)%edges%child_blk  = 0

      ELSE

        CALL setup_comm_cpy_interpolation(p_patch(jg), p_patch(jgp))
        CALL setup_comm_grf_interpolation(p_patch(jg), p_patch(jgp))
        CALL setup_comm_ubc_interpolation(p_patch(jg), p_patch(jgp))

        CALL set_comm_pat_bound_exch(p_patch_local_parent(jg))
        CALL set_comm_pat_gather(p_patch_local_parent(jg))

        CALL set_parent_child_relations(p_patch_local_parent(jg), p_patch(jg))

        CALL set_glb_loc_comm(p_patch(jgp), p_patch_local_parent(jg), &
          &                   p_patch(jg)%parent_child_index)

        CALL set_owner_mask(p_patch_local_parent(jg))
      ENDIF

    ENDDO

    CALL set_patch_communicators(p_patch)

  END SUBROUTINE complete_parallel_setup

  !-----------------------------------------------------------------------------

  SUBROUTINE finalize_decomposition

    implicit none
    integer jg

    ! Remap indices in patches and local parents

    DO jg = n_dom_start, n_dom

      CALL remap_patch_indices(p_patch(jg))

      IF(jg>n_dom_start) THEN
        CALL remap_patch_indices(p_patch_local_parent(jg))
      ENDIF

    ENDDO
                 
  END SUBROUTINE finalize_decomposition

  !-----------------------------------------------------------------------------
  !>
  ! Fills the in_patch%cells%owner_local using the in_patch%cells%owner_g
  ! Note: At the moment it uses the p_work_pe number which is not the same
  ! as the my_mpi_all_id. It requires 
!   SUBROUTINE fill_owner_local(in_patch)
! 
!     TYPE(t_patch), INTENT(inout) :: in_patch
! 
!     INTEGER :: local_cell_idx, global_cell_idx
!     INTEGER :: i, jb, jl, jb_e, jl_e, jb_v, jl_v, jv, je
!     INTEGER :: owner_id
! 
!     in_patch%edges%owner_local(:) = -1
!     in_patch%verts%owner_local(:) = -1
! 
!     DO local_cell_idx = 1, in_patch%n_patch_cells
!       global_cell_idx = in_patch%cells%glb_index(local_cell_idx)
!       owner_id = in_patch%cells%owner_g(global_cell_idx)
!       in_patch%cells%owner_local(local_cell_idx) = in_patch%cells%owner_g(global_cell_idx)
!       IF (owner_id < 0) CYCLE
! 
!       ! go around the cell edges mark the owner
!       jb = blk_no(local_cell_idx) ! block index
!       jl = idx_no(local_cell_idx) ! line index
!       DO i = 1,in_patch%cells%num_edges(jl,jb)
!         jl_e = in_patch%cells%edge_idx(jl,jb,i)
!         jb_e = in_patch%cells%edge_blk(jl,jb,i)
!         je = idx_1d(jl_e, jb_e)
!         jl_v = in_patch%cells%vertex_idx(jl,jb,i)
!         jb_v = in_patch%cells%vertex_blk(jl,jb,i)
!         jv = idx_1d(jl_v, jb_v)
!         IF (owner_id == p_pe_work) THEN
!           ! Assume we own the edges and vertices of the owned cells
!           ! No process can claim actual ownershipe, as these can be calculated
!           ! only using the halos. No reason to communicate them
!           in_patch%edges%owner_local(je) = p_pe_work
!           in_patch%verts%owner_local(jv) = p_pe_work
!         ELSEIF (in_patch%edges%owner_local(je) < 0) THEN
!           in_patch%edges%owner_local(je) = owner_id
!           in_patch%verts%owner_local(jv) = owner_id
!         ENDIF
! 
!       ENDDO
! 
!     ENDDO
! 
! 
!   END SUBROUTINE fill_owner_local
 
  !-----------------------------------------------------------------------------
  !>
  !! Sets the owner mask
  SUBROUTINE set_owner_mask(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch

    INTEGER :: j, jb, jl, jg

    p_patch%cells%owner_mask = .false.
    p_patch%edges%owner_mask = .false.
    p_patch%verts%owner_mask = .false.
      
    p_patch%cells%owner_local(:) = -1
    p_patch%edges%owner_local(:) = -1
    p_patch%verts%owner_local(:) = -1

    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      jg = p_patch%cells%glb_index(j) ! index in global patch

      p_patch%cells%owner_mask(jl,jb) = p_patch%cells%owner_g(jg)==p_pe_work
      ! fill local owner
      p_patch%cells%owner_local(j) = p_patch%cells%owner_g(jg)

    ENDDO

    DO j = 1, p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      jg = p_patch%edges%glb_index(j) ! index in global patch

      p_patch%edges%owner_mask(jl,jb) = p_patch%edges%owner_g(jg)==p_pe_work
      ! fill local owner
      p_patch%edges%owner_local(j) = p_patch%edges%owner_g(jg)
    
    ENDDO

    DO j = 1, p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      jg = p_patch%verts%glb_index(j) ! index in global patch

      p_patch%verts%owner_mask(jl,jb) = p_patch%verts%owner_g(jg)==p_pe_work
      ! fill local owner
      p_patch%verts%owner_local(j) = p_patch%verts%owner_g(jg)

    ENDDO

  END SUBROUTINE set_owner_mask
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !>
  !! Divides the cells of a global patch (in wrk_p_patch_g) for parallelization.
  !!
  !! Outputs the subdivsion in cell_owner(:) which must already be allocated
  !! to the correct size (wrk_p_patch_g%n_patch_cells)
  !!
  !! If  wrk_p_parent_patch_g is associated, this indicates that the patch has a parent
  !! which has consquences for subdivision (cell with the same parent must not
  !! get to different PEs).
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Split out as a separate routine, Rainer Johanni, Oct 2010

  SUBROUTINE divide_patch_cells(n_proc, proc0, cell_owner)

    INTEGER, INTENT(IN)  :: n_proc !> Number of processors for split
    INTEGER, INTENT(IN)  :: proc0  !> First processor of patch
    INTEGER, INTENT(OUT) :: cell_owner(:) !> Cell division

    INTEGER :: n, i, j, jl, jb, jl_p, jb_p
    INTEGER, ALLOCATABLE :: flag_c(:), tmp(:)
    CHARACTER(LEN=filename_max) :: use_division_file_name ! if div_from_file

    ! Please note: Unfortunatly we cannot use p_io for doing I/O,
    ! since this might be the test PE which is never calling this routine
    ! (this is the case in the actual setup).
    ! Thus we use the worker PE 0 for I/O and don't use message() for output.


    IF(division_method==div_from_file) THEN

      ! Area subdivision is read from file

      IF(p_pe_work == 0) THEN

        IF (division_file_name == "") THEN
          use_division_file_name = &
            & TRIM(get_filename_noext(wrk_p_patch_g%grid_filename))//'.cell_domain_ids'
        ELSE
          use_division_file_name = division_file_name
        ENDIF
        
        WRITE(0,*) "Read decomposition from file: ", TRIM(use_division_file_name)
        n = find_next_free_unit(10,99)

        OPEN(n,FILE=TRIM(use_division_file_name),STATUS='OLD',IOSTAT=i)
        IF(i /= 0) CALL finish('divide_patch',&
          & 'Unable to open input file: '//TRIM(use_division_file_name))

        DO j = 1, wrk_p_patch_g%n_patch_cells
          READ(n,*,IOSTAT=i) cell_owner(j)
          IF(i /= 0) CALL finish('divide_patch','Error reading: '//TRIM(use_division_file_name))
        ENDDO
        CLOSE(n)

        ! Quick check for correct values

        IF(MINVAL(cell_owner(:)) < 0 .or. MAXVAL(cell_owner(:)) >= n_proc) THEN
          WRITE(0,*) "n_porc=",n_proc, " MINAVAL=", MINVAL(cell_owner(:)), &
            " MAXVAL=", MAXVAL(cell_owner(:))
          CALL finish('divide_patch','Illegal subdivision in input file')
        ENDIF

      ENDIF

      CALL p_bcast(cell_owner, 0, comm=p_comm_work)

      IF(p_pe_work==0) PRINT *,'Successfully read: '//TRIM(division_file_name)

    ELSE

      ! Build in subdivison methods

      IF(ASSOCIATED(wrk_p_parent_patch_g)) THEN

        ! Cells with the same parent must not go to different PEs.
        ! Thus we have to divide in reality the subset of the parent cells
        ! which cover the actual patch and then assign the ownership
        ! of the cells of the actual patch according to the parent cells

        ! Determine the subset of the parent patch covering the actual patch
        ! by flagging the according cells

        ALLOCATE(flag_c(wrk_p_parent_patch_g%n_patch_cells))
        flag_c(:) = 0

        DO j = 1, wrk_p_patch_g%n_patch_cells

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          jl_p = wrk_p_patch_g%cells%parent_idx(jl,jb)
          jb_p = wrk_p_patch_g%cells%parent_blk(jl,jb)

          flag_c(idx_1d(jl_p, jb_p)) = 1
        ENDDO

        ! Divide subset of patch
        ! Receives the PE  numbers for every cell
        ALLOCATE(tmp(wrk_p_parent_patch_g%n_patch_cells))

        IF(division_method==div_geometric) THEN
          wrk_divide_patch => wrk_p_parent_patch_g
          CALL divide_subset_geometric( flag_c, n_proc, tmp)
#ifdef HAVE_METIS
        ELSE IF(division_method==div_metis) THEN
          wrk_divide_patch => wrk_p_parent_patch_g
          CALL divide_subset_metis( flag_c, n_proc, tmp)
#endif
        ELSE
          CALL finish('divide_patch','Illegal division_method setting')
        ENDIF

        ! Owners of the cells of the parent patch are now in tmp.
        ! Set owners in current patch from this

        DO j = 1, wrk_p_patch_g%n_patch_cells

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          jl_p = wrk_p_patch_g%cells%parent_idx(jl,jb)
          jb_p = wrk_p_patch_g%cells%parent_blk(jl,jb)

          cell_owner(j) = tmp(idx_1d(jl_p,jb_p))

        ENDDO

        DEALLOCATE(flag_c, tmp)

      ELSE

        ! No parent patch, simply divide current patch

        ! Set subset flags where the "subset" is the whole patch

        ALLOCATE(flag_c(wrk_p_patch_g%n_patch_cells))
        flag_c(:) = 1

        ! Divide complete patch

        IF(division_method==div_geometric) THEN
          wrk_divide_patch => wrk_p_patch_g
          CALL divide_subset_geometric(flag_c, n_proc, cell_owner)
#ifdef HAVE_METIS
        ELSE IF(division_method==div_metis) THEN
          wrk_divide_patch => wrk_p_patch_g
          CALL divide_subset_metis(flag_c, n_proc, cell_owner)
#endif
        ELSE
          CALL finish('divide_patch','Illegal division_method setting')
        ENDIF

        DEALLOCATE(flag_c)

      ENDIF

    ENDIF

    ! Set processor offset
    cell_owner(:) = cell_owner(:) + proc0

    ! Output how many cells go to every PE

  !  IF(p_pe_work==0) THEN
  !    PRINT '(a,i0,a,i0)','Patch: ',wrk_p_patch_g%id,&
  !      & ' Total number of cells: ',wrk_p_patch_g%n_patch_cells
  !    DO n = 0, p_n_work-1
  !      PRINT '(a,i5,a,i8)','PE',n,' # cells: ',COUNT(cell_owner(:) == n)
  !    ENDDO
  !  ENDIF

  END SUBROUTINE divide_patch_cells

  !-------------------------------------------------------------------------------------------------
  !>
  !! Sets the owner for the division of the cells of the parent patch
  !! with the same subdivision as the child

  SUBROUTINE divide_parent_cells(p_patch_g, cell_owner, cell_owner_p)

    TYPE(t_patch), INTENT(IN) :: p_patch_g  !> Global patch for which the parent should be divided
    INTEGER, INTENT(IN)  :: cell_owner(:)   !> Ownership of cells in p_patch_g
    INTEGER, INTENT(OUT) :: cell_owner_p(:) !> Output: Cell division for parent.
                                            !> Must be allocated to n_patch_cells of the global parent

    INTEGER :: j, jl, jb, jp
    INTEGER :: cnt(SIZE(cell_owner_p))

    cell_owner_p(:) = -1
    cnt(:) = 0

    DO j = 1, p_patch_g%n_patch_cells

      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index
      jp = idx_1d(p_patch_g%cells%parent_idx(jl,jb), p_patch_g%cells%parent_blk(jl,jb))

      IF(cell_owner_p(jp) < 0) THEN
        cell_owner_p(jp) = cell_owner(j)
      ELSEIF(cell_owner_p(jp) /= cell_owner(j)) THEN
        CALL finish('divide_parent_cells','Divided parent cell encountered')
      ENDIF
      cnt(jp) = cnt(jp) + 1 ! Only for safety check below
    ENDDO

    ! Safety check
    IF(ANY(cnt(:)/=0 .AND. cnt(:)/=4)) &
      & CALL finish('divide_parent_cells','Incomplete parent cell encountered')

  END SUBROUTINE divide_parent_cells

  !-------------------------------------------------------------------------------------------------
  !>
  !! Divides a global patch (in wrk_p_patch_g) for parallelization.
  !!
  !! Parameters:
  !! cell_owner          owner PE numbers of the global cells
  !! n_boundary_rows     number of boundary rows to be added
  !! l_compute_grid      if true, a "normal" grid for prognstic computations is processed,
  !!                     if false, a local parent grid is processed
  !!                     in the first case, a finer distinction between different halo
  !!                     cell/edge/vertex levels is made to optimize communication
  !!
  !! On exit, the entries of wrk_p_patch are set.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Changed for usage for parent patch division, Rainer Johanni, Oct 2010

  SUBROUTINE divide_patch(cell_owner, n_boundary_rows, l_compute_grid, my_proc)

    INTEGER, INTENT(IN) :: cell_owner(:)
    INTEGER, INTENT(IN) :: n_boundary_rows
    LOGICAL, INTENT(IN) :: l_compute_grid
    INTEGER, INTENT(IN) :: my_proc

    INTEGER :: n, i, j, jv, je, jl, jb, jl_g, jb_g, jl_e, jb_e, jl_v, jb_v, ilev, iown,   &
               ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, jl_c, jb_c, jc, irlev, ilev1, ilev_st, &
               jg, i_nchdom, irl0, irl1, irl2, irl3

    INTEGER, ALLOCATABLE :: flag_c(:), flag_e(:), flag_v(:)
    INTEGER, ALLOCATABLE :: flag2_c(:), flag2_e(:), flag2_v(:)
    LOGICAL, ALLOCATABLE :: lcount_c(:), lcount_e(:), lcount_v(:)

    !-----------------------------------------------------------------------------------------------
    ! Find inner cells/edges/verts and ghost rows for our patch:
    ! flag_c/e/v = 0 is set for inner cells/edges/verts
    ! flag_c/e/v > 0 is set for ghost rows (counting the level of displacement)
    ! flag_c/e/v = -1 for cells/edges/verts which don't belong to our patch at all
    !-----------------------------------------------------------------------------------------------

    ALLOCATE(flag_c(wrk_p_patch_g%n_patch_cells))
    ALLOCATE(flag_e(wrk_p_patch_g%n_patch_edges))
    ALLOCATE(flag_v(wrk_p_patch_g%n_patch_verts))

    flag_c(:) = -1
    flag_e(:) = -1
    flag_v(:) = -1

    !-----------------------------------------------------------------------------------------------
    ! The purpose of the second set of flags is to control moving the halo points
    ! to the end of the index vector for nearly all cells even if l_compute_grid = true
    ! flag2_c/e/v = 0 (-1) wherever flag_c/e/v = 0 (-1)
    ! flag2_c = 2*flag_c-1 if the cell has a common edge with a cell having flag_c-1
    ! flag2_c = 2*flag_c if the cell has no common edge with a cell having flag_c-1
    ! flag2_e/v = 1 for boundary edges/vertices not owned by the current PE,
    !             otherwise flag2_v = flag_v+1
    ! flag2_e is the sum of flag_c of the two neighboring cells and 2*flag_c+1 at outer boundary edges
    !-----------------------------------------------------------------------------------------------

    ALLOCATE(flag2_c(wrk_p_patch_g%n_patch_cells))
    ALLOCATE(flag2_e(wrk_p_patch_g%n_patch_edges))
    ALLOCATE(flag2_v(wrk_p_patch_g%n_patch_verts))

    ALLOCATE(lcount_c(wrk_p_patch_g%n_patch_cells))
    ALLOCATE(lcount_e(wrk_p_patch_g%n_patch_edges))
    ALLOCATE(lcount_v(wrk_p_patch_g%n_patch_verts))

    flag2_c(:) = -1
    flag2_e(:) = -1
    flag2_v(:) = -1

    lcount_c(:) = .FALSE.
    lcount_e(:) = .FALSE.
    lcount_v(:) = .FALSE.

    ! flag inner cells

    WHERE(cell_owner(:)==my_proc) flag_c(:) = 0
    WHERE(cell_owner(:)==my_proc) flag2_c(:) = 0

    jg = wrk_p_patch_g%id
    i_nchdom = MAX(1,wrk_p_patch_g%n_childdom)

    ! find inner edges/verts and ghost cells/edges/verts

    DO ilev = 0, n_boundary_rows

      ! Flag cells belonging to this level.
      ! Cells belonging to the kernel (ilev==0) are already flagged.
      ! The patch always needs a complete halo, even for local parents
      ! where some of the boundary cells have no global owner.

      IF(ilev>0) THEN

        DO j = 1, wrk_p_patch_g%n_patch_cells

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index

          IF(flag_c(j)<0) THEN

            ! Check if any vertex of this cell is already flagged.
            ! If this is the case, this cell goes to level ilev

            DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
              jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
              jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
              jv = idx_1d(jl_v, jb_v)
              IF(flag_v(jv)>=0) flag_c(j) = ilev
              ! TEST only, must never be triggered!
              ! IF(flag_v(jv)>=0 .and. flag_v(jv) /= ilev-1) &
              !  & PRINT *,'ERR: lev=',ilev,'flag_v=',flag_v(jv)
            ENDDO
          ENDIF
        ENDDO

      ENDIF

      ! Flag edges/verts belongig to this level.
      ! An edge/vert is flagged in this level if it belongs to a cell in this level
      ! and is not yet flagged (in which case it already belongs to a lower level).

      DO j = 1, wrk_p_patch_g%n_patch_cells

        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index

        IF(flag_c(j)==ilev) THEN
          DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
            jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
            jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
            je = idx_1d(jl_e, jb_e)
            IF(flag_e(je)<0) flag_e(je) = ilev
            jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
            jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
            jv = idx_1d(jl_v, jb_v)
            IF(flag_v(jv)<0) flag_v(jv) = ilev
          ENDDO
        ENDIF
      ENDDO

    ENDDO

    ! Now compute second set of flags
    DO ilev = 1, n_boundary_rows

      DO j = 1, wrk_p_patch_g%n_patch_cells

        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index

        IF (flag_c(j) == ilev) THEN

          flag2_c(j) = 2*flag_c(j)

          ! Check if any of the edges borders to a cell with flag_c(j) = ilev-1
          ! In this case, flag2_c = 2*flag_c-1

          DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
            jl_c = wrk_p_patch_g%cells%neighbor_idx(jl,jb,i)
            jb_c = wrk_p_patch_g%cells%neighbor_blk(jl,jb,i)
            jc = idx_1d(jl_c,jb_c)
            IF (jc < 1 .OR. jc > wrk_p_patch_g%n_patch_cells) CYCLE
            IF (flag_c(jc) == ilev-1) flag2_c(j) = 2*flag_c(j)-1
          ENDDO
        ENDIF
      ENDDO

    ENDDO

    ! Preset edges and vertices bordering to interior cells with 0
    DO j = 1, wrk_p_patch_g%n_patch_cells

      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index

      IF (flag_c(j)==0) THEN
        DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
          jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
          jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
          je = idx_1d(jl_e,jb_e)
          flag2_e(je) = 0

          jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
          jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
          jv = idx_1d(jl_v,jb_v)
          flag2_v(jv) = 0
        ENDDO
      ENDIF

    ENDDO

    DO ilev = 1, n_boundary_rows
      DO j = 1, wrk_p_patch_g%n_patch_cells

        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index

        IF (flag_c(j)==ilev) THEN
          DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
            jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
            jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
            je = idx_1d(jl_e,jb_e)
            jl_c = wrk_p_patch_g%cells%neighbor_idx(jl,jb,i)
            jb_c = wrk_p_patch_g%cells%neighbor_blk(jl,jb,i)
            jc = idx_1d(jl_c,jb_c)
            IF (jc < 1 .OR. jc > wrk_p_patch_g%n_patch_cells) THEN
              flag2_e(je) = 2*flag_c(j)+1
            ELSE IF (flag_c(jc) == -1) THEN
              flag2_e(je) = 2*flag_c(j)+1
            ELSE
              flag2_e(je) = flag_c(j)+flag_c(jc)
            ENDIF
            jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
            jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
            jv = idx_1d(jl_v,jb_v)
            flag2_v(jv) = flag_v(jv)+1
          ENDDO
        ENDIF

      ENDDO
    ENDDO

    !-----------------------------------------------------------------------------------------------
    ! Get the number of cells/edges/verts and other data for patch allocation
    !-----------------------------------------------------------------------------------------------

    wrk_p_patch%n_patch_cells = COUNT(flag_c(:)>=0)
    wrk_p_patch%n_patch_edges = COUNT(flag_e(:)>=0)
    wrk_p_patch%n_patch_verts = COUNT(flag_v(:)>=0)

    ! save the number of cells/edges/verts of the global patch
    wrk_p_patch%n_patch_cells_g = wrk_p_patch_g%n_patch_cells
    wrk_p_patch%n_patch_edges_g = wrk_p_patch_g%n_patch_edges
    wrk_p_patch%n_patch_verts_g = wrk_p_patch_g%n_patch_verts
    !
    ! calculate and save values for the blocking, these are needed for patch allocation.
    ! NB: Avoid the case nblks=0 for empty patches, this might cause troubles
    ! if a empty patch is used somewhere (and npromz gets wrong in the formulas below).
    !
    ! ... for the cells
    wrk_p_patch%nblks_c       = blk_no(wrk_p_patch%n_patch_cells)
    wrk_p_patch%nblks_int_c   = wrk_p_patch%nblks_c
    wrk_p_patch%npromz_c      = wrk_p_patch%n_patch_cells - (wrk_p_patch%nblks_c - 1)*nproma
    wrk_p_patch%npromz_int_c  = wrk_p_patch%npromz_c

    ! ... for the edges
    wrk_p_patch%nblks_e       = blk_no(wrk_p_patch%n_patch_edges)
    wrk_p_patch%nblks_int_e   = wrk_p_patch%nblks_e
    wrk_p_patch%npromz_e      = wrk_p_patch%n_patch_edges - (wrk_p_patch%nblks_e - 1)*nproma
    wrk_p_patch%npromz_int_e  = wrk_p_patch%npromz_e

    ! ... for the vertices
    wrk_p_patch%nblks_v       = blk_no(wrk_p_patch%n_patch_verts)
    wrk_p_patch%nblks_int_v   = wrk_p_patch%nblks_v
    wrk_p_patch%npromz_v      = wrk_p_patch%n_patch_verts - (wrk_p_patch%nblks_v - 1)*nproma
    wrk_p_patch%npromz_int_v  = wrk_p_patch%npromz_v

    ! Also needed for patch allocation
    wrk_p_patch%max_childdom  = wrk_p_patch_g%max_childdom

    ! Set other scalar members of patch here too ..
    wrk_p_patch%grid_filename = wrk_p_patch_g%grid_filename
    wrk_p_patch%level = wrk_p_patch_g%level
    wrk_p_patch%id    = wrk_p_patch_g%id
    wrk_p_patch%parent_id = wrk_p_patch_g%parent_id
    wrk_p_patch%parent_child_index = wrk_p_patch_g%parent_child_index
    wrk_p_patch%child_id(:) = wrk_p_patch_g%child_id(:)
    wrk_p_patch%child_id_list(:) = wrk_p_patch_g%child_id_list(:)
    wrk_p_patch%n_childdom = wrk_p_patch_g%n_childdom
    wrk_p_patch%n_chd_total = wrk_p_patch_g%n_chd_total
    wrk_p_patch%nlev   = wrk_p_patch_g%nlev
    wrk_p_patch%nlevp1 = wrk_p_patch_g%nlevp1
    wrk_p_patch%nshift = wrk_p_patch_g%nshift
    wrk_p_patch%nshift_total = wrk_p_patch_g%nshift_total
    wrk_p_patch%nshift_child = wrk_p_patch_g%nshift_child

    !-----------------------------------------------------------------------------------------------
    ! Allocate all data arrays in patch
    !-----------------------------------------------------------------------------------------------

    CALL allocate_patch(wrk_p_patch)

    !-----------------------------------------------------------------------------------------------
    ! Set the global ownership for cells, edges and verts (needed for boundary exchange).
    ! Please note that opposed to cells, the global owner for edges/verts is
    ! not unique at the boundaries of the divided patch.
    ! To minimize load imbalance at runtime, we set the owner to the PE with the highest processor number
    ! which is participating at this edge/vert if both PE numbers are even or odd, otherwise
    ! the lowest processor number is chosen
    !-----------------------------------------------------------------------------------------------

    wrk_p_patch%cells%owner_g(:) = cell_owner(:)

    wrk_p_patch%edges%owner_g(:) = -1
    wrk_p_patch%verts%owner_g(:) = -1

    DO j = 1, wrk_p_patch_g%n_patch_cells
      iown = wrk_p_patch%cells%owner_g(j)
      IF(iown<0) CYCLE

      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index
      DO i = 1,wrk_p_patch_g%cells%num_edges(jl,jb)
        jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
        jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
        je = idx_1d(jl_e, jb_e)
        IF (wrk_p_patch%edges%owner_g(je) < 0) THEN
          wrk_p_patch%edges%owner_g(je) = iown
        ELSE IF (MOD(iown,2)==0 .AND. MOD(wrk_p_patch%edges%owner_g(je),2)==0 .OR. &
                 MOD(iown,2)==1 .AND. MOD(wrk_p_patch%edges%owner_g(je),2)==1) THEN
          wrk_p_patch%edges%owner_g(je) = MAX(iown,wrk_p_patch%edges%owner_g(je))
        ELSE
          wrk_p_patch%edges%owner_g(je) = MIN(iown,wrk_p_patch%edges%owner_g(je))
        ENDIF
        jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
        jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
        jv = idx_1d(jl_v, jb_v)
        IF (wrk_p_patch%verts%owner_g(jv) < 0) THEN
          wrk_p_patch%verts%owner_g(jv) = iown
        ELSE IF (MOD(iown,2)==0 .AND. MOD(wrk_p_patch%verts%owner_g(jv),2)==0 .OR. &
                 MOD(iown,2)==1 .AND. MOD(wrk_p_patch%verts%owner_g(jv),2)==1) THEN
          wrk_p_patch%verts%owner_g(jv) = MAX(iown,wrk_p_patch%verts%owner_g(jv))
        ELSE
          wrk_p_patch%verts%owner_g(jv) = MIN(iown,wrk_p_patch%verts%owner_g(jv))
        ENDIF
      ENDDO
    ENDDO

    ! Reset flag2_e/v to 0 at boundary points owned by the current PE
    DO j = 1, wrk_p_patch_g%n_patch_cells
      IF (flag_c(j) == 0) THEN
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        DO i = 1,wrk_p_patch_g%cells%num_edges(jl,jb)
          jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
          jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
          je = idx_1d(jl_e,jb_e)
          IF (wrk_p_patch%edges%owner_g(je) == my_proc) flag2_e(je)=0
          IF (.NOT.l_compute_grid .AND. flag2_e(je)==1) flag2_e(je)=0

          jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
          jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
          jv = idx_1d(jl_v, jb_v)
          IF (wrk_p_patch%verts%owner_g(jv) == my_proc) flag2_v(jv)=0
          IF (.NOT.l_compute_grid .AND. flag2_v(jv)==1) flag2_v(jv)=0
        ENDDO
      ENDIF
    ENDDO

    !-----------------------------------------------------------------------------------------------
    ! Get the indices of local cells/edges/verts within global patch and vice versa.
    !-----------------------------------------------------------------------------------------------

    n = 0
    IF (.NOT. l_compute_grid) THEN
      DO j = 1, wrk_p_patch_g%n_patch_cells
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF (flag_c(j)==0) THEN
          n = n + 1
          wrk_p_patch%cells%glb_index(n) = j
          wrk_p_patch%cells%loc_index(j) = n
          lcount_c(j) = .TRUE.
        ELSE
          wrk_p_patch%cells%loc_index(j) = -(n+1)
        ENDIF
      ENDDO
    ELSE
      DO j = 1, wrk_p_patch_g%n_patch_cells
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        irl0 = wrk_p_patch_g%cells%refin_ctrl(jl,jb)
        IF (flag2_c(j)==0 .OR. flag2_c(j)>=1 .AND. irl0==1 .OR. &
          & flag2_c(j)==1 .AND. irl0==2) THEN
          n = n + 1
          wrk_p_patch%cells%glb_index(n) = j
          wrk_p_patch%cells%loc_index(j) = n
          lcount_c(j) = .TRUE.
        ELSE
          wrk_p_patch%cells%loc_index(j) = -(n+1)
        ENDIF
      ENDDO
    ENDIF

    ! Set start_idx/blk ... end_idx/blk for cells.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the loc_index array.

    DO j=1,wrk_p_patch%max_childdom
      ! Note: the segments between min_rlcell and min_rlcell_int-1 are reserved for
      ! halo cells; they are set below
      DO i=min_rlcell_int,max_rlcell
        CALL get_local_index(wrk_p_patch%cells%loc_index, &
          & wrk_p_patch_g%cells%start_idx(i,j),           &
          & wrk_p_patch_g%cells%start_blk(i,j),           &
          & wrk_p_patch%cells%start_idx(i,j),             &
          & wrk_p_patch%cells%start_blk(i,j),             &
          & +1 )
        CALL get_local_index(wrk_p_patch%cells%loc_index, &
          & wrk_p_patch_g%cells%end_idx(i,j),             &
          & wrk_p_patch_g%cells%end_blk(i,j),             &
          & wrk_p_patch%cells%end_idx(i,j),               &
          & wrk_p_patch%cells%end_blk(i,j),               &
          & -1 )
      ENDDO
      ! Preset remaining index sections with dummy values
      wrk_p_patch%cells%start_idx(min_rlcell:min_rlcell_int-1,j) = -9999
      wrk_p_patch%cells%start_blk(min_rlcell:min_rlcell_int-1,j) = -9999
      wrk_p_patch%cells%end_idx(min_rlcell:min_rlcell_int-1,j)   = -9999
      wrk_p_patch%cells%end_blk(min_rlcell:min_rlcell_int-1,j)   = -9999
    ENDDO

    IF(.NOT.l_compute_grid) THEN
      ! Gather all halo cells at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlcell_int-1
      irlev = min_rlcell_int-1
      DO j = 1, wrk_p_patch_g%n_patch_cells
        IF(flag2_c(j)>0) THEN
          n = n + 1
          wrk_p_patch%cells%glb_index(n) = j
          wrk_p_patch%cells%loc_index(j) = n
          jb_c = blk_no(n) ! block index
          jl_c = idx_no(n) ! line index
          ! This ensures that just the first point found at this irlev is saved as start point
          IF (wrk_p_patch%cells%start_idx(irlev,1)==-9999) THEN
            wrk_p_patch%cells%start_idx(irlev,:) = jl_c
            wrk_p_patch%cells%start_blk(irlev,:) = jb_c
          ENDIF
        ENDIF
      ENDDO
      ! Set end index when the last point is processed
      jb_c = blk_no(n) ! block index
      jl_c = idx_no(n) ! line index
      wrk_p_patch%cells%end_idx(min_rlcell:irlev,:) = jl_c
      wrk_p_patch%cells%end_blk(min_rlcell:irlev,:) = jb_c
      wrk_p_patch%cells%start_idx(min_rlcell:irlev-1,:) = wrk_p_patch%cells%start_idx(irlev,1)
      wrk_p_patch%cells%start_blk(min_rlcell:irlev-1,:) = wrk_p_patch%cells%start_blk(irlev,1)
    ELSE
      ! Gather all halo cells except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2_c value and will be placed at the
      ! index levels between min_rlcell_int-1 and min_rlcell
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here
      DO ilev = 1, 2*n_boundary_rows
        irlev = MAX(min_rlcell, min_rlcell_int - ilev)  ! index section into which the halo points are put
        DO j = 1, wrk_p_patch_g%n_patch_cells
          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          IF (flag2_c(j)==ilev .AND. .NOT. lcount_c(j)) THEN
            n = n + 1
            wrk_p_patch%cells%glb_index(n) = j
            wrk_p_patch%cells%loc_index(j) = n
            jb_c = blk_no(n) ! block index
            jl_c = idx_no(n) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (wrk_p_patch%cells%start_idx(irlev,1)==-9999) THEN
              wrk_p_patch%cells%start_idx(irlev,:) = jl_c
              wrk_p_patch%cells%start_blk(irlev,:) = jb_c
            ENDIF
            wrk_p_patch%cells%end_idx(irlev,:) = jl_c
            wrk_p_patch%cells%end_blk(irlev,:) = jb_c
          ENDIF
        ENDDO
      ENDDO
      ! Fill start and end indices for remaining index sections
      IF (jg == 0) THEN
        ilev1 = min_rlcell_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlcell, min_rlcell_int -2*n_boundary_rows)
        ilev_st = 2*n_boundary_rows+1
      ENDIF
      IF (wrk_p_patch%cell_type==6) THEN ! for hexagons, there are no even-order halo cells
        DO ilev = 2,2*max_hw,2
          irlev = MAX(min_rlcell, min_rlcell_int - ilev)  ! index section into which the halo points are put
          wrk_p_patch%cells%start_idx(irlev,:) = wrk_p_patch%cells%end_idx(irlev+1,1) + 1
          wrk_p_patch%cells%start_blk(irlev,:) = wrk_p_patch%cells%end_blk(irlev+1,1)
          wrk_p_patch%cells%end_idx(irlev,:)   = wrk_p_patch%cells%end_idx(irlev+1,1)
          wrk_p_patch%cells%end_blk(irlev,:)   = wrk_p_patch%cells%end_blk(irlev+1,1)
        ENDDO
      ENDIF
      DO ilev = ilev_st,2*max_hw
        irlev = MAX(min_rlcell, min_rlcell_int - ilev)  ! index section into which the halo points are put
        wrk_p_patch%cells%start_idx(irlev,:) = wrk_p_patch%cells%end_idx(ilev1,1) + 1
        wrk_p_patch%cells%start_blk(irlev,:) = wrk_p_patch%cells%end_blk(ilev1,1)
        wrk_p_patch%cells%end_idx(irlev,:)   = wrk_p_patch%cells%end_idx(ilev1,1)
        wrk_p_patch%cells%end_blk(irlev,:)   = wrk_p_patch%cells%end_blk(ilev1,1)
      ENDDO
    ENDIF

    ! If a PE owns only nest boundary points, it may happen that one or more index
    ! sections of halo cells are empty. These are filled here
    IF (wrk_p_patch%cells%start_blk(0,1)     > wrk_p_patch%cells%end_blk(min_rlcell_int,i_nchdom) &
      .OR. wrk_p_patch%cells%start_blk(0,1) == wrk_p_patch%cells%end_blk(min_rlcell_int,i_nchdom) &
      .AND. wrk_p_patch%cells%start_idx(0,1) > wrk_p_patch%cells%end_idx(min_rlcell_int,i_nchdom) &
      ) THEN
      DO i = min_rlcell_int-1, min_rlcell, -1
        IF (wrk_p_patch%cells%start_idx(i,1) == -9999 .OR. &
            wrk_p_patch%cells%start_blk(i,1) == -9999 .OR. &
            wrk_p_patch%cells%end_idx(i,1)   == -9999 .OR. &
            wrk_p_patch%cells%end_blk(i,1)   == -9999 ) THEN
          wrk_p_patch%cells%start_idx(i,:) = wrk_p_patch%cells%start_idx(i+1,1) 
          wrk_p_patch%cells%start_blk(i,:) = wrk_p_patch%cells%start_blk(i+1,1) 
          wrk_p_patch%cells%end_idx(i,:)   = wrk_p_patch%cells%end_idx(i+1,1) 
          wrk_p_patch%cells%end_blk(i,:)   = wrk_p_patch%cells%end_blk(i+1,1) 
        ENDIF
      ENDDO
    ENDIF

    ! Finally, fill start indices of halo rows with meaningful values for empty patches
    ! (which occur when processor splitting is applied)
    IF (wrk_p_patch%nblks_c <= 0 .OR. wrk_p_patch%npromz_c <= 0) THEN
      DO i=min_rlcell,min_rlcell_int-1
        wrk_p_patch%cells%start_idx(i,:) = wrk_p_patch%cells%start_idx(min_rlcell_int,1)
        wrk_p_patch%cells%start_blk(i,:) = wrk_p_patch%cells%start_blk(min_rlcell_int,1)
        wrk_p_patch%cells%end_idx(i,:)   = wrk_p_patch%cells%end_idx(min_rlcell_int,1)
        wrk_p_patch%cells%end_blk(i,:)   = wrk_p_patch%cells%end_blk(min_rlcell_int,1)
      ENDDO
    ENDIF

    !---------------------------------------------------------------------------------------

    n = 0
    IF (.NOT. l_compute_grid) THEN
      DO j = 1, wrk_p_patch_g%n_patch_edges
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF (flag2_e(j)==0) THEN
          n = n + 1
          wrk_p_patch%edges%glb_index(n) = j
          wrk_p_patch%edges%loc_index(j) = n
          lcount_e(j) = .TRUE.
        ELSE
          wrk_p_patch%edges%loc_index(j) = -(n+1)
        ENDIF
      ENDDO
    ELSE
      DO j = 1, wrk_p_patch_g%n_patch_edges
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        irl0 = wrk_p_patch_g%edges%refin_ctrl(jl,jb)
        IF (flag2_e(j)==0 .OR. ((irl0==1 .OR. irl0==2) .AND. flag2_e(j)>=1) .OR. &
          & (irl0==3 .OR. irl0==4) .AND. (flag2_e(j)==1 .OR. flag2_e(j)==2) .OR. &
          & (irl0==5 .OR. irl0==6 .OR. irl0<0) .AND. flag2_e(j)==1 ) THEN
          n = n + 1
          wrk_p_patch%edges%glb_index(n) = j
          wrk_p_patch%edges%loc_index(j) = n
          lcount_e(j) = .TRUE.
        ELSE
          wrk_p_patch%edges%loc_index(j) = -(n+1)
        ENDIF
      ENDDO
    ENDIF

    ! Set start_idx/blk ... end_idx/blk for edges.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the loc_index array.

    DO j=1,wrk_p_patch%max_childdom
      DO i=min_rledge_int,max_rledge
        CALL get_local_index(wrk_p_patch%edges%loc_index, &
          & wrk_p_patch_g%edges%start_idx(i,j),           &
          & wrk_p_patch_g%edges%start_blk(i,j),           &
          & wrk_p_patch%edges%start_idx(i,j),             &
          & wrk_p_patch%edges%start_blk(i,j),             &
          & +1 )
        CALL get_local_index(wrk_p_patch%edges%loc_index, &
          & wrk_p_patch_g%edges%end_idx(i,j),             &
          & wrk_p_patch_g%edges%end_blk(i,j),             &
          & wrk_p_patch%edges%end_idx(i,j),               &
          & wrk_p_patch%edges%end_blk(i,j),               &
          & -1 )
      ENDDO
      ! Preset remaining index sections with dummy values
      wrk_p_patch%edges%start_idx(min_rledge:min_rledge_int-1,j) = -9999
      wrk_p_patch%edges%start_blk(min_rledge:min_rledge_int-1,j) = -9999
      wrk_p_patch%edges%end_idx(min_rledge:min_rledge_int-1,j)   = -9999
      wrk_p_patch%edges%end_blk(min_rledge:min_rledge_int-1,j)   = -9999
    ENDDO

    IF(.NOT.l_compute_grid) THEN
      ! Gather all halo edges at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rledge_int-1
      irlev = min_rledge_int-1
      DO j = 1, wrk_p_patch_g%n_patch_edges
        IF(flag2_e(j)>1) THEN
          n = n + 1
          wrk_p_patch%edges%glb_index(n) = j
          wrk_p_patch%edges%loc_index(j) = n
          jb_e = blk_no(n) ! block index
          jl_e = idx_no(n) ! line index
          ! This ensures that just the first point found at this irlev is saved as start point
          IF (wrk_p_patch%edges%start_idx(irlev,1)==-9999) THEN
            wrk_p_patch%edges%start_idx(irlev,:) = jl_e
            wrk_p_patch%edges%start_blk(irlev,:) = jb_e
          ENDIF
        ENDIF
      ENDDO
      ! Set end index when the last point is processed
      jb_e = blk_no(n) ! block index
      jl_e = idx_no(n) ! line index
      wrk_p_patch%edges%end_idx(min_rledge:irlev,:) = jl_e
      wrk_p_patch%edges%end_blk(min_rledge:irlev,:) = jb_e
      wrk_p_patch%edges%start_idx(min_rledge:irlev-1,:) = wrk_p_patch%edges%start_idx(irlev,1)
      wrk_p_patch%edges%start_blk(min_rledge:irlev-1,:) = wrk_p_patch%edges%start_blk(irlev,1)
    ELSE
      ! Gather all halo edges except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2_e value and will be placed at the
      ! index levels between min_rledge_int-1 and min_rledge
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here
      DO ilev = 1, 2*n_boundary_rows+1
        irlev = MAX(min_rledge, min_rledge_int - ilev)  ! index section into which the halo points are put
        DO j = 1, wrk_p_patch_g%n_patch_edges
          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          IF (flag2_e(j)==ilev .AND. .NOT. lcount_e(j)) THEN
            n = n + 1
            wrk_p_patch%edges%glb_index(n) = j
            wrk_p_patch%edges%loc_index(j) = n
            jb_e = blk_no(n) ! block index
            jl_e = idx_no(n) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (wrk_p_patch%edges%start_idx(irlev,1)==-9999) THEN
              wrk_p_patch%edges%start_idx(irlev,:) = jl_e
              wrk_p_patch%edges%start_blk(irlev,:) = jb_e
            ENDIF
            wrk_p_patch%edges%end_idx(irlev,:) = jl_e
            wrk_p_patch%edges%end_blk(irlev,:) = jb_e
          ENDIF
        ENDDO
        ! Just in case that no grid point is found (may happen for ilev=1)
        IF (wrk_p_patch%edges%start_idx(irlev,1)==-9999) THEN
          wrk_p_patch%edges%start_idx(irlev,:) = wrk_p_patch%edges%end_idx(irlev+1,i_nchdom)+1
          wrk_p_patch%edges%start_blk(irlev,:) = wrk_p_patch%edges%end_blk(irlev+1,i_nchdom)
          wrk_p_patch%edges%end_idx(irlev,:)   = wrk_p_patch%edges%end_idx(irlev+1,i_nchdom)
          wrk_p_patch%edges%end_blk(irlev,:)   = wrk_p_patch%edges%end_blk(irlev+1,i_nchdom)
        ENDIF
      ENDDO
      ! Fill start and end indices for remaining index sections
      IF (jg == 0) THEN
        ilev1 = min_rledge_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rledge, min_rledge_int -(2*n_boundary_rows+1))
        ilev_st = 2*n_boundary_rows+2
      ENDIF
      DO ilev = ilev_st,2*max_hw+1
        irlev = MAX(min_rledge, min_rledge_int - ilev)  ! index section into which the halo points are put
        wrk_p_patch%edges%start_idx(irlev,:) = wrk_p_patch%edges%end_idx(ilev1,1) + 1
        wrk_p_patch%edges%start_blk(irlev,:) = wrk_p_patch%edges%end_blk(ilev1,1)
        wrk_p_patch%edges%end_idx(irlev,:)   = wrk_p_patch%edges%end_idx(ilev1,1)
        wrk_p_patch%edges%end_blk(irlev,:)   = wrk_p_patch%edges%end_blk(ilev1,1)
      ENDDO
    ENDIF

    ! If a PE owns only nest boundary points, it may happen that one or more index
    ! sections of halo edges are empty. These are filled here
    IF (wrk_p_patch%edges%start_blk(0,1)     > wrk_p_patch%edges%end_blk(min_rledge_int,i_nchdom) &
      .OR. wrk_p_patch%edges%start_blk(0,1) == wrk_p_patch%edges%end_blk(min_rledge_int,i_nchdom) &
      .AND. wrk_p_patch%edges%start_idx(0,1) > wrk_p_patch%edges%end_idx(min_rledge_int,i_nchdom) &
      ) THEN
      DO i = min_rledge_int-1, min_rledge, -1
        IF (wrk_p_patch%edges%start_idx(i,1) == -9999 .OR. &
            wrk_p_patch%edges%start_blk(i,1) == -9999 .OR. &
            wrk_p_patch%edges%end_idx(i,1)   == -9999 .OR. &
            wrk_p_patch%edges%end_blk(i,1)   == -9999 ) THEN
          wrk_p_patch%edges%start_idx(i,:) = wrk_p_patch%edges%start_idx(i+1,1) 
          wrk_p_patch%edges%start_blk(i,:) = wrk_p_patch%edges%start_blk(i+1,1) 
          wrk_p_patch%edges%end_idx(i,:)   = wrk_p_patch%edges%end_idx(i+1,1) 
          wrk_p_patch%edges%end_blk(i,:)   = wrk_p_patch%edges%end_blk(i+1,1) 
        ENDIF
      ENDDO
    ENDIF

    ! Finally, fill start indices of halo rows with meaningful values for empty patches
    ! (which occur when processor splitting is applied)
    IF (wrk_p_patch%nblks_e <= 0 .OR. wrk_p_patch%npromz_e <= 0) THEN
      DO i=min_rledge,min_rledge_int-1
        wrk_p_patch%edges%start_idx(i,:) = wrk_p_patch%edges%start_idx(min_rledge_int,1)
        wrk_p_patch%edges%start_blk(i,:) = wrk_p_patch%edges%start_blk(min_rledge_int,1)
        wrk_p_patch%edges%end_idx(i,:)   = wrk_p_patch%edges%end_idx(min_rledge_int,1)
        wrk_p_patch%edges%end_blk(i,:)   = wrk_p_patch%edges%end_blk(min_rledge_int,1)
      ENDDO
    ENDIF

    !---------------------------------------------------------------------------------------

    n = 0
    IF (.NOT. l_compute_grid) THEN
      DO j = 1, wrk_p_patch_g%n_patch_verts
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF (flag2_v(j)==0 ) THEN
          n = n + 1
          wrk_p_patch%verts%glb_index(n) = j
          wrk_p_patch%verts%loc_index(j) = n
          lcount_v(j) = .TRUE.
        ELSE
          wrk_p_patch%verts%loc_index(j) = -(n+1)
        ENDIF
      ENDDO
    ELSE
      DO j = 1, wrk_p_patch_g%n_patch_verts
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        irl0 = wrk_p_patch_g%verts%refin_ctrl(jl,jb)
        IF (flag2_v(j)==0 .OR. irl0==1 .AND. flag2_v(j)>0 .OR. &
          & (irl0==2 .OR. irl0<0) .AND. flag2_v(j)==1 ) THEN
          n = n + 1
          wrk_p_patch%verts%glb_index(n) = j
          wrk_p_patch%verts%loc_index(j) = n
          lcount_v(j) = .TRUE.
        ELSE
          wrk_p_patch%verts%loc_index(j) = -(n+1)
        ENDIF
      ENDDO
    ENDIF


    ! Set start_idx/blk ... end_idx/blk for verts.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the loc_index array.

    DO j=1,wrk_p_patch%max_childdom
      DO i=min_rlvert_int,max_rlvert
        CALL get_local_index(wrk_p_patch%verts%loc_index, &
          & wrk_p_patch_g%verts%start_idx(i,j),           &
          & wrk_p_patch_g%verts%start_blk(i,j),           &
          & wrk_p_patch%verts%start_idx(i,j),             &
          & wrk_p_patch%verts%start_blk(i,j),             &
          & +1 )
        CALL get_local_index(wrk_p_patch%verts%loc_index, &
          & wrk_p_patch_g%verts%end_idx(i,j),             &
          & wrk_p_patch_g%verts%end_blk(i,j),             &
          & wrk_p_patch%verts%end_idx(i,j),               &
          & wrk_p_patch%verts%end_blk(i,j),               &
          & -1 )
      ENDDO
      ! Preset remaining index sections with dummy values
      wrk_p_patch%verts%start_idx(min_rlvert:min_rlvert_int-1,j) = -9999
      wrk_p_patch%verts%start_blk(min_rlvert:min_rlvert_int-1,j) = -9999
      wrk_p_patch%verts%end_idx(min_rlvert:min_rlvert_int-1,j)   = -9999
      wrk_p_patch%verts%end_blk(min_rlvert:min_rlvert_int-1,j)   = -9999
    ENDDO

    IF(.NOT.l_compute_grid) THEN
      ! Gather all halo vertices at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlvert_int-1
      irlev = min_rlvert_int-1
      DO j = 1, wrk_p_patch_g%n_patch_verts
        IF(flag2_v(j)>1) THEN
          n = n + 1
          wrk_p_patch%verts%glb_index(n) = j
          wrk_p_patch%verts%loc_index(j) = n
          jb_v = blk_no(n) ! block index
          jl_v = idx_no(n) ! line index
          ! This ensures that just the first point found at this irlev is saved as start point
          IF (wrk_p_patch%verts%start_idx(irlev,1)==-9999) THEN
            wrk_p_patch%verts%start_idx(irlev,:) = jl_v
            wrk_p_patch%verts%start_blk(irlev,:) = jb_v
          ENDIF
        ENDIF
      ENDDO
      ! Set end index when the last point is processed
      jb_v = blk_no(n) ! block index
      jl_v = idx_no(n) ! line index
      wrk_p_patch%verts%end_idx(min_rlvert:irlev,:) = jl_v
      wrk_p_patch%verts%end_blk(min_rlvert:irlev,:) = jb_v
      wrk_p_patch%verts%start_idx(min_rlvert:irlev-1,:) = wrk_p_patch%verts%start_idx(irlev,1)
      wrk_p_patch%verts%start_blk(min_rlvert:irlev-1,:) = wrk_p_patch%verts%start_blk(irlev,1)
    ELSE
      ! Gather all halo vertices except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2_v value and will be placed at the
      ! index levels between min_rlvert_int-1 and min_rlvert
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here
      DO ilev = 1, n_boundary_rows+1
        irlev = MAX(min_rlvert, min_rlvert_int - ilev)  ! index section into which the halo points are put
        DO j = 1, wrk_p_patch_g%n_patch_verts
          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          IF (flag2_v(j)==ilev .AND. .NOT. lcount_v(j)) THEN
            n = n + 1
            wrk_p_patch%verts%glb_index(n) = j
            wrk_p_patch%verts%loc_index(j) = n
            jb_v = blk_no(n) ! block index
            jl_v = idx_no(n) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (wrk_p_patch%verts%start_idx(irlev,1)==-9999) THEN
              wrk_p_patch%verts%start_idx(irlev,:) = jl_v
              wrk_p_patch%verts%start_blk(irlev,:) = jb_v
            ENDIF
            wrk_p_patch%verts%end_idx(irlev,:) = jl_v
            wrk_p_patch%verts%end_blk(irlev,:) = jb_v
          ENDIF
        ENDDO
        ! Just in case that no grid point is found (may happen for ilev=1)
        IF (wrk_p_patch%verts%start_idx(irlev,1)==-9999) THEN
          wrk_p_patch%verts%start_idx(irlev,:) = wrk_p_patch%verts%end_idx(irlev+1,i_nchdom)+1
          wrk_p_patch%verts%start_blk(irlev,:) = wrk_p_patch%verts%end_blk(irlev+1,i_nchdom)
          wrk_p_patch%verts%end_idx(irlev,:)   = wrk_p_patch%verts%end_idx(irlev+1,i_nchdom)
          wrk_p_patch%verts%end_blk(irlev,:)   = wrk_p_patch%verts%end_blk(irlev+1,i_nchdom)
        ENDIF
      ENDDO
      ! Fill start and end indices for remaining index sections
      IF (jg == 0) THEN
        ilev1 = min_rlvert_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlvert, min_rlvert_int -(n_boundary_rows+1))
        ilev_st = n_boundary_rows+2
      ENDIF
      DO ilev = ilev_st,max_hw+1
        irlev = MAX(min_rlvert, min_rlvert_int - ilev)  ! index section into which the halo points are put
        wrk_p_patch%verts%start_idx(irlev,:) = wrk_p_patch%verts%end_idx(ilev1,1) + 1
        wrk_p_patch%verts%start_blk(irlev,:) = wrk_p_patch%verts%end_blk(ilev1,1)
        wrk_p_patch%verts%end_idx(irlev,:)   = wrk_p_patch%verts%end_idx(ilev1,1)
        wrk_p_patch%verts%end_blk(irlev,:)   = wrk_p_patch%verts%end_blk(ilev1,1)
      ENDDO
    ENDIF

    ! If a PE owns only nest boundary points, it may happen that one or more index
    ! sections of halo vertices are empty. These are filled here
    IF (wrk_p_patch%verts%start_blk(0,1)     > wrk_p_patch%verts%end_blk(min_rlvert_int,i_nchdom) &
      .OR. wrk_p_patch%verts%start_blk(0,1) == wrk_p_patch%verts%end_blk(min_rlvert_int,i_nchdom) &
      .AND. wrk_p_patch%verts%start_idx(0,1) > wrk_p_patch%verts%end_idx(min_rlvert_int,i_nchdom) &
      ) THEN
      DO i = min_rlvert_int-1, min_rlvert, -1
        IF (wrk_p_patch%verts%start_idx(i,1) == -9999 .OR. &
            wrk_p_patch%verts%start_blk(i,1) == -9999 .OR. &
            wrk_p_patch%verts%end_idx(i,1)   == -9999 .OR. &
            wrk_p_patch%verts%end_blk(i,1)   == -9999 ) THEN
          wrk_p_patch%verts%start_idx(i,:) = wrk_p_patch%verts%start_idx(i+1,1) 
          wrk_p_patch%verts%start_blk(i,:) = wrk_p_patch%verts%start_blk(i+1,1) 
          wrk_p_patch%verts%end_idx(i,:)   = wrk_p_patch%verts%end_idx(i+1,1) 
          wrk_p_patch%verts%end_blk(i,:)   = wrk_p_patch%verts%end_blk(i+1,1) 
        ENDIF
      ENDDO
    ENDIF

    ! Finally, fill start indices of halo rows with meaningful values for empty patches
    ! (which occur when processor splitting is applied)
    IF (wrk_p_patch%nblks_v <= 0 .OR. wrk_p_patch%npromz_v <= 0) THEN
      DO i=min_rlvert,min_rlvert_int-1
        wrk_p_patch%verts%start_idx(i,:) = wrk_p_patch%verts%start_idx(min_rlvert_int,1)
        wrk_p_patch%verts%start_blk(i,:) = wrk_p_patch%verts%start_blk(min_rlvert_int,1)
        wrk_p_patch%verts%end_idx(i,:)   = wrk_p_patch%verts%end_idx(min_rlvert_int,1)
        wrk_p_patch%verts%end_blk(i,:)   = wrk_p_patch%verts%end_blk(min_rlvert_int,1)
      ENDDO
    ENDIF

    ! Sanity checks: are there still elements of the index lists filled with dummy values?
    ! Note: the use of write(0,*) instead of CALL message is intended here
    !       because the debug output is wanted for all PEs where an error occurs
    IF (ANY(wrk_p_patch%cells%start_idx(:,:)<0) .OR. ANY(wrk_p_patch%cells%start_blk(:,:)<0) .OR.&
        ANY(wrk_p_patch%cells%end_idx(:,:)  <0) .OR. ANY(wrk_p_patch%cells%end_blk(:,:)  <0)) THEN
      DO j = 1, MAX(1,wrk_p_patch%n_childdom)
        DO i = min_rlcell, max_rlcell
          write(0,'(a,2i5,2i4,4i7)') 'cells',my_proc,wrk_p_patch%id,i,j,     &
            wrk_p_patch%cells%start_blk(i,j),wrk_p_patch%cells%start_idx(i,j), &
            wrk_p_patch%cells%end_blk(i,j),  wrk_p_patch%cells%end_idx(i,j)
        ENDDO
      ENDDO
      CALL finish('divide_patch','Error in cell start/end indices')
    ENDIF

    IF (ANY(wrk_p_patch%edges%start_idx(:,:)<0) .OR. ANY(wrk_p_patch%edges%start_blk(:,:)<0) .OR.&
        ANY(wrk_p_patch%edges%end_idx(:,:)  <0) .OR. ANY(wrk_p_patch%edges%end_blk(:,:)  <0)) THEN
      DO j = 1, MAX(1,wrk_p_patch%n_childdom)
        DO i = min_rledge, max_rledge
          write(0,'(a,2i5,2i4,4i7)') 'edges',my_proc,wrk_p_patch%id,i,j,     &
            wrk_p_patch%edges%start_blk(i,j),wrk_p_patch%edges%start_idx(i,j), &
            wrk_p_patch%edges%end_blk(i,j),  wrk_p_patch%edges%end_idx(i,j)
        ENDDO
      ENDDO
      CALL finish('divide_patch','Error in edge start/end indices')
    ENDIF

    IF (ANY(wrk_p_patch%verts%start_idx(:,:)<0) .OR. ANY(wrk_p_patch%verts%start_blk(:,:)<0) .OR.&
        ANY(wrk_p_patch%verts%end_idx(:,:)  <0) .OR. ANY(wrk_p_patch%verts%end_blk(:,:)  <0)) THEN
      DO j = 1, MAX(1,wrk_p_patch%n_childdom)
        DO i = min_rlvert, max_rlvert
          write(0,'(a,2i5,2i4,4i7)') 'verts',my_proc,wrk_p_patch%id,i,j,     &
            wrk_p_patch%verts%start_blk(i,j),wrk_p_patch%verts%start_idx(i,j), &
            wrk_p_patch%verts%end_blk(i,j),  wrk_p_patch%verts%end_idx(i,j)
        ENDDO
      ENDDO
      CALL finish('divide_patch','Error in vertex start/end indices')
    ENDIF

    !-----------------------------------------------------------------------------------------------
    ! Set arrays of divided patch
    !-----------------------------------------------------------------------------------------------

    ! decomp_domain is -1 for invalid locations (at the end of the last strip)

    wrk_p_patch%cells%decomp_domain = -1
    wrk_p_patch%edges%decomp_domain = -1
    wrk_p_patch%verts%decomp_domain = -1

    !---------------------------------------------------------------------------------------

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%cells%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%cells%glb_index(j)) ! Line  index in global patch

      ! parent_idx/parent_blk and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

      wrk_p_patch%cells%parent_idx(jl,jb)  = wrk_p_patch_g%cells%parent_idx(jl_g,jb_g)
      wrk_p_patch%cells%parent_blk(jl,jb)  = wrk_p_patch_g%cells%parent_blk(jl_g,jb_g)
      wrk_p_patch%cells%pc_idx(jl,jb)      = wrk_p_patch_g%cells%pc_idx(jl_g,jb_g)
      wrk_p_patch%cells%child_idx(jl,jb,:) = wrk_p_patch_g%cells%child_idx(jl_g,jb_g,:)
      wrk_p_patch%cells%child_blk(jl,jb,:) = wrk_p_patch_g%cells%child_blk(jl_g,jb_g,:)
      wrk_p_patch%cells%child_id (jl,jb)   = wrk_p_patch_g%cells%child_id(jl_g,jb_g)

      DO i=1,wrk_p_patch%cell_type

        CALL get_local_index(wrk_p_patch%cells%loc_index, &
          & wrk_p_patch_g%cells%neighbor_idx(jl_g,jb_g,i),&
          & wrk_p_patch_g%cells%neighbor_blk(jl_g,jb_g,i),&
          & wrk_p_patch%cells%neighbor_idx(jl,jb,i),      &
          & wrk_p_patch%cells%neighbor_blk(jl,jb,i))

        CALL get_local_index(wrk_p_patch%edges%loc_index, &
          & wrk_p_patch_g%cells%edge_idx(jl_g,jb_g,i),    &
          & wrk_p_patch_g%cells%edge_blk(jl_g,jb_g,i),    &
          & wrk_p_patch%cells%edge_idx(jl,jb,i),          &
          & wrk_p_patch%cells%edge_blk(jl,jb,i))

        CALL get_local_index(wrk_p_patch%verts%loc_index, &
          & wrk_p_patch_g%cells%vertex_idx(jl_g,jb_g,i),  &
          & wrk_p_patch_g%cells%vertex_blk(jl_g,jb_g,i),  &
          & wrk_p_patch%cells%vertex_idx(jl,jb,i),        &
          & wrk_p_patch%cells%vertex_blk(jl,jb,i))

      ENDDO

      ! Safety check only: edge_idx and vertex_idx must always be valid!
      if(any(wrk_p_patch%cells%edge_idx(jl,jb,:) <= 0) .or. &
       & any(wrk_p_patch%cells%edge_blk(jl,jb,:) <= 0) )    &
       & CALL finish('divide_patch','Illegal value for patch%cells%edge_idx')
      if(any(wrk_p_patch%cells%vertex_idx(jl,jb,:) <= 0) .or. &
       & any(wrk_p_patch%cells%vertex_blk(jl,jb,:) <= 0) )    &
       & CALL finish('divide_patch','Illegal value for patch%cells%vertex_idx')

      wrk_p_patch%cells%num_edges(jl,jb)          = wrk_p_patch_g%cells%num_edges(jl_g,jb_g)
      wrk_p_patch%cells%center(jl,jb)             = wrk_p_patch_g%cells%center(jl_g,jb_g)
      wrk_p_patch%cells%refin_ctrl(jl,jb)         = wrk_p_patch_g%cells%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%cells%child_id(jl,jb)           = wrk_p_patch_g%cells%child_id(jl_g,jb_g)
      wrk_p_patch%cells%decomp_domain(jl,jb)      = flag2_c(wrk_p_patch%cells%glb_index(j))

    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%edges%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%edges%glb_index(j)) ! Line  index in global patch

      ! parent_idx/parent_blk and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

      wrk_p_patch%edges%parent_idx(jl,jb)    = wrk_p_patch_g%edges%parent_idx(jl_g,jb_g)
      wrk_p_patch%edges%parent_blk(jl,jb)    = wrk_p_patch_g%edges%parent_blk(jl_g,jb_g)
      wrk_p_patch%edges%pc_idx(jl,jb)        = wrk_p_patch_g%edges%pc_idx(jl_g,jb_g)
      wrk_p_patch%edges%child_idx(jl,jb,:)   = wrk_p_patch_g%edges%child_idx(jl_g,jb_g,:)
      wrk_p_patch%edges%child_blk(jl,jb,:)   = wrk_p_patch_g%edges%child_blk(jl_g,jb_g,:)
      wrk_p_patch%edges%child_id (jl,jb)     = wrk_p_patch_g%edges%child_id(jl_g,jb_g)

      wrk_p_patch%edges%refin_ctrl(jl,jb)    = wrk_p_patch_g%edges%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%edges%decomp_domain(jl,jb) = flag2_e(wrk_p_patch%edges%glb_index(j))

    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%verts%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%verts%glb_index(j)) ! Line  index in global patch

      wrk_p_patch%verts%vertex(jl,jb)        = wrk_p_patch_g%verts%vertex(jl_g,jb_g)
      wrk_p_patch%verts%refin_ctrl(jl,jb)    = wrk_p_patch_g%verts%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%verts%decomp_domain(jl,jb) = flag2_v(wrk_p_patch%verts%glb_index(j))

    ENDDO
            
    DEALLOCATE(flag_c, flag_e, flag_v, flag2_c, flag2_e, flag2_v, lcount_c, lcount_e, lcount_v)

  END SUBROUTINE divide_patch
  !-------------------------------------------------------------------------------------------------
  !>
  !! Remaps negative index entries to valid ones
  !!

  SUBROUTINE remap_patch_indices(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch

    INTEGER :: i, j, jb, jl, ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, irl0, irl1, irl2, irl3

    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,p_patch%cell_type

        CALL remap_index(p_patch%cells%loc_index, &
          & p_patch%cells%neighbor_idx(jl,jb,i),      &
          & p_patch%cells%neighbor_blk(jl,jb,i))

        ! edge_idx and vertex_idx should not need a remap !!!
        ! This is only left here if there should change something in the decomposition
        CALL remap_index(p_patch%edges%loc_index, &
          & p_patch%cells%edge_idx(jl,jb,i),          &
          & p_patch%cells%edge_blk(jl,jb,i))

        CALL remap_index(p_patch%verts%loc_index, &
          & p_patch%cells%vertex_idx(jl,jb,i),        &
          & p_patch%cells%vertex_blk(jl,jb,i))

      ENDDO
    ENDDO

    ! ensure that cells%neighbor_idx lies in the correct grid row along the lateral boundary
    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (p_patch%cells%refin_ctrl(jl,jb) > 0 .AND.             &
          p_patch%cells%refin_ctrl(jl,jb) <= grf_bdywidth_c) THEN

        ilc1 = p_patch%cells%neighbor_idx(jl,jb,1)
        ibc1 = p_patch%cells%neighbor_blk(jl,jb,1)
        ilc2 = p_patch%cells%neighbor_idx(jl,jb,2)
        ibc2 = p_patch%cells%neighbor_blk(jl,jb,2)
        ilc3 = p_patch%cells%neighbor_idx(jl,jb,3)
        ibc3 = p_patch%cells%neighbor_blk(jl,jb,3)

        irl0 = p_patch%cells%refin_ctrl(jl,jb)

        IF (ilc1 > 0 .AND. ibc1 > 0) THEN
          irl1 = p_patch%cells%refin_ctrl(ilc1,ibc1)
          IF ( irl1 > 0 .AND. ABS(irl0 - irl1) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,1) = jl
            p_patch%cells%neighbor_blk(jl,jb,1) = jb
          ENDIF
        ENDIF

        IF (ilc2 > 0 .AND. ibc2 > 0) THEN
          irl2 = p_patch%cells%refin_ctrl(ilc2,ibc2)
          IF ( irl2 > 0 .AND. ABS(irl0 - irl2) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,2) = jl
            p_patch%cells%neighbor_blk(jl,jb,2) = jb
          ENDIF
        ENDIF

        IF (ilc3 > 0 .AND. ibc3 > 0) THEN
          irl3 = p_patch%cells%refin_ctrl(ilc3,ibc3)
          IF ( irl3 > 0 .AND. ABS(irl0 - irl3) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,3) = jl
            p_patch%cells%neighbor_blk(jl,jb,3) = jb
          ENDIF
        ENDIF

      ENDIF
    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,2
        CALL remap_index(p_patch%cells%loc_index, &
          & p_patch%edges%cell_idx(jl,jb,i),          &
          & p_patch%edges%cell_blk(jl,jb,i))
      ENDDO

      DO i=1,4
        CALL remap_index(p_patch%verts%loc_index, &
          & p_patch%edges%vertex_idx(jl,jb,i),        &
          & p_patch%edges%vertex_blk(jl,jb,i))
      ENDDO

      DO i=1,4
        CALL remap_index(p_patch%edges%loc_index, &
          & p_patch%edges%quad_idx(jl,jb,i),          &
          & p_patch%edges%quad_blk(jl,jb,i))
      ENDDO

      ! ensure that edges%cell_idx lies in the correct grid row along the lateral boundary
      IF (p_patch%edges%refin_ctrl(jl,jb) >= 2 .AND.            &
          p_patch%edges%refin_ctrl(jl,jb) <= grf_bdywidth_e) THEN

        ilc1 = p_patch%edges%cell_idx(jl,jb,1)
        ibc1 = p_patch%edges%cell_blk(jl,jb,1)
        ilc2 = p_patch%edges%cell_idx(jl,jb,2)
        ibc2 = p_patch%edges%cell_blk(jl,jb,2)
        irl1 = p_patch%cells%refin_ctrl(ilc1,ibc1)
        irl2 = p_patch%cells%refin_ctrl(ilc2,ibc2)

        IF (MOD(p_patch%edges%refin_ctrl(jl,jb),2)==0) THEN
          IF (irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            p_patch%edges%cell_idx(jl,jb,1) = p_patch%edges%cell_idx(jl,jb,2)
            p_patch%edges%cell_blk(jl,jb,1) = p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            p_patch%edges%cell_idx(jl,jb,2) = p_patch%edges%cell_idx(jl,jb,1)
            p_patch%edges%cell_blk(jl,jb,2) = p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ELSE
          IF (irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            p_patch%edges%cell_idx(jl,jb,1) = p_patch%edges%cell_idx(jl,jb,2)
            p_patch%edges%cell_blk(jl,jb,1) = p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            p_patch%edges%cell_idx(jl,jb,2) = p_patch%edges%cell_idx(jl,jb,1)
            p_patch%edges%cell_blk(jl,jb,2) = p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ENDIF

      ENDIF

    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,9-p_patch%cell_type

        CALL remap_index(p_patch%verts%loc_index, &
          & p_patch%verts%neighbor_idx(jl,jb,i),      &
          & p_patch%verts%neighbor_blk(jl,jb,i))

        CALL remap_index(p_patch%edges%loc_index, &
          & p_patch%verts%edge_idx(jl,jb,i),          &
          & p_patch%verts%edge_blk(jl,jb,i))

        CALL remap_index(p_patch%cells%loc_index, &
          & p_patch%verts%cell_idx(jl,jb,i),          &
          & p_patch%verts%cell_blk(jl,jb,i))
      ENDDO

    ENDDO

  END SUBROUTINE remap_patch_indices

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets parent_idx/blk in child and child_idx/blk in parent patches.

  SUBROUTINE set_parent_child_relations(p_pp, p_pc)

    TYPE(t_patch), INTENT(INOUT) :: p_pp   !> divided local parent patch
    TYPE(t_patch), INTENT(INOUT) :: p_pc   !> divided child patch

    INTEGER :: i, j, jl, jb, jc, jc_g, jp, jp_g

    ! Before this call, parent_idx/parent_blk and child_idx/child_blk still point to the global values.
    ! This is changed here.

    ! Attention:
    ! Only inner cells/edges get a valid child index,
    ! indexes for boundary cells/edges are not set.
    ! Therefore when these indexes are used, the code must assure that they are
    ! used only for inner cells/edges!
    ! The main reason for this is that - depending on the number of ghost rows -
    ! there are cells/edges in the parent boundary with missing childs (n_ghost_rows==1)

    ! Set child indices in parent ...

    ! ... cells

    DO j = 1, p_pp%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pp%cells%decomp_domain(jl,jb)>0) THEN
        p_pp%cells%child_idx(jl,jb,:) = 0
        p_pp%cells%child_blk(jl,jb,:) = 0
        CYCLE
      ENDIF

      DO i= 1, 4
        jc_g = idx_1d(p_pp%cells%child_idx(jl,jb,i),p_pp%cells%child_blk(jl,jb,i))
        IF(jc_g<1 .OR. jc_g>p_pc%n_patch_cells_g) &
          & CALL finish('set_parent_child_relations','Invalid cell child index in global parent')
        jc = p_pc%cells%loc_index(jc_g)
        IF(jc <= 0) &
          & CALL finish('set_parent_child_relations','cell child index outside child domain')
        p_pp%cells%child_blk(jl,jb,i) = blk_no(jc)
        p_pp%cells%child_idx(jl,jb,i) = idx_no(jc)
      ENDDO

    ENDDO

    ! ... edges

    DO j = 1, p_pp%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pp%edges%decomp_domain(jl,jb)>1) THEN
        p_pp%edges%child_idx(jl,jb,:) = 0
        p_pp%edges%child_blk(jl,jb,:) = 0
        CYCLE ! only inner edges get a valid parent index
      ENDIF

      DO i= 1, 4

        IF(i==4 .AND. p_pp%edges%refin_ctrl(jl,jb) == -1) THEN
          p_pp%edges%child_blk(jl,jb,i) = blk_no(0)
          p_pp%edges%child_idx(jl,jb,i) = idx_no(0)
          CYCLE
        ENDIF

        jc_g = idx_1d(p_pp%edges%child_idx(jl,jb,i),p_pp%edges%child_blk(jl,jb,i))
        IF(jc_g<1 .OR. jc_g>p_pc%n_patch_edges_g) &
          & CALL finish('set_parent_child_relations','Inv. edge child index in global parent')
        jc = p_pc%edges%loc_index(ABS(jc_g))
        IF(jc <= 0) &
          & CALL finish('set_parent_child_relations','edge child index outside child domain')
        p_pp%edges%child_blk(jl,jb,i) = blk_no(jc)
        p_pp%edges%child_idx(jl,jb,i) = SIGN(idx_no(jc),jc_g)
      ENDDO

      p_pp%edges%child_id(jl,jb) = p_pc%id

    ENDDO

    ! Set parent indices in child ...

    ! ... cells

    DO j = 1, p_pc%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jp_g = idx_1d(p_pc%cells%parent_idx(jl,jb),p_pc%cells%parent_blk(jl,jb))
      IF(jp_g<1 .OR. jp_g>p_pp%n_patch_cells_g) &
        & CALL finish('set_parent_child_relations','Inv. cell parent index in global child')

      jp = p_pp%cells%loc_index(jp_g)
      IF(jp <= 0) THEN
        p_pc%cells%parent_blk(jl,jb) = 0
        p_pc%cells%parent_idx(jl,jb) = 0
      ELSE
        p_pc%cells%parent_blk(jl,jb) = blk_no(jp)
        p_pc%cells%parent_idx(jl,jb) = idx_no(jp)
      ENDIF

    ENDDO

    ! ... edges

    DO j = 1, p_pc%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jp_g = idx_1d(p_pc%edges%parent_idx(jl,jb),p_pc%edges%parent_blk(jl,jb))
      IF(jp_g<1 .OR. jp_g>p_pp%n_patch_edges_g) &
        & CALL finish('set_parent_child_relations','Inv. edge parent index in global child')

      jp = p_pp%edges%loc_index(jp_g)
      IF(jp <= 0) THEN
        p_pc%edges%parent_blk(jl,jb) = 0
        p_pc%edges%parent_idx(jl,jb) = 0
      ELSE
        p_pc%edges%parent_blk(jl,jb) = blk_no(jp)
        p_pc%edges%parent_idx(jl,jb) = idx_no(jp)
      ENDIF

    ENDDO

    ! Although this is not really necessary, we set the child index in child
    ! and the parent index in parent to 0 since these have no significance
    ! in the parallel code (and must not be used as they are).

    p_pc%cells%child_idx  = 0
    p_pc%cells%child_blk  = 0
    p_pp%cells%parent_idx = 0
    p_pp%cells%parent_blk = 0

    p_pc%edges%child_idx  = 0
    p_pc%edges%child_blk  = 0
    p_pp%edges%parent_idx = 0
    p_pp%edges%parent_blk = 0

  END SUBROUTINE set_parent_child_relations


  !-------------------------------------------------------------------------------------------------
  !> Sets the communication patterns for boundary exchange

  SUBROUTINE set_comm_pat_bound_exch(p)

    TYPE(t_patch), INTENT(INOUT):: p

    INTEGER, ALLOCATABLE :: owner_c(:), owner_e(:), owner_v(:)
    INTEGER :: jc

    ALLOCATE(owner_c(p%n_patch_cells))
    ALLOCATE(owner_e(p%n_patch_edges))
    ALLOCATE(owner_v(p%n_patch_verts))

    ! Set the owner arrays for cells/edges/verts which have to be transferred.
    ! The following setting always transfers edges/verts if they are owned
    ! (according to the MIN/MAX PE setting above) by another PE, which implies
    ! that boundary edges/verts (with flag2_e/flag2_v = 1) can be excluded from
    ! prognostic computations if they are followed by a synchronization call.

    owner_c(:) = p%cells%owner_g(p%cells%glb_index(:))
    WHERE(owner_c(:) == p_pe_work) owner_c(:) = -1

    owner_e(:) = p%edges%owner_g(p%edges%glb_index(:))
    WHERE(owner_e(:) == p_pe_work) owner_e(:) = -1

    owner_v(:) = p%verts%owner_g(p%verts%glb_index(:))
    WHERE(owner_v(:) == p_pe_work) owner_v(:) = -1

    ! Set communication patterns for boundary exchange
    CALL setup_comm_pattern(p%n_patch_cells, owner_c, p%cells%glb_index, &
      & p%cells%loc_index, p%comm_pat_c)

    CALL setup_comm_pattern(p%n_patch_edges, owner_e, p%edges%glb_index, &
      & p%edges%loc_index, p%comm_pat_e)

    CALL setup_comm_pattern(p%n_patch_verts, owner_v, p%verts%glb_index, &
      & p%verts%loc_index, p%comm_pat_v)

    DEALLOCATE(owner_c, owner_e, owner_v)

    ! Set reduced communication pattern containing only level-1 halo cells (immediate neighbors)
    jc = idx_1d(p%cells%end_idx(min_rlcell_int-1,1), &
                p%cells%end_blk(min_rlcell_int-1,1))
    ALLOCATE(owner_c(jc))
    owner_c(1:jc) = p%cells%owner_g(p%cells%glb_index(1:jc))
    WHERE(owner_c(:) == p_pe_work) owner_c(:) = -1

    CALL setup_comm_pattern(jc, owner_c, p%cells%glb_index, &
      & p%cells%loc_index, p%comm_pat_c1)

    DEALLOCATE(owner_c)

  END SUBROUTINE set_comm_pat_bound_exch

  !-------------------------------------------------------------------------------------------------
  !> Sets the gather communication patterns of a patch

  SUBROUTINE set_comm_pat_gather(p)

    TYPE(t_patch), INTENT(INOUT):: p

    INTEGER, ALLOCATABLE :: tmp(:)
    INTEGER :: j

    ! For gathering the global fields on p_pe_work==0
    ALLOCATE(tmp(MAX(p%n_patch_cells_g, p%n_patch_edges_g, p%n_patch_verts_g)))

    DO j = 1, SIZE(tmp)
      tmp(j) = j ! Global/local index in global array, i.e. identity!
    ENDDO

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(p%n_patch_cells_g, p%cells%owner_g, tmp, &
        & p%cells%loc_index, p%comm_pat_gather_c)
    ELSE
      ! We don't want to receive any data, i.e. the number of cells is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, p%cells%owner_g, tmp, &
        & p%cells%loc_index, p%comm_pat_gather_c)
    ENDIF

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(p%n_patch_edges_g, p%edges%owner_g, tmp, &
        & p%edges%loc_index, p%comm_pat_gather_e)
    ELSE
      ! We don't want to receive any data, i.e. the number of edges is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, p%edges%owner_g, tmp, &
        & p%edges%loc_index, p%comm_pat_gather_e)
    ENDIF

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(p%n_patch_verts_g, p%verts%owner_g, tmp, &
        & p%verts%loc_index, p%comm_pat_gather_v)
    ELSE
      ! We don't want to receive any data, i.e. the number of edges is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, p%verts%owner_g, tmp, &
        & p%verts%loc_index, p%comm_pat_gather_v)
    ENDIF

    DEALLOCATE(tmp)

  END SUBROUTINE set_comm_pat_gather

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets up communication patterns between global and local parent patches

  SUBROUTINE set_glb_loc_comm(p_pglb, p_ploc, i_chidx)
    TYPE(t_patch), INTENT(IN)    :: p_pglb !> global parent
    TYPE(t_patch), INTENT(INOUT) :: p_ploc !> local parent
    INTEGER, INTENT(IN) :: i_chidx

    INTEGER, ALLOCATABLE :: owner(:)
    INTEGER :: j, js, je, icid, jb, jl

    ! Please note:
    ! For creating communication patterns for different amount of data to be transferred
    ! (e.g. only the boundary interpolation zone), create a new communication pattern
    ! by copying the code below and adjusting the limits in the calculation of js/je.

    !-----------------------------------------------------------------------------------------------

    ! child ID
    icid = p_pglb%child_id(i_chidx)

    ! Communication global -> local
    ! Only one pattern is set which can be used everywhere since it doesn't matter
    ! if more cells/edges than needed are set in the local parent.
    ! Please note that only the inner area of the local patch is set and not the boundary cells
    ! (which are needed only for gradient calculation in the moment).

    ! ... cells

    ALLOCATE(owner(p_ploc%n_patch_cells))

    js = idx_1d(p_ploc%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
      &         p_ploc%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    je = idx_1d(p_ploc%cells%end_idx(min_rlcell_int,i_chidx), &
      &         p_ploc%cells%end_blk(min_rlcell_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_pglb%cells%owner_g(p_ploc%cells%glb_index(j))
    ENDDO

    CALL setup_comm_pattern(p_ploc%n_patch_cells, owner, p_ploc%cells%glb_index,  &
      & p_pglb%cells%loc_index, p_ploc%comm_pat_glb_to_loc_c)

    DEALLOCATE(owner)

    ! ... edges

    ALLOCATE(owner(p_ploc%n_patch_edges))

    js = idx_1d(p_ploc%edges%start_idx(grf_bdyintp_start_e,i_chidx), &
      &         p_ploc%edges%start_blk(grf_bdyintp_start_e,i_chidx))
    je = idx_1d(p_ploc%edges%end_idx(min_rledge_int,i_chidx), &
      &         p_ploc%edges%end_blk(min_rledge_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_pglb%edges%owner_g(p_ploc%edges%glb_index(j))
    ENDDO

    CALL setup_comm_pattern(p_ploc%n_patch_edges, owner, p_ploc%edges%glb_index,  &
      & p_pglb%edges%loc_index, p_ploc%comm_pat_glb_to_loc_e)

    DEALLOCATE(owner)

    !-----------------------------------------------------------------------------------------------

    ! Communication local -> global
    ! Here it might get necessary to have different patterns for different start levels
    ! of the copy to the global parent since we may not overwrite arbitrary values there.
    ! Currently only one pattern for feedback is needed (starting at grf_fbk_start_c/e).

    ! ... cells

    ALLOCATE(owner(p_pglb%n_patch_cells))

    js = idx_1d(p_pglb%cells%start_idx(grf_fbk_start_c,i_chidx), &
      &         p_pglb%cells%start_blk(grf_fbk_start_c,i_chidx))
    je = idx_1d(p_pglb%cells%end_idx(min_rlcell_int,i_chidx), &
      &         p_pglb%cells%end_blk(min_rlcell_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_ploc%cells%owner_g(p_pglb%cells%glb_index(j))
    ENDDO

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      js = idx_1d(p_pglb%cells%start_idx(min_rlcell_int-1,1), &
        &         p_pglb%cells%start_blk(min_rlcell_int-1,1))
      je = idx_1d(p_pglb%cells%end_idx(min_rlcell,MAX(1,p_pglb%n_childdom)), &
        &         p_pglb%cells%end_blk(min_rlcell,MAX(1,p_pglb%n_childdom)))

      DO j = js, je

        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index
        IF (p_pglb%cells%child_id(jl,jb)   == icid            .AND. &
            p_pglb%cells%refin_ctrl(jl,jb) <= grf_fbk_start_c .AND. &
            p_pglb%cells%refin_ctrl(jl,jb) >= min_rlcell_int )   THEN

          owner(j) = p_ploc%cells%owner_g(p_pglb%cells%glb_index(j))
        ENDIF
      ENDDO

    ENDIF

    CALL setup_comm_pattern(p_pglb%n_patch_cells, owner, p_pglb%cells%glb_index, &
      & p_ploc%cells%loc_index, p_ploc%comm_pat_loc_to_glb_c_fbk)

    DEALLOCATE(owner)

    ! ... edges

    ALLOCATE(owner(p_pglb%n_patch_edges))

    js = idx_1d(p_pglb%edges%start_idx(grf_fbk_start_e,i_chidx), &
      &         p_pglb%edges%start_blk(grf_fbk_start_e,i_chidx))
    je = idx_1d(p_pglb%edges%end_idx(min_rledge_int,i_chidx), &
      &         p_pglb%edges%end_blk(min_rledge_int,i_chidx))

    owner(:) = -1 ! By default don't include into comm pattern
    DO j = js, je
      owner(j) = p_ploc%edges%owner_g(p_pglb%edges%glb_index(j))
    ENDDO

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      js = idx_1d(p_pglb%edges%start_idx(min_rledge_int-1,1), &
        &         p_pglb%edges%start_blk(min_rledge_int-1,1))
      je = idx_1d(p_pglb%edges%end_idx(min_rledge,MAX(1,p_pglb%n_childdom)), &
        &         p_pglb%edges%end_blk(min_rledge,MAX(1,p_pglb%n_childdom)))

      DO j = js, je

        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index
        IF (p_pglb%edges%child_id(jl,jb)   == icid            .AND. &
            p_pglb%edges%refin_ctrl(jl,jb) <= grf_fbk_start_e .AND. &
            p_pglb%edges%refin_ctrl(jl,jb) >= min_rledge_int )   THEN

          owner(j) = p_ploc%edges%owner_g(p_pglb%edges%glb_index(j))
        ENDIF
      ENDDO

    ENDIF

    CALL setup_comm_pattern(p_pglb%n_patch_edges, owner, p_pglb%edges%glb_index, &
      & p_ploc%edges%loc_index, p_ploc%comm_pat_loc_to_glb_e_fbk)

    DEALLOCATE(owner)

   END SUBROUTINE set_glb_loc_comm

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for interpolation by direct copying.
  !!
  !! This routine sets up a communication pattern for interpolation by direct copying.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_cpy_interpolation(p_patch, p_parent_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch, p_parent_patch

    INTEGER :: j, jc, jb, jp, p_index_s, p_index_e, i_chidx
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(my_process_is_mpi_seq()) &
      & CALL finish('setup_comm_cpy_interpolation','must not be called in a single CPU run')

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Get start and end index of the GLOBAL parent cells as used in the interpolation

    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(grf_bdyintp_end_c,i_chidx), &
                       p_parent_patch%cells%end_blk(grf_bdyintp_end_c,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    glb_index(:) = -1
    owner(:)     = -1

    DO j = 1, p_patch%n_patch_cells
      jc = idx_no(j)
      jb = blk_no(j)
      jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
      IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
      glb_index(j) = jp
      owner(j) = p_parent_patch%cells%owner_g(jp)
    ENDDO

    ! Set up communication pattern

    CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
      & p_parent_patch%cells%loc_index, &
      & p_patch%comm_pat_interpolation_c)

    DEALLOCATE(owner, glb_index)

  END SUBROUTINE setup_comm_cpy_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for grf interpolation.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_grf_interpolation(p_patch, p_parent_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch, p_parent_patch

    INTEGER :: j, n, jc, je, jb, jp, p_index_s, p_index_e, i_chidx
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(my_process_is_mpi_seq()) &
      & CALL finish('setup_comm_grf_interpolation','must not be called in a single CPU run')

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Start and end index of the GLOBAL parent cells as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_bdyintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_bdyintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(grf_bdyintp_end_c,i_chidx), &
                       p_parent_patch%cells%end_blk(grf_bdyintp_end_c,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_cells
        jc = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%cells%pc_idx(jc,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%cells%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
        & p_parent_patch%cells%loc_index, &
        & p_patch%comm_pat_interpol_scal_grf(n))

    ENDDO

    DEALLOCATE(owner, glb_index)

    !--------------------------------------------------------------------
    ! Edges

    ! Start and end index of the GLOBAL parent edges as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%edges%start_idx(grf_bdyintp_start_e,i_chidx), &
                       p_parent_patch%edges%start_blk(grf_bdyintp_start_e,i_chidx))
    p_index_e = idx_1d(p_parent_patch%edges%end_idx(grf_bdyintp_end_e,i_chidx), &
                       p_parent_patch%edges%end_blk(grf_bdyintp_end_e,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%edges%glb_index(p_index_s)
      p_index_e = p_parent_patch%edges%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(p_patch%n_patch_edges))
    ALLOCATE(owner(p_patch%n_patch_edges))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_edges
        je = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%edges%parent_idx(je,jb),p_patch%edges%parent_blk(je,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%edges%pc_idx(je,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%edges%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_edges, owner, glb_index,  &
        & p_parent_patch%edges%loc_index, &
        & p_patch%comm_pat_interpol_vec_grf(n))

    ENDDO

    DEALLOCATE(owner, glb_index)

  END SUBROUTINE setup_comm_grf_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for ubc interpolation.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_ubc_interpolation(p_patch, p_parent_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch, p_parent_patch

    INTEGER :: j, n, jc, je, jb, jp, p_index_s, p_index_e, i_chidx
    INTEGER, ALLOCATABLE :: owner(:), glb_index(:)

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(my_process_is_mpi_seq()) &
      & CALL finish('setup_comm_ubc_interpolation','must not be called in a single CPU run')

    i_chidx = p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Start and end index of the GLOBAL parent cells as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%cells%start_idx(grf_nudgintp_start_c,i_chidx), &
                       p_parent_patch%cells%start_blk(grf_nudgintp_start_c,i_chidx))
    p_index_e = idx_1d(p_parent_patch%cells%end_idx(min_rlcell_int,i_chidx), &
                       p_parent_patch%cells%end_blk(min_rlcell_int,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%cells%glb_index(p_index_s)
      p_index_e = p_parent_patch%cells%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(p_patch%n_patch_cells))
    ALLOCATE(owner(p_patch%n_patch_cells))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_cells
        jc = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%cells%parent_idx(jc,jb),p_patch%cells%parent_blk(jc,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%cells%pc_idx(jc,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%cells%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_cells, owner, glb_index,  &
        & p_parent_patch%cells%loc_index, &
        & p_patch%comm_pat_interpol_scal_ubc(n))

    ENDDO

    DEALLOCATE(owner, glb_index)

    !--------------------------------------------------------------------
    ! Edges

    ! Start and end index of the GLOBAL parent edges as used in the interpolation
    p_index_s = idx_1d(p_parent_patch%edges%start_idx(grf_nudgintp_start_e,i_chidx), &
                       p_parent_patch%edges%start_blk(grf_nudgintp_start_e,i_chidx))
    p_index_e = idx_1d(p_parent_patch%edges%end_idx(min_rledge_int,i_chidx), &
                       p_parent_patch%edges%end_blk(min_rledge_int,i_chidx))
    IF(p_index_s <= p_index_e) THEN
      p_index_s = p_parent_patch%edges%glb_index(p_index_s)
      p_index_e = p_parent_patch%edges%glb_index(p_index_e)
    ELSE
      p_index_s =  HUGE(0)
      p_index_e = -HUGE(0)
    ENDIF
    p_index_s = p_min(p_index_s, p_comm_work)
    p_index_e = p_max(p_index_e, p_comm_work)

    ! For our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(p_patch%n_patch_edges))
    ALLOCATE(owner(p_patch%n_patch_edges))

    DO n = 1, 4

      glb_index(:) = -1
      owner(:)     = -1

      DO j = 1,p_patch%n_patch_edges
        je = idx_no(j)
        jb = blk_no(j)
        jp = idx_1d(p_patch%edges%parent_idx(je,jb),p_patch%edges%parent_blk(je,jb))
        IF(jp<p_index_s .OR. jp>p_index_e) CYCLE
        IF(p_patch%edges%pc_idx(je,jb) /= n) CYCLE
        glb_index(j) = jp
        owner(j) = p_parent_patch%edges%owner_g(jp)
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(p_patch%n_patch_edges, owner, glb_index,  &
        & p_parent_patch%edges%loc_index, &
        & p_patch%comm_pat_interpol_vec_ubc(n))

    ENDDO

    DEALLOCATE(owner, glb_index)

  END SUBROUTINE setup_comm_ubc_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a the physical patches
  !! It works on the divided patch state, so it can be used after a restore also
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2011

  SUBROUTINE setup_phys_patches

    INTEGER :: jp, jg, n, i, j, jb, jl
    INTEGER, ALLOCATABLE :: glb_phys_id_c(:), glb_phys_id_e(:), glb_phys_id_v(:)
    INTEGER, ALLOCATABLE :: owner(:), glbidx(:)
    CHARACTER(LEN=*), PARAMETER :: routine = 'setup_phys_patches'

    p_phys_patch(:)%logical_id = -1

    DO jg = 1, n_dom

      ! Get global arrays for phys_id

      ! Allocate and set to 0
      ALLOCATE(glb_phys_id_c(p_patch(jg)%n_patch_cells_g))
      ALLOCATE(glb_phys_id_e(p_patch(jg)%n_patch_edges_g))
      ALLOCATE(glb_phys_id_v(p_patch(jg)%n_patch_verts_g))
      glb_phys_id_c(:) = 0
      glb_phys_id_e(:) = 0
      glb_phys_id_v(:) = 0

      ! Fill with own values of phys_id

      DO j = 1, p_patch(jg)%n_patch_cells
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF(.NOT.p_patch(jg)%cells%owner_mask(jl,jb)) CYCLE
        glb_phys_id_c(p_patch(jg)%cells%glb_index(j)) = p_patch(jg)%cells%phys_id(jl,jb)
      ENDDO

      DO j = 1, p_patch(jg)%n_patch_edges
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF(.NOT.p_patch(jg)%edges%owner_mask(jl,jb)) CYCLE
        glb_phys_id_e(p_patch(jg)%edges%glb_index(j)) = p_patch(jg)%edges%phys_id(jl,jb)
      ENDDO

      DO j = 1, p_patch(jg)%n_patch_verts
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        IF(.NOT.p_patch(jg)%verts%owner_mask(jl,jb)) CYCLE
        glb_phys_id_v(p_patch(jg)%verts%glb_index(j)) = p_patch(jg)%verts%phys_id(jl,jb)
      ENDDO

      ! Get global arrays by obtaining the global maximum
      ! Since p_max works on real valued arrays only, we have to convert to real and back

      glb_phys_id_c = INT( p_max(REAL(glb_phys_id_c,wp), p_comm_work) )
      glb_phys_id_e = INT( p_max(REAL(glb_phys_id_e,wp), p_comm_work) )
      glb_phys_id_v = INT( p_max(REAL(glb_phys_id_v,wp), p_comm_work) )

      ! Get the physical patches contained within current patch

      ! Set logical id of p_phys_patch from cells

      DO j = 1, p_patch(jg)%n_patch_cells_g
        jp = glb_phys_id_c(j)
        IF(jp<1 .OR. jp>max_phys_dom) THEN
          WRITE(message_text,'(a,i4,a,i12)') &
            & 'Patch ',jg,' contains illegal value for cells phys_id: ',jp
          CALL finish  (routine, TRIM(message_text))
        ENDIF
        IF(p_phys_patch(jp)%logical_id < 0) THEN
          p_phys_patch(jp)%logical_id = jg
        ELSE
          ! Check if no other patch uses the same phys_id
          IF(p_phys_patch(jp)%logical_id /= jg) THEN
            WRITE(message_text,'(a,i4,a,i12,a,i12)') &
             & 'Patch ',jg,' contains cells phys_id: ',jp, &
             & ', already used by patch: ',p_phys_patch(jp)%logical_id
            CALL finish  (routine, TRIM(message_text))
          ENDIF
        ENDIF
      ENDDO

      ! Make a check for egdes and verts

      DO j = 1, p_patch(jg)%n_patch_edges_g
        jp = glb_phys_id_e(j)
        IF(jp<1 .OR. jp>max_phys_dom) THEN
          WRITE(message_text,'(a,i4,a,i12)') &
            & 'Patch ',jg,' contains illegal value for edges phys_id: ',jp
          CALL finish  (routine, TRIM(message_text))
        ENDIF
        ! Check if no other patch uses the same phys_id
        IF(p_phys_patch(jp)%logical_id /= jg) THEN
          WRITE(message_text,'(a,i4,a,i12,a,i12)') &
           & 'Patch ',jg,' contains edges phys_id: ',jp, &
           & ', already used by patch: ',p_phys_patch(jp)%logical_id
          CALL finish  (routine, TRIM(message_text))
        ENDIF
      ENDDO

      DO j = 1, p_patch(jg)%n_patch_verts_g
        jp = glb_phys_id_v(j)
        IF(jp<1 .OR. jp>max_phys_dom) THEN
          WRITE(message_text,'(a,i4,a,i12)') &
            & 'Patch ',jg,' contains illegal value for verts phys_id: ',jp
          CALL finish  (routine, TRIM(message_text))
        ENDIF
        ! Check if no other patch uses the same phys_id
        IF(p_phys_patch(jp)%logical_id /= jg) THEN
          WRITE(message_text,'(a,i4,a,i12,a,i12)') &
           & 'Patch ',jg,' contains verts phys_id: ',jp, &
           & ', already used by patch: ',p_phys_patch(jp)%logical_id
          CALL finish  (routine, TRIM(message_text))
        ENDIF
      ENDDO

      ! Set up communication patterns in physical patches

      n = MAX(p_patch(jg)%n_patch_cells_g,p_patch(jg)%n_patch_edges_g,p_patch(jg)%n_patch_verts_g)
      ALLOCATE(owner (n))
      ALLOCATE(glbidx(n))

      DO jp = 1, max_phys_dom ! Loop over physical patches

        IF(p_phys_patch(jp)%logical_id /= jg) CYCLE ! do only for physical patches belonging to jg

        ! cells

        n = 0
        DO i = 1, p_patch(jg)%n_patch_cells_g
          IF(glb_phys_id_c(i) == jp) THEN
            n = n+1
            owner (n) = p_patch(jg)%cells%owner_g(i)
            glbidx(n) = i
          ENDIF
        ENDDO

        p_phys_patch(jp)%n_patch_cells = n

        IF(.NOT.my_process_is_mpi_seq()) THEN
          IF(p_pe_work == 0) THEN
            CALL setup_comm_pattern(n, owner, glbidx, &
              & p_patch(jg)%cells%loc_index, p_phys_patch(jp)%comm_pat_gather_c)
          ELSE
            ! We don't want to receive any data, i.e. the number of cells is 0
            ! and owner/global index are dummies!
            CALL setup_comm_pattern(0, owner, glbidx, &
              & p_patch(jg)%cells%loc_index, p_phys_patch(jp)%comm_pat_gather_c)
          ENDIF
        ENDIF

        ! edges

        n = 0
        DO i = 1, p_patch(jg)%n_patch_edges_g
          IF(glb_phys_id_e(i) == jp) THEN
            n = n+1
            owner (n) = p_patch(jg)%edges%owner_g(i)
            glbidx(n) = i
          ENDIF
        ENDDO

        p_phys_patch(jp)%n_patch_edges = n

        IF(.NOT.my_process_is_mpi_seq()) THEN
          IF(p_pe_work == 0) THEN
            CALL setup_comm_pattern(n, owner, glbidx, &
              & p_patch(jg)%edges%loc_index, p_phys_patch(jp)%comm_pat_gather_e)
          ELSE
            ! We don't want to receive any data, i.e. the number of edges is 0
            ! and owner/global index are dummies!
            CALL setup_comm_pattern(0, owner, glbidx, &
              & p_patch(jg)%edges%loc_index, p_phys_patch(jp)%comm_pat_gather_e)
          ENDIF
        ENDIF

        ! verts

        n = 0
        DO i = 1, p_patch(jg)%n_patch_verts_g
          IF(glb_phys_id_v(i) == jp) THEN
            n = n+1
            owner (n) = p_patch(jg)%verts%owner_g(i)
            glbidx(n) = i
          ENDIF
        ENDDO

        p_phys_patch(jp)%n_patch_verts = n

        IF(.NOT.my_process_is_mpi_seq()) THEN
          IF(p_pe_work == 0) THEN
            CALL setup_comm_pattern(n, owner, glbidx, &
              & p_patch(jg)%verts%loc_index, p_phys_patch(jp)%comm_pat_gather_v)
          ELSE
            ! We don't want to receive any data, i.e. the number of verts is 0
            ! and owner/global index are dummies!
            CALL setup_comm_pattern(0, owner, glbidx, &
              & p_patch(jg)%verts%loc_index, p_phys_patch(jp)%comm_pat_gather_v)
          ENDIF
        ENDIF

      ENDDO

      DEALLOCATE(owner)
      DEALLOCATE(glbidx)
      DEALLOCATE(glb_phys_id_c)
      DEALLOCATE(glb_phys_id_e)
      DEALLOCATE(glb_phys_id_v)

    ENDDO

    ! Get total number of physical patches

    DO jp = max_phys_dom, 1, -1
      IF(p_phys_patch(jp)%logical_id > 0) EXIT
    ENDDO

    n_phys_dom = jp

    ! Print info and check if there are no unused physical patches

    DO jp = 1, n_phys_dom
      WRITE(message_text,'(a,i4,a,i4)') &
        & 'Physical domain ',jp,' belongs to logical domain ',p_phys_patch(jp)%logical_id
      CALL message (routine, message_text)
      IF(p_phys_patch(jp)%logical_id < 0) THEN
        WRITE(message_text,'(a,i4)') 'Missing physical domain # ',jp
        CALL finish (routine, message_text)
      ENDIF
    ENDDO

  END SUBROUTINE setup_phys_patches

  !-------------------------------------------------------------------------
  !>
  !!               Calculates local line/block indices l_idx, l_blk
  !!               from global line/block indices g_idx, g_blk
  !!               using the mapping in loc_index
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE get_local_index(loc_index, g_idx, g_blk, l_idx, l_blk, opt_mode)

    !

    INTEGER, INTENT(in) :: loc_index(:), g_idx, g_blk
    INTEGER, INTENT(in), OPTIONAL :: opt_mode

    INTEGER, INTENT(out) :: l_idx, l_blk

    INTEGER mode, j_g, j_l

    !-----------------------------------------------------------------------

    ! mode controls the behaviour if the global index is valid, but outside the local domain:
    ! mode > 0 : Map it to the next greater index (even if this is outside the local domain)
    ! mode < 0 : Map it to the next smaller index (even if this is 0)
    ! mode = 0 : Map it to a negative value (negative of global index)

    IF(PRESENT(opt_mode)) THEN
      mode = opt_mode
    ELSE
      mode = 0
    ENDIF

    ! Get 1D global index from g_idx, g_blk

    j_g = idx_1d(g_idx, g_blk)

    ! If j_g is invalid, return 0

    IF(j_g < 1 .OR. j_g > UBOUND(loc_index,1)) THEN
      l_idx = idx_no(0)
      l_blk = blk_no(0)
      RETURN
    ENDIF

    ! Get local index corresponding to j_g; if this is outside the local domain,
    ! do what was requested by the mode setting

    j_l = loc_index(j_g)

    IF(j_l < 0) THEN
      ! Please note: ABS(j_l) is the (last valid one + 1) or n_local+1 if after the last one
      IF(mode > 0) THEN
        j_l = ABS(j_l)
      ELSE IF(mode < 0) THEN
        j_l = ABS(j_l)-1
      ELSE
        j_l = -j_g
      ENDIF
    ENDIF

    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE get_local_index

  !-------------------------------------------------------------------------
  !>
  !! Maps indices which point outside domain (returned as <0 by get_local_index)
  !! to valid ones which are from the same set of neighbors for the given
  !! cell/edge/vert
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Oct 2011
  !!
  SUBROUTINE remap_index(loc_index, l_idx, l_blk)

    INTEGER, INTENT(in) :: loc_index(:)
    INTEGER, INTENT(inout) :: l_idx, l_blk

    INTEGER :: j_l, j_g

    IF(l_idx>=0) RETURN ! Nothing to do

    j_g = idx_1d(-l_idx, l_blk) ! Global index

    ! Safety check only, global index should be in correct range
    IF(j_g < 1 .OR. j_g > UBOUND(loc_index,1)) CALL finish('remap_index','Invalid global index')

    j_l = loc_index(j_g)

    ! Safety check only, local index should be negative
    IF(j_l >= 0) CALL finish('remap_index','Invalid local index')

    ! Remap in the same way as previously in get_local_index (for compatibility only)
    j_l = MAX(ABS(j_l)-1,1)

    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE remap_index

  !-------------------------------------------------------------------------
  !>
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE divide_subset_geometric(subset_flag, n_proc, owner)

    INTEGER, INTENT(in)    :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    INTEGER, INTENT(out)   :: owner(:) ! receives the owner PE for every cell
    ! (-1 for cells not in subset)

    INTEGER :: i, j, jl, jb, jl_v, jb_v, nc
    REAL(wp), ALLOCATABLE :: cell_desc(:,:)
    REAL(wp)              :: cclat, cclon

    !-----------------------------------------------------------------------

    IF(p_pe_work==0) THEN
      IF(divide_for_radiation) THEN
        PRINT *,'divide_patch: Using geometric area subdivision for radiation'
      ELSE
        PRINT *,'divide_patch: Using geometric area subdivision (normal)'
      ENDIF
    ENDIF

    ! Fill the cell_desc array, it must contain:
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner

    ALLOCATE(cell_desc(4,wrk_divide_patch%n_patch_cells))

    cell_desc(1:2,:) = 1.d99 ! for fining min lat/lon

    nc = 0
    DO j = 1, wrk_divide_patch%n_patch_cells

      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index

      IF(subset_flag(j)<=0) CYCLE ! Cell not in subset

      nc = nc+1 ! Cell counter

      IF(.NOT. divide_for_radiation) THEN

        ! Using the center of the cells for geometric subdivision leads
        ! to "toothed" edges of the subdivision area
        ! Thus we use the minimum lat/lon as subdision criterion.

        DO i=1,wrk_divide_patch%cells%num_edges(jl,jb)
          jl_v = wrk_divide_patch%cells%vertex_idx(jl,jb,i)
          jb_v = wrk_divide_patch%cells%vertex_blk(jl,jb,i)
          cell_desc(1,nc) = MIN(cell_desc(1,nc),wrk_divide_patch%verts%vertex(jl_v,jb_v)%lat)
          cell_desc(2,nc) = MIN(cell_desc(2,nc),wrk_divide_patch%verts%vertex(jl_v,jb_v)%lon)
        ENDDO

      ELSE

        ! Patch division for radiation calculations:
        ! To minimize load imbalance, every patch contains 10 areas
        ! distributed in a way similar as the "diamonds" in GME
        ! This is accomplished by mapping all cells to one section  
        ! lying in the NH and having a width of 0.4*pi (72 deg) 

        cclat = wrk_divide_patch%cells%center(jl,jb)%lat
        cclon = wrk_divide_patch%cells%center(jl,jb)%lon

        IF (cclat>=0._wp .AND. cclon>=-0.2_wp*pi .AND. cclon<=0.2_wp*pi) THEN
          cell_desc(1,nc) = cclat
          cell_desc(2,nc) = cclon
        ELSE IF (cclat>=0._wp .AND. cclon>=0.2_wp*pi .AND. cclon<=0.6_wp*pi) THEN
          cell_desc(1,nc) = cclat
          cell_desc(2,nc) = cclon - 0.4_wp*pi
        ELSE IF (cclat>=0._wp .AND. cclon>=0.6_wp*pi) THEN
          cell_desc(1,nc) = cclat
          cell_desc(2,nc) = cclon - 0.8_wp*pi
        ELSE IF (cclat>=0._wp .AND. cclon<=-0.6_wp*pi) THEN
          cell_desc(1,nc) = cclat
          cell_desc(2,nc) = cclon + 0.8_wp*pi
        ELSE IF (cclat>=0._wp .AND. cclon>=-0.6_wp*pi .AND. cclon<=-0.2_wp*pi) THEN
          cell_desc(1,nc) = cclat
          cell_desc(2,nc) = cclon + 0.4_wp*pi
        ELSE IF (cclat<0._wp .AND. (cclon<=-0.8_wp*pi .OR. cclon>=0.8_wp*pi)) THEN
          cell_desc(1,nc) = -cclat
          cell_desc(2,nc) = cclon + pi
        ELSE IF (cclat<0._wp .AND. cclon>=-0.8_wp*pi .AND. cclon<=-0.4_wp*pi) THEN
          cell_desc(1,nc) = -cclat
          cell_desc(2,nc) = cclon + 0.6_wp*pi
        ELSE IF (cclat<0._wp .AND. cclon>=-0.4_wp*pi .AND. cclon<=0.0_wp*pi) THEN
          cell_desc(1,nc) = -cclat
          cell_desc(2,nc) = cclon + 0.2_wp*pi
        ELSE IF (cclat<0._wp .AND. cclon>=0.0_wp*pi .AND. cclon<=0.4_wp*pi) THEN
          cell_desc(1,nc) = -cclat
          cell_desc(2,nc) = cclon - 0.2_wp*pi
        ELSE IF (cclat<0._wp .AND. cclon>=0.4_wp*pi .AND. cclon<=0.8_wp*pi) THEN
          cell_desc(1,nc) = -cclat
          cell_desc(2,nc) = cclon - 0.6_wp*pi
        ENDIF

        IF (cell_desc(2,nc)>pi) THEN
          cell_desc(2,nc) = cell_desc(2,nc) - 2._wp*pi
        ELSE IF (cell_desc(2,nc)<-pi) THEN
          cell_desc(2,nc) = cell_desc(2,nc) + 2._wp*pi
        ENDIF

      ENDIF

      cell_desc(3,nc) = REAL(nc,wp)
      cell_desc(4,nc) = 0.0_wp

    ENDDO

    CALL divide_cells_by_location(nc, cell_desc, 0, n_proc-1)

    ! After divide_cells_by_location the cells are sorted by owner,
    ! order them by original cell numbers again

    CALL sort_array_by_row(cell_desc(:,1:nc), 3)

    ! Set owner list (of complete patch)

    owner(:) = -1
    nc = 0 ! Counts cells in subset

    DO j = 1, wrk_divide_patch%n_patch_cells
      IF(subset_flag(j)>0) THEN
        nc = nc+1
        owner(j) = INT(cell_desc(4,nc))
      ENDIF
    ENDDO

    DEALLOCATE(cell_desc)

  END SUBROUTINE divide_subset_geometric

  !-------------------------------------------------------------------------
  !>
  !! Actually divides geometrically by location on cpu_a .. cpu_b
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  RECURSIVE SUBROUTINE divide_cells_by_location(n_cells,cell_desc,cpu_a,cpu_b)

    INTEGER, INTENT(in) :: n_cells, cpu_a, cpu_b

    REAL(wp), INTENT(inout) :: cell_desc(4,n_cells)
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner

    INTEGER cpu_m, n_cells_m
    REAL(wp) :: xmax(2), xmin(2)
    !-----------------------------------------------------------------------

    ! If there is only 1 CPU for distribution, we are done

    IF(cpu_a==cpu_b) THEN
      cell_desc(4,:) = REAL(cpu_a,wp)
      RETURN
    ENDIF

    ! Get geometric extensions and total number of points of all batches

    xmin(1) = MINVAL(cell_desc(1,:))
    xmin(2) = MINVAL(cell_desc(2,:))
    xmax(1) = MAXVAL(cell_desc(1,:))
    xmax(2) = MAXVAL(cell_desc(2,:))

    ! Get dimension with biggest distance from min to max
    ! and sort cells in this dimension

    IF(xmax(1)-xmin(1) >= xmax(2)-xmin(2)) THEN
      CALL sort_array_by_row(cell_desc, 1)
    ELSE
      CALL sort_array_by_row(cell_desc, 2)
    ENDIF

    ! CPU number where to split CPU set

    cpu_m = (cpu_a+cpu_b-1)/2

    ! If the number of CPUs is not even, we have to split the set
    ! of cells accordingly into to differntly sized halfes
    ! in order to end with an equal number of points on every CPU.
    ! Note that DOUBLE arithmetic is used since the integer size
    ! may be exceeded in this calculation!

    n_cells_m = INT(DBLE(n_cells)*DBLE(cpu_m-cpu_a+1)/DBLE(cpu_b-cpu_a+1))

    ! If there are only two CPUs, we are done

    IF(cpu_b == cpu_a+1) THEN
      cell_desc(4,1:n_cells_m)         = REAL(cpu_a,wp)
      cell_desc(4,n_cells_m+1:n_cells) = REAL(cpu_b,wp)
      RETURN
    ENDIF

    ! Further divide both halves recursively

    CALL divide_cells_by_location(n_cells_m,cell_desc(:,1:n_cells_m),cpu_a,cpu_m)
    CALL divide_cells_by_location(n_cells-n_cells_m,&
      & cell_desc(:,n_cells_m+1:n_cells),cpu_m+1,cpu_b)

  END SUBROUTINE divide_cells_by_location

  !-------------------------------------------------------------------------
  !>
  !! Special quicksort implementation for sorting a 2D array by one selected row.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  RECURSIVE SUBROUTINE sort_array_by_row(x,row)

    !

    INTEGER, INTENT(in) :: row ! number of row for sorting

    REAL(wp), INTENT(inout) :: x(:,:) ! array to be sorted

    REAL(wp) :: p
    REAL(wp), ALLOCATABLE :: y(:,:)
    INTEGER :: n, ipiv, ix, iy, i
    !-----------------------------------------------------------------------

    n = SIZE(x,2)

    IF(n<=1) RETURN

    ALLOCATE(y(SIZE(x,1),n))

    ipiv = (n+1)/2
    p = x(row,ipiv)
    ix = 0
    iy = 1
    y(:,1) = x(:,ipiv) ! Store pivot

    DO i=1,n
      IF(i==ipiv) CYCLE
      IF(x(row,i) < p) THEN
        ix = ix+1
        x(:,ix) = x(:,i)
      ELSE
        iy = iy+1
        y(:,iy) = x(:,i)
      ENDIF
    ENDDO

    x(:,ix+1:ix+iy) = y(:,1:iy)

    ipiv = ix+1 ! New pivot location

    DEALLOCATE(y)

    IF(ipiv>2)   CALL sort_array_by_row(x(:,:ipiv-1),row)
    IF(ipiv<n-1) CALL sort_array_by_row(x(:,ipiv+1:),row)

  END SUBROUTINE sort_array_by_row

  !-----------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Utility function: Compute number of points in nproma-based array
  !!
  !! @par Revision History
  !! F. Prill, Nov 2011
  !!
  FUNCTION count_entries(start_blk, start_idx, &
    &                    end_blk, end_idx)
    INTEGER :: count_entries
    INTEGER, INTENT(IN) :: start_blk, start_idx, end_blk, end_idx

    count_entries = nproma * (end_blk - start_blk)  &
      &             + end_idx - start_idx + 1
  END FUNCTION count_entries

  !-----------------------------------------------------------------------


#ifdef HAVE_METIS
  !-------------------------------------------------------------------------
  !>
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE divide_subset_metis(subset_flag, n_proc, owner)

    INTEGER, INTENT(in)    :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    INTEGER, INTENT(out)   :: owner(:) ! receives the owner PE for every cell
    ! (-1 for cells not in subset)


    INTEGER :: i, j, jl, jb, jl_n, jb_n, na, nc
    INTEGER, ALLOCATABLE :: local_index(:), tmp(:)
    INTEGER, ALLOCATABLE :: metis_xadj(:), metis_adjncy(:)
    INTEGER :: metis_options(5), metis_edgecut

    !-----------------------------------------------------------------------

    IF(p_pe_work==0) PRINT *,'divide_patch: Using METIS for area subdivision'

    ! Get the local index (i.e. contiguous numbers) of the cells in the subset.
    ! Since the METIS rountine is called with the option for C-Style numbering,
    ! we let local_index start with 0 here!

    ALLOCATE(local_index(wrk_divide_patch%n_patch_cells))
    local_index(:) = -1

    nc = 0
    DO j = 1, wrk_divide_patch%n_patch_cells
      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index
      IF(subset_flag(j)>0) THEN
        local_index(j) = nc
        nc = nc+1
      ENDIF
    ENDDO

    ! Construct adjacency structure of graph
    ! Please note that the Metis vertices are our grid cells!

    ALLOCATE(metis_xadj(0:wrk_divide_patch%n_patch_cells))
    ALLOCATE(metis_adjncy(wrk_divide_patch%n_patch_cells*wrk_divide_patch%cell_type))
    metis_options(:) = 0

    metis_xadj(0) = 0
    na = 0 ! Counts adjacency entries
    nc = 0 ! Counts cells in subset

    DO j = 1, wrk_divide_patch%n_patch_cells
      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index

      IF(subset_flag(j)<=0) CYCLE ! Cell not in subset

      ! Loop over all neighbors of the cell and include neighbors
      ! which are also in the subset in the adjacency list

      nc = nc+1
      DO i = 1,wrk_divide_patch%cells%num_edges(jl,jb)
        jl_n = wrk_divide_patch%cells%neighbor_idx(jl,jb,i)
        jb_n = wrk_divide_patch%cells%neighbor_blk(jl,jb,i)
        jn = idx_1d(jl_n,jb_n)

        ! Neighbor not existing
        IF(jl_n<1 .or. jb_n<1 .or. jn>wrk_divide_patch%n_patch_cells) CYCLE

        ! Neighbor not in subset
        IF(subset_flag(jn)<0) CYCLE

        na = na+1
        metis_adjncy(na) = local_index(jn)
      ENDDO
      metis_xadj(nc) = na
    ENDDO

    ! Divide graph with Metis - currently no weights are assigned to vertices/edges

    ALLOCATE(tmp(nc))
    CALL metis_partgraphrecursive(nc, metis_xadj, metis_adjncy, 0, 0, 0, 0, &
      & n_proc, metis_options, metis_edgecut, tmp)

    ! The owner list delivered by mets (in tmp) contains only owners for
    ! cells in subset, scatter it into the global owner list

    owner(:) = -1
    nc = 0 ! Counts cells in subset

    DO j = 1, wrk_divide_patch%n_patch_cells
      IF(subset_flag(j)>0) THEN
        nc = nc+1
        owner(j) = tmp(nc)
      ENDIF
    ENDDO

    DEALLOCATE(metis_xadj, metis_adjncy, local_index, tmp)

  END SUBROUTINE divide_subset_metis
#endif

END MODULE mo_subdivision

