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
MODULE mo_setup_subdivision
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
    & min_rledge, max_rledge, min_rlvert, max_rlvert, max_phys_dom,  &
    & min_rlcell_int, min_rledge_int, min_rlvert_int, max_hw, max_dom
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish, message, message_text,    &
    &                              get_filename_noext

  USE mo_run_config,         ONLY: msg_level
  USE mo_io_units,           ONLY: find_next_free_unit, filename_max
  USE mo_model_domain,       ONLY: t_patch, p_patch_local_parent
  USE mo_mpi,                ONLY: p_bcast, p_sum, proc_split, get_my_global_mpi_id, global_mpi_barrier
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_UNDEFINED, MPI_COMM_NULL, p_int, &
       mpi_in_place, mpi_max, mpi_success, process_mpi_all_comm
#endif
  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_test, &
    & my_process_is_mpi_seq, process_mpi_all_test_id, process_mpi_all_workroot_id, &
    & p_pe_work, p_n_work, get_my_mpi_all_id, my_process_is_mpi_parallel

  USE mo_parallel_config,       ONLY:  nproma, ldiv_phys_dom, &
    & division_method, division_file_name, n_ghost_rows, div_from_file,   &
    & div_geometric, ext_div_medial, ext_div_medial_cluster, ext_div_medial_redrad, &
    & ext_div_medial_redrad_cluster, ext_div_from_file, redrad_split_factor

#ifdef HAVE_METIS
  USE mo_parallel_config,    ONLY: div_metis
#endif
  USE mo_communication,      ONLY: blk_no, idx_no, idx_1d
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, patch_weight
  USE mo_alloc_patches,ONLY: allocate_basic_patch, allocate_remaining_patch, &
                             deallocate_basic_patch, deallocate_patch
  USE mo_decomposition_tools, ONLY: t_decomposition_structure, divide_geometric_medial, &
    & read_ascii_decomposition
#ifndef __ICON_OCEAN_ONLY__
  USE mo_dump_restore,        ONLY: dump_all_domain_decompositions
#endif
  USE mo_math_utilities,      ONLY: geographical_to_cartesian
  USE mo_master_control,      ONLY: get_my_process_type, ocean_process
  USE mo_grid_config,         ONLY: use_dummy_cell_closure
  USE mo_util_sort, ONLY: quicksort, insertion_sort

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: decompose_domain

  ! Private flag if patch should be divided for radiation calculation
  LOGICAL :: divide_for_radiation = .FALSE.

  TYPE nb_flag_list_elem
    INTEGER, ALLOCATABLE :: idx(:)
  END TYPE nb_flag_list_elem

CONTAINS

  !------------------------------------------------------------------
  !>
  !!  Divide patches and interpolation states.
  !!
  !!  If the parameter @p opt_n_procs is given, a domain decomposition
  !!  for @p opt_n_procsprocessors is done (called from a single
  !!  processor run).
  !!
  !!  @note Despite its name, this subroutine is also called in
  !!        sequential runs where it simply copies (and initializes)
  !!        the patch data structure.
  !!
  SUBROUTINE decompose_domain( p_patch, p_patch_global, opt_n_procs )

    TYPE(t_patch), INTENT(INOUT), TARGET :: p_patch(n_dom_start:)
    TYPE(t_patch), INTENT(INOUT), TARGET :: p_patch_global(n_dom_start:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_n_procs

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_subdivision:decompose_domain")

    ! Local variables:
    INTEGER :: jg, jgp, jc, jgc, n, &
      &        n_procs_decomp
    INTEGER :: nprocs(p_patch_global(1)%n_childdom)
    INTEGER, POINTER :: cell_owner(:)
    INTEGER, POINTER :: cell_owner_p(:)
    REAL(wp) :: weight(p_patch_global(1)%n_childdom)
    TYPE(t_patch), ALLOCATABLE, TARGET :: p_patch_out(:), p_patch_lp_out(:)
    TYPE(t_patch), POINTER :: wrk_p_parent_patch_g
    ! LOGICAL :: l_compute_grid
    INTEGER :: order_type_of_halos

    CALL message(routine, 'start of domain decomposition')

    ! This subroutine has 3 different operating modes:
    !
    ! - If opt_n_procs is set, this routine should be called from a
    !   single processor run and it splits the domain into opt_n_procs
    !   parts.
    !
    ! - Otherwise, if the subroutine is called from a parallel run
    !   then the domain is split according to the number of processors
    !
    ! - Third, when called from a single processor run, the domain is
    !   not split but simply copied (and the p_patch_local_parent is
    !   initialized properly).
    !
    n_procs_decomp = p_n_work
    IF (PRESENT(opt_n_procs)) THEN
      n_procs_decomp = opt_n_procs
      IF (my_process_is_mpi_parallel()) &
        CALL finish(routine, 'Call with opt_n_procs in parallel mode')

      ! Allocate p_patch_out, p_patch_lp_out so that they can keep
      ! all domain decompositions for one patch
      ALLOCATE(p_patch_out(opt_n_procs))
      ALLOCATE(p_patch_lp_out(opt_n_procs))
    ENDIF

    ! -----------------------------------------------------------------------------
    ! Check for processor splitting

    ! Check if the processor set should be split for patches of the 1st generation.
    ! This is the case if patch_weight > 0 for at least one of the root's childs.
    ! For clarity, we require that patch_weight=0 for any other patch

    ! Get weights for 1st generation patches
    proc_split = .FALSE.
    DO jc = 1, p_patch_global(1)%n_childdom
      jgc = p_patch_global(1)%child_id(jc)
      weight(jc) = patch_weight(jgc)
      IF(weight(jc) > 0._wp) proc_split = .TRUE.
    ENDDO

    ! Check if weights for other patches are 0 (for clarity only)
    IF(patch_weight(1) /= 0._wp) &
      CALL finish(routine,'Weight for root patch must be 0')
    DO jg = 2, n_dom
      jgp = p_patch_global(jg)%parent_id
      IF(jgp /= 1 .AND. patch_weight(jg) > 0._wp) &
        CALL finish(routine,'Weight for higher level patch must be 0')
    ENDDO
    ! -----------------------------------------------------------------------------

    IF(proc_split) THEN

      ! -----------------------------------------------------------------------------
      ! CASE 1: DECOMPOSE WITH PROCESSOR SPLITTING
      ! -----------------------------------------------------------------------------

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

      ! -----------------------------------------------------------------------------
      ! CASE 2: STANDARD DECOMPOSITION, NO SPLITTING
      ! -----------------------------------------------------------------------------
      ! No splitting, proc0, n_proc are identical for all patches

      IF(p_pe_work==0) WRITE(0,*) 'No splitting of processor grid'
      p_patch(:)%n_proc = n_procs_decomp
      p_patch(:)%proc0  = 0

    ENDIF

#ifdef NOMPI
    p_patch(:)%comm = 0
#else
    p_patch(:)%comm = MPI_COMM_NULL
#endif
    p_patch(:)%rank = -1

    ! -----------------------------------------------------------------------------
    ! Divide patches
    ! -----------------------------------------------------------------------------

    DO jg = n_dom_start, n_dom

      jgp = p_patch_global(jg)%parent_id

      IF(jg == n_dom_start) THEN
        NULLIFY(wrk_p_parent_patch_g) ! Must be NULL for global patch
      ELSE
        wrk_p_parent_patch_g => p_patch_global(jgp)
      ENDIF

      ! Set division method, divide_for_radiation is only used for patch 0

      divide_for_radiation = (jg == 0)

      ! Here comes the actual domain decomposition:
      ! Every cells gets assigned an owner.

      ALLOCATE(cell_owner(p_patch_global(jg)%n_patch_cells))
      IF(jg > n_dom_start) THEN
        ALLOCATE(cell_owner_p(p_patch_global(jgp)%n_patch_cells))
      END IF

      IF (my_process_is_mpi_parallel()) THEN
        CALL divide_patch_cells(p_patch_global(jg), jg, p_patch(jg)%n_proc, &
             p_patch(jg)%proc0, cell_owner, wrk_p_parent_patch_g, &
             p_patch(jg)%cells%radiation_owner)
      ELSE
        cell_owner(:) = 0 ! trivial "decomposition"
      END IF

      IF(jg > n_dom_start) THEN
        ! Assign the cell owners of the current patch to the parent cells
        ! for the construction of the local parent, set "-1" elsewhere.
        CALL divide_parent_cells(p_patch_global(jg),cell_owner,cell_owner_p)
      END IF

      DEALLOCATE(p_patch_global(jg)%cells%phys_id)

      ! Please note: Previously, for jg==0 no ghost rows were set.
      ! Currently, we need ghost rows for jg==0 also for dividing the int state and grf state
      ! Have still to check if int state/grf state is needed at all for jg==0,
      ! if this is not the case, the ghost rows can be dropped again.

      IF(PRESENT(opt_n_procs)) THEN

        ! ----------------------------------------------
        ! operation mode 1: "opt_n_proc" present
        ! ----------------------------------------------


        ! Calculate all basic patches and output them at once.
        ! This is done in this way (calculating all before output) since
        ! we need the biggest dimension of all patches before real NetCDF output starts.
        ! All arrays in the patch which are not scaling (i.e. having global dimensions)
        ! are discarded so that the stored patches shouldn't need much more space
        ! than the global patch descriptions.

        ! If the storage for the patches should get a problem nontheless, there has to be found
        ! a way to output them before all are calculated.

        ! Do all domain decompositions

#ifndef __xlC__
!$OMP PARALLEL DO PRIVATE(n)
#endif
        DO n = 1, opt_n_procs

          WRITE(0,'(2(a,i0))') 'Dividing patch ',jg,' for proc ',n-1
          p_patch_out(n)%n_proc = p_patch(jg)%n_proc
          p_patch_out(n)%proc0  = p_patch(jg)%proc0
          order_type_of_halos = 1
          ! CALL divide_patch(p_patch_out(n), p_patch_global(jg), cell_owner, n_ghost_rows, .TRUE., n-1)
          CALL divide_patch(p_patch_out(n), p_patch_global(jg), cell_owner, n_ghost_rows, order_type_of_halos, n-1)
          CALL discard_large_arrays(p_patch_out(n), n)

        ENDDO
#ifndef __xlC__
!$OMP END PARALLEL DO
#endif

        IF(jg > n_dom_start) THEN
#ifndef __xlC__
!$OMP PARALLEL DO PRIVATE(n)
#endif
          DO n = 1, opt_n_procs

            ! Divide local parent
            order_type_of_halos = 0
            ! CALL divide_patch(p_patch_lp_out(n), p_patch_global(jgp), cell_owner_p, 1, .FALSE., n-1)
            CALL divide_patch(p_patch_lp_out(n), p_patch_global(jgp), cell_owner_p, 1, order_type_of_halos, n-1)
            CALL discard_large_arrays(p_patch_lp_out(n), n)

          ENDDO
#ifndef __xlC__
!$OMP END PARALLEL DO
#endif
        ENDIF

#ifndef __ICON_OCEAN_ONLY__
        ! Dump domain decompositions to NetCDF
        IF(jg > n_dom_start) THEN
          CALL dump_all_domain_decompositions(p_patch_out, p_patch_lp_out)
        ELSE
          CALL dump_all_domain_decompositions(p_patch_out)
        ENDIF
#endif

        ! Deallocate patch arrays
        DO n = 1, opt_n_procs
          CALL deallocate_patch(p_patch_out(n),lddmode=.TRUE.)
          IF(jg > n_dom_start) CALL deallocate_patch(p_patch_lp_out(n),lddmode=.TRUE.)
        ENDDO

      ELSE

        ! ----------------------------------------------
        ! operation modes 2,3: "opt_n_procs" not present
        ! ----------------------------------------------
        ! order_type_of_halos = 0 order for parent (l_compute_grid = false)
        !                       1 order root grid  (l_compute_grid = true)
        !                       2 all halos go to the end, for ocean
        IF (get_my_process_type() == ocean_process) THEN
          order_type_of_halos = 2
        ELSE
          order_type_of_halos = 1
        ENDIF

        ! CALL divide_patch(p_patch(jg), p_patch_global(jg), cell_owner, n_ghost_rows, l_compute_grid, p_pe_work)
        CALL divide_patch(p_patch(jg), p_patch_global(jg), cell_owner, n_ghost_rows, order_type_of_halos, p_pe_work)

        IF(jg > n_dom_start) THEN
          order_type_of_halos = 0
          ! CALL divide_patch(p_patch_local_parent(jg), p_patch_global(jgp), cell_owner_p, 1, .FALSE., p_pe_work)
          CALL divide_patch(p_patch_local_parent(jg), p_patch_global(jgp), cell_owner_p, 1, &
            & order_type_of_halos,  p_pe_work)
        ENDIF

      ENDIF ! IF (PRESENT(opt_n_procs))

      DEALLOCATE(cell_owner)
      IF(jg > n_dom_start) DEALLOCATE(cell_owner_p)

    ENDDO

    ! The global patches may be discarded now
    DO jg = n_dom_start, n_dom
      CALL deallocate_basic_patch( p_patch_global(jg) )
    ENDDO

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

  SUBROUTINE divide_patch_cells(wrk_p_patch_g, patch_no, n_proc, proc0, cell_owner, &
                                wrk_p_parent_patch_g, radiation_owner)

    TYPE(t_patch), INTENT(INOUT) :: wrk_p_patch_g
    INTEGER, INTENT(IN)  :: patch_no !> The patch number,
                                     ! used to identify patch specific decomposition
    INTEGER, INTENT(IN)  :: n_proc !> Number of processors for split
    INTEGER, INTENT(IN)  :: proc0  !> First processor of patch
    INTEGER, POINTER :: cell_owner(:) !> Cell division
    TYPE(t_patch), POINTER :: wrk_p_parent_patch_g
    INTEGER, POINTER :: radiation_owner(:)

    TYPE(t_decomposition_structure)  :: decomposition_struct

    INTEGER :: n, i, j, jl, jb, jl_p, jb_p
    INTEGER, ALLOCATABLE :: flag_c(:), tmp(:)
    CHARACTER(LEN=filename_max) :: use_division_file_name ! if div_from_file

    ! Please note: Unfortunatly we cannot use p_io for doing I/O,
    ! since this might be the test PE which is never calling this routine
    ! (this is the case in the actual setup).
    ! Thus we use the worker PE 0 for I/O and don't use message() for output.

    NULLIFY(radiation_owner)

    IF (division_method(patch_no) == div_from_file     .OR. &
        division_method(patch_no) == ext_div_from_file .OR. &
        division_method(patch_no) > 100) THEN

      ! only procs 0 will decompose and then broadcast
      IF (p_pe_work == 0) THEN

        IF (division_method(patch_no) == div_from_file) THEN

          ! Area subdivision is read from file

          IF (division_file_name(patch_no) == "") THEN
            use_division_file_name = &
              & TRIM(get_filename_noext(wrk_p_patch_g%grid_filename))//'.cell_domain_ids'
          ELSE
            use_division_file_name = division_file_name(patch_no)
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

        ELSEIF (division_method(patch_no) == ext_div_from_file) THEN

          IF (division_file_name(patch_no) == "") THEN
            use_division_file_name = &
              & TRIM(get_filename_noext(wrk_p_patch_g%grid_filename))//'.cell_domain_ids'
          ELSE
            use_division_file_name = division_file_name(patch_no)
          ENDIF

          CALL read_ascii_decomposition(use_division_file_name, cell_owner, wrk_p_patch_g%n_patch_cells)

        ELSE

          !----------------------------------------------------------
          ! external decompositions
          ! just to make sure that the radiation onwer is not acitve

          ! fill decomposition_structure
          CALL fill_wrk_decomposition_struct(decomposition_struct, wrk_p_patch_g)

          SELECT CASE (division_method(patch_no))

          CASE (ext_div_medial)
            CALL divide_geometric_medial(decomposition_struct, &
              & decomposition_size = n_proc, &
              & cluster = .false.,           &
              & cells_owner = cell_owner)

          CASE (ext_div_medial_cluster)
            CALL divide_geometric_medial(decomposition_struct, &
              & decomposition_size = n_proc, &
              & cluster = .true.,            &
              & cells_owner = cell_owner)

          CASE (ext_div_medial_redrad)
            CALL divide_geometric_medial(decomposition_struct, &
              & decomposition_size = n_proc, &
              & cluster = .false.,           &
              & cells_owner = cell_owner,    &
              & radiation_onwer = radiation_owner,      &
              & radiation_split_factor = redrad_split_factor)

          CASE (ext_div_medial_redrad_cluster)
            CALL divide_geometric_medial(decomposition_struct, &
              & decomposition_size = n_proc, &
              & cluster = .true.,            &
              & cells_owner = cell_owner,    &
              & radiation_onwer = radiation_owner,      &
              & radiation_split_factor = redrad_split_factor)


          CASE DEFAULT
            CALL finish('divide_patch_cells', 'Unkown division_method')

          END SELECT

          ! clean decomposition_struct
          CALL clean_wrk_decomposition_struct(decomposition_struct)

        ENDIF ! subddivision_method(patch_no)

        ! Quick check for correct values
        IF(MINVAL(cell_owner(:)) < 0 .or. &
           MAXVAL(cell_owner(:)) >= n_proc) THEN
          WRITE(0,*) "n_porc=",n_proc, " MINAVAL=", MINVAL(cell_owner(:)), &
            " MAXVAL=", MAXVAL(cell_owner(:))
          CALL finish('divide_patch','Illegal subdivision in input file')
        ENDIF

      ENDIF ! IF (p_pe_work == 0)

      ! broadcast cell_owner array to other processes
      CALL p_bcast(cell_owner, 0, comm=p_comm_work)

      IF (division_method(patch_no) == ext_div_medial_redrad_cluster .OR. &
        & division_method(patch_no) == ext_div_medial_redrad) THEN
        ! distribute the radiation owner
        IF (p_pe_work /= 0) &
            ALLOCATE(radiation_owner(wrk_p_patch_g%n_patch_cells))

        CALL p_bcast(radiation_owner, 0, comm=p_comm_work)
        wrk_p_patch_g%cells%radiation_owner => radiation_owner

      ENDIF

    ELSE ! IF (division_method(patch_no) == div_from_file     .OR. &
         !     division_method(patch_no) == ext_div_from_file .OR. &
         !     division_method(patch_no) > 100)


      ! Built-in subdivison methods

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

          flag_c(idx_1d(jl_p, jb_p)) = MAX(1,wrk_p_patch_g%cells%phys_id(jl,jb))
        ENDDO

        ! Divide subset of patch
        ! Receives the PE  numbers for every cell
        ALLOCATE(tmp(wrk_p_parent_patch_g%n_patch_cells))

        IF(division_method(patch_no) == div_geometric) THEN
          CALL divide_subset_geometric( flag_c, n_proc, wrk_p_parent_patch_g, tmp)
#ifdef HAVE_METIS
        ELSE IF(division_method(patch_no) == div_metis) THEN
          CALL divide_subset_metis( flag_c, n_proc, wrk_p_parent_patch_g, tmp)
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

        IF(division_method(patch_no) == div_geometric) THEN
          CALL divide_subset_geometric(flag_c, n_proc, wrk_p_patch_g, cell_owner)
#ifdef HAVE_METIS
        ELSE IF(division_method(patch_no) == div_metis) THEN
          CALL divide_subset_metis(flag_c, n_proc, wrk_p_patch_g, cell_owner)
#endif
        ELSE
          CALL finish('divide_patch','Illegal division_method setting')
        ENDIF

        DEALLOCATE(flag_c)

      ENDIF

    ENDIF ! division_method

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
    INTEGER, POINTER  :: cell_owner(:)   !> Ownership of cells in p_patch_g
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
  !! order_type_of_halos 1= a "normal" grid for prognostic computations is processed,
  !!                     0= a local parent grid is processed
  !!                       in the first case, a finer distinction between different halo
  !!                       cell/edge/vertex levels is made to optimize communication
  !!                     2=move all halos to the end (for ocean)


  !!
  !! On exit, the entries of wrk_p_patch are set.
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !! Changed for usage for parent patch division, Rainer Johanni, Oct 2010

  SUBROUTINE divide_patch(wrk_p_patch, wrk_p_patch_g, cell_owner, n_boundary_rows, order_type_of_halos, my_proc)

    TYPE(t_patch), INTENT(INOUT) :: wrk_p_patch ! output patch, designated as INOUT because
                                                ! a few attributes are already set
    TYPE(t_patch), INTENT(in) :: wrk_p_patch_g

    INTEGER, POINTER :: cell_owner(:)
    INTEGER, INTENT(IN) :: n_boundary_rows
    INTEGER, INTENT(IN) :: order_type_of_halos
    INTEGER, INTENT(IN) :: my_proc

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_setup_subdivision::divide_patch'
    CHARACTER(len=132) :: msg

    INTEGER :: n, i, j, jv, je, jl, jb, jl_g, jb_g, jl_e, jb_e, jl_v, jb_v, ilev, iown,   &
               jl_c, jb_c, jc, irlev, ilev1, ilev_st, &
               jg, i_nchdom, irl0, mod_iown, k_c, k_e, k_v, &
               a_iown(2), a_mod_iown(2)

    INTEGER, ALLOCATABLE :: flag_c(:), flag_v(:)
    INTEGER, ALLOCATABLE :: flag2_c(:), flag2_e(:), flag2_v(:)
    LOGICAL :: is_outer, swap

    TYPE(nb_flag_list_elem), DIMENSION(-1:n_boundary_rows) :: &
         flag_c_list, flag_v_list, flag_e_list
    TYPE(nb_flag_list_elem) :: flag2_c_list(0:2*n_boundary_rows), &
         flag2_v_list(0:n_boundary_rows+1), flag2_e_list(0:2*n_boundary_rows+1)
    INTEGER, DIMENSION(-1:n_boundary_rows) :: n_ilev_c, n_ilev_v, n_ilev_e
    INTEGER :: n2_ilev_c(0:2*n_boundary_rows), n2_ilev_v(0:n_boundary_rows+1), &
         n2_ilev_e(0:2*n_boundary_rows+1), a_idx(2), loc_index_fill
    LOGICAL, ALLOCATABLE :: promote(:)
    INTEGER :: t_cells(1:wrk_p_patch_g%verts%max_connectivity), &
         t_cell_owner(1:wrk_p_patch_g%verts%max_connectivity)
    INTEGER, ALLOCATABLE :: owned_edges(:), owned_verts(:)
    INTEGER :: ierror

    IF (msg_level >= 10)  CALL message(routine, 'dividing patch')

    !---------------------------------------------------------------------
    ! flag_c_list(-1)%idx empty dummy list
    ! flag_c_list(0)%idx  all cells jg where cell_owner(jg) == my_proc
    ! flag_c_list(j)%idx j in 1..n_boundary_cells all cells bordering
    !                    on cells in flag_c_list(j - 1)
    ! flag_e_list(-1)%idx empty dummy list
    ! flag_e_list(j)%idx all edges where an adjacent cell is in
    !                    flag_c_list(j)%idx for minimal j,
    !                    i.e. no other cell in list flag_c_list(k)%idx
    !                    with k < j is adjacent to the edge
    ! flag_v_list(-1)%idx empty dummy list
    ! flag_v_list(j)%idx all vertices where an adjacent cell is in
    !                    flag_c_list(j)%idx for minimal j,
    !                    i.e. no other cell in list flag_c_list(k)%idx
    !                    with k < j is adjacent to the vertex
    !---------------------------------------------------------------------

    ALLOCATE(flag_c(wrk_p_patch_g%n_patch_cells))
    ALLOCATE(flag_v(wrk_p_patch_g%n_patch_verts))

    flag_c(:) = -1
    flag_v(:) = -1

    n2_ilev_c(:) = 0
    n2_ilev_e(:) = 0
    n2_ilev_v(:) = 0

    ! flag inner cells
    WHERE(cell_owner(:)==my_proc) flag_c(:) = 0
    n_ilev_c(0) = COUNT(cell_owner(:)==my_proc)
    n2_ilev_c(0) = n_ilev_c(0)
    n_ilev_c(1:n_boundary_rows) = 0
    ALLOCATE(flag_c_list(-1)%idx(0), flag_v_list(-1)%idx(0), &
         flag_e_list(-1)%idx(0))
    n_ilev_c(-1) = 0
    n_ilev_e(-1) = 0
    n_ilev_v(-1) = 0
    n = wrk_p_patch_g%n_patch_cells
    ALLOCATE(flag_c_list(0)%idx(n), flag2_c_list(0)%idx(n))
    j = 0
    DO i = 1, n
      IF (cell_owner(i) == my_proc) THEN
        j = j + 1
        flag_c_list(0)%idx(j) = i
        flag2_c_list(0)%idx(j) = i
      END IF
    END DO

    !-----------------------------------------------------------------------------------------------
    ! The purpose of the second set of flags is to control moving the halo points
    ! to the end of the index vector for nearly all cells even if order_type_of_halos = 1
    ! flag2_c/e/v = 0 (-1) wherever flag_c/e/v = 0 (-1)
    ! flag2_c = 2*flag_c-1 if the cell has a common edge with a cell having flag_c-1
    ! flag2_c = 2*flag_c if the cell has no common edge with a cell having flag_c-1
    ! flag2_e/v = 1 for boundary edges/vertices not owned by the current PE,
    !             otherwise flag2_v = flag_v+1
    ! flag2_e is the sum of flag_c of the two neighboring cells and 2*flag_c+1 at outer boundary edges
    !-----------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! flag2_e_list(0)%idx all edges which are owned by the local process
    ! flag2_e_list(1)%idx all edges adjacent to cells in flag_c_list(0)
    !                     but owned by other tasks
    ! flag2_e_list(2*i)%idx where i > 0, edges in idx are ONLY adjacent
    !                       to cells in flag_c_list(i)%idx
    ! flag2_e_list(2*i+1)%idx where i > 0, edges adjacent to one cell in
    !                         flag_c_list(i) and either
    !                            one cell in flag_c_list(i+1)
    !                         or edge is an outer edge,
    !                            i.e. adjacent to only one cell
    ! flag2_c_list(0)%idx == flag_c_list(0)%idx
    ! flag2_c_list(2*j-1)%idx where j > 0, contains all cells from
    !                         flag_c_list(j) which are
    !                         adjacent to cells in flag_c_list(j-1)
    ! flag2_c_list(2*j)%idx   where j > 0, contains all cells from
    !                         flag_c_list(j) which are NOT
    !                         adjacent to cells in flag_c_list(j-1)
    ! flag2_v_list(0)%idx contains all vertices from flag_v_list(0)%idx
    !                         owned by the current task
    ! flag2_v_list(1)%idx contains all vertices from flag_v_list(0)%idx
    !                     NOT owned by the current task
    ! flag2_v_list(j)%idx for j > 1 == flag_v_list(j-1)%idx
    !------------------------------------------------------------------------

    ALLOCATE(flag2_c(wrk_p_patch_g%n_patch_cells))
    ALLOCATE(flag2_e(wrk_p_patch_g%n_patch_edges))
    ALLOCATE(flag2_v(wrk_p_patch_g%n_patch_verts))

    flag2_c(:) = flag_c(:)
    flag2_e(:) = -1
    flag2_v(:) = -1

    i_nchdom = MAX(1,wrk_p_patch_g%n_childdom)

    ! find inner edges/verts and ghost cells/edges/verts

    DO ilev = 0, n_boundary_rows

      ! Flag cells belonging to this level.
      ! Cells belonging to the kernel (ilev==0) are already flagged.
      ! The patch always needs a complete halo, even for local parents
      ! where some of the boundary cells have no global owner.

      IF(ilev>0) THEN

#ifdef __SX__
        DO i = 1, wrk_p_patch_g%cell_type
          DO j = 1, wrk_p_patch_g%n_patch_cells

            IF(flag_c(j)<0) THEN

              jb = blk_no(j) ! block index
              jl = idx_no(j) ! line index

              IF (i > wrk_p_patch_g%cells%num_edges(jl,jb)) CYCLE

              ! Check if any vertex of this cell is already flagged.
              ! If this is the case, this cell goes to level ilev

              jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
              jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
              jv = idx_1d(jl_v, jb_v)

              IF(flag_v(jv)>=0) flag_c(j) = ilev
            ENDIF
          ENDDO
        ENDDO
#else
        DO j = 1, wrk_p_patch_g%n_patch_cells

          IF(flag_c(j)<0) THEN

            jb = blk_no(j) ! block index
            jl = idx_no(j) ! line index

            ! Check if any vertex of this cell is already flagged.
            ! If this is the case, this cell goes to level ilev

            DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
              jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
              jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
              jv = idx_1d(jl_v, jb_v)

              IF (flag_v(jv)>=0)  flag_c(j) = ilev
              ! the last line may be replaced by the following:
              !
              ! do not set "flag_c" for local_parent_patches, if the
              ! refin_ctrl flag indicates that we are outside of the
              ! nesting region.
              !  IF (flag_v(jv)>=0) THEN
              !    IF ( (jg > 1) .AND. .NOT. lcompute_grid ) THEN
              !      IF (wrk_p_patch_g%cells%refin_ctrl(jl,jb) < 0)  flag_c(j) = ilev
              !    ELSE
              !      flag_c(j) = ilev
              !    END IF
              !  END IF
            ENDDO
          ENDIF
        ENDDO
#endif
        n_ilev_c(ilev) = 2 * n_ilev_e(ilev - 1)
        ALLOCATE(flag_c_list(ilev)%idx(n_ilev_c(ilev)))
        k_c = 0
        k_v = n_ilev_v(ilev - 1)
        DO j = 1, k_v
          jv = flag_v_list(ilev - 1)%idx(j)
          jl_v = idx_no(jv)
          jb_v = blk_no(jv)
          DO i = 1, wrk_p_patch_g%verts%num_edges(jl_v, jb_v)
            jc = idx_1d(wrk_p_patch_g%verts%cell_idx(jl_v, jb_v, i), &
                 wrk_p_patch_g%verts%cell_blk(jl_v, jb_v, i))
            IF (.NOT. is_in_lists(flag_c_list(ilev-1:ilev), &
                 n_ilev_c(ilev-1), k_c, jc)) THEN
              k_c = k_c + 1
              flag_c_list(ilev)%idx(k_c) = jc
            END IF
          END DO
        END DO
        n_ilev_c(ilev) = k_c
      ENDIF

      ! Flag edges/verts belongig to this level.
      ! An edge/vert is flagged in this level if it belongs to a cell in this level
      ! and is not yet flagged (in which case it already belongs to a lower level).

      DO j = 1, wrk_p_patch_g%n_patch_cells

        IF(flag_c(j)==ilev) THEN

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index

          DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
            jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
            jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
            jv = idx_1d(jl_v, jb_v)
            IF(flag_v(jv)<0) flag_v(jv) = ilev
          ENDDO
        ENDIF
      ENDDO
      n = n_ilev_c(ilev)
      n_ilev_v(ilev) = n * wrk_p_patch_g%cell_type
      n_ilev_e(ilev) = n * wrk_p_patch_g%cell_type
      ALLOCATE(flag_v_list(ilev)%idx(n_ilev_v(ilev)), &
           flag_e_list(ilev)%idx(n_ilev_e(ilev)))
      k_e = 0
      k_v = 0
      DO j = 1, n
        jb = blk_no(flag_c_list(ilev)%idx(j))
        jl = idx_no(flag_c_list(ilev)%idx(j))
        DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
          jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
          jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
          je = idx_1d(jl_e, jb_e)
          IF (.NOT. is_in_lists(flag_e_list(ilev-1:ilev), &
               &                n_ilev_e(ilev-1), k_e, je)) THEN
            k_e = k_e + 1
            flag_e_list(ilev)%idx(k_e) = je
          END IF
          jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
          jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
          jv = idx_1d(jl_v, jb_v)
          IF (.NOT. is_in_lists(flag_v_list(ilev-1:ilev), &
               &                n_ilev_v(ilev-1), k_v, jv)) THEN
            k_v = k_v + 1
            flag_v_list(ilev)%idx(k_v) = jv
          END IF
        ENDDO
      END DO
      n_ilev_e(ilev) = k_e
      n_ilev_v(ilev) = k_v
      ! ilev single pass application of flag2_[ev] in
      ! build_patch_start_end requires them to be sorted
      CALL insertion_sort(flag_c_list(ilev)%idx(1:n_ilev_c(ilev)))
      CALL insertion_sort(flag_e_list(ilev)%idx(1:n_ilev_e(ilev)))
      CALL insertion_sort(flag_v_list(ilev)%idx(1:n_ilev_v(ilev)))
    ENDDO

    ! Now compute second set of flags
    DO ilev = 1, n_boundary_rows

      DO j = 1, wrk_p_patch_g%n_patch_cells

        IF (flag_c(j) == ilev) THEN

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index

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

      ALLOCATE(flag2_c_list(2*ilev)%idx(n_ilev_c(ilev)), &
           flag2_c_list(2*ilev - 1)%idx(n_ilev_c(ilev)))
      DO j = 1, n_ilev_c(ilev)
        jg = 2 * ilev
        jb = blk_no(flag_c_list(ilev)%idx(j)) ! block index
        jl = idx_no(flag_c_list(ilev)%idx(j)) ! line index
        DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
          jl_c = wrk_p_patch_g%cells%neighbor_idx(jl,jb,i)
          jb_c = wrk_p_patch_g%cells%neighbor_blk(jl,jb,i)
          jc = idx_1d(jl_c,jb_c)
          IF (jc < 1 .OR. jc > wrk_p_patch_g%n_patch_cells) CYCLE
          IF (ANY(flag_c_list(ilev - 1)%idx(1:n_ilev_c(ilev-1)) == jc)) THEN
            jg = jg - 1
            EXIT
          END IF
        ENDDO
        n2_ilev_c(jg) = n2_ilev_c(jg) + 1
        flag2_c_list(jg)%idx(n2_ilev_c(jg)) = flag_c_list(ilev)%idx(j)
      END DO

    ENDDO

    ! Preset edges and vertices bordering to interior cells with 0
    DO j = 1, wrk_p_patch_g%n_patch_cells

      IF (flag_c(j)==0) THEN

        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index

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

    ALLOCATE(promote(MAX(MAXVAL(n_ilev_e), n_ilev_v(0))))
    promote(1:n_ilev_e(0)) = .FALSE.
    k_e = 0
    DO j = 1, n_ilev_e(0)
      je = flag_e_list(0)%idx(j)
      jl_e = idx_no(je)
      jb_e = blk_no(je)
      a_idx(1:2) = idx_1d(wrk_p_patch_g%edges%cell_idx(jl_e, jb_e, 1:2), &
           wrk_p_patch_g%edges%cell_blk(jl_e, jb_e, 1:2))
      ! outer boundary edges always belong to single adjacent cell owner
      is_outer = ANY(a_idx == 0) .OR. a_idx(1) == a_idx(2)
      a_idx = MERGE(a_idx, 1, a_idx /= 0)
      is_outer = is_outer .OR. ANY(cell_owner(a_idx(1:2)) < 0)
      a_iown(1:2) = cell_owner(a_idx(1:2))
      a_mod_iown(1:2) = MOD(a_iown(1:2), 2)
      ! 50% chance of aquiring an edge, unless inner to flag_c == 0 or is_outer
      promote(j) = .NOT. is_outer .AND. a_iown(1) /= a_iown(2) &
           .AND. (MAXVAL(a_iown) == my_proc &
           &      .NEQV. a_mod_iown(1) == a_mod_iown(2))
      k_e = k_e + MERGE(1, 0, promote(j))
    END DO
    IF (n_boundary_rows > 0 .AND. order_type_of_halos /= 0) THEN
      n2_ilev_e(1) = k_e
      n2_ilev_e(0) = n_ilev_e(0) - n2_ilev_e(1)
      ALLOCATE(flag2_e_list(1)%idx(n2_ilev_e(1)))
      flag2_e_list(1)%idx &
           = PACK(flag_e_list(0)%idx(1:n_ilev_e(0)), &
           &      promote(1:n_ilev_e(0)))
      ALLOCATE(flag2_e_list(0)%idx(n2_ilev_e(0)))
      flag2_e_list(0)%idx &
           = PACK(flag_e_list(0)%idx(1:n_ilev_e(0)), &
           &      .NOT. promote(1:n_ilev_e(0)))
      ALLOCATE(owned_edges(1:n_ilev_e(0) - k_e))
      owned_edges = flag2_e_list(0)%idx
    ELSE
      n2_ilev_e(0) = n_ilev_e(0)
      ALLOCATE(flag2_e_list(0)%idx(n2_ilev_e(0)))
      flag2_e_list(0)%idx = flag_e_list(0)%idx(1:n_ilev_e(0))
      IF (n_boundary_rows > 0) THEN
        ALLOCATE(flag2_e_list(1)%idx(0))
        n2_ilev_e(1) = 0
        ALLOCATE(owned_edges(1:n_ilev_e(0) - k_e))
        owned_edges &
             = PACK(flag_e_list(0)%idx(1:n_ilev_e(0)), &
             &      .NOT. promote(1:n_ilev_e(0)))
      ELSE
        ALLOCATE(owned_edges(1:n2_ilev_e(0)))
        owned_edges = flag2_e_list(0)%idx
      END IF
    END IF

    DO ilev = 1, n_boundary_rows
      DO j = 1, wrk_p_patch_g%n_patch_cells

        IF (flag_c(j)==ilev) THEN

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index

          DO i = 1, wrk_p_patch_g%cells%num_edges(jl,jb)
            jl_e = wrk_p_patch_g%cells%edge_idx(jl,jb,i)
            jb_e = wrk_p_patch_g%cells%edge_blk(jl,jb,i)
            je = idx_1d(jl_e,jb_e)
            jl_c = wrk_p_patch_g%cells%neighbor_idx(jl,jb,i)
            jb_c = wrk_p_patch_g%cells%neighbor_blk(jl,jb,i)
            jc = idx_1d(jl_c,jb_c)
            IF (jc < 1 .OR. jc > wrk_p_patch_g%n_patch_cells) THEN
              flag2_e(je) = 2*ilev+1
            ELSE IF (flag_c(jc) == -1) THEN
              flag2_e(je) = 2*ilev+1
            ELSE
              flag2_e(je) = ilev + flag_c(jc)
            ENDIF
            jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
            jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
            jv = idx_1d(jl_v,jb_v)
            flag2_v(jv) = flag_v(jv)+1
          ENDDO
        ENDIF

      ENDDO

    ENDDO

    IF (n_boundary_rows > 0) THEN
      promote(1:n_ilev_v(0)) = .FALSE.
      k_v = 0
      DO j = 1, n_ilev_v(0)
        jv = flag_v_list(0)%idx(j)
        jb_v = blk_no(jv) ! block index
        jl_v = idx_no(jv) ! line index
        DO i = 1, wrk_p_patch_g%verts%num_edges(jl_v, jb_v)
          jl = wrk_p_patch_g%verts%cell_idx(jl_v, jb_v, i)
          jb = wrk_p_patch_g%verts%cell_blk(jl_v, jb_v, i)
          jc = idx_1d(jl,jb)
          promote(j) = promote(j) .OR. ANY(flag_c_list(1)%idx == jc)
        ENDDO
        k_v = k_v + MERGE(1, 0, promote(j))
      END DO


      !===================================================================
      DO j = 1, n_ilev_v(0)
        IF (.NOT. promote(j)) CYCLE
        jv = flag_v_list(0)%idx(j)
        jb_v = blk_no(jv) ! block index
        jl_v = idx_no(jv) ! line index
        n = wrk_p_patch_g%verts%num_edges(jl_v, jb_v)
        DO i = 1, n
          jl = wrk_p_patch_g%verts%cell_idx(jl_v, jb_v, i)
          jb = wrk_p_patch_g%verts%cell_blk(jl_v, jb_v, i)
          jc = idx_1d(jl,jb)
          t_cells(i) = jc
        END DO
        CALL insertion_sort(t_cells(1:n))
        t_cell_owner(1:n) = cell_owner(t_cells(1:n))
        jc = COUNT(t_cell_owner >= 0)
        t_cell_owner(1:jc) = PACK(t_cell_owner(1:n), t_cell_owner(1:n) >= 0)
        n = jc
        a_iown(1) = t_cell_owner(1)
        a_mod_iown(1) = MOD(a_iown(1), 2)
        DO i = 2, n
          a_iown(2) = t_cell_owner(i)
          a_mod_iown(2) = MOD(a_iown(2), 2)
          ! 50% chance of aquiring a vertex
          swap = a_iown(1) > a_iown(2) &
               .NEQV. a_mod_iown(1) == a_mod_iown(2)
          a_iown(1) = MERGE(a_iown(2), a_iown(1), swap)
          a_mod_iown(1) = MERGE(a_mod_iown(2), a_mod_iown(1), swap)
        END DO
        promote(j) = a_iown(1) /= my_proc
        k_v = k_v - MERGE(0, 1, promote(j))
      END DO

      ALLOCATE(owned_verts(1:n_ilev_v(0) - k_v))
      owned_verts(1:n_ilev_v(0) - k_v) &
           = PACK(flag_v_list(0)%idx(1:n_ilev_v(0)), &
           .NOT. promote(1:n_ilev_v(0)))

      IF (order_type_of_halos == 0) THEN
        promote = .FALSE.
        k_v = 0
      END IF

      n2_ilev_v(1) = k_v
      n2_ilev_v(0) = n_ilev_v(0) - n2_ilev_v(1)
      ALLOCATE(flag2_v_list(1)%idx(n2_ilev_v(1)), &
           flag2_v_list(0)%idx(n2_ilev_v(0)))
      flag2_v_list(0)%idx(:) = PACK(flag_v_list(0)%idx(1:n_ilev_v(0)), &
           .NOT. promote(1:n_ilev_v(0)))
      flag2_v_list(1)%idx(:) = PACK(flag_v_list(0)%idx(1:n_ilev_v(0)), &
           promote(1:n_ilev_v(0)))

      DO ilev = 1, n_boundary_rows
        ALLOCATE(flag2_v_list(ilev + 1)%idx(n_ilev_v(ilev)))
        flag2_v_list(ilev + 1)%idx = flag_v_list(ilev)%idx(1:n_ilev_v(ilev))
        n2_ilev_v(ilev + 1) = n_ilev_v(ilev)
      END DO

      DO ilev = 1, n_boundary_rows
        promote(1:n_ilev_e(ilev)) = .FALSE.
        DO j = 1, n_ilev_e(ilev)
          je = flag_e_list(ilev)%idx(j)
          jl_e = idx_no(je)
          jb_e = blk_no(je)
          a_idx(1:2) = idx_1d(wrk_p_patch_g%edges%cell_idx(jl_e, jb_e, 1:2), &
               wrk_p_patch_g%edges%cell_blk(jl_e, jb_e, 1:2))
          promote(j) = .NOT. (ANY(a_idx(1) == flag_c_list(ilev)%idx(1:n_ilev_c(ilev))) &
               .AND. ANY(a_idx(2) == flag_c_list(ilev)%idx(1:n_ilev_c(ilev))))
        END DO
        n2_ilev_e(2*ilev) = COUNT(.NOT. promote(1:n_ilev_e(ilev)))
        n2_ilev_e(2*ilev+1) = n_ilev_e(ilev) - n2_ilev_e(2*ilev)
        ALLOCATE(flag2_e_list(2*ilev)%idx(n2_ilev_e(2*ilev)), &
             flag2_e_list(2*ilev+1)%idx(n2_ilev_e(2*ilev+1)))
        flag2_e_list(2*ilev)%idx &
             = PACK(flag_e_list(ilev)%idx(1:n_ilev_e(ilev)), &
             &      .NOT. promote(1:n_ilev_e(ilev)))
        flag2_e_list(2*ilev + 1)%idx &
             = PACK(flag_e_list(ilev)%idx(1:n_ilev_e(ilev)), &
             &      promote(1:n_ilev_e(ilev)))

      END DO
    ELSE
      n2_ilev_v(0) = n_ilev_v(0)
      n2_ilev_v(1) = 0
      ALLOCATE(flag2_v_list(1)%idx(n2_ilev_v(1)), &
           flag2_v_list(0)%idx(n2_ilev_v(0)), owned_verts(n_ilev_v(0)))
      flag2_v_list(0)%idx(:) = flag_v_list(0)%idx(1:n_ilev_v(0))
      owned_verts(:) = flag2_v_list(0)%idx(:)
    END IF

    !-----------------------------------------------------------------------------------------------
    ! Get the number of cells/edges/verts and other data for patch allocation
    !-----------------------------------------------------------------------------------------------

    CALL prepare_patch(wrk_p_patch_g, wrk_p_patch, &
         SUM(n_ilev_c(:)), SUM(n_ilev_e(:)), SUM(n_ilev_v(:)))
    !-----------------------------------------------------------------------------------------------
    ! Set the global ownership for cells, edges and verts (needed for boundary exchange).
    ! Please note that opposed to cells, the global owner for edges/verts is
    ! not unique at the boundaries of the divided patch.
    ! To minimize load imbalance at runtime, we set the owner to the PE with the highest processor number
    ! which is participating at this edge/vert if both PE numbers are even or odd, otherwise
    ! the lowest processor number is chosen
    !-----------------------------------------------------------------------------------------------
    wrk_p_patch%cells%owner_g(:) = -1
    wrk_p_patch%edges%owner_g(:) = -1
    wrk_p_patch%verts%owner_g(:) = -1
    wrk_p_patch%cells%owner_g(flag_c_list(0)%idx(1:n_ilev_c(0))) = my_proc
    wrk_p_patch%edges%owner_g(owned_edges) = my_proc
    wrk_p_patch%verts%owner_g(owned_verts) = my_proc
#ifndef NOMPI
    CALL mpi_allreduce(mpi_in_place, wrk_p_patch%cells%owner_g, &
         SIZE(wrk_p_patch%cells%owner_g), &
         p_int, mpi_max, p_comm_work, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
    CALL mpi_allreduce(mpi_in_place, wrk_p_patch%edges%owner_g, &
         SIZE(wrk_p_patch%edges%owner_g), &
         p_int, mpi_max, p_comm_work, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
    CALL mpi_allreduce(mpi_in_place, wrk_p_patch%verts%owner_g, &
         SIZE(wrk_p_patch%verts%owner_g), &
         p_int, mpi_max, p_comm_work, ierror)
    IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
#endif


    IF (SUM(n_ilev_c(:)) /= wrk_p_patch%n_patch_cells) &
      CALL finish(routine, 'number of cells mismatch')
    IF (SUM(n_ilev_e(:)) /= wrk_p_patch%n_patch_edges) &
      CALL finish(routine, 'number of edges mismatch')
    IF (SUM(n_ilev_v(:)) /= wrk_p_patch%n_patch_verts) &
      CALL finish(routine, 'number of verts mismatch')

    DO ilev = 0, n_boundary_rows
      IF (.NOT. (ALL(flag_c(flag_c_list(ilev)%idx(1:n_ilev_c(ilev))) == ilev)) &
           .OR. COUNT(flag_c == ilev) /= n_ilev_c(ilev)) THEN
        WRITE (msg, '(a,i0)') 'cells ilev mismatch at ', ilev
        CALL finish(routine, msg)
      END IF
      IF (.NOT. (ALL(flag_v(flag_v_list(ilev)%idx(1:n_ilev_v(ilev))) == ilev)) &
           .OR. COUNT(flag_v == ilev) /= n_ilev_v(ilev)) THEN
        WRITE (msg, '(a,i0)') 'vertices ilev mismatch at ', ilev
        CALL finish(routine, msg)
      END IF
    END DO

    DEALLOCATE(flag_c, flag_v)

    !-----------------------------------------------------------------------------------------------
    ! Get the indices of local cells/edges/verts within global patch and vice versa.

    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_cells, &
         n_patch_cve_g = wrk_p_patch_g%n_patch_cells, &
         max_childdom = wrk_p_patch%max_childdom, &
         nchilddom = MAX(1,wrk_p_patch%n_childdom), &
         patch_id = wrk_p_patch_g%id, &
         nblks = wrk_p_patch%nblks_c, &
         npromz = wrk_p_patch%npromz_c, &
         cell_type = wrk_p_patch%cell_type, &
         min_rlcve = min_rlcell, &
         min_rlcve_int = min_rlcell_int, &
         max_rlcve = max_rlcell, &
         max_ilev = 2*n_boundary_rows, &
         max_hw_cve = 2*max_hw, &
         flag2_list = flag2_c_list, &
         n2_ilev = n2_ilev_c, &
         glb_index = wrk_p_patch%cells%glb_index, &
         loc_index = wrk_p_patch%cells%loc_index, &
         start_idx = wrk_p_patch%cells%start_idx, &
         start_blk = wrk_p_patch%cells%start_blk, &
         end_idx = wrk_p_patch%cells%end_idx, &
         end_blk = wrk_p_patch%cells%end_blk, &
         start_idx_g = wrk_p_patch_g%cells%start_idx, &
         start_blk_g = wrk_p_patch_g%cells%start_blk, &
         end_idx_g = wrk_p_patch_g%cells%end_idx, &
         end_blk_g = wrk_p_patch_g%cells%end_blk, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .TRUE., &
         refin_ctrl = wrk_p_patch_g%cells%refin_ctrl, &
         refinement_predicate = refine_cells)
    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_edges, &
         n_patch_cve_g = wrk_p_patch_g%n_patch_edges, &
         max_childdom = wrk_p_patch%max_childdom, &
         nchilddom = MAX(1,wrk_p_patch%n_childdom), &
         patch_id = wrk_p_patch_g%id, &
         nblks = wrk_p_patch%nblks_e, &
         npromz = wrk_p_patch%npromz_e, &
         cell_type = wrk_p_patch%cell_type, &
         min_rlcve = min_rledge, &
         min_rlcve_int = min_rledge_int, &
         max_rlcve = max_rledge, &
         max_ilev = 2*n_boundary_rows+1, &
         max_hw_cve = 2*max_hw + 1, &
         flag2_list = flag2_e_list, &
         n2_ilev = n2_ilev_e, &
         glb_index = wrk_p_patch%edges%glb_index, &
         loc_index = wrk_p_patch%edges%loc_index, &
         start_idx = wrk_p_patch%edges%start_idx, &
         start_blk = wrk_p_patch%edges%start_blk, &
         end_idx = wrk_p_patch%edges%end_idx, &
         end_blk = wrk_p_patch%edges%end_blk, &
         start_idx_g = wrk_p_patch_g%edges%start_idx, &
         start_blk_g = wrk_p_patch_g%edges%start_blk, &
         end_idx_g = wrk_p_patch_g%edges%end_idx, &
         end_blk_g = wrk_p_patch_g%edges%end_blk, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .FALSE., &
         refin_ctrl = wrk_p_patch_g%edges%refin_ctrl, &
         refinement_predicate = refine_edges)
    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_verts, &
         n_patch_cve_g = wrk_p_patch_g%n_patch_verts, &
         max_childdom = wrk_p_patch%max_childdom, &
         nchilddom = MAX(1,wrk_p_patch%n_childdom), &
         patch_id = wrk_p_patch_g%id, &
         nblks = wrk_p_patch%nblks_v, &
         npromz = wrk_p_patch%npromz_v, &
         cell_type = wrk_p_patch%cell_type, &
         min_rlcve = min_rlvert, &
         min_rlcve_int = min_rlvert_int, &
         max_rlcve = max_rlvert, &
         max_ilev = n_boundary_rows+1, &
         max_hw_cve = max_hw + 1, &
         flag2_list = flag2_v_list, &
         n2_ilev = n2_ilev_v, &
         glb_index = wrk_p_patch%verts%glb_index, &
         loc_index = wrk_p_patch%verts%loc_index, &
         start_idx = wrk_p_patch%verts%start_idx, &
         start_blk = wrk_p_patch%verts%start_blk, &
         end_idx = wrk_p_patch%verts%end_idx, &
         end_blk = wrk_p_patch%verts%end_blk, &
         start_idx_g = wrk_p_patch_g%verts%start_idx, &
         start_blk_g = wrk_p_patch_g%verts%start_blk, &
         end_idx_g = wrk_p_patch_g%verts%end_idx, &
         end_blk_g = wrk_p_patch_g%verts%end_blk, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .FALSE., &
         refin_ctrl = wrk_p_patch_g%verts%refin_ctrl, &
         refinement_predicate = refine_verts)

    ! Sanity checks: are there still elements of the index lists filled with dummy values?
    ! Note: the use of write(0,*) instead of CALL message is intended here
    !       because the debug output is wanted for all PEs where an error occurs
    IF (ANY(wrk_p_patch%cells%start_idx(:,:)<0) .OR. ANY(wrk_p_patch%cells%start_blk(:,:)<0) .OR.&
        ANY(wrk_p_patch%cells%end_idx(:,:)  <0) .OR. ANY(wrk_p_patch%cells%end_blk(:,:)  <0)) THEN
      DO j = 1, MAX(1,wrk_p_patch%n_childdom)
        DO i = min_rlcell, max_rlcell
          WRITE(0,'(a,2i5,2i4,4i7)') 'cells',my_proc,wrk_p_patch%id,i,j,     &
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

      DO i=1,wrk_p_patch%cell_type
!CDIR IEXPAND
        CALL get_local_index(wrk_p_patch%cells%loc_index, &
          & wrk_p_patch_g%cells%neighbor_idx(jl_g,jb_g,i),&
          & wrk_p_patch_g%cells%neighbor_blk(jl_g,jb_g,i),&
          & wrk_p_patch%cells%neighbor_idx(jl,jb,i),      &
          & wrk_p_patch%cells%neighbor_blk(jl,jb,i))

!CDIR IEXPAND
        CALL get_local_index(wrk_p_patch%edges%loc_index, &
          & wrk_p_patch_g%cells%edge_idx(jl_g,jb_g,i),    &
          & wrk_p_patch_g%cells%edge_blk(jl_g,jb_g,i),    &
          & wrk_p_patch%cells%edge_idx(jl,jb,i),          &
          & wrk_p_patch%cells%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL get_local_index(wrk_p_patch%verts%loc_index, &
          & wrk_p_patch_g%cells%vertex_idx(jl_g,jb_g,i),  &
          & wrk_p_patch_g%cells%vertex_blk(jl_g,jb_g,i),  &
          & wrk_p_patch%cells%vertex_idx(jl,jb,i),        &
          & wrk_p_patch%cells%vertex_blk(jl,jb,i))

      ENDDO
    ENDDO

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
      wrk_p_patch%cells%child_idx(jl,jb,1:4) = wrk_p_patch_g%cells%child_idx(jl_g,jb_g,1:4)
      wrk_p_patch%cells%child_blk(jl,jb,1:4) = wrk_p_patch_g%cells%child_blk(jl_g,jb_g,1:4)
      wrk_p_patch%cells%child_id (jl,jb)   = wrk_p_patch_g%cells%child_id(jl_g,jb_g)

      ! Safety check only: edge_idx and vertex_idx must always be valid!
 !     if(any(wrk_p_patch%cells%edge_idx(jl,jb,:) <= 0) .or. &
 !      & any(wrk_p_patch%cells%edge_blk(jl,jb,:) <= 0) )    &
 !      & CALL finish('divide_patch','Illegal value for patch%cells%edge_idx')
 !     if(any(wrk_p_patch%cells%vertex_idx(jl,jb,:) <= 0) .or. &
 !      & any(wrk_p_patch%cells%vertex_blk(jl,jb,:) <= 0) )    &
 !      & CALL finish('divide_patch','Illegal value for patch%cells%vertex_idx')

      wrk_p_patch%cells%num_edges(jl,jb)          = wrk_p_patch_g%cells%num_edges(jl_g,jb_g)
      wrk_p_patch%cells%center(jl,jb)%lat         = wrk_p_patch_g%cells%center(jl_g,jb_g)%lat
      wrk_p_patch%cells%center(jl,jb)%lon         = wrk_p_patch_g%cells%center(jl_g,jb_g)%lon
      wrk_p_patch%cells%refin_ctrl(jl,jb)         = wrk_p_patch_g%cells%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%cells%child_id(jl,jb)           = wrk_p_patch_g%cells%child_id(jl_g,jb_g)
    ENDDO

    DO j = 0, 2 * n_boundary_rows
      DO i = 1, n2_ilev_c(j)
        jg = flag2_c_list(j)%idx(i)
        jc = wrk_p_patch%cells%loc_index(jg)
        jb = blk_no(jc)
        jl = idx_no(jc)
        wrk_p_patch%cells%decomp_domain(jl,jb) = j
      END DO
    END DO

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
      wrk_p_patch%edges%child_idx(jl,jb,1:4) = wrk_p_patch_g%edges%child_idx(jl_g,jb_g,1:4)
      wrk_p_patch%edges%child_blk(jl,jb,1:4) = wrk_p_patch_g%edges%child_blk(jl_g,jb_g,1:4)
      wrk_p_patch%edges%child_id (jl,jb)     = wrk_p_patch_g%edges%child_id(jl_g,jb_g)

      wrk_p_patch%edges%refin_ctrl(jl,jb)    = wrk_p_patch_g%edges%refin_ctrl(jl_g,jb_g)
    ENDDO

    DO j = 0, 2 * n_boundary_rows + 1
      DO i = 1, n2_ilev_e(j)
        jg = flag2_e_list(j)%idx(i)
        je = wrk_p_patch%edges%loc_index(jg)
        jb = blk_no(je)
        jl = idx_no(je)
        wrk_p_patch%edges%decomp_domain(jl,jb) = j
      END DO
    END DO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%verts%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%verts%glb_index(j)) ! Line  index in global patch

      wrk_p_patch%verts%vertex(jl,jb)%lat    = wrk_p_patch_g%verts%vertex(jl_g,jb_g)%lat
      wrk_p_patch%verts%vertex(jl,jb)%lon    = wrk_p_patch_g%verts%vertex(jl_g,jb_g)%lon
      wrk_p_patch%verts%refin_ctrl(jl,jb)    = wrk_p_patch_g%verts%refin_ctrl(jl_g,jb_g)
    ENDDO

    DO j = 0, n_boundary_rows + 1
      DO i = 1, n2_ilev_v(j)
        jg = flag2_v_list(j)%idx(i)
        jv = wrk_p_patch%verts%loc_index(jg)
        jb = blk_no(jv)
        jl = idx_no(jv)
        wrk_p_patch%verts%decomp_domain(jl,jb) = j
      END DO
    END DO

    DEALLOCATE(flag2_c, flag2_e, flag2_v)

  CONTAINS

    SUBROUTINE prepare_patch(wrk_p_patch_g, wrk_p_patch, &
         n_patch_cells, n_patch_edges, n_patch_verts)
      !> output patch, designated as INOUT because
      !! a few attributes are already set
      TYPE(t_patch), INTENT(inout) :: wrk_p_patch
      TYPE(t_patch), INTENT(in) :: wrk_p_patch_g
      INTEGER, INTENT(in) :: n_patch_cells, n_patch_edges, n_patch_verts

      wrk_p_patch%n_patch_cells = n_patch_cells
      wrk_p_patch%n_patch_edges = n_patch_edges
      wrk_p_patch%n_patch_verts = n_patch_verts

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
      wrk_p_patch%npromz_c      = wrk_p_patch%n_patch_cells - (wrk_p_patch%nblks_c - 1)*nproma
      wrk_p_patch%alloc_cell_blocks = wrk_p_patch%nblks_c
      IF (use_dummy_cell_closure) THEN
        IF (wrk_p_patch%npromz_c == nproma) &
             wrk_p_patch%alloc_cell_blocks = wrk_p_patch%nblks_c + 1
      ENDIF

      ! ... for the edges
      wrk_p_patch%nblks_e       = blk_no(wrk_p_patch%n_patch_edges)
      wrk_p_patch%npromz_e      = wrk_p_patch%n_patch_edges - (wrk_p_patch%nblks_e - 1)*nproma

      ! ... for the vertices
      wrk_p_patch%nblks_v       = blk_no(wrk_p_patch%n_patch_verts)
      wrk_p_patch%npromz_v      = wrk_p_patch%n_patch_verts - (wrk_p_patch%nblks_v - 1)*nproma

      ! Also needed for patch allocation
      wrk_p_patch%max_childdom  = wrk_p_patch_g%max_childdom

      ! Set other scalar members of patch here too ..
      wrk_p_patch%grid_filename      = wrk_p_patch_g%grid_filename
      wrk_p_patch%level              = wrk_p_patch_g%level
      wrk_p_patch%id                 = wrk_p_patch_g%id
      wrk_p_patch%cells%max_connectivity = wrk_p_patch_g%cells%max_connectivity
      wrk_p_patch%parent_id          = wrk_p_patch_g%parent_id
      wrk_p_patch%parent_child_index = wrk_p_patch_g%parent_child_index
      wrk_p_patch%child_id(:)        = wrk_p_patch_g%child_id(:)
      wrk_p_patch%child_id_list(:)   = wrk_p_patch_g%child_id_list(:)
      wrk_p_patch%n_childdom         = wrk_p_patch_g%n_childdom
      wrk_p_patch%n_chd_total        = wrk_p_patch_g%n_chd_total
      wrk_p_patch%nlev               = wrk_p_patch_g%nlev
      wrk_p_patch%nlevp1             = wrk_p_patch_g%nlevp1
      wrk_p_patch%nshift             = wrk_p_patch_g%nshift
      wrk_p_patch%nshift_total       = wrk_p_patch_g%nshift_total
      wrk_p_patch%nshift_child       = wrk_p_patch_g%nshift_child
      wrk_p_patch%grid_uuid          = wrk_p_patch_g%grid_uuid

      !-----------------------------------------------------------------------------------------------
      ! Allocate the required data arrays in patch
      !-----------------------------------------------------------------------------------------------

      CALL allocate_basic_patch(wrk_p_patch)
      CALL allocate_remaining_patch(wrk_p_patch,2) ! 2 = only those needed for parallelization control


    END SUBROUTINE prepare_patch

    FUNCTION is_in_lists(list, lim1, lim2, idx) RESULT(p)
      TYPE(nb_flag_list_elem) :: list(1:2)
      INTEGER, INTENT(in) :: lim1, lim2, idx
      LOGICAL :: p
      p = ANY(list(1)%idx(1:lim1) == idx) .OR. ANY(list(2)%idx(1:lim2) == idx)
    END FUNCTION is_in_lists

  END SUBROUTINE divide_patch

  FUNCTION refine_cells(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag == 0 &
         .OR. ( flag >= 1 .AND. irl0 == 1 ) &
         .OR. ( flag == 1 .AND. irl0 == 2 )
  END FUNCTION refine_cells

  FUNCTION refine_edges(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag==0 &
         .OR. ((irl0==1 .OR. irl0==2) .AND. flag>=1) &
         .OR. ((irl0==3 .OR. irl0==4) .AND. (flag==1 .OR. flag==2)) &
         .OR. ((irl0==5 .OR. irl0==6 .OR. irl0<0) .AND. flag==1)
  END FUNCTION refine_edges

  FUNCTION refine_verts(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag == 0 &
         .OR.  (irl0 == 1                .AND. flag  > 0) &
         .OR. ((irl0 == 2 .OR. irl0 < 0) .AND. flag == 1)
  END FUNCTION refine_verts

  SUBROUTINE build_patch_start_end(n_patch_cve, n_patch_cve_g, &
       max_childdom, nchilddom, patch_id, nblks, npromz, cell_type, &
       min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
       flag2_list, n2_ilev, glb_index, loc_index, &
       start_idx, start_blk, end_idx, end_blk, &
       start_idx_g, start_blk_g, end_idx_g, end_blk_g, &
       order_type_of_halos, l_cell_correction, refin_ctrl, refinement_predicate)
    INTEGER, INTENT(in) :: n_patch_cve, n_patch_cve_g, max_childdom, &
         nchilddom, patch_id, nblks, npromz, cell_type, &
         min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve
    INTEGER, INTENT(in) :: order_type_of_halos
    LOGICAL, INTENT(in) :: l_cell_correction
    INTEGER, INTENT(in) :: refin_ctrl(:, :), n2_ilev(0:)
    TYPE(nb_flag_list_elem), INTENT(in) :: flag2_list(0:)
    INTEGER, INTENT(inout) :: glb_index(:), loc_index(:)
    INTEGER, DIMENSION(min_rlcve:, :), INTENT(inout) :: &
         start_idx, start_blk, end_idx, end_blk
    INTEGER, DIMENSION(min_rlcve:, :), INTENT(in) :: &
         start_idx_g, start_blk_g, end_idx_g, end_blk_g
    INTERFACE
      FUNCTION refinement_predicate(flag, irl0) RESULT(p)
        INTEGER, INTENT(in) :: flag, irl0
        LOGICAL :: p
      END FUNCTION refinement_predicate
    END INTERFACE

    CHARACTER(LEN=*), PARAMETER :: routine = 'mo_setup_subdivision::build_patch_start_end'

    INTEGER :: i, ilev, ilev1, ilev_st, irl0, irlev, j, jb, jf, jl, &
         exec_sequence, ref_flag, k, k_start, loc_index_fill

    INTEGER, ALLOCATABLE :: flag2_union(:,:)

    ALLOCATE(flag2_union(n_patch_cve, 2))
    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      glb_index(1:n2_ilev(0)) = flag2_list(0)%idx(1:n2_ilev(0))
      j = 1
      DO k = 1, n2_ilev(0)
        DO j = j, flag2_list(0)%idx(k) - 1
          loc_index(j) = -(k)
        END DO
        loc_index(j) = k
        j = j + 1
      END DO
      DO j = j, n_patch_cve_g
        loc_index(j) = -k
      END DO
      loc_index_fill = n2_ilev(0)

    CASE (1)
      k = 1
      DO ilev = 0, max_ilev
        flag2_union(k:k+n2_ilev(ilev)-1, 1) &
             = flag2_list(ilev)%idx(1:n2_ilev(ilev))
        flag2_union(k:k+n2_ilev(ilev)-1, 2) &
             = ilev
        k = k + n2_ilev(ilev)
      END DO
      CALL quicksort(flag2_union(:, 1), flag2_union(:, 2))
      k = 1
      jf = 1
      j = 0
      DO WHILE (j < n_patch_cve_g)
        j = j + 1
        DO WHILE (j < flag2_union(jf, 1))
          loc_index(j) = -k
          j = j + 1
        END DO
        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index
        irl0 = refin_ctrl(jl,jb)
        IF (refinement_predicate(flag2_union(jf, 2), irl0)) THEN
          loc_index(j) = k
          glb_index(k) = j
          k = k + 1
        ELSE
          loc_index(j) = -k
        END IF
        jf = jf + 1
        IF (jf > n_patch_cve) THEN
          DO WHILE (j < n_patch_cve_g)
            j = j + 1
            loc_index(j) = -k
          END DO
        END IF
      END DO
      loc_index_fill = k - 1

    CASE default
      CALL finish("", "Uknown order_type_of_halos")
    END SELECT

    ! Set start_idx/blk ... end_idx/blk for cells/verts/edges.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the loc_index array.

    DO j=1,max_childdom
      ! Note: the segments between min_rlcell and min_rlcell_int-1 are reserved for
      ! halo cells; they are set below
      DO i=min_rlcve_int,max_rlcve
!CDIR IEXPAND
        CALL get_local_index(loc_index, &
          & start_idx_g(i,j),           &
          & start_blk_g(i,j),           &
          & start_idx(i,j),             &
          & start_blk(i,j),             &
          & +1 )
!CDIR IEXPAND
        CALL get_local_index(loc_index, &
          & end_idx_g(i,j),             &
          & end_blk_g(i,j),             &
          & end_idx(i,j),               &
          & end_blk(i,j),               &
          & -1 )
      ENDDO
    ENDDO

    ! Preset remaining index sections with dummy values
    start_idx(min_rlcve:min_rlcve_int-1, :) = -9999
    start_blk(min_rlcve:min_rlcve_int-1, :) = -9999
    end_idx(min_rlcve:min_rlcve_int-1, :)   = -9999
    end_blk(min_rlcve:min_rlcve_int-1, :)   = -9999

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      ! processing local parent grid
      !
      ! Gather all halo cells/edges/verts at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlcve_int-1
      irlev = min_rlcve_int-1
      ref_flag = MERGE(0, 1, l_cell_correction .OR. order_type_of_halos == 2)

      j = 1
      DO ilev = ref_flag + 1, max_ilev
        flag2_union(j:j+n2_ilev(ilev)-1, 1) &
             = flag2_list(ilev)%idx(1:n2_ilev(ilev))
        j = j + n2_ilev(ilev)
      END DO
      CALL quicksort(flag2_union(1:j-1, 1))
      k_start = SUM(n2_ilev(0:ref_flag))
      IF (start_idx(irlev,1)==-9999 .AND. k_start + 1 <= n_patch_cve) THEN
        start_idx(irlev,:) = idx_no(k_start + 1) ! line index
        start_blk(irlev,:) = blk_no(k_start + 1) ! block index
      ENDIF

      DO k = k_start + 1, n_patch_cve
        j = flag2_union(k - k_start, 1)
        glb_index(k) = j
        loc_index(j) = k
      END DO

      ! Set end index when the last point is processed
      end_idx(min_rlcve:irlev,:) = idx_no(n_patch_cve)
      end_blk(min_rlcve:irlev,:) = blk_no(n_patch_cve)
      start_idx(min_rlcve:irlev-1,:) = start_idx(irlev,1)
      start_blk(min_rlcve:irlev-1,:) = start_blk(irlev,1)

    CASE(1)
      k = loc_index_fill
      ! Gather all halo cells/verts/edges except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2 value and will be placed at the
      ! index levels between min_rlcve_int-1 and min_rlcve
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here
      DO ilev = 1, max_ilev
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        DO j = 1, n2_ilev(ilev)
          jf = flag2_list(ilev)%idx(j)
          IF (loc_index(jf) < 0) THEN
            k = k + 1
            glb_index(k) = jf
            loc_index(jf) = k
            jb = blk_no(k) ! block index
            jl = idx_no(k) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (start_idx(irlev,1)==-9999) THEN
              start_idx(irlev,:) = jl
              start_blk(irlev,:) = jb
            ENDIF
            end_idx(irlev,:) = jl
            end_blk(irlev,:) = jb
          END IF
        END DO

        ! Just in case that no grid point is found (may happen for ilev=1)
        IF (.NOT. l_cell_correction .AND. start_idx(irlev,1)==-9999) THEN
          start_idx(irlev,:) = end_idx(irlev+1,nchilddom)+1
          start_blk(irlev,:) = end_blk(irlev+1,nchilddom)
          end_idx(irlev,:)   = end_idx(irlev+1,nchilddom)
          end_blk(irlev,:)   = end_blk(irlev+1,nchilddom)
        ENDIF

      ENDDO
      loc_index_fill = k

      ! Fill start and end indices for remaining index sections
      IF (patch_id == 0) THEN
        ilev1 = min_rlcve_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlcve, min_rlcve_int - max_ilev)
        ilev_st = max_ilev + 1
      ENDIF
      IF (l_cell_correction .AND. cell_type==6) THEN ! for hexagons, there are no even-order halo cells
        DO ilev = 2, max_hw_cve, 2
          irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
          start_idx(irlev,:) = end_idx(irlev+1,1) + 1
          start_blk(irlev,:) = end_blk(irlev+1,1)
          end_idx(irlev,:)   = end_idx(irlev+1,1)
          end_blk(irlev,:)   = end_blk(irlev+1,1)
        ENDDO
      ENDIF
      DO ilev = ilev_st,max_hw_cve
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        start_idx(irlev,:) = end_idx(ilev1,1) + 1
        start_blk(irlev,:) = end_blk(ilev1,1)
        end_idx(irlev,:)   = end_idx(ilev1,1)
        end_blk(irlev,:)   = end_blk(ilev1,1)
      ENDDO
    END SELECT

    ! exec_sequence only serves the purpose to retain the execution
    ! sequence for these three correction as it was in the previous
    ! code version, it might be irrelevant.
    ! see redmine issue #3924
    DO exec_sequence = 0, 1
      IF (.NOT. l_cell_correction .AND. exec_sequence == 0 &
           .OR. l_cell_correction .AND. exec_sequence == 1) THEN
        ! If a PE owns only nest boundary points, it may happen that
        ! one or more index sections of halo cells/verts/edges are
        ! empty. These are filled here
        IF (start_blk(0,1)     > end_blk(min_rlcve_int,nchilddom) &
             .OR. start_blk(0,1) == end_blk(min_rlcve_int,nchilddom) &
             .AND. start_idx(0,1) > end_idx(min_rlcve_int,nchilddom) ) THEN
          DO i = min_rlcve_int-1, min_rlcve, -1
            IF (start_idx(i,1) == -9999 .OR. &
                 start_blk(i,1) == -9999 .OR. &
                 end_idx(i,1)   == -9999 .OR. &
                 end_blk(i,1)   == -9999 ) THEN
              start_idx(i,:) = start_idx(i+1,1)
              start_blk(i,:) = start_blk(i+1,1)
              end_idx(i,:)   = end_idx(i+1,1)
              end_blk(i,:)   = end_blk(i+1,1)
            ENDIF
          ENDDO
        ENDIF

        ! Finally, fill start indices of halo rows with meaningful
        ! values for empty patches (which occur when processor
        ! splitting is applied)
        IF (nblks <= 0 .OR. npromz <= 0) THEN
          DO i=min_rlcve,min_rlcve_int-1
            start_idx(i,:) = start_idx(min_rlcve_int,1)
            start_blk(i,:) = start_blk(min_rlcve_int,1)
            end_idx(i,:)   = end_idx(min_rlcve_int,1)
            end_blk(i,:)   = end_blk(min_rlcve_int,1)
          ENDDO
        ENDIF
      END IF
      ! special treatment for sequential runs with trivial decomposition
      IF (exec_sequence == 0 .AND. SUM(n2_ilev(1:)) == 0) THEN
        end_idx(min_rlcve:min_rlcve_int-1,:) &
             = end_idx(min_rlcve_int, nchilddom)
        end_blk(min_rlcve:min_rlcve_int-1,:) &
             = end_blk(min_rlcve_int, nchilddom)
        IF (end_idx(min_rlcve_int, nchilddom) == nproma) THEN
          start_blk(min_rlcve:min_rlcve_int-1,:) = &
               &   end_blk(min_rlcve_int, nchilddom) + 1
          start_idx(min_rlcve:min_rlcve_int-1,:) = 1
        ELSE
          start_idx(min_rlcve:min_rlcve_int-1,:) = &
               &    end_idx(min_rlcve_int, nchilddom) + 1
          start_blk(min_rlcve:min_rlcve_int-1,:) = &
               &   end_blk(min_rlcve_int, nchilddom)
        END IF
      END IF

    END DO

  END SUBROUTINE build_patch_start_end

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
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE divide_subset_geometric(subset_flag, n_proc, wrk_divide_patch, &
                                     owner)

    INTEGER, INTENT(in)    :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    TYPE(t_patch), INTENT(in) :: wrk_divide_patch
    INTEGER, INTENT(out)   :: owner(:) ! receives the owner PE for every cell
    ! (-1 for cells not in subset)

    INTEGER :: i, ii, j, jl, jb, jn, jl_v, jb_v, nc, nn, npt, jd, idp, ncs, nce, jm(1)
    INTEGER :: count_physdom(max_phys_dom), count_total, id_physdom(max_phys_dom), &
               num_physdom, proc_count(max_phys_dom), proc_offset(max_phys_dom), checksum, &
               ncell_offset(0:max_phys_dom)
    REAL(wp), ALLOCATABLE :: cell_desc(:,:), workspace(:,:)
    REAL(wp) :: cclat, cclon, corr_ratio(max_phys_dom)
    LOGICAL  :: lsplit_merged_domains

    !-----------------------------------------------------------------------

    IF(p_pe_work==0) THEN
      IF(divide_for_radiation) THEN
        WRITE(0,*) 'divide_patch: Using geometric area subdivision for radiation'
      ELSE
        WRITE(0,*) 'divide_patch: Using geometric area subdivision (normal)'
      ENDIF
    ENDIF

    ! Initialization of cell owner field
    owner(:) = -1

    lsplit_merged_domains = .FALSE.

    ! Check if domain merging has been applied; for large PE numbers, calculating the DD
    ! for each physical domain separately tends to yield a more balanced distribution
    ! of the halo points
    IF (ldiv_phys_dom .AND. .NOT. divide_for_radiation) THEN
      count_physdom(:) = 0
      count_total      = 0
      num_physdom      = 0
      id_physdom(:)    = 0
      proc_count(:)    = 0
      ncell_offset(:)  = 0
      proc_offset(:)   = 0
      corr_ratio(:)    = 1._wp

      DO j = 1, wrk_divide_patch%n_patch_cells
        IF (subset_flag(j) > 0) THEN
          count_physdom(subset_flag(j)) = count_physdom(subset_flag(j)) + 1
          count_total = count_total + 1
        ENDIF
      ENDDO

      IF (MAXVAL(count_physdom) < count_total) THEN
        lsplit_merged_domains = .TRUE.
        DO j = 1, max_phys_dom
          IF (count_physdom(j) > 0) THEN
            num_physdom = num_physdom + 1
            id_physdom(num_physdom) = j
            proc_count(num_physdom) = NINT(REAL(count_physdom(j),wp)/REAL(count_total,wp)*REAL(n_proc,wp))
          ENDIF
        ENDDO

        ! Ensure that the sum of partial processor counts matches n_proc
        checksum = SUM(proc_count(1:num_physdom)) - n_proc
        IF (checksum /= 0) THEN
          DO j = 1, num_physdom
            corr_ratio(j) = REAL(proc_count(j),wp) / &
              (REAL(count_physdom(id_physdom(j)),wp)/REAL(count_total,wp)*REAL(n_proc,wp))
          ENDDO
        ENDIF
        IF (checksum > 0) THEN
          DO WHILE (checksum > 0)
            jm = MAXLOC(corr_ratio)
            j = jm(1)
            corr_ratio(j) = 1._wp
            proc_count(j) = proc_count(j) - 1
            checksum = checksum - 1
          ENDDO
        ELSE IF (checksum < 0) THEN
          DO WHILE (checksum < 0)
            jm = MINLOC(corr_ratio)
            j = jm(1)
            corr_ratio(j) = 1._wp
            proc_count(j) = proc_count(j) + 1
            checksum = checksum + 1
          ENDDO
        ENDIF

        ! Compute offset for processor IDs
        DO j = 2, num_physdom
          proc_offset(j) = proc_offset(j-1) + proc_count(j-1)
        ENDDO

        IF(p_pe_work==0) THEN
          WRITE(0,*) 'divide_patch: partial processor counts used for decomposition of merged domain:'
          WRITE(0,*) proc_count(1:num_physdom), 'total: ', n_proc
        ENDIF

      ENDIF

    ENDIF

    IF (.NOT. lsplit_merged_domains) THEN
      num_physdom       = 1
      proc_count(1)     = n_proc
      proc_offset(1)    = 0
      id_physdom(1)     = MAXVAL(subset_flag)
    ENDIF


    ! Fill the cell_desc array, it must contain:
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner

    ALLOCATE(cell_desc(4,wrk_divide_patch%n_patch_cells))

    nc = 0
    nn = 0

    IF(divide_for_radiation) THEN

      cell_desc(1:2,:) = 1.d99 ! for finding min lat/lon

      DO j = 1, wrk_divide_patch%n_patch_cells

        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset

        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index

        nc = nc+1 ! Cell counter

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

        cell_desc(3,nc) = REAL(nc,wp)
        cell_desc(4,nc) = 0.0_wp

      ENDDO

      ncell_offset(0) = 0
      ncell_offset(1) = nc

    ELSE IF (wrk_divide_patch%cell_type == 6) THEN

      DO j = 1, wrk_divide_patch%n_patch_cells

        IF (subset_flag(j) <= 0) CYCLE

        jb = blk_no(j) ! block index
        jl = idx_no(j) ! line index

        nc = nc+1 ! Cell counter

        cell_desc(1,nc) = wrk_divide_patch%cells%center(jl,jb)%lat
        cell_desc(2,nc) = wrk_divide_patch%cells%center(jl,jb)%lon
        cell_desc(3,nc) = REAL(nc,wp)
        cell_desc(4,nc) = 0.0_wp

      ENDDO

      ncell_offset(1) = nc

    ELSE ! domain decomposition for triangular cells with optional splitting into physical domains

      npt = wrk_divide_patch%n_patch_cells+1

      DO jd = 1, num_physdom

        idp = id_physdom(jd)

!CDIR NODEP
        DO j = 1, wrk_divide_patch%n_patch_cells

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index

          ! Disregard outer nest boundary points for the time being. They do very little
          ! computational work, so they can be added to the closest PEs afterwards
          IF (wrk_divide_patch%cells%refin_ctrl(jl,jb) == -1) THEN
            nn = nn+1
            cell_desc(3,npt-nn) = REAL(j,wp)
            CYCLE
          ELSE
            nc = nc+1 ! Cell counter
          ENDIF

          cell_desc(1,nc) = wrk_divide_patch%cells%center(jl,jb)%lat
          cell_desc(2,nc) = wrk_divide_patch%cells%center(jl,jb)%lon

          ! Using the center of the cells for geometric subdivision leads
          ! to "toothed" edges of the subdivision area
          ! Thus we use the minimum lat/lon as subdision criterion.

          IF (cell_desc(1,nc) >= 0._wp) THEN
            DO i = 1, 3
              jl_v = wrk_divide_patch%cells%vertex_idx(jl,jb,i)
              jb_v = wrk_divide_patch%cells%vertex_blk(jl,jb,i)
              cell_desc(1,nc) = MAX(cell_desc(1,nc),wrk_divide_patch%verts%vertex(jl_v,jb_v)%lat)
              cell_desc(2,nc) = MAX(cell_desc(2,nc),wrk_divide_patch%verts%vertex(jl_v,jb_v)%lon)
            ENDDO
          ELSE
            DO i = 1, 3
              jl_v = wrk_divide_patch%cells%vertex_idx(jl,jb,i)
              jb_v = wrk_divide_patch%cells%vertex_blk(jl,jb,i)
              cell_desc(1,nc) = MIN(cell_desc(1,nc),wrk_divide_patch%verts%vertex(jl_v,jb_v)%lat)
              cell_desc(2,nc) = MAX(cell_desc(2,nc),wrk_divide_patch%verts%vertex(jl_v,jb_v)%lon)
            ENDDO
          ENDIF

          cell_desc(3,nc) = REAL(nc,wp)
          cell_desc(4,nc) = 0.0_wp

        ENDDO

        ncell_offset(jd) = nc

      ENDDO

    ENDIF

    DO j = 1, num_physdom

      ncs = ncell_offset(j-1)+1
      nce = ncell_offset(j)
      nc  = ncell_offset(j) - ncell_offset(j-1)

      ALLOCATE(workspace(4,nc))

      CALL divide_cells_by_location(nc, cell_desc(:,ncs:nce), workspace, 0, proc_count(j)-1)

      ! After divide_cells_by_location the cells are sorted by owner,
      ! order them by original cell numbers again

      CALL sort_array_by_row(cell_desc(:,ncs:nce), workspace, 3)

      DEALLOCATE(workspace)

      ! Apply shift of processor IDs
      IF (j > 1) cell_desc(4,ncs:nce) = cell_desc(4,ncs:nce) + proc_offset(j)

    ENDDO

    ! Set owner list (of complete patch)

    nc = 0 ! Counts cells in subset

    DO jd = 1, num_physdom
      idp = id_physdom(jd)

      DO j = 1, wrk_divide_patch%n_patch_cells
        IF(subset_flag(j) == idp .OR. .NOT. lsplit_merged_domains .AND. subset_flag(j)> 0) THEN
          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          IF (wrk_divide_patch%cells%refin_ctrl(jl,jb) /= -1) THEN
            nc = nc+1
            owner(j) = NINT(cell_desc(4,nc))
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! Add outer nest boundary points that have been disregarded so far
    IF (nn > 0) THEN
      nc = 0
      DO WHILE (nc < nn) ! Iterations are needed because outer nest boundary row contains indirect neighbors
        DO i = 1, nn
          j = NINT(cell_desc(3,npt-i))
          IF (owner(j) >= 0) CYCLE
          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index
          DO ii = 1, wrk_divide_patch%cells%num_edges(jl,jb)
            jn = idx_1d(wrk_divide_patch%cells%neighbor_idx(jl,jb,ii),wrk_divide_patch%cells%neighbor_blk(jl,jb,ii))
            IF (owner(jn) >= 0) THEN
              owner(j) = owner(jn)
              nc = nc + 1
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DEALLOCATE(cell_desc)

  END SUBROUTINE divide_subset_geometric

  !-------------------------------------------------------------------------
  !>
  !! Actually divides geometrically by location on cpu_a .. cpu_b
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  RECURSIVE SUBROUTINE divide_cells_by_location(n_cells,cell_desc,work,cpu_a,cpu_b)

    INTEGER, INTENT(in) :: n_cells, cpu_a, cpu_b

    REAL(wp), INTENT(inout) :: cell_desc(4,n_cells),work(4,n_cells)
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner
    !
    ! work contains workspace for sort_array_by_row to avoid local allocation there

    INTEGER cpu_m, n_cells_m
    REAL(wp) :: xmax(2), xmin(2), avglat, scalexp
    !-----------------------------------------------------------------------

    ! If there is only 1 CPU for distribution, we are done

    IF(cpu_a==cpu_b) THEN
      cell_desc(4,:) = REAL(cpu_a,wp)
      RETURN
    ENDIF

    ! Get geometric extensions and total number of points of all patches

    xmin(1) = MINVAL(cell_desc(1,:))
    xmin(2) = MINVAL(cell_desc(2,:))
    xmax(1) = MAXVAL(cell_desc(1,:))
    xmax(2) = MAXVAL(cell_desc(2,:))

    ! average latitude in patch
    avglat  = SUM(cell_desc(1,:))/REAL(n_cells,wp)

    ! account somehow for convergence of meridians - this formula is just empiric
    scalexp = 1._wp - MAX(0._wp,ABS(xmax(1))-1._wp,ABS(xmin(1))-1._wp)
    xmin(2) = xmin(2)*(COS(avglat))**scalexp
    xmax(2) = xmax(2)*(COS(avglat))**scalexp

    ! Get dimension with biggest distance from min to max
    ! and sort cells in this dimension

    IF(xmax(1)-xmin(1) >= xmax(2)-xmin(2)) THEN
      CALL sort_array_by_row(cell_desc, work, 1)
    ELSE
      CALL sort_array_by_row(cell_desc, work, 2)
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

    CALL divide_cells_by_location(n_cells_m,cell_desc(:,1:n_cells_m),work(:,1:n_cells_m),&
      cpu_a,cpu_m)
    CALL divide_cells_by_location(n_cells-n_cells_m,cell_desc(:,n_cells_m+1:n_cells),&
      work(:,n_cells_m+1:n_cells),cpu_m+1,cpu_b)

  END SUBROUTINE divide_cells_by_location

  !-------------------------------------------------------------------------
  !>
  !! Special quicksort implementation for sorting a 2D array by one selected row.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  RECURSIVE SUBROUTINE sort_array_by_row(x,y,row)

    !

    INTEGER, INTENT(in) :: row ! number of row for sorting

    REAL(wp), INTENT(inout) :: x(:,:) ! array to be sorted
    REAL(wp), INTENT(inout) :: y(:,:) ! workspace

    REAL(wp) :: p
    INTEGER :: n, ipiv, ix, iy, i
    !-----------------------------------------------------------------------

    n = SIZE(x,2)

    IF(n<=1) RETURN

    ipiv = (n+1)/2
    p = x(row,ipiv)
    ix = 0
    iy = 1

#ifdef __SX__
    IF (n >= 12) THEN ! Use vectorized version
      DO i=1,n
        IF(i==ipiv) CYCLE
        IF(x(row,i) < p) THEN
          ix = ix+1
          y(1:4,ix) = x(1:4,i)
        ENDIF
      ENDDO

      y(1:4,ix+1) = x(1:4,ipiv) ! Store pivot

      DO i=1,n
        IF(i==ipiv) CYCLE
        IF(x(row,i) >= p) THEN
          iy = iy+1
          y(1:4,ix+iy) = x(1:4,i)
        ENDIF
      ENDDO
!CDIR COLLAPSE
      x(:,:) = y(:,:)
    ELSE  ! use non-vectorized version
      y(1:4,1) = x(1:4,ipiv) ! Store pivot

      DO i=1,n
        IF(i==ipiv) CYCLE
        IF(x(row,i) < p) THEN
          ix = ix+1
          x(1:4,ix) = x(1:4,i)
        ELSE
          iy = iy+1
          y(1:4,iy) = x(1:4,i)
        ENDIF
      ENDDO

      x(1:4,ix+1:ix+iy) = y(1:4,1:iy)
    ENDIF
#else
    y(1:4,1) = x(1:4,ipiv) ! Store pivot

    DO i=1,n
      IF(i==ipiv) CYCLE
      IF(x(row,i) < p) THEN
        ix = ix+1
        x(1:4,ix) = x(1:4,i)
      ELSE
        iy = iy+1
        y(1:4,iy) = x(1:4,i)
      ENDIF
    ENDDO

    x(1:4,ix+1:ix+iy) = y(1:4,1:iy)
#endif

    ipiv = ix+1 ! New pivot location

    IF(ipiv>2)   CALL sort_array_by_row(x(:,:ipiv-1),y(:,:ipiv-1),row)
    IF(ipiv<n-1) CALL sort_array_by_row(x(:,ipiv+1:),y(:,ipiv+1:),row)

  END SUBROUTINE sort_array_by_row


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


#ifdef HAVE_METIS
  !-------------------------------------------------------------------------
  !>
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  !!
  SUBROUTINE divide_subset_metis(subset_flag, n_proc, wrk_divide_patch, &
                                 owner)

    INTEGER, INTENT(in)    :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    TYPE(t_patch), INTENT(in) :: wrk_divide_patch
    INTEGER, INTENT(out)   :: owner(:) ! receives the owner PE for every cell
    ! (-1 for cells not in subset)


    INTEGER :: i, j, jl, jb, jl_n, jb_n, na, nc
    INTEGER, ALLOCATABLE :: local_index(:), tmp(:)
    INTEGER, ALLOCATABLE :: metis_xadj(:), metis_adjncy(:)
    INTEGER :: metis_options(5), metis_edgecut

    !-----------------------------------------------------------------------

    IF(p_pe_work==0) WRITE(0,*) 'divide_patch: Using METIS for area subdivision'

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

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE fill_wrk_decomposition_struct(decomposition_struct, patch)
    TYPE(t_decomposition_structure) :: decomposition_struct
    TYPE(t_patch) :: patch

    INTEGER :: no_of_cells, no_of_verts, cell, vertex
    INTEGER :: jb, jl, i, jl_v, jb_v, return_status

    CHARACTER(*), PARAMETER :: method_name = "fill_wrk_decomposition_struct"

    decomposition_struct%no_of_cells = patch%n_patch_cells
    decomposition_struct%no_of_edges = patch%n_patch_edges
    decomposition_struct%no_of_verts = patch%n_patch_verts
    no_of_cells = decomposition_struct%no_of_cells
    no_of_verts = decomposition_struct%no_of_verts

    ALLOCATE( decomposition_struct%cell_geo_center(no_of_cells), &
      &  decomposition_struct%cells_vertex(3,no_of_cells), &
      &  decomposition_struct%vertex_geo_coord(no_of_verts),  &
      &  stat=return_status)
    IF (return_status > 0) &
      & CALL finish (method_name, "ALLOCATE(decomposition_struct")

    DO cell = 1, no_of_cells
      jb = blk_no(cell) ! block index
      jl = idx_no(cell) ! line index
      decomposition_struct%cell_geo_center(cell)%lat = patch%cells%center(jl,jb)%lat
      decomposition_struct%cell_geo_center(cell)%lon = patch%cells%center(jl,jb)%lon
      DO i = 1, 3
        jl_v = patch%cells%vertex_idx(jl,jb,i)
        jb_v = patch%cells%vertex_blk(jl,jb,i)
        vertex = idx_1d(jl_v, jb_v)
        decomposition_struct%cells_vertex(i, cell) = vertex
      ENDDO
    ENDDO

    DO vertex = 1, no_of_verts
      jb = blk_no(vertex) ! block index
      jl = idx_no(vertex) ! line index
      decomposition_struct%vertex_geo_coord(vertex)%lat = patch%verts%vertex(jl,jb)%lat
      decomposition_struct%vertex_geo_coord(vertex)%lon = patch%verts%vertex(jl,jb)%lon
    ENDDO

    NULLIFY(decomposition_struct%cell_cartesian_center)
    CALL geographical_to_cartesian(decomposition_struct%cell_geo_center, no_of_cells, &
      & decomposition_struct%cell_cartesian_center)

  END SUBROUTINE fill_wrk_decomposition_struct
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE clean_wrk_decomposition_struct(decomposition_struct)
    TYPE(t_decomposition_structure) :: decomposition_struct

    DEALLOCATE( decomposition_struct%cell_geo_center, &
      &  decomposition_struct%cells_vertex,           &
      &  decomposition_struct%vertex_geo_coord,       &
      &  decomposition_struct%cell_cartesian_center)

  END SUBROUTINE clean_wrk_decomposition_struct
  !-------------------------------------------------------------------------

END MODULE mo_setup_subdivision
