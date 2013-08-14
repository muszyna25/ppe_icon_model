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
MODULE mo_ext_decompose_patches
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
  USE mo_model_domain,       ONLY: t_patch
  USE mo_mpi,                ONLY: p_bcast, p_sum, proc_split
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_UNDEFINED, MPI_COMM_NULL
#endif
  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_test, &
    & my_process_is_mpi_seq, process_mpi_all_test_id, process_mpi_all_workroot_id, &
    & p_pe_work, p_n_work,  get_my_mpi_all_id, my_process_is_mpi_parallel

  USE mo_parallel_config,       ONLY:  nproma, p_test_run, ldiv_phys_dom, &
    & division_method, division_file_name, n_ghost_rows, div_from_file,   &
    & div_geometric, ext_div_medial, ext_div_medial_cluster, ext_div_medial_redrad, &
    & ext_div_medial_redrad_cluster, ext_div_from_file, redrad_split_factor

#ifdef HAVE_METIS
  USE mo_parallel_config,    ONLY: div_metis
#endif
  USE mo_communication,      ONLY: blk_no, idx_no, idx_1d
  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
    & grf_bdyintp_end_c, grf_bdyintp_end_e, grf_fbk_start_c, grf_fbk_start_e, &
    & grf_bdywidth_c, grf_bdywidth_e, grf_nudgintp_start_c, grf_nudgintp_start_e
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, patch_weight
  USE mo_alloc_patches,ONLY: allocate_basic_patch, allocate_remaining_patch, &
                             deallocate_basic_patch, deallocate_patch
  USE mo_decomposition_tools, ONLY: t_decomposition_structure, divide_geometric_medial, &
    & read_ascii_decomposition
  USE mo_math_utilities,      ONLY: geographical_to_cartesian
  USE mo_ocean_config,        ONLY: ignore_land_points

  USE mo_setup_subdivision, ONLY: discard_large_arrays,     &
                                  divide_parent_cells,      &
#ifdef HAVE_METIS
                                  divide_subset_metis,      &
#endif
                                  get_local_index,          &
                                  divide_subset_geometric,  &
                                  sort_array_by_row,        &
                                  count_entries,            &
                                  divide_cells_by_location, &
                                  divide_patch

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: ext_decompose_patches, discard_large_arrays, &
            divide_parent_cells, get_local_index,        &
            divide_subset_geometric, sort_array_by_row,  &
            count_entries, divide_cells_by_location,     &
            divide_patch
#ifdef HAVE_METIS
  PUBLIC :: divide_subset_metis
#endif

  ! pointers to the work patches
  TYPE(t_patch), POINTER :: wrk_p_parent_patch_g

  ! Private flag if patch should be divided for radiation calculation
  LOGICAL :: divide_for_radiation = .FALSE.

CONTAINS

  !------------------------------------------------------------------
  !>
  !!  Divide patches and interpolation states for mpi parallel runs.
  !!
  !!  @note This implementation uses an ext decomposition library driver.
  !!
  SUBROUTINE ext_decompose_patches( patch_2D, p_patch_global )

    TYPE(t_patch), INTENT(INOUT), TARGET :: patch_2D(n_dom_start:)
    TYPE(t_patch), INTENT(INOUT), TARGET :: p_patch_global(n_dom_start:)

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_subdivision:ext_decompose_patches")

    ! Local variables:
    INTEGER :: jg, jgp, jc, jgc, n, &
      &        l1, i_nchdom, n_procs_decomp
    INTEGER :: nprocs(p_patch_global(1)%n_childdom)
    INTEGER, POINTER :: cell_owner(:)
    INTEGER, POINTER :: cell_owner_p(:)
    REAL(wp) :: weight(p_patch_global(1)%n_childdom)
    INTEGER(i8) :: npts_global(4)
    ! (Optional:) Print a detailed summary on model grid points
    LOGICAL :: is_compute_grid
    TYPE(t_patch), ALLOCATABLE, TARGET :: p_patch_out(:), p_patch_lp_out(:)

    CALL message(routine, 'start of ext domain decomposition')

    n_procs_decomp = p_n_work


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

      patch_2D(n_dom_start:1)%n_proc = n_procs_decomp
      patch_2D(n_dom_start:1)%proc0  = 0

      ! ... for 1st generation childs

      n = 0
      DO jc = 1, p_patch_global(1)%n_childdom
        jgc = p_patch_global(1)%child_id(jc)
        patch_2D(jgc)%proc0  = n
        patch_2D(jgc)%n_proc = nprocs(jc)
        n = n + nprocs(jc)
      ENDDO

      ! ... for deeper level descandants

      DO jg = 2, n_dom

        jgp = p_patch_global(jg)%parent_id

        IF(jgp /= 1) THEN
          patch_2D(jg)%n_proc = patch_2D(jgp)%n_proc
          patch_2D(jg)%proc0  = patch_2D(jgp)%proc0
        ENDIF

      ENDDO

    ELSE

      ! No splitting, proc0, n_proc are identical for all patches

      IF(p_pe_work==0) WRITE(0,*) 'No splitting of processor grid'
      patch_2D(:)%n_proc = n_procs_decomp
      patch_2D(:)%proc0  = 0

    ENDIF

#ifdef NOMPI
    patch_2D(:)%comm = 0
#else
    patch_2D(:)%comm = MPI_COMM_NULL
#endif
    patch_2D(:)%rank = -1

    ! Divide patches

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
      CALL divide_patch_cells(p_patch_global(jg), jg, patch_2D(jg)%n_proc, &
                              patch_2D(jg)%proc0, cell_owner,              &
                              patch_2D(jg)%cells%radiation_owner )

      DEALLOCATE(p_patch_global(jg)%cells%phys_id)

!      IF(jg > n_dom_start) THEN
!        ! Assign the cell owners of the current patch to the parent cells
!        ! for the construction of the local parent:
!        ALLOCATE(cell_owner_p(p_patch_global(jgp)%n_patch_cells))
!        CALL divide_parent_cells(p_patch_global(jg),cell_owner,cell_owner_p)
!      ENDIF

      ! Please note: Previously, for jg==0 no ghost rows were set.
      ! Currently, we need ghost rows for jg==0 also for dividing the int state and grf state
      ! Have still to check if int state/grf state is needed at all for jg==0,
      ! if this is not the case, the ghost rows can be dropped again.

        is_compute_grid = .true.
!        IF (ignore_land_points) &
!          is_compute_grid = .false.

        CALL divide_patch(patch_2D(jg), p_patch_global(jg), cell_owner, n_ghost_rows, is_compute_grid, p_pe_work)

!        IF(jg > n_dom_start) THEN
!          CALL divide_patch(p_patch_local_parent(jg), p_patch_global(jgp), cell_owner_p, 1, .FALSE., p_pe_work)
!        ENDIF

      DEALLOCATE(cell_owner)
      IF(jg > n_dom_start) DEALLOCATE(cell_owner_p)

    ENDDO


    ! The global patches may be discarded now
    DO jg = n_dom_start, n_dom
      CALL deallocate_basic_patch( p_patch_global(jg) )
    ENDDO

  END SUBROUTINE ext_decompose_patches

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

  SUBROUTINE divide_patch_cells(wrk_p_patch_g, patch_no, n_proc, proc0, cell_owner, radiation_owner)

    TYPE(t_patch), INTENT(INOUT) :: wrk_p_patch_g
    INTEGER, INTENT(IN)  :: patch_no !> Patch number, identifies how to decompose it
    INTEGER, INTENT(IN)  :: n_proc !> Number of processors for split
    INTEGER, INTENT(IN)  :: proc0  !> First processor of patch
    INTEGER, POINTER :: cell_owner(:) !> Cell division
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
    ! only procs 0 will decompoise and then will broadcast
    IF(p_pe_work == 0) THEN
     write(0,*) "Divide into ", n_proc

     IF (division_method(patch_no) == ext_div_from_file) THEN

        IF (division_file_name(patch_no) == "") THEN
          use_division_file_name = &
            & TRIM(get_filename_noext(wrk_p_patch_g%grid_filename))//'.cell_domain_ids'
        ELSE
          use_division_file_name = division_file_name(patch_no)
        ENDIF

        CALL read_ascii_decomposition(use_division_file_name, cell_owner, wrk_p_patch_g%n_patch_cells)

      ELSE IF (division_method(patch_no) > 100) THEN
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

        !----------------------------------------------------------
      ENDIF ! subdivision method

      ! Quick check for correct values
      IF(MINVAL(cell_owner(:)) < 0 .or. MAXVAL(cell_owner(:)) >= n_proc) THEN
        WRITE(0,*) "n_porc=",n_proc, " MINAVAL=", MINVAL(cell_owner(:)), &
          " MAXVAL=", MAXVAL(cell_owner(:))
        CALL finish('divide_patch','Illegal subdivision in input file')
      ENDIF

    ENDIF !p_pe_work == 0

    CALL p_bcast(cell_owner, 0, comm=p_comm_work)

    IF (division_method(patch_no) == ext_div_medial_redrad_cluster .OR. &
      & division_method(patch_no) == ext_div_medial_redrad) THEN
      ! distribute the radiation owner
      IF (p_pe_work /= 0) &
          ALLOCATE(radiation_owner(wrk_p_patch_g%n_patch_cells))

      CALL p_bcast(radiation_owner, 0, comm=p_comm_work)
      wrk_p_patch_g%cells%radiation_owner => radiation_owner

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

END MODULE mo_ext_decompose_patches

