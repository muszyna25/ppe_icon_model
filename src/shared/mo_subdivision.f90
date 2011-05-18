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
  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: success, min_rlcell, max_rlcell,  &
    & min_rledge, max_rledge, min_rlvert, max_rlvert,                &
    & min_rlcell_int, min_rledge_int, min_rlvert_int, max_hw
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish, message

  USE mo_run_nml,            ONLY: nproma, i_cell_type, ltransport, &
    & iequations
  USE mo_io_units,           ONLY: find_next_free_unit, filename_max
  USE mo_model_domain,       ONLY: t_patch, t_grid_cells
  USE mo_interpolation,      ONLY: t_int_state, rbf_vec_dim_c, rbf_vec_dim_e, &
    & rbf_vec_dim_v, rbf_c2grad_dim
  USE mo_grf_interpolation,  ONLY: t_gridref_state, t_gridref_single_state
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_mpi,                ONLY: p_pe, p_nprocs, p_bcast, p_send, p_recv
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_UNDEFINED, MPI_COMM_NULL
#endif
  USE mo_parallel_nml,       ONLY: p_test_pe, p_test_run, p_n_work, p_pe_work, &
    & p_comm_work, p_work_pe0, division_method, n_ghost_rows, div_from_file,   &
    & div_geometric
#ifdef HAVE_METIS
  USE mo_parallel_nml        ONLY: div_metis
#endif
  USE mo_communication,      ONLY: setup_comm_pattern, blk_no, idx_no, idx_1d
  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
    & grf_bdyintp_end_c, grf_bdyintp_end_e, grf_fbk_start_c, grf_fbk_start_e, &
    & grf_bdywidth_c, grf_bdywidth_e, grf_nudgintp_start_c, grf_nudgintp_start_e
  USE mo_grid_nml,           ONLY: n_dom, n_dom_start, nroot, patch_weight
  USE mo_model_domimp_patches,ONLY: destruct_patches, allocate_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp_state,          ONLY: allocate_int_state, destruct_2d_interpol_state
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_grf_intp_state,      ONLY: allocate_grf_state, destruct_2d_gridref_state
  USE mo_interpol_nml,        ONLY: i_cori_method, lsq_lin_set, lsq_high_set
  USE mo_atmo_control,        ONLY: p_patch_global, p_patch_subdiv, p_patch, &
    & p_int_state_global, p_int_state_subdiv, p_int_state, &
    & p_grf_state_global, p_grf_state_subdiv, p_grf_state


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: version = '$Id$'

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: decompose_atmo_domain, copy_processor_splitting
  PUBLIC :: set_patch_communicators

  ! pointers to the work patches
  TYPE(t_patch), POINTER :: wrk_p_patch, wrk_p_parent_patch
  TYPE(t_patch), POINTER :: wrk_p_patch_g, wrk_p_parent_patch_g
  TYPE(t_patch), POINTER :: wrk_divide_patch

  ! pointers to the work states
  TYPE(t_int_state), POINTER :: wrk_int_state_in, wrk_int_state_out
  TYPE(t_gridref_state), POINTER :: wrk_gridref_state_in, wrk_gridref_state_out

  !-------------------------------------------------------------------------
  ! Definition of local parent patches
  ! For any given patch p_patch(jg) and jgp = p_patch(jg)%parent_id,
  ! p_patch_local_parent(jg) has the same resolution as p_patch(jgp)
  ! but it covers only the area of p_patch(jgp) which is covered by its child p_patch(jg)
  ! and it is divided in the same manner as p_patch(jg).
  ! Please note that p_patch_local_parent(1) is undefined if n_dom_start = 1

  TYPE(t_patch),         ALLOCATABLE, TARGET :: p_patch_local_parent(:)
  TYPE(t_int_state),     ALLOCATABLE, TARGET :: p_int_state_local_parent(:)
  TYPE(t_gridref_state), ALLOCATABLE, TARGET :: p_grf_state_local_parent(:)
  PUBLIC :: p_patch_local_parent, p_int_state_local_parent, p_grf_state_local_parent
  !-------------------------------------------------------------------------

  ! Flag if processor splitting is active
  LOGICAL, PUBLIC :: proc_split = .FALSE.

  ! Private flag if patch should be divided for radiation calculation
  LOGICAL :: divide_for_radiation = .FALSE.

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!  In case of a test run: Copies processor splitting to test PE
  SUBROUTINE copy_processor_splitting(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER :: ibuf(n_dom_start:n_dom,2)

    IF(.NOT. p_test_run) RETURN ! Nothing to do

    IF(p_pe == p_work_pe0) THEN
      CALL p_send(proc_split, p_test_pe, 1)
      ibuf(:,1) = p_patch(:)%n_proc
      ibuf(:,2) = p_patch(:)%proc0
      CALL p_send(ibuf, p_test_pe, 2)
    ENDIF

    IF(p_pe == p_test_pe) THEN
      CALL p_recv(proc_split, p_work_pe0, 1)
      CALL p_recv(ibuf, p_work_pe0, 2)
      p_patch(:)%n_proc = ibuf(:,1)
      p_patch(:)%proc0  = ibuf(:,2)
    ENDIF

  END SUBROUTINE copy_processor_splitting
  !-------------------------------------------------------------------------
  !>
  !!  Sets the communicators in the patches if these have been read.
  SUBROUTINE set_patch_communicators(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

#ifdef NOMPI
    CALL finish('mo_subdivision','set_patch_communicators must only be called in parallel runs')
#else
    INTEGER jc, jgc, jg, jgp, n_proc_total, comm, mpierr
    INTEGER, ALLOCATABLE :: patch_no(:)


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

      DO jc = 1, p_patch_global(1)%n_childdom
        jgc = p_patch_global(1)%child_id(jc)
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

        jgp = p_patch_global(jg)%parent_id

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
  SUBROUTINE decompose_atmo_domain()

#ifdef NOMPI
    CALL finish('mo_subdivision','decompose_atmo_domain must only be called in parallel runs')
#else
    INTEGER :: mpierr

    INTEGER :: ist, jg, jgp, jc, jgc, comm, my_color, n
    INTEGER :: nprocs(p_patch_global(1)%n_childdom)
    INTEGER, ALLOCATABLE :: cell_owner(:)
    REAL(wp) :: weight(p_patch_global(1)%n_childdom)

    CALL message('mo_subdivision:decompose_atmo_domain',   &
                 'start of domain decomposition')

    ALLOCATE (p_patch_subdiv(n_dom_start:n_dom), &
      & p_int_state_subdiv(n_dom_start:n_dom),   &
      & p_grf_state_subdiv(n_dom_start:n_dom),stat=ist)
    IF (ist /= success) THEN
      CALL finish('decompose_atmo_domain','allocation for subdivided patches/states failed')
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
      CALL finish('decompose_atmo_domain','Weight for root patch must be 0')
    DO jg = 2, n_dom
      jgp = p_patch_global(jg)%parent_id
      IF(jgp /= 1 .AND. patch_weight(jg) > 0._wp) &
        CALL finish('decompose_atmo_domain','Weight for higher level patch must be 0')
    ENDDO

    IF(proc_split) THEN

      IF(p_pe_work==0) PRINT *,'Splitting processor grid for first level patches'
      IF(p_pe_work==0) PRINT '(a,10f12.3)','Weights for first level patches:',weight(:)

      ! In this case, the working processor set must be at least as big
      ! as the number of childs of the root patch
      IF(p_patch_global(1)%n_childdom > p_n_work) &
        CALL finish('decompose_atmo_domain','Too few procs for processor grid splitting')

      ! Get the number of procs per patch according to weight(:).
      ! Every patch gets at least 1 proc (of course!):
      nprocs(:) = 1

      ! The remaining procs are divided among patches similar to
      ! the d'Hondt method for elections

      DO n = p_patch_global(1)%n_childdom+1, p_n_work
        jg = MAXLOC(weight(:)/REAL(nprocs(:)+1,wp),1)
        nprocs(jg) = nprocs(jg)+1
      ENDDO

      IF(p_pe_work==0) THEN
        PRINT *,'Processor splitting:'
        DO jc = 1, p_patch_global(1)%n_childdom
          jgc =  p_patch_global(1)%child_id(jc)
          PRINT '(a,i0,a,f10.3,a,i0,a,i0)',                                             &
            &   'Patch ',jgc,' weight ',weight(jc),' gets ',nprocs(jc),' of ',p_n_work
        ENDDO
      ENDIF

      ! Set proc0, n_proc, comm, rank for all patches ...

      ! ... for the root patch and patch 0 if it exists

      p_patch_subdiv(n_dom_start:1)%comm   = p_comm_work
      p_patch_subdiv(n_dom_start:1)%rank   = p_pe_work
      p_patch_subdiv(n_dom_start:1)%n_proc = p_n_work
      p_patch_subdiv(n_dom_start:1)%proc0  = 0

      ! ... for 1st generation childs

      my_color = MPI_UNDEFINED
      n = 0
      DO jc = 1, p_patch_global(1)%n_childdom
        jgc = p_patch_global(1)%child_id(jc)
        p_patch_subdiv(jgc)%proc0  = n
        p_patch_subdiv(jgc)%n_proc = nprocs(jc)
        IF(p_pe_work >= n .AND. p_pe_work < n+nprocs(jc)) my_color = jc
        n = n + nprocs(jc)
      ENDDO

      ! Safety only
      IF(my_color == MPI_UNDEFINED) &
        CALL finish('decompose_atmo_domain','Internal error: my_color == MPI_UNDEFINED')

      ! Split communicator among childs of root patch

      CALL MPI_Comm_split(p_comm_work, my_color, p_pe_work, comm, mpierr)

      ! Set comm and rank for childs of root patch

      DO jc = 1, p_patch_global(1)%n_childdom
        jgc = p_patch_global(1)%child_id(jc)
        IF(my_color == jc) THEN
          p_patch_subdiv(jgc)%comm = comm
          CALL MPI_Comm_rank(comm, p_patch_subdiv(jgc)%rank, mpierr)
        ELSE
          p_patch_subdiv(jgc)%comm = MPI_COMM_NULL ! We should never use this comm
          p_patch_subdiv(jgc)%rank = -1
        ENDIF
      ENDDO

      ! ... for deeper level descandants

      DO jg = 2, n_dom

        jgp = p_patch_global(jg)%parent_id

        IF(jgp /= 1) THEN
          p_patch_subdiv(jg)%comm   = p_patch_subdiv(jgp)%comm
          p_patch_subdiv(jg)%rank   = p_patch_subdiv(jgp)%rank
          p_patch_subdiv(jg)%n_proc = p_patch_subdiv(jgp)%n_proc
          p_patch_subdiv(jg)%proc0  = p_patch_subdiv(jgp)%proc0
        ENDIF

      ENDDO

    ELSE

      ! No splitting, proc0, n_proc, comm, rank are identical for all patches

      IF(p_pe_work==0) PRINT *,'No splitting of processor grid'
      p_patch_subdiv(:)%comm   = p_comm_work
      p_patch_subdiv(:)%rank   = p_pe_work
      p_patch_subdiv(:)%n_proc = p_n_work
      p_patch_subdiv(:)%proc0  = 0

    ENDIF

    ! Divide patches and set up communication patterns

    DO jg = n_dom_start, n_dom

      jgp = p_patch_global(jg)%parent_id

      NULLIFY(wrk_p_parent_patch) ! Safety only, not needed for patch division

      IF(jg == n_dom_start) THEN
        NULLIFY(wrk_p_parent_patch_g) ! Must be NULL for global patch
      ELSE
        wrk_p_parent_patch_g => p_patch_global(jgp)
      ENDIF

      wrk_p_patch_g => p_patch_global(jg)
      wrk_p_patch   => p_patch_subdiv(jg)

      ! Set division method, divide_for_radiation is only used for patch 0

      divide_for_radiation = (jg == 0)

      ALLOCATE(cell_owner(p_patch_global(jg)%n_patch_cells))
      CALL divide_patch_cells(p_patch_subdiv(jg)%n_proc, p_patch_subdiv(jg)%proc0, cell_owner)

      ! Please note: For jg==0 no ghost rows are set
      CALL divide_patch(cell_owner, MERGE(0, n_ghost_rows, jg==0), .TRUE.)
      DEALLOCATE(cell_owner)

      IF(jg == n_dom_start) CYCLE

      wrk_p_parent_patch   => p_patch_subdiv(jgp)

      CALL setup_comm_cpy_interpolation()
      CALL setup_comm_grf_interpolation()
      CALL setup_comm_ubc_interpolation()

    ENDDO

    ! Create local parents for patches

    ALLOCATE(p_patch_local_parent(n_dom_start+1:n_dom),     &
             p_int_state_local_parent(n_dom_start+1:n_dom), &
             p_grf_state_local_parent(n_dom_start+1:n_dom) )

    DO jg = n_dom_start+1, n_dom

      jgp = p_patch_global(jg)%parent_id

      ALLOCATE(cell_owner(p_patch_global(jgp)%n_patch_cells))
      CALL divide_parent_cells(p_patch_subdiv(jg),p_patch_global(jg),cell_owner)

      wrk_p_patch_g => p_patch_global(jgp)
      wrk_p_patch   => p_patch_local_parent(jg)
      CALL divide_patch(cell_owner, 1, .FALSE.)
      DEALLOCATE(cell_owner)

      CALL set_parent_child_relations(p_patch_local_parent(jg), p_patch_subdiv(jg), &
        &                             p_patch_global(jgp), p_patch_global(jg))
      CALL set_glb_loc_comm(p_patch_subdiv(jgp), p_patch_local_parent(jg), &
        &                   p_patch_subdiv(jg)%parent_child_index)

    ENDDO

    ! Divide interpolation states (for regular patches and local parents)

    DO jg = n_dom_start, n_dom

      CALL allocate_int_state(p_patch_subdiv(jg), p_int_state_subdiv(jg))

      wrk_p_patch       => p_patch_subdiv(jg)
      wrk_int_state_in  => p_int_state_global(jg)
      wrk_int_state_out => p_int_state_subdiv(jg)
      CALL divide_int_state()

      IF(n_dom_start==0 .OR. n_dom > 1) THEN
        CALL allocate_grf_state(p_patch_subdiv(jg), p_grf_state_subdiv(jg))
        wrk_gridref_state_in  => p_grf_state_global(jg)
        wrk_gridref_state_out => p_grf_state_subdiv(jg)
        CALL divide_grf_state()
      ENDIF

      IF(jg == n_dom_start) CYCLE

      jgp = p_patch_global(jg)%parent_id

      CALL allocate_int_state(p_patch_local_parent(jg), p_int_state_local_parent(jg))
      CALL allocate_grf_state(p_patch_local_parent(jg), p_grf_state_local_parent(jg))

      wrk_p_patch       => p_patch_local_parent(jg)
      wrk_int_state_in  => p_int_state_global(jgp)
      wrk_int_state_out => p_int_state_local_parent(jg)
      CALL divide_int_state()

      wrk_gridref_state_in  => p_grf_state_global(jgp)
      wrk_gridref_state_out => p_grf_state_local_parent(jg)
      CALL divide_grf_state()
    ENDDO

    ! Set pointers to subdivided patches/states

    p_patch => p_patch_subdiv
    p_int_state => p_int_state_subdiv
    p_grf_state => p_grf_state_subdiv

    ! The global patches/states may be discarded now

    IF (n_dom_start==0 .OR. n_dom > 1) THEN
      CALL destruct_2d_gridref_state( p_patch_global, p_grf_state_global )
      DEALLOCATE (p_grf_state_global, STAT = ist)
      IF (ist/=SUCCESS)THEN
        CALL message('mo_subdivision:decompose_atmo_domain',   &
                     'deallocation of p_grf_state_global failed')
      ENDIF
    ENDIF

    CALL destruct_2d_interpol_state( p_int_state_global )
    DEALLOCATE (p_int_state_global, STAT = ist)
    IF (ist/=SUCCESS)THEN
      CALL message('mo_subdivision:decompose_atmo_domain',   &
                   'deallocation of p_int_state_global failed')
    ENDIF

    CALL destruct_patches( p_patch_global )
    DEALLOCATE( p_patch_global, STAT = ist)
    IF (ist/=SUCCESS)THEN
      CALL message('mo_subdivision:decompose_atmo_domain',   &
                   'deallocation of p_patch_global failed')
    ENDIF

    CALL message('mo_subdivision:decompose_atmo_domain',   &
                 'end of domain decomposition')
#endif

  END SUBROUTINE decompose_atmo_domain


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

    CHARACTER(filename_max) :: div_file

    !-----------------------------------------------------------------------
    ! This routine must not be called in a single CPU run
    IF(p_nprocs == 1 .or. p_pe == p_test_pe) &
      & CALL finish('divide_patch','must not be called in a single CPU run')


    ! Please note: Unfortunatly we cannot use p_io for doing I/O,
    ! since this might be the test PE which is never calling this routine
    ! (this is the case in the actual setup).
    ! Thus we use the worker PE 0 for I/O and don't use message() for output.


    IF(division_method==div_from_file) THEN

      ! Area subdivision is read from file

      IF(p_pe_work == 0) THEN

        WRITE (div_file,'(a,i0,2(a,i2.2),a)') 'iconR',nroot,'B',wrk_p_patch_g%level, &
          & '_DOM',wrk_p_patch_g%id,'-div.txt'

        n = find_next_free_unit(10,99)

        OPEN(n,FILE=div_file,STATUS='OLD',IOSTAT=i)
        IF(i /= 0) CALL finish('divide_patch','Unable to open input file: '//TRIM(div_file))

        DO j = 1, wrk_p_patch_g%n_patch_cells
          READ(n,*,IOSTAT=i) cell_owner(j)
          IF(i /= 0) CALL finish('divide_patch','Error reading: '//TRIM(div_file))
        ENDDO
        CLOSE(n)

        ! Quick check for correct values

        IF(MINVAL(cell_owner(:)) < 0 .or. MAXVAL(cell_owner(:)) >= n_proc) &
          & CALL finish('divide_patch','Illegal subdivision in input file')

      ENDIF

      CALL p_bcast(cell_owner, 0, comm=p_comm_work)

      IF(p_pe_work==0) PRINT *,'Successfully read: '//TRIM(div_file)

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

    IF(p_pe_work==0) THEN
      PRINT '(a,i0,a,i0)','Patch: ',wrk_p_patch_g%id,&
        & ' Total number of cells: ',wrk_p_patch_g%n_patch_cells
      DO n = 0, p_n_work-1
        PRINT '(a,i5,a,i8)','PE',n,' # cells: ',COUNT(cell_owner(:) == n)
      ENDDO
    ENDIF

  END SUBROUTINE divide_patch_cells

  !-------------------------------------------------------------------------------------------------
  !>
  !! Sets the owner for the division of the cells of the parent patch
  !! with the same subdivision as the child

  SUBROUTINE divide_parent_cells(p_patch, p_patch_g, cell_owner)

    TYPE(t_patch), INTENT(IN) :: p_patch    !> Divided patch for which the parent should be divided
    TYPE(t_patch), INTENT(IN) :: p_patch_g  !> Global patch to p_patch
    INTEGER, INTENT(OUT) :: cell_owner(:) !> Output: Cell division for parent.
                                          !> Must be allocated to n_patch_cells of the global parent

    INTEGER :: j, jl, jb, jl_p, jb_p
    INTEGER :: cnt(SIZE(cell_owner))

    cell_owner(:) = -1
    cnt(:) = 0

    DO j = 1, p_patch_g%n_patch_cells

      jb = blk_no(j) ! block index
      jl = idx_no(j) ! line index
      jl_p = p_patch_g%cells%parent_idx(jl,jb)
      jb_p = p_patch_g%cells%parent_blk(jl,jb)

      IF(cell_owner(idx_1d(jl_p,jb_p)) < 0) THEN
        cell_owner(idx_1d(jl_p,jb_p)) = p_patch%cells%owner_g(j)
      ELSEIF(cell_owner(idx_1d(jl_p,jb_p)) /= p_patch%cells%owner_g(j)) THEN
        CALL finish('divide_parent_cells','Divided parent cell encountered')
      ENDIF
      cnt(idx_1d(jl_p,jb_p)) = cnt(idx_1d(jl_p,jb_p)) + 1 ! Only for safety check below
    ENDDO

    ! Safety check
    IF(ANY(cnt(:)/=0 .AND. cnt(:)/=4)) &
      & CALL finish('divide_parent_cells','Incomplete parent cell encountered')

    IF(p_pe_work==0) THEN
      PRINT '(a,i0,a,i0)','Patch: ',p_patch%id,&
        & ' Total number of parent cells: ',COUNT(cell_owner(:) >= 0)
      DO j = 0, p_n_work-1
        PRINT '(a,i5,a,i8)','PE',j,' # parent cells: ',COUNT(cell_owner(:) == j)
      ENDDO
    ENDIF

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

  SUBROUTINE divide_patch(cell_owner, n_boundary_rows, l_compute_grid)

    INTEGER, INTENT(IN) :: cell_owner(:)
    INTEGER, INTENT(IN) :: n_boundary_rows
    LOGICAL, INTENT(IN) :: l_compute_grid

    INTEGER :: n, i, j, jv, je, jl, jb, jl_g, jb_g, jl_e, jb_e, jl_v, jb_v, ilev, iown,   &
               ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, jl_c, jb_c, jc, irlev, ilev1, ilev_st, &
               jg, i_nchdom, irl0, irl1, irl2, irl3

    INTEGER, ALLOCATABLE :: owner_c(:), owner_e(:), owner_v(:), tmp(:)
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

    WHERE(cell_owner(:)==p_pe_work) flag_c(:) = 0
    WHERE(cell_owner(:)==p_pe_work) flag2_c(:) = 0

    jg = wrk_p_patch_g%id
    i_nchdom = MAX(1,wrk_p_patch_g%n_childdom)

    ! find inner edges/verts and ghost cells/edges/verts

    DO ilev = 0, n_boundary_rows

      ! Flag cells belonging to this level.
      ! Cells belonging to the kernel (ilev==0) are already flagged.
      ! Cells which have no global owner are never included to the ghost cells!

      IF(ilev>0) THEN

        DO j = 1, wrk_p_patch_g%n_patch_cells

          jb = blk_no(j) ! block index
          jl = idx_no(j) ! line index

          IF(cell_owner(j)>=0 .AND. flag_c(j)<0) THEN

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
    wrk_p_patch%level = wrk_p_patch_g%level
    wrk_p_patch%id    = wrk_p_patch_g%id
    wrk_p_patch%parent_id = wrk_p_patch_g%parent_id
    wrk_p_patch%parent_child_index = wrk_p_patch_g%parent_child_index
    wrk_p_patch%child_id(:) = wrk_p_patch_g%child_id(:)
    wrk_p_patch%n_childdom = wrk_p_patch_g%n_childdom
    wrk_p_patch%nlev   = wrk_p_patch_g%nlev
    wrk_p_patch%nlevp1 = wrk_p_patch_g%nlevp1
    wrk_p_patch%nshift = wrk_p_patch_g%nshift
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
          IF (wrk_p_patch%edges%owner_g(je) == p_pe_work) flag2_e(je)=0
          IF (.NOT.l_compute_grid .AND. flag2_e(je)==1) flag2_e(je)=0

          jl_v = wrk_p_patch_g%cells%vertex_idx(jl,jb,i)
          jb_v = wrk_p_patch_g%cells%vertex_blk(jl,jb,i)
          jv = idx_1d(jl_v, jb_v)
          IF (wrk_p_patch%verts%owner_g(jv) == p_pe_work) flag2_v(jv)=0
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
      IF (i_cell_type==6) THEN ! for hexagons, there are no even-order halo cells
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
          write(0,'(a,2i5,2i4,4i7)') 'cells',p_pe_work,wrk_p_patch%id,i,j,     &
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
          write(0,'(a,2i5,2i4,4i7)') 'edges',p_pe_work,wrk_p_patch%id,i,j,     &
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
          write(0,'(a,2i5,2i4,4i7)') 'verts',p_pe_work,wrk_p_patch%id,i,j,     &
            wrk_p_patch%verts%start_blk(i,j),wrk_p_patch%verts%start_idx(i,j), &
            wrk_p_patch%verts%end_blk(i,j),  wrk_p_patch%verts%end_idx(i,j)
        ENDDO
      ENDDO
      CALL finish('divide_patch','Error in vertex start/end indices')
    ENDIF


    !-----------------------------------------------------------------------------------------------
    ! Set up communication patterns
    !-----------------------------------------------------------------------------------------------

    ALLOCATE(owner_c(wrk_p_patch%n_patch_cells))
    ALLOCATE(owner_e(wrk_p_patch%n_patch_edges))
    ALLOCATE(owner_v(wrk_p_patch%n_patch_verts))

    ! Set the owner arrays for cells/edges/verts which have to be transferred.
    ! The following setting always transfers edges/verts if they are owned
    ! (according to the MIN/MAX PE setting above) by another PE, which implies
    ! that boundary edges/verts (with flag2_e/flag2_v = 1) can be excluded from
    ! prognostic computations if they are followed by a synchronization call.

    owner_c(:) = wrk_p_patch%cells%owner_g(wrk_p_patch%cells%glb_index(:))
    WHERE(owner_c(:) == p_pe_work) owner_c(:) = -1

    owner_e(:) = wrk_p_patch%edges%owner_g(wrk_p_patch%edges%glb_index(:))
    WHERE(owner_e(:) == p_pe_work) owner_e(:) = -1

    owner_v(:) = wrk_p_patch%verts%owner_g(wrk_p_patch%verts%glb_index(:))
    WHERE(owner_v(:) == p_pe_work) owner_v(:) = -1

    ! Set communication patterns for boundary exchange
    CALL setup_comm_pattern(wrk_p_patch%n_patch_cells, owner_c, wrk_p_patch%cells%glb_index, &
      & wrk_p_patch%cells%loc_index, wrk_p_patch%comm_pat_c)

    CALL setup_comm_pattern(wrk_p_patch%n_patch_edges, owner_e, wrk_p_patch%edges%glb_index, &
      & wrk_p_patch%edges%loc_index, wrk_p_patch%comm_pat_e)

    CALL setup_comm_pattern(wrk_p_patch%n_patch_verts, owner_v, wrk_p_patch%verts%glb_index, &
      & wrk_p_patch%verts%loc_index, wrk_p_patch%comm_pat_v)

    DEALLOCATE(owner_c, owner_e, owner_v)

    ! Set reduced communication pattern containing only level-1 halo cells (immediate neighbors)
    jc = idx_1d(wrk_p_patch%cells%end_idx(min_rlcell_int-1,1), &
                wrk_p_patch%cells%end_blk(min_rlcell_int-1,1))
    ALLOCATE(owner_c(jc))
    owner_c(1:jc) = wrk_p_patch%cells%owner_g(wrk_p_patch%cells%glb_index(1:jc))
    WHERE(owner_c(:) == p_pe_work) owner_c(:) = -1

    CALL setup_comm_pattern(jc, owner_c, wrk_p_patch%cells%glb_index, &
      & wrk_p_patch%cells%loc_index, wrk_p_patch%comm_pat_c1)

    DEALLOCATE(owner_c)

    ! For gathering the global fields on p_pe_work==0
    ALLOCATE(tmp(MAX(wrk_p_patch_g%n_patch_cells, wrk_p_patch_g%n_patch_edges)))

    DO j = 1, SIZE(tmp)
      tmp(j) = j ! Global/local index in global array, i.e. identity!
    ENDDO

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(wrk_p_patch_g%n_patch_cells, wrk_p_patch%cells%owner_g, tmp, &
        & wrk_p_patch%cells%loc_index, wrk_p_patch%comm_pat_gather_c)
    ELSE
      ! We don't want to receive any data, i.e. the number of cells is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, wrk_p_patch%cells%owner_g, tmp, &
        & wrk_p_patch%cells%loc_index, wrk_p_patch%comm_pat_gather_c)
    ENDIF

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(wrk_p_patch_g%n_patch_edges, wrk_p_patch%edges%owner_g, tmp, &
        & wrk_p_patch%edges%loc_index, wrk_p_patch%comm_pat_gather_e)
    ELSE
      ! We don't want to receive any data, i.e. the number of edges is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, wrk_p_patch%edges%owner_g, tmp, &
        & wrk_p_patch%edges%loc_index, wrk_p_patch%comm_pat_gather_e)
    ENDIF

    IF(p_pe_work == 0) THEN
      CALL setup_comm_pattern(wrk_p_patch_g%n_patch_verts, wrk_p_patch%verts%owner_g, tmp, &
        & wrk_p_patch%verts%loc_index, wrk_p_patch%comm_pat_gather_v)
    ELSE
      ! We don't want to receive any data, i.e. the number of edges is 0
      ! and owner/global index are dummies!
      CALL setup_comm_pattern(0, wrk_p_patch%verts%owner_g, tmp, &
        & wrk_p_patch%verts%loc_index, wrk_p_patch%comm_pat_gather_v)
    ENDIF


    ! For scattering the global fields from PE 0, i.e. owner is 0

    CALL setup_comm_pattern(wrk_p_patch%n_patch_cells, (/ (0,j=1,wrk_p_patch%n_patch_cells) /), &
      & wrk_p_patch%cells%glb_index, tmp, wrk_p_patch%comm_pat_scatter_c)

    CALL setup_comm_pattern(wrk_p_patch%n_patch_edges, (/ (0,j=1,wrk_p_patch%n_patch_edges) /), &
      & wrk_p_patch%edges%glb_index, tmp, wrk_p_patch%comm_pat_scatter_e)

    !-----------------------------------------------------------------------------------------------
    ! Set arrays of divided patch
    !-----------------------------------------------------------------------------------------------

    ! decomp_domain is -1 for invalid locations (at the end of the last strip)
    ! owner_mask is set to .FALSE.

    wrk_p_patch%cells%decomp_domain = -1
    wrk_p_patch%edges%decomp_domain = -1
    wrk_p_patch%verts%decomp_domain = -1
    wrk_p_patch%cells%owner_mask = .false.
    wrk_p_patch%edges%owner_mask = .false.
    wrk_p_patch%verts%owner_mask = .false.

    ! parent_idx/blk, child_idx/blk/id are set later (where it makes sense)
    ! just set all to 0 here

    wrk_p_patch%cells%parent_idx = 0
    wrk_p_patch%cells%parent_blk = 0
    wrk_p_patch%cells%child_idx  = 0
    wrk_p_patch%cells%child_blk  = 0
    wrk_p_patch%cells%child_id   = 0

    !---------------------------------------------------------------------------------------

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%cells%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%cells%glb_index(j)) ! Line  index in global patch

      wrk_p_patch%cells%idx(jl,jb) = jl
      wrk_p_patch%cells%blk(jl,jb) = jb

      DO i=1,i_cell_type

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

      wrk_p_patch%cells%edge_orientation(jl,jb,:) = &
        & wrk_p_patch_g%cells%edge_orientation(jl_g,jb_g,:)
      wrk_p_patch%cells%num_edges(jl,jb)          = wrk_p_patch_g%cells%num_edges(jl_g,jb_g)
      wrk_p_patch%cells%center(jl,jb)             = wrk_p_patch_g%cells%center(jl_g,jb_g)
      wrk_p_patch%cells%area(jl,jb)               = wrk_p_patch_g%cells%area(jl_g,jb_g)
      wrk_p_patch%cells%f_c(jl,jb)                = wrk_p_patch_g%cells%f_c(jl_g,jb_g)
      wrk_p_patch%cells%refin_ctrl(jl,jb)         = wrk_p_patch_g%cells%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%cells%child_id(jl,jb)           = wrk_p_patch_g%cells%child_id(jl_g,jb_g)
      wrk_p_patch%cells%phys_id(jl,jb)            = wrk_p_patch_g%cells%phys_id(jl_g,jb_g)
      wrk_p_patch%cells%decomp_domain(jl,jb)      = flag2_c(wrk_p_patch%cells%glb_index(j))
      wrk_p_patch%cells%owner_mask(jl,jb)         = &
        & wrk_p_patch%cells%owner_g(idx_1d(jl_g,jb_g)) == p_pe_work

    ENDDO

    ! ensure that cells%neighbor_idx lies in the correct grid row along the lateral boundary
    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (wrk_p_patch%cells%refin_ctrl(jl,jb) > 0 .AND.             &
          wrk_p_patch%cells%refin_ctrl(jl,jb) <= grf_bdywidth_c) THEN

        ilc1 = wrk_p_patch%cells%neighbor_idx(jl,jb,1)
        ibc1 = wrk_p_patch%cells%neighbor_blk(jl,jb,1)
        ilc2 = wrk_p_patch%cells%neighbor_idx(jl,jb,2)
        ibc2 = wrk_p_patch%cells%neighbor_blk(jl,jb,2)
        ilc3 = wrk_p_patch%cells%neighbor_idx(jl,jb,3)
        ibc3 = wrk_p_patch%cells%neighbor_blk(jl,jb,3)

        irl0 = wrk_p_patch%cells%refin_ctrl(jl,jb)

        IF (ilc1 > 0 .AND. ibc1 > 0) THEN
          irl1 = wrk_p_patch%cells%refin_ctrl(ilc1,ibc1)
          IF ( irl1 > 0 .AND. ABS(irl0 - irl1) >1 ) THEN
            wrk_p_patch%cells%neighbor_idx(jl,jb,1) = jl
            wrk_p_patch%cells%neighbor_blk(jl,jb,1) = jb
          ENDIF
        ENDIF

        IF (ilc2 > 0 .AND. ibc2 > 0) THEN
          irl2 = wrk_p_patch%cells%refin_ctrl(ilc2,ibc2)
          IF ( irl2 > 0 .AND. ABS(irl0 - irl2) >1 ) THEN
            wrk_p_patch%cells%neighbor_idx(jl,jb,2) = jl
            wrk_p_patch%cells%neighbor_blk(jl,jb,2) = jb
          ENDIF
        ENDIF

        IF (ilc3 > 0 .AND. ibc3 > 0) THEN
          irl3 = wrk_p_patch%cells%refin_ctrl(ilc3,ibc3)
          IF ( irl3 > 0 .AND. ABS(irl0 - irl3) >1 ) THEN
            wrk_p_patch%cells%neighbor_idx(jl,jb,3) = jl
            wrk_p_patch%cells%neighbor_blk(jl,jb,3) = jb
          ENDIF
        ENDIF

      ENDIF
    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%edges%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%edges%glb_index(j)) ! Line  index in global patch

      wrk_p_patch%edges%idx(jl,jb) = jl
      wrk_p_patch%edges%blk(jl,jb) = jb

      DO i=1,2
        CALL get_local_index(wrk_p_patch%cells%loc_index, &
          & wrk_p_patch_g%edges%cell_idx(jl_g,jb_g,i),    &
          & wrk_p_patch_g%edges%cell_blk(jl_g,jb_g,i),    &
          & wrk_p_patch%edges%cell_idx(jl,jb,i),          &
          & wrk_p_patch%edges%cell_blk(jl,jb,i))
      ENDDO

      DO i=1,4
        CALL get_local_index(wrk_p_patch%verts%loc_index, &
          & wrk_p_patch_g%edges%vertex_idx(jl_g,jb_g,i),  &
          & wrk_p_patch_g%edges%vertex_blk(jl_g,jb_g,i),  &
          & wrk_p_patch%edges%vertex_idx(jl,jb,i),        &
          & wrk_p_patch%edges%vertex_blk(jl,jb,i))
      ENDDO

      DO i=1,4

        CALL get_local_index(wrk_p_patch%edges%loc_index, &
          & wrk_p_patch_g%edges%quad_idx(jl_g,jb_g,i),    &
          & wrk_p_patch_g%edges%quad_blk(jl_g,jb_g,i),    &
          & wrk_p_patch%edges%quad_idx(jl,jb,i),          &
          & wrk_p_patch%edges%quad_blk(jl,jb,i))
      ENDDO

      wrk_p_patch%edges%system_orientation(jl,jb)    =&
        & wrk_p_patch_g%edges%system_orientation(jl_g,jb_g)
      wrk_p_patch%edges%quad_orientation(jl,jb,:)    =&
        & wrk_p_patch_g%edges%quad_orientation(jl_g,jb_g,:)

      wrk_p_patch%edges%center(jl,jb)                =&
        & wrk_p_patch_g%edges%center(jl_g,jb_g)
      wrk_p_patch%edges%primal_normal(jl,jb)         =&
        & wrk_p_patch_g%edges%primal_normal(jl_g,jb_g)
      wrk_p_patch%edges%primal_cart_normal(jl,jb)    =&
        & wrk_p_patch_g%edges%primal_cart_normal(jl_g,jb_g)
      wrk_p_patch%edges%dual_normal(jl,jb)           =&
        & wrk_p_patch_g%edges%dual_normal(jl_g,jb_g)
      wrk_p_patch%edges%dual_cart_normal(jl,jb)      =&
        & wrk_p_patch_g%edges%dual_cart_normal(jl_g,jb_g)
      wrk_p_patch%edges%primal_normal_cell(jl,jb,:)  =&
        & wrk_p_patch_g%edges%primal_normal_cell(jl_g,jb_g,:)
      wrk_p_patch%edges%dual_normal_cell(jl,jb,:)    =&
        & wrk_p_patch_g%edges%dual_normal_cell(jl_g,jb_g,:)
      wrk_p_patch%edges%primal_normal_vert(jl,jb,:)  =&
        & wrk_p_patch_g%edges%primal_normal_vert(jl_g,jb_g,:)
      wrk_p_patch%edges%dual_normal_vert(jl,jb,:)    =&
        & wrk_p_patch_g%edges%dual_normal_vert(jl_g,jb_g,:)

      wrk_p_patch%edges%primal_edge_length(jl,jb)    =&
        & wrk_p_patch_g%edges%primal_edge_length(jl_g,jb_g)
      wrk_p_patch%edges%inv_primal_edge_length(jl,jb)=&
        & wrk_p_patch_g%edges%inv_primal_edge_length(jl_g,jb_g)
      wrk_p_patch%edges%dual_edge_length(jl,jb)      =&
        & wrk_p_patch_g%edges%dual_edge_length(jl_g,jb_g)
      wrk_p_patch%edges%inv_dual_edge_length(jl,jb)  =&
        & wrk_p_patch_g%edges%inv_dual_edge_length(jl_g,jb_g)

      wrk_p_patch%edges%edge_vert_length(jl,jb,:)    =&
        & wrk_p_patch_g%edges%edge_vert_length(jl_g,jb_g,:)
      wrk_p_patch%edges%edge_cell_length(jl,jb,:)    =&
        & wrk_p_patch_g%edges%edge_cell_length(jl_g,jb_g,:)
      wrk_p_patch%edges%area_edge(jl,jb)             =&
        & wrk_p_patch_g%edges%area_edge(jl_g,jb_g)
      wrk_p_patch%edges%quad_area(jl,jb)             =&
        & wrk_p_patch_g%edges%quad_area(jl_g,jb_g)
      wrk_p_patch%edges%f_e(jl,jb)                   =&
        & wrk_p_patch_g%edges%f_e(jl_g,jb_g)
      wrk_p_patch%edges%refin_ctrl(jl,jb)            =&
        & wrk_p_patch_g%edges%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%edges%child_id(jl,jb)              =&
        & wrk_p_patch_g%edges%child_id(jl_g,jb_g)
      wrk_p_patch%edges%phys_id(jl,jb)               =wrk_p_patch_g%edges%phys_id(jl_g,jb_g)
      wrk_p_patch%edges%decomp_domain(jl,jb)         =flag2_e(wrk_p_patch%edges%glb_index(j))
      wrk_p_patch%edges%owner_mask(jl,jb)         = &
        & wrk_p_patch%edges%owner_g(idx_1d(jl_g,jb_g))==p_pe_work

      IF (i_cell_type==3) THEN
        wrk_p_patch%edges%inv_vert_vert_length(jl,jb)   =&
          & wrk_p_patch_g%edges%inv_vert_vert_length(jl_g,jb_g)
      ENDIF

      ! ensure that edges%cell_idx lies in the correct grid row along the lateral boundary
      IF (wrk_p_patch%edges%refin_ctrl(jl,jb) >= 2 .AND.            &
          wrk_p_patch%edges%refin_ctrl(jl,jb) <= grf_bdywidth_e) THEN

        ilc1 = wrk_p_patch%edges%cell_idx(jl,jb,1)
        ibc1 = wrk_p_patch%edges%cell_blk(jl,jb,1)
        ilc2 = wrk_p_patch%edges%cell_idx(jl,jb,2)
        ibc2 = wrk_p_patch%edges%cell_blk(jl,jb,2)
        irl1 = wrk_p_patch%cells%refin_ctrl(ilc1,ibc1)
        irl2 = wrk_p_patch%cells%refin_ctrl(ilc2,ibc2)

        IF (MOD(wrk_p_patch%edges%refin_ctrl(jl,jb),2)==0) THEN
          IF (irl1 /= wrk_p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            wrk_p_patch%edges%cell_idx(jl,jb,1) = wrk_p_patch%edges%cell_idx(jl,jb,2)
            wrk_p_patch%edges%cell_blk(jl,jb,1) = wrk_p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= wrk_p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            wrk_p_patch%edges%cell_idx(jl,jb,2) = wrk_p_patch%edges%cell_idx(jl,jb,1)
            wrk_p_patch%edges%cell_blk(jl,jb,2) = wrk_p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ELSE
          IF (irl1 /= wrk_p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl1 /= wrk_p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            wrk_p_patch%edges%cell_idx(jl,jb,1) = wrk_p_patch%edges%cell_idx(jl,jb,2)
            wrk_p_patch%edges%cell_blk(jl,jb,1) = wrk_p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= wrk_p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl2 /= wrk_p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            wrk_p_patch%edges%cell_idx(jl,jb,2) = wrk_p_patch%edges%cell_idx(jl,jb,1)
            wrk_p_patch%edges%cell_blk(jl,jb,2) = wrk_p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ENDIF

      ENDIF

    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%verts%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%verts%glb_index(j)) ! Line  index in global patch

      wrk_p_patch%verts%idx(jl,jb) = jl
      wrk_p_patch%verts%blk(jl,jb) = jb

      DO i=1,9-i_cell_type

        CALL get_local_index(wrk_p_patch%verts%loc_index, &
          & wrk_p_patch_g%verts%neighbor_idx(jl_g,jb_g,i),&
          & wrk_p_patch_g%verts%neighbor_blk(jl_g,jb_g,i),&
          & wrk_p_patch%verts%neighbor_idx(jl,jb,i),      &
          & wrk_p_patch%verts%neighbor_blk(jl,jb,i))

        CALL get_local_index(wrk_p_patch%edges%loc_index, &
          & wrk_p_patch_g%verts%edge_idx(jl_g,jb_g,i),    &
          & wrk_p_patch_g%verts%edge_blk(jl_g,jb_g,i),    &
          & wrk_p_patch%verts%edge_idx(jl,jb,i),          &
          & wrk_p_patch%verts%edge_blk(jl,jb,i))

        CALL get_local_index(wrk_p_patch%cells%loc_index, &
          & wrk_p_patch_g%verts%cell_idx(jl_g,jb_g,i),    &
          & wrk_p_patch_g%verts%cell_blk(jl_g,jb_g,i),    &
          & wrk_p_patch%verts%cell_idx(jl,jb,i),          &
          & wrk_p_patch%verts%cell_blk(jl,jb,i))
      ENDDO

      wrk_p_patch%verts%edge_orientation(jl,jb,:) = &
        & wrk_p_patch_g%verts%edge_orientation(jl_g,jb_g,:)
      wrk_p_patch%verts%num_edges(jl,jb)          = wrk_p_patch_g%verts%num_edges(jl_g,jb_g)
      wrk_p_patch%verts%vertex(jl,jb)             = wrk_p_patch_g%verts%vertex(jl_g,jb_g)
      wrk_p_patch%verts%dual_area(jl,jb)          = wrk_p_patch_g%verts%dual_area(jl_g,jb_g)
      wrk_p_patch%verts%f_v(jl,jb)                = wrk_p_patch_g%verts%f_v(jl_g,jb_g)
      wrk_p_patch%verts%refin_ctrl(jl,jb)         = wrk_p_patch_g%verts%refin_ctrl(jl_g,jb_g)
      wrk_p_patch%verts%decomp_domain(jl,jb)      = flag2_v(wrk_p_patch%verts%glb_index(j))
      wrk_p_patch%verts%owner_mask(jl,jb)         = &
        & wrk_p_patch%verts%owner_g(idx_1d(jl_g,jb_g)) == p_pe_work

    ENDDO

    DEALLOCATE(flag_c, flag_e, flag_v, flag2_c, flag2_e, flag2_v, lcount_c, lcount_e, lcount_v)

  END SUBROUTINE divide_patch

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets parent_idx/blk in child and child_idx/blk in parent patches.

  SUBROUTINE set_parent_child_relations(p_pp, p_pc, p_pp_g, p_pc_g)

    TYPE(t_patch), INTENT(INOUT) :: p_pp   !> divided parent patch
    TYPE(t_patch), INTENT(INOUT) :: p_pc   !> divided child patch
    TYPE(t_patch), INTENT(INOUT) :: p_pp_g !> global parent patch
    TYPE(t_patch), INTENT(INOUT) :: p_pc_g !> global child patch

    INTEGER :: i, j, jl, jb, jl_g, jb_g, jc, jc_g, jp, jp_g

    ! Attention:
    ! Only inner cells/edges get a valid parent or child index,
    ! indexes for boundary cells/edges are not set.
    ! Therefore when these indexes are used, the code must assure that they are
    ! used only for inner cells/edges!
    ! The main reason for this is that - depending on the number of ghost rows -
    ! there are cells/edges in the parent boundary with missing childs (n_ghost_rows==1)
    ! or cells/edges in the child boundary with missing parents (n_ghost_rows>2).

    ! Set child indices in parent ...

    ! ... cells

    DO j = 1, p_pp%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pp%cells%decomp_domain(jl,jb)>0) CYCLE ! only inner cells get a valid parent index

      jb_g = blk_no(p_pp%cells%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(p_pp%cells%glb_index(j)) ! Line  index in global patch

      DO i= 1, 4
        jc_g = idx_1d(p_pp_g%cells%child_idx(jl_g,jb_g,i),p_pp_g%cells%child_blk(jl_g,jb_g,i))
        IF(jc_g<1 .OR. jc_g>p_pc%n_patch_cells_g) &
          & CALL finish('set_parent_child_relations','Invalid cell child index in global parent')
        jc = p_pc%cells%loc_index(jc_g)
        IF(jc <= 0) &
          & CALL finish('set_parent_child_relations','cell child index outside child domain')
        p_pp%cells%child_blk(jl,jb,i) = blk_no(jc)
        p_pp%cells%child_idx(jl,jb,i) = idx_no(jc)
      ENDDO

      ! Check child_id, it must be p_pc%id

      IF(p_pp_g%cells%child_id(jl_g,jb_g) /= p_pc%id) &
        & CALL finish('set_parent_child_relations','Invalid cell child ID in global parent')

      p_pp%cells%child_id(jl,jb) = p_pc%id

    ENDDO

    ! ... edges

    DO j = 1, p_pp%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pp%edges%decomp_domain(jl,jb)>1) CYCLE ! only inner edges get a valid parent index

      jb_g = blk_no(p_pp%edges%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(p_pp%edges%glb_index(j)) ! Line  index in global patch

      DO i= 1, 4

        IF(i==4 .AND. p_pp_g%edges%refin_ctrl(jl_g,jb_g) == -1) THEN
          p_pp%edges%child_blk(jl,jb,i) = blk_no(0)
          p_pp%edges%child_idx(jl,jb,i) = idx_no(0)
          CYCLE
        ENDIF

        jc_g = idx_1d(p_pp_g%edges%child_idx(jl_g,jb_g,i),p_pp_g%edges%child_blk(jl_g,jb_g,i))
        IF(ABS(jc_g)<1 .OR. ABS(jc_g)>p_pc%n_patch_edges_g) &
          & CALL finish('set_parent_child_relations','Inv. edge child index in global parent')
        jc = p_pc%edges%loc_index(ABS(jc_g))
        IF(jc <= 0) &
          & CALL finish('set_parent_child_relations','edge child index outside child domain')
        p_pp%edges%child_blk(jl,jb,i) = blk_no(jc)
        p_pp%edges%child_idx(jl,jb,i) = SIGN(idx_no(jc),jc_g)
      ENDDO

      ! Check child_id, it must be p_pc%id

      IF(p_pp_g%edges%child_id(jl_g,jb_g) /= p_pc%id) &
        & CALL finish('set_parent_child_relations','Invalid edge child ID in global parent')

      p_pp%edges%child_id(jl,jb) = p_pc%id

    ENDDO

    ! Set parent indices in child ...

    ! ... cells

    DO j = 1, p_pc%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pc%cells%decomp_domain(jl,jb)>0) CYCLE ! only inner cells get a valid parent index

      jb_g = blk_no(p_pc%cells%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(p_pc%cells%glb_index(j)) ! Line  index in global patch

      jp_g = idx_1d(p_pc_g%cells%parent_idx(jl_g,jb_g),p_pc_g%cells%parent_blk(jl_g,jb_g))
      IF(jp_g<1 .OR. jp_g>p_pp%n_patch_cells_g) &
        & CALL finish('set_parent_child_relations','Inv. cell parent index in global child')

      jp = p_pp%cells%loc_index(jp_g)
      IF(jp <= 0) &
        & CALL finish('set_parent_child_relations','cell parent index outside parent domain')
      p_pc%cells%parent_blk(jl,jb) = blk_no(jp)
      p_pc%cells%parent_idx(jl,jb) = idx_no(jp)

    ENDDO

    ! ... edges

    DO j = 1, p_pc%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF(p_pc%edges%decomp_domain(jl,jb)>1) CYCLE ! only inner edges get a valid parent index

      jb_g = blk_no(p_pc%edges%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(p_pc%edges%glb_index(j)) ! Line  index in global patch

      jp_g = idx_1d(p_pc_g%edges%parent_idx(jl_g,jb_g),p_pc_g%edges%parent_blk(jl_g,jb_g))
      IF(jp_g<1 .OR. jp_g>p_pp%n_patch_edges_g) &
        & CALL finish('set_parent_child_relations','Inv. edge parent index in global child')

      jp = p_pp%edges%loc_index(jp_g)
      IF(jp <= 0) &
        & CALL finish('set_parent_child_relations','edge parent index outside parent domain')
      p_pc%edges%parent_blk(jl,jb) = blk_no(jp)
      p_pc%edges%parent_idx(jl,jb) = idx_no(jp)

    ENDDO

  END SUBROUTINE set_parent_child_relations

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets up communication patterns between global and local parent patches

  SUBROUTINE set_glb_loc_comm(p_pglb, p_ploc, i_chidx)
    TYPE(t_patch), INTENT(IN)    :: p_pglb !> global parent
    TYPE(t_patch), INTENT(INOUT) :: p_ploc !> local parent
    INTEGER, INTENT(IN) :: i_chidx

    INTEGER, ALLOCATABLE :: owner(:)
    INTEGER :: j, js, je

    ! Please note:
    ! For creating communication patterns for different amount of data to be transferred
    ! (e.g. only the boundary interpolation zone), create a new communication pattern
    ! by copying the code below and adjusting the limits in the calculation of js/je.

    !-----------------------------------------------------------------------------------------------

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
  SUBROUTINE setup_comm_cpy_interpolation()

    INTEGER :: j, jc, jb

    INTEGER, ALLOCATABLE :: parent_index(:,:), owner(:), glb_index(:)
    INTEGER ici1, icb1, ici2, icb2, ici3, icb3, ici4, icb4, jl_g, jb_g
    INTEGER i_chidx, i_startblk, i_endblk, i_startidx, i_endidx

    TYPE(t_grid_cells), POINTER::p_gcp => NULL()

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(p_nprocs == 1 .or. p_pe == p_test_pe) &
      & CALL finish('setup_comm_cpy_interpolation','must not be called in a single CPU run')

    i_chidx = wrk_p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Assign the global parent index (1D) to every cell in the interpolation zone
    ! of the global patch (in parent_index).
    ! This is done in the same way as in mo_hierarchy_management/interpolate_tendencies.

    ALLOCATE(parent_index(nproma,wrk_p_patch%n_patch_cells_g/nproma+1)) ! spans GLOBAL patch
    parent_index = 0

    p_gcp => wrk_p_parent_patch_g%cells

    ! Start and end blocks for which interpolation is needed
    i_startblk = p_gcp%start_blk(grf_bdyintp_start_c,i_chidx)
    i_endblk   = p_gcp%end_blk(grf_bdyintp_end_c,i_chidx)

    DO jb =  i_startblk, i_endblk

      CALL get_indices_c(wrk_p_parent_patch_g, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)
      DO jc = i_startidx, i_endidx

        ici1 = p_gcp%child_idx(jc,jb,1)
        icb1 = p_gcp%child_blk(jc,jb,1)
        ici2 = p_gcp%child_idx(jc,jb,2)
        icb2 = p_gcp%child_blk(jc,jb,2)
        ici3 = p_gcp%child_idx(jc,jb,3)
        icb3 = p_gcp%child_blk(jc,jb,3)
        ici4 = p_gcp%child_idx(jc,jb,4)
        icb4 = p_gcp%child_blk(jc,jb,4)

        parent_index(ici1,icb1) = idx_1d(jc,jb)
        parent_index(ici2,icb2) = idx_1d(jc,jb)
        parent_index(ici3,icb3) = idx_1d(jc,jb)
        parent_index(ici4,icb4) = idx_1d(jc,jb)

      ENDDO
    ENDDO

    ! Now, for our local child patch, gather which cells receive values from which parent cell

    ALLOCATE(glb_index(wrk_p_patch%n_patch_cells))
    ALLOCATE(owner(wrk_p_patch%n_patch_cells))

    DO j = 1,wrk_p_patch%n_patch_cells
      jl_g = idx_no(wrk_p_patch%cells%glb_index(j))
      jb_g = blk_no(wrk_p_patch%cells%glb_index(j))
      IF(parent_index(jl_g,jb_g) /= 0) THEN
        owner(j) = wrk_p_parent_patch%cells%owner_g(ABS(parent_index(jl_g,jb_g)))
        glb_index(j) = parent_index(jl_g,jb_g)
      ELSE
        owner(j) = -1
        glb_index(j) = -1
      ENDIF
    ENDDO

    ! Set up communication pattern

    CALL setup_comm_pattern(wrk_p_patch%n_patch_cells, owner, glb_index,  &
      & wrk_p_parent_patch%cells%loc_index, &
      & wrk_p_patch%comm_pat_interpolation_c)

    DEALLOCATE(parent_index, owner, glb_index)

  END SUBROUTINE setup_comm_cpy_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for grf interpolation.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_grf_interpolation()

    INTEGER :: j, n, jc, je, jb, jl_g, jb_g

    INTEGER, ALLOCATABLE :: parent_index(:,:), owner(:), glb_index(:)
    INTEGER ici1, icb1, ici2, icb2, ici3, icb3, ici4, icb4
    INTEGER i_chidx, i_startblk, i_endblk, i_startidx, i_endidx

    TYPE(t_patch), POINTER:: ptr_pp

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(p_nprocs == 1 .or. p_pe == p_test_pe) &
      & CALL finish('setup_comm_grf_interpolation','must not be called in a single CPU run')

    i_chidx = wrk_p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Assign the global parent index (1D) to every cell in the interpolation zone
    ! of the global patch (in parent_index).
    ! This is done in the same way as in mo_grf_interpolation/interpol_scal_grf.
    ! parent_index also gets the info about the number of the child in the parent.

    ALLOCATE(parent_index(nproma,wrk_p_patch%n_patch_cells_g/nproma+1)) ! spans GLOBAL patch
    parent_index = 0

    ptr_pp => wrk_p_parent_patch_g

    i_startblk = ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(grf_bdyintp_end_c,i_chidx)

    DO jb =  i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)

      DO jc = i_startidx, i_endidx

        ici1 = ptr_pp%cells%child_idx(jc,jb,1)
        icb1 = ptr_pp%cells%child_blk(jc,jb,1)
        ici2 = ptr_pp%cells%child_idx(jc,jb,2)
        icb2 = ptr_pp%cells%child_blk(jc,jb,2)
        ici3 = ptr_pp%cells%child_idx(jc,jb,3)
        icb3 = ptr_pp%cells%child_blk(jc,jb,3)
        ici4 = ptr_pp%cells%child_idx(jc,jb,4)
        icb4 = ptr_pp%cells%child_blk(jc,jb,4)

        parent_index(ici1,icb1) = idx_1d(jc,jb)*4 + 0
        parent_index(ici2,icb2) = idx_1d(jc,jb)*4 + 1
        parent_index(ici3,icb3) = idx_1d(jc,jb)*4 + 2
        parent_index(ici4,icb4) = idx_1d(jc,jb)*4 + 3

      ENDDO
    ENDDO

    ! Now, for our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(wrk_p_patch%n_patch_cells))
    ALLOCATE(owner(wrk_p_patch%n_patch_cells))

    DO n = 1, 4

      DO j = 1,wrk_p_patch%n_patch_cells
        jl_g = idx_no(wrk_p_patch%cells%glb_index(j))
        jb_g = blk_no(wrk_p_patch%cells%glb_index(j))
        IF(parent_index(jl_g,jb_g) /= 0 .and.     &
          & MOD(parent_index(jl_g,jb_g),4) == n-1) THEN
          owner(j) = wrk_p_parent_patch%cells%owner_g(parent_index(jl_g,jb_g)/4)
          glb_index(j) = parent_index(jl_g,jb_g)/4
        ELSE
          owner(j) = -1
          glb_index(j) = -1
        ENDIF
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(wrk_p_patch%n_patch_cells, owner, glb_index,  &
        & wrk_p_parent_patch%cells%loc_index, &
        & wrk_p_patch%comm_pat_interpol_scal_grf(n))

    ENDDO

    DEALLOCATE(parent_index, owner, glb_index)


    !--------------------------------------------------------------------
    ! Edges

    ! Assign the global parent index (1D) to every edge in the interpolation zone
    ! of the global patch (in parent_index).
    ! This is done in the same way as in mo_grf_interpolation/interpol_vec_grf.
    ! parent_index also gets the info about the number of the child in the parent.

    ALLOCATE(parent_index(nproma,wrk_p_patch%n_patch_edges_g/nproma+1))
    parent_index = 0

    ptr_pp => wrk_p_parent_patch_g

    ! Start and end blocks for which vector interpolation is needed
    i_startblk = ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx)
    i_endblk   = ptr_pp%edges%end_blk(grf_bdyintp_end_e,i_chidx)

    DO jb =  i_startblk, i_endblk

      CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_bdyintp_start_e, grf_bdyintp_end_e, i_chidx)

      DO je = i_startidx, i_endidx

        ici1 = ptr_pp%edges%child_idx(je,jb,1)
        icb1 = ptr_pp%edges%child_blk(je,jb,1)
        ici2 = ptr_pp%edges%child_idx(je,jb,2)
        icb2 = ptr_pp%edges%child_blk(je,jb,2)
        ici3 = ABS(ptr_pp%edges%child_idx(je,jb,3))
        icb3 = ptr_pp%edges%child_blk(je,jb,3)

        parent_index(ici1,icb1) = idx_1d(je,jb)*4 + 0
        parent_index(ici2,icb2) = idx_1d(je,jb)*4 + 1
        parent_index(ici3,icb3) = idx_1d(je,jb)*4 + 2

        IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
          ici4 = ABS(ptr_pp%edges%child_idx(je,jb,4))
          icb4 = ptr_pp%edges%child_blk(je,jb,4)
          parent_index(ici4,icb4) = idx_1d(je,jb)*4 + 3
        ENDIF

      ENDDO
    ENDDO

    ! Now, for our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(wrk_p_patch%n_patch_edges))
    ALLOCATE(owner(wrk_p_patch%n_patch_edges))

    DO n = 1, 4

      DO j = 1,wrk_p_patch%n_patch_edges
        jl_g = idx_no(wrk_p_patch%edges%glb_index(j))
        jb_g = blk_no(wrk_p_patch%edges%glb_index(j))
        IF(parent_index(jl_g,jb_g) /= 0 .and.     &
          & MOD(parent_index(jl_g,jb_g),4) == n-1) THEN
          owner(j) = wrk_p_parent_patch%edges%owner_g(parent_index(jl_g,jb_g)/4)
          glb_index(j) = parent_index(jl_g,jb_g)/4
        ELSE
          owner(j) = -1
          glb_index(j) = -1
        ENDIF
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(wrk_p_patch%n_patch_edges, owner, glb_index,  &
        & wrk_p_parent_patch%edges%loc_index, &
        & wrk_p_patch%comm_pat_interpol_vec_grf(n))

    ENDDO

    DEALLOCATE(parent_index, owner, glb_index)

  END SUBROUTINE setup_comm_grf_interpolation

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a communication pattern for ubc interpolation.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE setup_comm_ubc_interpolation()

    INTEGER :: j, n, jc, je, jb, jl_g, jb_g

    INTEGER, ALLOCATABLE :: parent_index(:,:), owner(:), glb_index(:)
    INTEGER ici1, icb1, ici2, icb2, ici3, icb3, ici4, icb4
    INTEGER i_chidx, i_startblk, i_endblk, i_startidx, i_endidx

    TYPE(t_patch), POINTER:: ptr_pp

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(p_nprocs == 1 .or. p_pe == p_test_pe) &
      & CALL finish('setup_comm_ubc_interpolation','must not be called in a single CPU run')

    i_chidx = wrk_p_patch%parent_child_index

    !--------------------------------------------------------------------
    ! Cells

    ! Assign the global parent index (1D) to every cell in the interpolation zone
    ! of the global patch (in parent_index).
    ! This is done in the same way as in mo_grf_interpolation/interpol_scal_grf.
    ! parent_index also gets the info about the number of the child in the parent.

    ALLOCATE(parent_index(nproma,wrk_p_patch%n_patch_cells_g/nproma+1)) ! spans GLOBAL patch
    parent_index = 0

    ptr_pp => wrk_p_parent_patch_g

    i_startblk = ptr_pp%cells%start_blk(grf_nudgintp_start_c,i_chidx)
    i_endblk   = ptr_pp%cells%end_blk(min_rlcell_int,i_chidx)

    DO jb =  i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_nudgintp_start_c, min_rlcell_int, i_chidx)

      DO jc = i_startidx, i_endidx

        ici1 = ptr_pp%cells%child_idx(jc,jb,1)
        icb1 = ptr_pp%cells%child_blk(jc,jb,1)
        ici2 = ptr_pp%cells%child_idx(jc,jb,2)
        icb2 = ptr_pp%cells%child_blk(jc,jb,2)
        ici3 = ptr_pp%cells%child_idx(jc,jb,3)
        icb3 = ptr_pp%cells%child_blk(jc,jb,3)
        ici4 = ptr_pp%cells%child_idx(jc,jb,4)
        icb4 = ptr_pp%cells%child_blk(jc,jb,4)

        parent_index(ici1,icb1) = idx_1d(jc,jb)*4 + 0
        parent_index(ici2,icb2) = idx_1d(jc,jb)*4 + 1
        parent_index(ici3,icb3) = idx_1d(jc,jb)*4 + 2
        parent_index(ici4,icb4) = idx_1d(jc,jb)*4 + 3

      ENDDO
    ENDDO

    ! Now, for our local child patch, gather which cells receive values from which parent cell
    ! This is done once for every of the four child cells

    ALLOCATE(glb_index(wrk_p_patch%n_patch_cells))
    ALLOCATE(owner(wrk_p_patch%n_patch_cells))

    DO n = 1, 4

      DO j = 1,wrk_p_patch%n_patch_cells
        jl_g = idx_no(wrk_p_patch%cells%glb_index(j))
        jb_g = blk_no(wrk_p_patch%cells%glb_index(j))
        IF(parent_index(jl_g,jb_g) /= 0 .and.     &
          & MOD(parent_index(jl_g,jb_g),4) == n-1) THEN
          owner(j) = wrk_p_parent_patch%cells%owner_g(parent_index(jl_g,jb_g)/4)
          glb_index(j) = parent_index(jl_g,jb_g)/4
        ELSE
          owner(j) = -1
          glb_index(j) = -1
        ENDIF
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(wrk_p_patch%n_patch_cells, owner, glb_index,  &
        & wrk_p_parent_patch%cells%loc_index, &
        & wrk_p_patch%comm_pat_interpol_scal_ubc(n))

    ENDDO

    DEALLOCATE(parent_index, owner, glb_index)


    !--------------------------------------------------------------------
    ! Edges

    ! Assign the global parent index (1D) to every edge in the interpolation zone
    ! of the global patch (in parent_index).
    ! This is done in the same way as in mo_grf_interpolation/interpol_vec_grf.
    ! parent_index also gets the info about the number of the child in the parent.

    ALLOCATE(parent_index(nproma,wrk_p_patch%n_patch_edges_g/nproma+1))
    parent_index = 0

    ptr_pp => wrk_p_parent_patch_g

    ! Start and end blocks for which vector interpolation is needed
    i_startblk = ptr_pp%edges%start_blk(grf_nudgintp_start_e,i_chidx)
    i_endblk   = ptr_pp%edges%end_blk(min_rledge_int,i_chidx)

    DO jb =  i_startblk, i_endblk

      CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, grf_nudgintp_start_e, min_rledge_int, i_chidx)

      DO je = i_startidx, i_endidx

        ici1 = ptr_pp%edges%child_idx(je,jb,1)
        icb1 = ptr_pp%edges%child_blk(je,jb,1)
        ici2 = ptr_pp%edges%child_idx(je,jb,2)
        icb2 = ptr_pp%edges%child_blk(je,jb,2)
        ici3 = ABS(ptr_pp%edges%child_idx(je,jb,3))
        icb3 = ptr_pp%edges%child_blk(je,jb,3)

        parent_index(ici1,icb1) = idx_1d(je,jb)*4 + 0
        parent_index(ici2,icb2) = idx_1d(je,jb)*4 + 1
        parent_index(ici3,icb3) = idx_1d(je,jb)*4 + 2

        IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
          ici4 = ABS(ptr_pp%edges%child_idx(je,jb,4))
          icb4 = ptr_pp%edges%child_blk(je,jb,4)
          parent_index(ici4,icb4) = idx_1d(je,jb)*4 + 3
        ENDIF

      ENDDO
    ENDDO

    ! Now, for our local child patch, gather which edges receive values from which parent edge
    ! This is done once for every of the four child edges

    ALLOCATE(glb_index(wrk_p_patch%n_patch_edges))
    ALLOCATE(owner(wrk_p_patch%n_patch_edges))

    DO n = 1, 4

      DO j = 1,wrk_p_patch%n_patch_edges
        jl_g = idx_no(wrk_p_patch%edges%glb_index(j))
        jb_g = blk_no(wrk_p_patch%edges%glb_index(j))
        IF(parent_index(jl_g,jb_g) /= 0 .and.     &
          & MOD(parent_index(jl_g,jb_g),4) == n-1) THEN
          owner(j) = wrk_p_parent_patch%edges%owner_g(parent_index(jl_g,jb_g)/4)
          glb_index(j) = parent_index(jl_g,jb_g)/4
        ELSE
          owner(j) = -1
          glb_index(j) = -1
        ENDIF
      ENDDO

      ! Set up communication pattern

      CALL setup_comm_pattern(wrk_p_patch%n_patch_edges, owner, glb_index,  &
        & wrk_p_parent_patch%edges%loc_index, &
        & wrk_p_patch%comm_pat_interpol_vec_ubc(n))

    ENDDO

    DEALLOCATE(parent_index, owner, glb_index)

  END SUBROUTINE setup_comm_ubc_interpolation

  !-------------------------------------------------------------------------
  !>
  !! Divides all variable in int_state.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE divide_int_state()

    ! Local scalars:

    INTEGER :: j, jb, jl, jb_g, jl_g, i, nincr

    !-----------------------------------------------------------------------

    ! This routine must not be called in a single CPU run
    IF(p_nprocs == 1 .or. p_pe == p_test_pe) &
      & CALL finish('divide_int_state','must not be called in a single CPU run')


    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%cells%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%cells%glb_index(j)) ! Line  index in global patch

      wrk_int_state_out%e_inn_c(jl,:,jb)          = wrk_int_state_in%e_inn_c(jl_g,:,jb_g)
      wrk_int_state_out%verts_aw_cells(jl,:,jb)   = wrk_int_state_in%verts_aw_cells(jl_g,:,jb_g)

      IF (i_cell_type == 6) THEN
        wrk_int_state_out%e_aw_c(jl,:,jb)         = wrk_int_state_in%e_aw_c(jl_g,:,jb_g)
        wrk_int_state_out%r_aw_c(jl,:,jb)         = wrk_int_state_in%r_aw_c(jl_g,:,jb_g)
        wrk_int_state_out%hex_north(jl,:,jb)      = wrk_int_state_in%hex_north(jl_g,:,jb_g)
        wrk_int_state_out%hex_east(jl,:,jb)       = wrk_int_state_in%hex_east(jl_g,:,jb_g)
      ENDIF

      IF (i_cell_type == 3) THEN

        wrk_int_state_out%e_bln_c_s(jl,:,jb)      = wrk_int_state_in%e_bln_c_s(jl_g,:,jb_g)
        wrk_int_state_out%e_bln_c_u(jl,:,jb)      = wrk_int_state_in%e_bln_c_u(jl_g,:,jb_g)
        wrk_int_state_out%e_bln_c_v(jl,:,jb)      = wrk_int_state_in%e_bln_c_v(jl_g,:,jb_g)
        wrk_int_state_out%c_bln_avg(jl,:,jb)      = wrk_int_state_in%c_bln_avg(jl_g,:,jb_g)
        wrk_int_state_out%nudgecoeff_c(jl,jb)     = wrk_int_state_in%nudgecoeff_c(jl_g,jb_g)

        DO i=1,rbf_vec_dim_c
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%rbf_vec_idx_c(i,jl_g,jb_g), &
            & wrk_int_state_in%rbf_vec_blk_c(i,jl_g,jb_g), &
            & wrk_int_state_out%rbf_vec_idx_c(i,jl,jb), &
            & wrk_int_state_out%rbf_vec_blk_c(i,jl,jb))
        ENDDO
        wrk_int_state_out%rbf_vec_stencil_c(jl, jb)  =&
          & wrk_int_state_in%rbf_vec_stencil_c(jl_g, jb_g)
        wrk_int_state_out%rbf_vec_coeff_c(:,:,jl,jb) = &
          & wrk_int_state_in%rbf_vec_coeff_c(:,:,jl_g,jb_g)

        DO i=1,rbf_c2grad_dim
          CALL get_local_index(wrk_p_patch%cells%loc_index, &
            & wrk_int_state_in%rbf_c2grad_idx(i,jl_g,jb_g), &
            & wrk_int_state_in%rbf_c2grad_blk(i,jl_g,jb_g), &
            & wrk_int_state_out%rbf_c2grad_idx(i,jl,jb), &
            & wrk_int_state_out%rbf_c2grad_blk(i,jl,jb))
        ENDDO
        wrk_int_state_out%rbf_c2grad_coeff(:,:,jl,jb) = &
          & wrk_int_state_in%rbf_c2grad_coeff(:,:,jl_g,jb_g)

        !
        ! quadrature points on triangles
        !
        wrk_int_state_out%gquad%qpts_tri_l(jl,jb)   =          &
          & wrk_int_state_in%gquad%qpts_tri_l(jl_g,jb_g)
        wrk_int_state_out%gquad%qpts_tri_q(jl,jb,:)   =        &
          & wrk_int_state_in%gquad%qpts_tri_q(jl_g,jb_g,:)
        wrk_int_state_out%gquad%qpts_tri_c(jl,jb,:)   =        &
          & wrk_int_state_in%gquad%qpts_tri_c(jl_g,jb_g,:)
        wrk_int_state_out%gquad%weights_tri_q(:)   =           &
          & wrk_int_state_in%gquad%weights_tri_q(:)
        wrk_int_state_out%gquad%weights_tri_c(:)   =           &
          & wrk_int_state_in%gquad%weights_tri_c(:)
      ENDIF

      IF( ltransport .OR. iequations == 3) THEN

        wrk_int_state_out%lsq_lin%lsq_dim_stencil(jl, jb)   = &
          & wrk_int_state_in%lsq_lin%lsq_dim_stencil(jl_g, jb_g)

        DO i=1,lsq_lin_set%dim_c
          CALL get_local_index(wrk_p_patch%cells%loc_index, &
            & wrk_int_state_in%lsq_lin%lsq_idx_c(jl_g,jb_g,i), &
            & wrk_int_state_in%lsq_lin%lsq_blk_c(jl_g,jb_g,i), &
            & wrk_int_state_out%lsq_lin%lsq_idx_c(jl,jb,i), &
            & wrk_int_state_out%lsq_lin%lsq_blk_c(jl,jb,i))
        ENDDO

        wrk_int_state_out%lsq_lin%lsq_weights_c(jl, :, jb)   = &
          & wrk_int_state_in%lsq_lin%lsq_weights_c(jl_g, :, jb_g)
        wrk_int_state_out%lsq_lin%lsq_qtmat_c(jl, :, :, jb)  = &
          & wrk_int_state_in%lsq_lin%lsq_qtmat_c(jl_g, :, :, jb_g)
        wrk_int_state_out%lsq_lin%lsq_rmat_rdiag_c(jl, :, jb)= &
          & wrk_int_state_in%lsq_lin%lsq_rmat_rdiag_c(jl_g, :, jb_g)
        wrk_int_state_out%lsq_lin%lsq_rmat_utri_c(jl, :, jb) = &
          & wrk_int_state_in%lsq_lin%lsq_rmat_utri_c(jl_g, :, jb_g)

        wrk_int_state_out%lsq_lin%lsq_moments(jl, jb, :)     = &
          & wrk_int_state_in%lsq_lin%lsq_moments(jl_g, jb_g, :)
        wrk_int_state_out%lsq_lin%lsq_moments_hat(jl,jb,:,:) = &
          & wrk_int_state_in%lsq_lin%lsq_moments_hat(jl_g,jb_g,:,:)

        !
        ! high order lsq
        !
        wrk_int_state_out%lsq_high%lsq_dim_stencil(jl, jb)   = &
          & wrk_int_state_in%lsq_high%lsq_dim_stencil(jl_g, jb_g)

        DO i=1,lsq_high_set%dim_c
          CALL get_local_index(wrk_p_patch%cells%loc_index, &
            & wrk_int_state_in%lsq_high%lsq_idx_c(jl_g,jb_g,i), &
            & wrk_int_state_in%lsq_high%lsq_blk_c(jl_g,jb_g,i), &
            & wrk_int_state_out%lsq_high%lsq_idx_c(jl,jb,i), &
            & wrk_int_state_out%lsq_high%lsq_blk_c(jl,jb,i))
        ENDDO

        wrk_int_state_out%lsq_high%lsq_weights_c(jl, :, jb)   = &
          & wrk_int_state_in%lsq_high%lsq_weights_c(jl_g, :, jb_g)
        wrk_int_state_out%lsq_high%lsq_qtmat_c(jl, :, :, jb)  = &
          & wrk_int_state_in%lsq_high%lsq_qtmat_c(jl_g, :, :, jb_g)
        wrk_int_state_out%lsq_high%lsq_rmat_rdiag_c(jl, :, jb)= &
          & wrk_int_state_in%lsq_high%lsq_rmat_rdiag_c(jl_g, :, jb_g)
        wrk_int_state_out%lsq_high%lsq_rmat_utri_c(jl, :, jb) = &
          & wrk_int_state_in%lsq_high%lsq_rmat_utri_c(jl_g, :, jb_g)

        wrk_int_state_out%lsq_high%lsq_moments(jl, jb, :)     = &
          & wrk_int_state_in%lsq_high%lsq_moments(jl_g, jb_g, :)
        wrk_int_state_out%lsq_high%lsq_moments_hat(jl,jb,:,:) = &
          & wrk_int_state_in%lsq_high%lsq_moments_hat(jl_g,jb_g,:,:)
      ENDIF

      wrk_int_state_out%geofac_div(jl, :, jb)             = &
        & wrk_int_state_in%geofac_div(jl_g, :, jb_g)
      wrk_int_state_out%geofac_n2s(jl, :, jb)             = &
        & wrk_int_state_in%geofac_n2s(jl_g, :, jb_g)
      wrk_int_state_out%geofac_grg(jl, :, jb, :)          = &
        & wrk_int_state_in%geofac_grg(jl_g, :, jb_g, :)

      wrk_int_state_out%cart_cell_coord(jl, jb, :)        = &
        & wrk_int_state_in%cart_cell_coord(jl_g, jb_g, :)
      wrk_int_state_out%primal_normal_ec(jl, jb, :, :)    = &
        & wrk_int_state_in%primal_normal_ec(jl_g, jb_g, :, :)
      wrk_int_state_out%edge_cell_length(jl, jb, :)       = &
        & wrk_int_state_in%edge_cell_length(jl_g, jb_g, :)

    ENDDO

    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%edges%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%edges%glb_index(j)) ! Line  index in global patch

      wrk_int_state_out%c_lin_e(jl,:,jb)         = wrk_int_state_in%c_lin_e(jl_g,:,jb_g)
      wrk_int_state_out%v_1o2_e(jl,:,jb)         = wrk_int_state_in%v_1o2_e(jl_g,:,jb_g)

      IF( i_cell_type==6 ) THEN

        wrk_int_state_out%tria_aw_rhom(jl,:,jb)  = wrk_int_state_in%tria_aw_rhom  (jl_g,:,jb_g)
        wrk_int_state_out%heli_coeff  (:,jl,jb)  = wrk_int_state_in%heli_coeff    (:,jl_g,jb_g)
        wrk_int_state_out%dir_gradhux_c1(:,jl,jb)= wrk_int_state_in%dir_gradhux_c1(:,jl_g,jb_g)
        wrk_int_state_out%dir_gradhux_c2(:,jl,jb)= wrk_int_state_in%dir_gradhux_c2(:,jl_g,jb_g)
        wrk_int_state_out%strain_def_c1(:,jl,jb) = wrk_int_state_in%strain_def_c1(:,jl_g,jb_g)
        wrk_int_state_out%strain_def_c2(:,jl,jb) = wrk_int_state_in%strain_def_c2(:,jl_g,jb_g)
        wrk_int_state_out%dir_gradtxy_v1(:,jl,jb)= wrk_int_state_in%dir_gradtxy_v1(:,jl_g,jb_g)
        wrk_int_state_out%dir_gradtxy_v2(:,jl,jb)= wrk_int_state_in%dir_gradtxy_v2(:,jl_g,jb_g)
        wrk_int_state_out%dir_gradtyx_v1(:,jl,jb)= wrk_int_state_in%dir_gradtyx_v1(:,jl_g,jb_g)
        wrk_int_state_out%dir_gradtyx_v2(:,jl,jb)= wrk_int_state_in%dir_gradtyx_v2(:,jl_g,jb_g)
        wrk_int_state_out%shear_def_v1(:,jl,jb)= wrk_int_state_in%shear_def_v1(:,jl_g,jb_g)
        wrk_int_state_out%shear_def_v2(:,jl,jb)= wrk_int_state_in%shear_def_v2(:,jl_g,jb_g)

        IF(i_cori_method>=3)THEN
          wrk_int_state_out%quad_east   (:,jl,jb)= wrk_int_state_in%quad_east     (:,jl_g,jb_g)
          wrk_int_state_out%quad_north  (:,jl,jb)= wrk_int_state_in%quad_north    (:,jl_g,jb_g)
        ENDIF

        wrk_int_state_out%cno_en(jl,:,jb)= wrk_int_state_in%cno_en(jl_g,:,jb_g)
        wrk_int_state_out%cea_en(jl,:,jb)= wrk_int_state_in%cea_en(jl_g,:,jb_g)

        SELECT CASE (i_cori_method)
        CASE(1,3,4)
          nincr = 14
        CASE(2)
          nincr = 10
        END SELECT
        IF (i_cori_method < 3) THEN
          DO i=1,nincr
            CALL get_local_index(wrk_p_patch%edges%loc_index, &
              & wrk_int_state_in%heli_vn_idx(i,jl_g,jb_g), &
              & wrk_int_state_in%heli_vn_blk(i,jl_g,jb_g), &
              & wrk_int_state_out%heli_vn_idx(i,jl,jb), &
              & wrk_int_state_out%heli_vn_blk(i,jl,jb))
          ENDDO
        ENDIF

        DO i=1,6
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%dir_gradh_i1(i,jl_g,jb_g), &
            & wrk_int_state_in%dir_gradh_b1(i,jl_g,jb_g), &
            & wrk_int_state_out%dir_gradh_i1(i,jl,jb), &
            & wrk_int_state_out%dir_gradh_b1(i,jl,jb))
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%dir_gradh_i2(i,jl_g,jb_g), &
            & wrk_int_state_in%dir_gradh_b2(i,jl_g,jb_g), &
            & wrk_int_state_out%dir_gradh_i2(i,jl,jb), &
            & wrk_int_state_out%dir_gradh_b2(i,jl,jb))
        ENDDO

        DO i=1,9
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%dir_gradt_i1(i,jl_g,jb_g), &
            & wrk_int_state_in%dir_gradt_b1(i,jl_g,jb_g), &
            & wrk_int_state_out%dir_gradt_i1(i,jl,jb), &
            & wrk_int_state_out%dir_gradt_b1(i,jl,jb))
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%dir_gradt_i2(i,jl_g,jb_g), &
            & wrk_int_state_in%dir_gradt_b2(i,jl_g,jb_g), &
            & wrk_int_state_out%dir_gradt_i2(i,jl,jb), &
            & wrk_int_state_out%dir_gradt_b2(i,jl,jb))
        ENDDO

      ENDIF

      IF (i_cell_type == 3) THEN

        wrk_int_state_out%e_flx_avg(jl,:,jb) = wrk_int_state_in%e_flx_avg(jl_g,:,jb_g)

        DO i=1,rbf_vec_dim_e
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%rbf_vec_idx_e(i,jl_g,jb_g), &
            & wrk_int_state_in%rbf_vec_blk_e(i,jl_g,jb_g), &
            & wrk_int_state_out%rbf_vec_idx_e(i,jl,jb), &
            & wrk_int_state_out%rbf_vec_blk_e(i,jl,jb))
        ENDDO
        wrk_int_state_out%rbf_vec_stencil_e(jl, jb)  = &
          & wrk_int_state_in%rbf_vec_stencil_e(jl_g, jb_g)
        wrk_int_state_out%rbf_vec_coeff_e(:, jl, jb) = &
          & wrk_int_state_in%rbf_vec_coeff_e(:, jl_g, jb_g)
      ENDIF

      IF( i_cell_type == 3 .AND. (ltransport .OR. iequations == 3) ) THEN
        wrk_int_state_out%pos_on_tplane_e(jl, jb, :, :)     = &
          & wrk_int_state_in%pos_on_tplane_e(jl_g, jb_g, :, :)
        wrk_int_state_out%tplane_e_dotprod(jl, jb, :, :)    = &
          & wrk_int_state_in%tplane_e_dotprod(jl_g, jb_g, :, :)
      ENDIF

      IF (i_cell_type == 3 ) THEN
        wrk_int_state_out%geofac_qdiv(jl, :, jb)        = &
          & wrk_int_state_in%geofac_qdiv(jl_g, :, jb_g)
        wrk_int_state_out%nudgecoeff_e(jl, jb)        = &
          & wrk_int_state_in%nudgecoeff_e(jl_g, jb_g)
      ENDIF
      wrk_int_state_out%cart_edge_coord(jl, jb, :)    = &
        & wrk_int_state_in%cart_edge_coord(jl_g, jb_g, :)

    ENDDO


    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%verts%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%verts%glb_index(j)) ! Line  index in global patch

      IF (i_cell_type == 6) THEN
        wrk_int_state_out%tria_north(:,jl,jb) = wrk_int_state_in%tria_north(:,jl_g,jb_g)
        wrk_int_state_out%tria_east(:,jl,jb)  = wrk_int_state_in%tria_east(:,jl_g,jb_g)
        wrk_int_state_out%e_1o3_v(jl,:,jb)    = wrk_int_state_in%e_1o3_v(jl_g,:,jb_g)
        wrk_int_state_out%e_aw_v(jl,:,jb)     = wrk_int_state_in%e_aw_v(jl_g,:,jb_g)
        wrk_int_state_out%e_inn_v(jl,:,jb)    = wrk_int_state_in%e_inn_v(jl_g,:,jb_g)
      ENDIF

      wrk_int_state_out%cells_aw_verts(jl,:,jb) = wrk_int_state_in%cells_aw_verts(jl_g,:,jb_g)

      IF (i_cell_type == 3) THEN
        DO i=1,rbf_vec_dim_v
          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%rbf_vec_idx_v(i,jl_g,jb_g), &
            & wrk_int_state_in%rbf_vec_blk_v(i,jl_g,jb_g), &
            & wrk_int_state_out%rbf_vec_idx_v(i,jl,jb), &
            & wrk_int_state_out%rbf_vec_blk_v(i,jl,jb))
        ENDDO
        wrk_int_state_out%rbf_vec_stencil_v(jl, jb)      = &
          & wrk_int_state_in%rbf_vec_stencil_v(jl_g, jb_g)
        wrk_int_state_out%rbf_vec_coeff_v(:, :, jl, jb) = &
          & wrk_int_state_in%rbf_vec_coeff_v(:, :, jl_g, jb_g)
      ENDIF

      wrk_int_state_out%geofac_rot(jl, :, jb) = wrk_int_state_in%geofac_rot(jl_g, :, jb_g)

    ENDDO

  END SUBROUTINE divide_int_state

  !-------------------------------------------------------------------------
  !>
  !! Divides all variables in gridref_single_state.
  !!
  !!
  !! @par Revision History
  !! Initial version by Rainer Johanni, Nov 2009
  SUBROUTINE divide_grf_state()

    ! Local scalars:
    INTEGER :: jcd, j, js, je, jb, jl, jb_g, jl_g, i

    TYPE(t_gridref_single_state), POINTER :: wrk_int_state_out   => NULL()
    TYPE(t_gridref_single_state), POINTER :: wrk_int_state_in => NULL()

    INTEGER :: grf_vec_dim_1, grf_vec_dim_2

    !-----------------------------------------------------------------------
    grf_vec_dim_1 = 6
    grf_vec_dim_2 = 5

    ! This routine must not be called in a single CPU run
    IF(p_nprocs == 1 .or. p_pe == p_test_pe) &
      & CALL finish('divide_grf_state','must not be called in a single CPU run')

    ! Loop over all child domains of the present domain
    ! If there are no childs (wrk_p_patch%n_childdom==0) nothing needs to be done

    DO jcd = 1, wrk_p_patch%n_childdom

      wrk_int_state_out   => wrk_gridref_state_out%p_dom(jcd)
      wrk_int_state_in => wrk_gridref_state_in%p_dom(jcd)

      ! Copy area of feedback domain
      wrk_gridref_state_out%fbk_dom_area(jcd) = wrk_gridref_state_in%fbk_dom_area(jcd)

      ! Set cell related entries

      js = idx_1d(wrk_p_patch%cells%start_idx(grf_bdyintp_start_c,jcd), &
        &         wrk_p_patch%cells%start_blk(grf_bdyintp_start_c,jcd))

      je = idx_1d(wrk_p_patch%cells%end_idx(min_rlcell_int,jcd), &
        &         wrk_p_patch%cells%end_blk(min_rlcell_int,jcd))

      DO j=js,je

        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index

        jb_g = blk_no(wrk_p_patch%cells%glb_index(j)) ! Block index in global patch
        jl_g = idx_no(wrk_p_patch%cells%glb_index(j)) ! Line  index in global patch

        wrk_int_state_out%grf_dist_pc2cc (jl,:,:,jb) = &
          & wrk_int_state_in%grf_dist_pc2cc (jl_g,:,:,jb_g)

      ENDDO


      ! Set edge related entries

      js = idx_1d(wrk_p_patch%edges%start_idx(grf_bdyintp_start_e,jcd), &
        &         wrk_p_patch%edges%start_blk(grf_bdyintp_start_e,jcd))

      je = idx_1d(wrk_p_patch%edges%end_idx(min_rledge_int,jcd), &
        &         wrk_p_patch%edges%end_blk(min_rledge_int,jcd))

      DO j=js,je

        jb = blk_no(j) ! Block index
        jl = idx_no(j) ! Line  index

        jb_g = blk_no(wrk_p_patch%edges%glb_index(j)) ! Block index in global patch
        jl_g = idx_no(wrk_p_patch%edges%glb_index(j)) ! Line  index in global patch

        wrk_int_state_out%grf_dist_pe2ce (jl,:,jb) = &
          & wrk_int_state_in%grf_dist_pe2ce (jl_g,:,jb_g)

        DO i = 1, grf_vec_dim_1

          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%grf_vec_ind_1a(jl_g,i,jb_g), &
            & wrk_int_state_in%grf_vec_blk_1a(jl_g,i,jb_g), &
            & wrk_int_state_out%grf_vec_ind_1a(jl,i,jb), &
            & wrk_int_state_out%grf_vec_blk_1a(jl,i,jb))

          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%grf_vec_ind_1b(jl_g,i,jb_g), &
            & wrk_int_state_in%grf_vec_blk_1b(jl_g,i,jb_g), &
            & wrk_int_state_out%grf_vec_ind_1b(jl,i,jb), &
            & wrk_int_state_out%grf_vec_blk_1b(jl,i,jb))
        ENDDO

        DO i = 1, grf_vec_dim_2

          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%grf_vec_ind_2a(jl_g,i,jb_g), &
            & wrk_int_state_in%grf_vec_blk_2a(jl_g,i,jb_g), &
            & wrk_int_state_out%grf_vec_ind_2a(jl,i,jb), &
            & wrk_int_state_out%grf_vec_blk_2a(jl,i,jb))

          CALL get_local_index(wrk_p_patch%edges%loc_index, &
            & wrk_int_state_in%grf_vec_ind_2b(jl_g,i,jb_g), &
            & wrk_int_state_in%grf_vec_blk_2b(jl_g,i,jb_g), &
            & wrk_int_state_out%grf_vec_ind_2b(jl,i,jb), &
            & wrk_int_state_out%grf_vec_blk_2b(jl,i,jb))
        ENDDO


        wrk_int_state_out%grf_vec_stencil_1a(jl,jb) = &
          & wrk_int_state_in%grf_vec_stencil_1a(jl_g,jb_g)
        wrk_int_state_out%grf_vec_stencil_1b(jl,jb) = &
          & wrk_int_state_in%grf_vec_stencil_1b(jl_g,jb_g)
        wrk_int_state_out%grf_vec_stencil_2a(jl,jb) = &
          & wrk_int_state_in%grf_vec_stencil_2a(jl_g,jb_g)
        wrk_int_state_out%grf_vec_stencil_2b(jl,jb) = &
          & wrk_int_state_in%grf_vec_stencil_2b(jl_g,jb_g)

        wrk_int_state_out%grf_vec_coeff_1a(:,jl,jb) = &
          & wrk_int_state_in%grf_vec_coeff_1a(:,jl_g,jb_g)
        wrk_int_state_out%grf_vec_coeff_1b(:,jl,jb) = &
          & wrk_int_state_in%grf_vec_coeff_1b(:,jl_g,jb_g)
        wrk_int_state_out%grf_vec_coeff_2a(:,jl,jb) = &
          & wrk_int_state_in%grf_vec_coeff_2a(:,jl_g,jb_g)
        wrk_int_state_out%grf_vec_coeff_2b(:,jl,jb) = &
          & wrk_int_state_in%grf_vec_coeff_2b(:,jl_g,jb_g)

      ENDDO

    ENDDO

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%cells%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%cells%glb_index(j)) ! Line  index in global patch

      wrk_gridref_state_out%fbk_wgt_c(jl,jb,:) = wrk_gridref_state_in%fbk_wgt_c(jl_g,jb_g,:)
      wrk_gridref_state_out%fbk_wgt_ct(jl,jb,:) = wrk_gridref_state_in%fbk_wgt_ct(jl_g,jb_g,:)
      wrk_gridref_state_out%pc_idx_c(jl,jb) = wrk_gridref_state_in%pc_idx_c(jl_g,jb_g)
    ENDDO

    DO j = 1, wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jb_g = blk_no(wrk_p_patch%edges%glb_index(j)) ! Block index in global patch
      jl_g = idx_no(wrk_p_patch%edges%glb_index(j)) ! Line  index in global patch

      wrk_gridref_state_out%fbk_wgt_e(jl,jb,:) = wrk_gridref_state_in%fbk_wgt_e(jl_g,jb_g,:)
      wrk_gridref_state_out%pc_idx_e(jl,jb) = wrk_gridref_state_in%pc_idx_e(jl_g,jb_g)
    ENDDO

    ! fbk_dom_volume is not yet set, we copy it nevertheless
    wrk_gridref_state_out%fbk_dom_area(:) = wrk_gridref_state_in%fbk_dom_area(:)
    wrk_gridref_state_out%fbk_dom_volume(:,:) = wrk_gridref_state_in%fbk_dom_volume(:,:)

    ! tracer_bdyflx is not yet set and might go away eventually - no need to care

  END SUBROUTINE divide_grf_state

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
    ! mode = 0 : Mapping doesn't matter, but it must point to a valid local index

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
        ! It doesn't matter which one we choose, but it must be valid
        j_l = MAX(ABS(j_l)-1,1)
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
  SUBROUTINE divide_subset_geometric(subset_flag, n_proc, owner)

    INTEGER, INTENT(in)    :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    INTEGER, INTENT(out)   :: owner(:) ! receives the owner PE for every cell
    ! (-1 for cells not in subset)

    INTEGER :: i, j, jl, jb, jl_v, jb_v, nc
    REAL(wp), ALLOCATABLE :: cell_desc(:,:)

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

        ! Patch should be divided for radiation calculations,
        ! i.e. every patch should contain 2 distinct areas from
        ! opposite regions of the earth.
        ! We do that by mapping cells with negative lattitudes
        ! to the opposed position on earth.
        ! Since this is for radiation calculations, "toothed" edges
        ! don't matter and so we use the cell centers for subdivision.

        IF(wrk_divide_patch%cells%center(jl,jb)%lat>=0._wp) THEN
          cell_desc(1,nc) = wrk_divide_patch%cells%center(jl,jb)%lat
          cell_desc(2,nc) = wrk_divide_patch%cells%center(jl,jb)%lon
        ELSE
          cell_desc(1,nc) = -wrk_divide_patch%cells%center(jl,jb)%lat
          cell_desc(2,nc) =  wrk_divide_patch%cells%center(jl,jb)%lon + pi
          IF(cell_desc(2,nc)>pi) cell_desc(2,nc) = cell_desc(2,nc) - 2._wp*pi
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
    ALLOCATE(metis_adjncy(wrk_divide_patch%n_patch_cells*i_cell_type))
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

