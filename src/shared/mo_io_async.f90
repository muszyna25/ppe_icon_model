!>
!! Contains routines for asynchronous I/O
!!
!!
!! @par Revision History
!! Initial implementation by Rainer Johanni (2010-11-11)
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
!!

! Define USE_CRAY_POINTER for platforms having problems with ISO_C_BINDING
! BUT understand CRAY pointers

!!! #define USE_CRAY_POINTER


MODULE mo_io_async

#ifndef NOMPI
  ! Please note: Since we use MPI_Alloc_mem with the wrong signature below (TYPE(c_ptr)
  ! where INTEGER(KIND=MPI_ADDRESS_KIND) is required), we are getting here only the
  ! required constants from module mpi. Otherwise we get problems on the NEC
  USE mpi, ONLY: MPI_ADDRESS_KIND, MPI_ROOT, MPI_PROC_NULL, MPI_ANY_SOURCE, MPI_INFO_NULL, &
   &             MPI_LOCK_SHARED, MPI_LOCK_EXCLUSIVE, MPI_MODE_NOCHECK
#endif

#ifndef USE_CRAY_POINTER
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_intptr_t, c_f_pointer
#endif
  USE mo_mpi,  ONLY: p_comm_work, p_comm_work_io, p_comm_work_2_io

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: max_dom, ihs_atm_temp, ihs_atm_theta, &
                                    inh_atmosphere, ishallow_water, inwp
  USE mo_datetime,            ONLY: t_datetime
  USE mo_mpi,                 ONLY: p_pe, p_bcast, p_barrier, p_stop, p_real_dp, p_send, &
    & p_recv, my_process_is_mpi_test
  USE mo_parallel_configuration,  ONLY: p_pe_work, p_work_pe0, p_io_pe0,     &
   &                                num_work_procs, num_io_procs, pio_type
  USE mo_global_variables,    ONLY: setup_physics
  USE mo_nonhydrostatic_nml,  ONLY: ivctype, nonhydrostatic_nml_setup
  USE mo_dynamics_nml,        ONLY: dynamics_nml_setup
  USE mo_diffusion_nml,       ONLY: diffusion_nml_setup
!  USE mo_io_nml,              ONLY: io_nml_setup
  USE mo_io_config            
  USE mo_dynamics_config,     ONLY: iequations 
  USE mo_run_config,          ONLY: ldump_states, ltransport, lforcing, num_lev, iforcing, nlev
 ! USE mo_atm_phy_nwp_nml,     ONLY: setup_nwp_phy, inwp_surface
  USE mo_atm_phy_nwp_config, ONLY: configure_atm_phy_nwp
  USE mo_io_units,            ONLY: filename_max
  USE mo_communication,       ONLY: idx_no, blk_no
  USE mo_io_vlist,            ONLY: GATHER_C, GATHER_E, GATHER_V,                                &
   &                                setup_vlist, destruct_vlist,                                 &
   &                                open_output_vlist, close_output_vlist,                       &
   &                                vlist_set_date_time, vlist_start_step, vlist_write_var,      &
   &                                num_output_vars, outvar_desc,                                &
   &                                get_outvar_ptr_ha, get_outvar_ptr_nh
  USE mo_grid_configuration,  ONLY: n_dom, parent_id
  USE mo_vertical_coord_table,ONLY: init_vertical_coord_table
  USE mo_vertical_grid,       ONLY: init_hybrid_coord, init_sleve_coord
  USE mo_advection_nml,       ONLY: transport_nml_setup
  USE mo_namelist,            ONLY: close_nml

  !-------------------------------------------------------------------------------------------------
  ! Needed only for compute PEs, patches are NOT set on I/O PEs

  USE mo_atmo_control,        ONLY: p_patch

  ! End of needed only for compute PEs
  !-------------------------------------------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'cdi.inc' ! for cdiEncodeDate/cdiEncodeTime

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_io_async'

  ! Tags for communication between compte PEs and I/O PEs

  INTEGER, PARAMETER :: msg_io_start    = 123456
  INTEGER, PARAMETER :: msg_io_done     = 654321
  INTEGER, PARAMETER :: msg_io_shutdown = 999999

  ! Variables used for MPI communication between Compute and I/O PEs

  INTEGER :: mpi_win ! Memory window

  REAL(wp), POINTER :: mem_ptr(:) ! Pointer to memory window

  ! Info about own cells/edges/verts

  TYPE t_owner_info

    INTEGER :: n_own_cells
    INTEGER :: n_own_edges
    INTEGER :: n_own_verts

    INTEGER, ALLOCATABLE :: own_cell_idx(:), own_cell_blk(:)
    INTEGER, ALLOCATABLE :: own_edge_idx(:), own_edge_blk(:)
    INTEGER, ALLOCATABLE :: own_vert_idx(:), own_vert_blk(:)
  END TYPE t_owner_info

  TYPE(t_owner_info), ALLOCATABLE, TARGET :: patch_owner_info(:)

  ! Output file names

  CHARACTER(LEN=filename_max) :: output_file_name(max_dom)
  LOGICAL :: new_output_files = .FALSE.

  !-------------------------------------------------------------------------------------------------
  ! Work distribution

  INTEGER :: io_task_no(max_dom) ! I/O task number caring about specific patch
  INTEGER :: num_io_tasks
  INTEGER :: my_io_task_no

  INTEGER, ALLOCATABLE :: var_pe_no(:,:) ! Absoulte PE number for every variable in output

  !-------------------------------------------------------------------------------------------------
  ! PIO

  LOGICAL :: use_pio = .FALSE. ! Flag if to use parallel I/O

  INTEGER n_pes_io_task    ! Number of PEs in my I/O task
  INTEGER my_pe_io_task    ! My rank within my I/O task
  INTEGER pio_comm         ! Communicator within I/O task

  INTEGER, ALLOCATABLE :: io_task_pe0(:) ! PE 0 of every I/O task (stored only on I/O PE 0)
  
  !-------------------------------------------------------------------------------------------------
  ! The following definitions are only used on I/O PEs

  TYPE t_patch_desc ! Contains only data necessary on I/O PEs

    ! total number of (global) cells, edges and vertices

    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g

    ! Global arrays of owners

    INTEGER, ALLOCATABLE :: cell_owner_g(:)
    INTEGER, ALLOCATABLE :: edge_owner_g(:)
    INTEGER, ALLOCATABLE :: vert_owner_g(:)

    ! Number of cells/edges/verts owned by every worker PE

    INTEGER, ALLOCATABLE :: n_own_cells(:)
    INTEGER, ALLOCATABLE :: n_own_edges(:)
    INTEGER, ALLOCATABLE :: n_own_verts(:)

    ! Sort index for getting the received arrays into native order

    INTEGER, ALLOCATABLE :: sort_index_cells(:)
    INTEGER, ALLOCATABLE :: sort_index_edges(:)
    INTEGER, ALLOCATABLE :: sort_index_verts(:)

  END TYPE t_patch_desc

  TYPE(t_patch_desc), ALLOCATABLE, TARGET :: p_patch_desc(:)

  !-------------------------------------------------------------------------------------------------

  ! Public routines:

  PUBLIC :: io_main_proc, setup_io_procs, shutdown_io_procs, set_output_file, output_async

CONTAINS

#ifdef NOMPI
  ! The whole purpose of this module is to run under MPI, so just define
  ! entry points for an error-free link here

  SUBROUTINE io_main_proc
  END SUBROUTINE io_main_proc

  SUBROUTINE setup_io_procs()
!     CHARACTER(LEN=*), INTENT(INOUT) :: gridfile(:)
  END SUBROUTINE setup_io_procs

  SUBROUTINE shutdown_io_procs
  END SUBROUTINE shutdown_io_procs

  SUBROUTINE set_output_file(file, jg)
    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: jg
  END SUBROUTINE set_output_file

  SUBROUTINE output_async(datetime, z_sim_time)
    TYPE(t_datetime), INTENT(in) :: datetime
    REAL(wp), OPTIONAL, INTENT(in) :: z_sim_time(n_dom)
  END SUBROUTINE output_async

#else
  !-------------------------------------------------------------------------------------------------
  !>
  !! Returns the root for the sender side of intercommunicator broadcasts (always sent from PE 0).
  !! Please note the special root setting for intercommunicators:
  !! The PE really sending must use MPI_ROOT, the others MPI_PROC_NULL

  INTEGER FUNCTION bcast_root()

    IF(p_pe_work == 0) THEN
      bcast_root = MPI_ROOT
    ELSE
      bcast_root = MPI_PROC_NULL
    ENDIF

  END FUNCTION bcast_root

  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  ! The following routines are called only from I/O PEs
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------
  !>
  !! Main routine for I/O PEs.
  !! Please note that this routine never returns.

  SUBROUTINE io_main_proc

    INTEGER jg, idate, itime, iostep
    LOGICAL done

    print '(a,i0)','============================================ Hello from I/O PE ',p_pe

    ! If ldump_states is set, the compute PEs will exit after dumping,
    ! there is nothing to do at all for I/O PEs

    IF(ldump_states) THEN
      CALL p_stop
      STOP
    ENDIF

    ! ----------------------------------------------------------------------------------------------
    ! This routine is called from control_model immediatly before reading patches.
    ! Unfortunately the setup (reading the namelists) is not completed at that point
    ! since it is interspersed with calculations (like setting up the interploation).
    ! Therefore we have to complete the setup here, at least for all relevant data
    ! which is needed in this module and for setting up the CDI vlists.
    !
    ! Please note that the setup sequence of control_model is duplicated here.
    ! !!! This implies that changes in control_model have to be repeated here !!!

    CALL diffusion_nml_setup(n_dom,parent_id,nlev)
    CALL dynamics_nml_setup(n_dom)  

    IF (ltransport) CALL transport_nml_setup

    SELECT CASE (iequations)
      CASE (ishallow_water)
        CALL init_vertical_coord_table(iequations, num_lev(1))
        !
      CASE (ihs_atm_temp, ihs_atm_theta)
        CALL init_vertical_coord_table(iequations, num_lev(1))
        !
      CASE (inh_atmosphere)
        CALL nonhydrostatic_nml_setup
        IF (ivctype == 1) THEN
          CALL init_hybrid_coord (iequations, num_lev(1))
        ELSE IF (ivctype == 2) THEN
          CALL init_sleve_coord (num_lev(1))
        ENDIF
        !
    END SELECT

    IF ( lforcing ) CALL setup_physics

    IF ( iforcing == inwp) THEN 
!      CALL setup_nwp_phy
      CALL configure_atm_phy_nwp
!KF temp
!      IF (inwp_surface > 0) CALL setup_nwp_lnd
    ENDIF

!    CALL io_nml_setup

    CALL close_nml

    ! ----------------------------------------------------------------------------------------------
    ! Done with setup, receive relevant parts of patch configuration from compute PEs

    CALL receive_patch_configuration

    ! Distribute work among the (possibly) several I/O PEs

    CALL distribute_work

    ! Open first set of output files

    DO jg = 1, n_dom
      IF(io_task_no(jg) == my_io_task_no) THEN
        CALL open_output_vlist(output_file_name(jg), jg)
      ENDIF
    ENDDO
    new_output_files = .FALSE.

    ! Tell the compute PEs that we are ready to work

    CALL io_send_ready_message

    ! Enter I/O loop

    iostep = 0

    DO

      ! Wait for a message from the compute PEs to start
      CALL io_wait_for_start_message(done, idate, itime)

      IF(done) EXIT ! leave loop, we are done

      IF(new_output_files) THEN
        ! Open new set of output files
        DO jg = 1, n_dom
          IF(io_task_no(jg) == my_io_task_no) THEN
            CALL close_output_vlist(jg)
            CALL open_output_vlist(output_file_name(jg), jg)
          ENDIF
        ENDDO
        new_output_files = .FALSE.
        iostep = 0
      ENDIF

      ! Initialize step
      DO jg = 1, n_dom
        CALL vlist_set_date_time(jg, idate, itime)
        IF(io_task_no(jg) == my_io_task_no) THEN
          CALL vlist_start_step(jg, iostep)
        ENDIF
      ENDDO

      ! perform I/O
      CALL do_io

      iostep = iostep + 1

      ! Inform compute PEs that we are done
      CALL io_send_ready_message

    ENDDO

    ! Finalization sequence:

    print '(a,i0,a)','============================================ I/O PE ',p_pe,' shutting down'
    DO jg = 1, n_dom
      IF(io_task_no(jg) == my_io_task_no) THEN
        CALL close_output_vlist(jg)
      ENDIF
      CALL destruct_vlist(jg)
    ENDDO

    IF(use_pio) CALL pioFinalize

    ! Shut down MPI
    !
    CALL p_stop

    STOP

  END SUBROUTINE io_main_proc

  !-------------------------------------------------------------------------------------------------
  !>
  !! Distributes work among the I/O PEs

  SUBROUTINE distribute_work

    INTEGER jg, res, i, n, nv1, nv2, mpierr
    CHARACTER(LEN=256) text

    ALLOCATE(var_pe_no(MAXVAL(num_output_vars(1:n_dom)), n_dom))
    var_pe_no(:,:) = -1

    IF(num_io_procs <= n_dom) THEN

      ! Every I/O PE gets at least 1 domain, no need to do parallel I/O
      ! We map the domains just cyclically to the I/O procs, maybe there exist better ways ...

      use_pio = .FALSE.

      my_io_task_no = p_pe_work
      num_io_tasks = num_io_procs

      DO jg = 1, n_dom
        io_task_no(jg) = MOD(jg-1,num_io_tasks)
        var_pe_no(:,jg) = io_task_no(jg)
      ENDDO

    ELSE

      ! Use parallel I/O

      use_pio = .TRUE.

      res = pioInit ( pio_type, p_comm_work, my_io_task_no, num_io_tasks, pio_comm )

      IF(res==0) THEN
        ! This is the return of the writer PEs at the very end
        CALL pioFinalize
        CALL p_stop
        STOP
      ENDIF

      IF(num_io_tasks > n_dom) THEN
        WRITE(text,'(2(a,i0))') 'ERROR: Number of PIO tasks = ',num_io_tasks, &
         &                      ' is greater than # of domains = ', n_dom
        CALL finish(modname, text)
      ENDIF

      ! my_io_task_no seems to be 1 based - make it 0 based

      my_io_task_no = my_io_task_no - 1

      CALL MPI_Comm_size(pio_comm, n_pes_io_task, mpierr)
      CALL MPI_Comm_rank(pio_comm, my_pe_io_task, mpierr)

      ! Map the domains cyclically to the I/O tasks
      ! As above: Maybe there exist better ways ...

      DO jg = 1, n_dom
        io_task_no(jg) = MOD(jg-1,num_io_tasks)
      ENDDO

      ! Map variables to PEs
      ! Note that var_pe_no is set a bit different than above:
      ! It contains either -1 or p_pe_work for current PE.
      ! This shouldn't matter during output.

      nv1 = 0
      nv2 = 0

      DO jg = 1, n_dom
        IF(io_task_no(jg) /= my_io_task_no) CYCLE
        DO n = 1, num_output_vars(jg)
          IF(outvar_desc(n,jg)%nlev == 1) THEN
            ! Map 1D array
            IF(MOD(nv1,n_pes_io_task) == my_pe_io_task) var_pe_no(n,jg) = p_pe_work
            nv1 = nv1 + 1
          ELSE
            ! Map 2D array
            IF(MOD(nv2,n_pes_io_task) == my_pe_io_task) var_pe_no(n,jg) = p_pe_work
            nv2 = nv2 + 1
          ENDIF
        ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      ! !!! ATTENTION !!!
      ! The following code assumes that I/O PE 0 will get PE 0 of io_task 0
      ! Otherwise the code will go horribly wrong !!!

      IF(my_pe_io_task == 0) THEN
        IF(my_io_task_no > 0) THEN
          ! The group leaders of every I/O task send their PE numbers to PE 0
          n = p_pe_work
          CALL p_send(n,0,100+my_io_task_no,comm=p_comm_work)
        ELSE
          ! PE 0 of io_task 0 receives PE numbers
          ALLOCATE(io_task_pe0(num_io_tasks-1))
          DO i=1,num_io_tasks-1
            CALL p_recv(n,MPI_ANY_SOURCE,100+i,comm=p_comm_work)
            io_task_pe0(i) = n
          ENDDO
        ENDIF
      ENDIF

    ENDIF

  END SUBROUTINE distribute_work

  !-------------------------------------------------------------------------------------------------
  !>
  !! Receive patch configuration from compute PEs, participate in memory window creation,
  !! set up and check CDI vlists

  SUBROUTINE receive_patch_configuration

    INTEGER :: jg, np, i, n, n_dims(3), numv, nbytes_real, mpierr
    INTEGER :: idx(0:num_work_procs-1)
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size
    INTEGER, ALLOCATABLE :: ncheck(:,:)
    REAL(wp), TARGET :: dummy(1)
    CHARACTER(LEN=filename_max) :: gridfile(max_dom)
    CHARACTER(LEN=256) text


    ! Please note: The broacasts below use an intercommunicator, i.e. we are getting all
    ! data from compute PE 0

    ! Allocate and fill patch descriptions

    ALLOCATE(p_patch_desc(n_dom))

    DO jg = 1, n_dom

      CALL p_bcast(n_dims, 0, p_comm_work_2_io)

      p_patch_desc(jg)%n_patch_cells_g = n_dims(1)
      p_patch_desc(jg)%n_patch_edges_g = n_dims(2)
      p_patch_desc(jg)%n_patch_verts_g = n_dims(3)

      ALLOCATE(p_patch_desc(jg)%cell_owner_g(n_dims(1)))
      ALLOCATE(p_patch_desc(jg)%edge_owner_g(n_dims(2)))
      ALLOCATE(p_patch_desc(jg)%vert_owner_g(n_dims(3)))

      CALL p_bcast(p_patch_desc(jg)%cell_owner_g, 0, p_comm_work_2_io)
      CALL p_bcast(p_patch_desc(jg)%edge_owner_g, 0, p_comm_work_2_io)
      CALL p_bcast(p_patch_desc(jg)%vert_owner_g, 0, p_comm_work_2_io)

      ! Get grid file

      CALL p_bcast(gridfile(jg), 0, p_comm_work_2_io)

      ! Get number of data points on every worker PE

      ALLOCATE(p_patch_desc(jg)%n_own_cells(0:num_work_procs-1))
      ALLOCATE(p_patch_desc(jg)%n_own_edges(0:num_work_procs-1))
      ALLOCATE(p_patch_desc(jg)%n_own_verts(0:num_work_procs-1))

      DO np = 0, num_work_procs-1
        p_patch_desc(jg)%n_own_cells(np) = COUNT(p_patch_desc(jg)%cell_owner_g(:) == np)
        p_patch_desc(jg)%n_own_edges(np) = COUNT(p_patch_desc(jg)%edge_owner_g(:) == np)
        p_patch_desc(jg)%n_own_verts(np) = COUNT(p_patch_desc(jg)%vert_owner_g(:) == np)
      ENDDO

      ! Set the sort_index arrays

      ALLOCATE(p_patch_desc(jg)%sort_index_cells(n_dims(1)))
      ALLOCATE(p_patch_desc(jg)%sort_index_edges(n_dims(2)))
      ALLOCATE(p_patch_desc(jg)%sort_index_verts(n_dims(3)))

      idx(0) = 0
      DO i=1, num_work_procs-1
        idx(i) = idx(i-1) + p_patch_desc(jg)%n_own_cells(i-1)
      ENDDO
      DO i = 1, p_patch_desc(jg)%n_patch_cells_g
        np = p_patch_desc(jg)%cell_owner_g(i)
        idx(np) = idx(np)+1
        p_patch_desc(jg)%sort_index_cells(i) = idx(np)
      ENDDO

      idx(0) = 0
      DO i=1, num_work_procs-1
        idx(i) = idx(i-1) + p_patch_desc(jg)%n_own_edges(i-1)
      ENDDO
      DO i = 1, p_patch_desc(jg)%n_patch_edges_g
        np = p_patch_desc(jg)%edge_owner_g(i)
        idx(np) = idx(np)+1
        p_patch_desc(jg)%sort_index_edges(i) = idx(np)
      ENDDO

      idx(0) = 0
      DO i=1, num_work_procs-1
        idx(i) = idx(i-1) + p_patch_desc(jg)%n_own_verts(i-1)
      ENDDO
      DO i = 1, p_patch_desc(jg)%n_patch_verts_g
        np = p_patch_desc(jg)%vert_owner_g(i)
        idx(np) = idx(np)+1
        p_patch_desc(jg)%sort_index_verts(i) = idx(np)
      ENDDO

    ENDDO

    ! ----------------------------------------------------------------------------------------------
    ! Memory window creation

    ! Get the amount of bytes per default REAL variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpierr)

    ! We use a NULL window for the I/O PEs
    mem_ptr => dummy
    mem_size = 0_MPI_ADDRESS_KIND

    CALL MPI_Win_create(mem_ptr,mem_size,nbytes_real,MPI_INFO_NULL,p_comm_work_io,mpi_win,mpierr)

    ! ----------------------------------------------------------------------------------------------
    ! All variables needed for output are initialized now, we can set up the vlists

    DO jg = 1, n_dom
      CALL setup_vlist( TRIM(gridfile(jg)), jg )
    ENDDO

    ! ----------------------------------------------------------------------------------------------
    ! Just for safety since it might be possible that the I/O PEs missed some initialization:
    ! Check if the number and types of output variables are the same on I/O and compute PEs:

    DO jg = 1, n_dom

      CALL p_bcast(numv, 0, p_comm_work_2_io)
      IF(numv /= num_output_vars(jg)) THEN
        WRITE(text,'(a,i0,a,i0)') 'Mismatch of number of output variables, I/O: ', &
         &                        num_output_vars(jg),' Compute: ',numv
        CALL finish(modname, TRIM(text))
      ENDIF

      ALLOCATE(ncheck(2,numv))
      CALL p_bcast(ncheck, 0, p_comm_work_2_io)
      DO n=1, numv
        IF(ncheck(1,n) /= outvar_desc(n,jg)%type .OR. ncheck(2,n) /= outvar_desc(n,jg)%nlev) THEN
          WRITE(text,'(3(a,2i5))') 'Mismatch of type/levels output variable ',n, jg,        &
           &                       ' I/O: ',outvar_desc(n,jg)%type, outvar_desc(n,jg)%nlev, &
           &                       ' Compute: ',ncheck(1,n),ncheck(2,n)
          CALL finish(modname, TRIM(text))
        ENDIF
      ENDDO
      DEALLOCATE(ncheck)

    ENDDO

    ! ----------------------------------------------------------------------------------------------
    ! Receive first set of files to be opened.

    CALL p_bcast(output_file_name(1:n_dom), 0, p_comm_work_2_io)

  END SUBROUTINE receive_patch_configuration

  !-------------------------------------------------------------------------------------------------
  !>
  !! Actually does the I/O on the I/O PE side

  SUBROUTINE do_io

    INTEGER jg, n, jk, np, i, nv, nv_off, nlev_max, n_glb, mpierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: ioff(0:num_work_procs-1)
    INTEGER, POINTER :: n_own(:), iidx(:)
    REAL(wp), ALLOCATABLE :: var1(:), var2(:)
    CHARACTER*10 ctime
#if defined (__SX__) && !defined (NOMPI)
! It may be necessary that var1 is in global memory on NEC
! (Note: this is only allowed when we compile with MPI.)
!CDIR GM_ARRAY(var1)
#endif

    CALL date_and_time(TIME=ctime)
    print '(a,i0,a)','#################### I/O PE ',p_pe,' starting I/O at '//ctime

    ! Get maximum number of data points in a slice and allocate tmp variables

    nv = 0
    nlev_max = 0
    DO jg = 1, n_dom
      nv = MAX(nv, p_patch_desc(jg)%n_patch_cells_g)
      nv = MAX(nv, p_patch_desc(jg)%n_patch_edges_g)
      nv = MAX(nv, p_patch_desc(jg)%n_patch_verts_g)
      DO n = 1, num_output_vars(jg)
        nlev_max = MAX(nlev_max, outvar_desc(n,jg)%nlev)
      ENDDO
    ENDDO

    ALLOCATE(var1(nv), var2(nv*nlev_max))

    ioff(:) = 0_MPI_ADDRESS_KIND

    DO jg = 1, n_dom

      DO n = 1, num_output_vars(jg)
!Re-enable to check performance problems!!!
!CALL date_and_time(TIME=ctime)
!print '(a,i0,a,2i5,a)','#################### I/O PE ',p_pe,' jg,n= ',jg,n,' at '//ctime

        ! Set n_own, n_glb, iidx according to the type of the variable

        SELECT CASE(outvar_desc(n,jg)%type)
          CASE (GATHER_C)
            n_own => p_patch_desc(jg)%n_own_cells
            n_glb =  p_patch_desc(jg)%n_patch_cells_g
            iidx  => p_patch_desc(jg)%sort_index_cells
          CASE (GATHER_E)
            n_own => p_patch_desc(jg)%n_own_edges
            n_glb =  p_patch_desc(jg)%n_patch_edges_g
            iidx  => p_patch_desc(jg)%sort_index_edges
          CASE (GATHER_V)
            n_own => p_patch_desc(jg)%n_own_verts
            n_glb =  p_patch_desc(jg)%n_patch_verts_g
            iidx  => p_patch_desc(jg)%sort_index_verts
          CASE DEFAULT
            CALL finish(modname, 'Illegal type in outvar_desc')
        END SELECT

        IF(var_pe_no(n, jg) /= p_pe_work) THEN
          ! Update the offset in the memory window and cycle
          ioff(:) = ioff(:) + n_own(:)*outvar_desc(n,jg)%nlev
          CYCLE
        ENDIF

        ! Loop over all levels in the output variable

        nv_off = 0

        DO jk = 1, outvar_desc(n,jg)%nlev

          ! Retrieve part of level from every worker PE using MPI_Get

          nv = 0

          DO np = 0, num_work_procs-1

            IF(n_own(np) == 0) CYCLE

            CALL MPI_Win_lock(MPI_LOCK_SHARED, np, MPI_MODE_NOCHECK, mpi_win, mpierr)

            CALL MPI_Get(var1(nv+1), n_own(np), p_real_dp, np, ioff(np), &
             &           n_own(np), p_real_dp, mpi_win, mpierr)

            CALL MPI_Win_unlock(np, mpi_win, mpierr)

            nv = nv + n_own(np)

          ENDDO

          ! Update the offset in the memory window on compute PEs

          ioff(:) = ioff(:) + n_own(:)

          ! var1 is stored in the order in which the variable was stored on compute PEs,
          ! get it back into the global storage order

          DO i = 1, n_glb
            var2(nv_off + i) = var1(iidx(i))
          ENDDO

          nv_off = nv_off + n_glb

        ENDDO ! Loop over levels

        CALL vlist_write_var(n, jg, var2)

      ENDDO ! Loop over output variables

    ENDDO ! Loop over patches

    DEALLOCATE(var1, var2)

    CALL date_and_time(TIME=ctime)
    print '(a,i0,a)','#################### I/O PE ',p_pe,' done at '//ctime

  END SUBROUTINE do_io

  !-------------------------------------------------------------------------------------------------
  !>
  !! io_send_ready_message: Send a message to the compute PEs that the I/O is ready
  !! The counterpart on the compute side is compute_wait_for_io_ready

  SUBROUTINE io_send_ready_message

    INTEGER msg

    ! make sure all are done
    IF(use_pio) THEN
      CALL pio_barrier
    ELSE
      CALL p_barrier(comm=p_comm_work)
    ENDIF

    ! Simply send a message from I/O PE 0 to compute PE 0
    IF(p_pe_work == 0) THEN
      msg = msg_io_done
      CALL p_send(msg, p_work_pe0, 0)
    ENDIF

  END SUBROUTINE io_send_ready_message

  !-------------------------------------------------------------------------------------------------
  !>
  !! io_wait_for_start_message: Wait for a message from I/O PEs that we should start I/O or finish
  !! The counterpart on the compute side is compute_start_io/compute_shutdown_io

  SUBROUTINE io_wait_for_start_message(done, idate, itime)

    LOGICAL, INTENT(OUT) :: done ! flag if we should shut down
    INTEGER, INTENT(OUT) :: idate, itime ! Date and time for I/O

    INTEGER msg(4), jg

    ! Receive message that we may start I/O (or should finish)
    ! If necessary, there may be more information transferred with this message
    ! in the future (e.g. what should be output in this turn)
    ! Currently, msg contains the following 4 values
    ! 1: Tag for I/O or shutdown
    ! 2: Date
    ! 3: Time
    ! 4: Flag if new output files follow

    IF(p_pe_work == 0) CALL p_recv(msg, p_work_pe0, 0)
    IF(use_pio) THEN
      CALL pio_bcast(msg)
    ELSE
      CALL p_bcast(msg, 0, comm=p_comm_work)
    ENDIF

    done  = .FALSE.
    idate = 0
    itime = 0

    SELECT CASE(msg(1))

    CASE(msg_io_start)

      idate = msg(2)
      itime = msg(3)

      IF(msg(4) /= 0) THEN
        ! The names of the new output files follow
        ! For PIO, it is currently not possible to switch output files
        ! so there is no need to bother with broadcasting them!
        DO jg = 1, n_dom
          IF(p_pe_work == 0) CALL p_recv(output_file_name(jg), p_work_pe0, 0)
          IF(.NOT.use_pio) &
           & CALL p_bcast(output_file_name(jg), 0, comm=p_comm_work)
        ENDDO
        IF(.NOT.use_pio) new_output_files = .TRUE.
      ENDIF

    CASE(msg_io_shutdown)
      done = .TRUE.

    CASE DEFAULT
      ! Anything else is an error
      CALL finish(modname,'I/O PE: Got illegal I/O tag')

    END SELECT

  END SUBROUTINE io_wait_for_start_message

  !-------------------------------------------------------------------------------------------------
  ! Simulate a barrier on all active PEs doing PIO

  SUBROUTINE pio_barrier

    INTEGER :: dummy(1)

    ! Do a broadcast, which is effectivly also a barrier
    CALL pio_bcast(dummy)

  END SUBROUTINE pio_barrier

  !-------------------------------------------------------------------------------------------------
  ! Simulate a broadcast on all active PEs doing PIO, root is I/O PE 0

  SUBROUTINE pio_bcast(iarray)

    INTEGER, INTENT(INOUT) :: iarray(:)
    INTEGER i

    ! First I/O PE 0 sends to all group leaders
    IF(my_pe_io_task == 0) THEN
      IF(my_io_task_no > 0) THEN
        CALL p_recv(iarray,0,0,comm=p_comm_work)
      ELSE
        DO i=1,num_io_tasks-1
          CALL p_send(iarray,io_task_pe0(i),0,comm=p_comm_work)
        ENDDO
      ENDIF
    ENDIF

    ! Then every group leader broadcasts within its group
    CALL p_bcast(iarray, 0, comm=pio_comm)

  END SUBROUTINE pio_bcast

  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  ! The following routines are called only from compute PEs
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------------------
  !>
  !! Send setup from compute PEs to I/O PEs

  SUBROUTINE setup_io_procs()

!     CHARACTER(LEN=*), INTENT(INOUT) :: gridfile(:)

    INTEGER :: jg, i, n, n_dims(3), nbytes_real, mpierr, nlev, type
    INTEGER (KIND=MPI_ADDRESS_KIND) :: mem_size, mem_bytes
    INTEGER, ALLOCATABLE :: ncheck(:,:)
#ifdef USE_CRAY_POINTER
    INTEGER (KIND=MPI_ADDRESS_KIND) :: iptr
    REAL(wp) :: tmp
    POINTER(tmp_ptr,tmp(*))
#else
    TYPE(c_ptr) :: c_mem_ptr
#endif

    IF(my_process_is_mpi_test()) THEN
      DO jg = 1, n_dom
        CALL setup_vlist( p_patch(jg)%grid_filename, jg )
      ENDDO
      RETURN
    ENDIF

    ! Send dimensions and owner distributions to the I/O PEs

    DO jg = 1, n_dom

      n_dims(1) = p_patch(jg)%n_patch_cells_g
      n_dims(2) = p_patch(jg)%n_patch_edges_g
      n_dims(3) = p_patch(jg)%n_patch_verts_g

      CALL p_bcast(n_dims, bcast_root(), p_comm_work_2_io)

      CALL p_bcast(p_patch(jg)%cells%owner_g, bcast_root(), p_comm_work_2_io)
      CALL p_bcast(p_patch(jg)%edges%owner_g, bcast_root(), p_comm_work_2_io)
      CALL p_bcast(p_patch(jg)%verts%owner_g, bcast_root(), p_comm_work_2_io)

      ! Send name of grid file for building vlist
      CALL p_bcast(p_patch(jg)%grid_filename, bcast_root(), p_comm_work_2_io)

    ENDDO

    ! Set up owner info for my patches

    ALLOCATE(patch_owner_info(n_dom))

    DO jg = 1, n_dom

      ! Get number of owned cells/edges/verts (without halos)
      ! and set index arrays to own cells

      ! cells

      n = COUNT(p_patch(jg)%cells%owner_g(:) == p_pe_work)
      patch_owner_info(jg)%n_own_cells = n
      ALLOCATE(patch_owner_info(jg)%own_cell_idx(n))
      ALLOCATE(patch_owner_info(jg)%own_cell_blk(n))

      n = 0
      DO i = 1, p_patch(jg)%n_patch_cells
        IF(p_patch(jg)%cells%owner_mask(idx_no(i),blk_no(i))) THEN
          n = n+1
          patch_owner_info(jg)%own_cell_idx(n) = idx_no(i)
          patch_owner_info(jg)%own_cell_blk(n) = blk_no(i)
        ENDIF
      ENDDO
      IF(n/=patch_owner_info(jg)%n_own_cells) CALL finish(modname,'Internal error cell owners')

      ! edges

      n = COUNT(p_patch(jg)%edges%owner_g(:) == p_pe_work)
      patch_owner_info(jg)%n_own_edges = n
      ALLOCATE(patch_owner_info(jg)%own_edge_idx(n))
      ALLOCATE(patch_owner_info(jg)%own_edge_blk(n))

      n = 0
      DO i = 1, p_patch(jg)%n_patch_edges
        IF(p_patch(jg)%edges%owner_mask(idx_no(i),blk_no(i))) THEN
          n = n+1
          patch_owner_info(jg)%own_edge_idx(n) = idx_no(i)
          patch_owner_info(jg)%own_edge_blk(n) = blk_no(i)
        ENDIF
      ENDDO
      IF(n/=patch_owner_info(jg)%n_own_edges) CALL finish(modname,'Internal error edge owners')

      ! verts

      n = COUNT(p_patch(jg)%verts%owner_g(:) == p_pe_work)
      patch_owner_info(jg)%n_own_verts = n
      ALLOCATE(patch_owner_info(jg)%own_vert_idx(n))
      ALLOCATE(patch_owner_info(jg)%own_vert_blk(n))

      n = 0
      DO i = 1, p_patch(jg)%n_patch_verts
        IF(p_patch(jg)%verts%owner_mask(idx_no(i),blk_no(i))) THEN
          n = n+1
          patch_owner_info(jg)%own_vert_idx(n) = idx_no(i)
          patch_owner_info(jg)%own_vert_blk(n) = blk_no(i)
        ENDIF
      ENDDO
      IF(n/=patch_owner_info(jg)%n_own_verts) CALL finish(modname,'Internal error vert owners')

    ENDDO

    ! Calculate the amount of memory needed to store all variables to be output

    mem_size = 0

    DO jg = 1, n_dom

      DO n = 1, num_output_vars(jg)

        nlev = outvar_desc(n,jg)%nlev
        type = outvar_desc(n,jg)%type

        ! Calculate total memory size
        SELECT CASE(type)
          CASE (GATHER_C); mem_size = mem_size + patch_owner_info(jg)%n_own_cells*nlev
          CASE (GATHER_E); mem_size = mem_size + patch_owner_info(jg)%n_own_edges*nlev
          CASE (GATHER_V); mem_size = mem_size + patch_owner_info(jg)%n_own_verts*nlev
          CASE DEFAULT
            CALL finish(modname, 'Illegal type from vlist_get_VarGrid')
        END SELECT

      ENDDO
    ENDDO

    ! mem_size is calculated as number of variables above, get number of bytes

    ! Get the amount of bytes per default REAL variable (as used in MPI communication)
    CALL MPI_Type_extent(p_real_dp, nbytes_real, mpierr)

    mem_bytes = mem_size*nbytes_real

    ! allocate amount of memory needed with MPI_Alloc_mem

#ifdef USE_CRAY_POINTER
    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, iptr, mpierr)

    tmp_ptr = iptr
    CALL set_mem_ptr(tmp, INT(mem_size))
#else
    ! TYPE(c_ptr) and INTEGER(KIND=MPI_ADDRESS_KIND) do NOT necessarily have the same size!!!
    ! So check if at least c_intptr_t and MPI_ADDRESS_KIND are the same, else we may get
    ! into deep, deep troubles!
    ! There is still a slight probability that TYPE(c_ptr) does not have the size indicated
    ! by c_intptr_t since the standard only requires c_intptr_t is big enough to hold pointers
    ! (so it may be bigger than a pointer), but I hope no vendor screws up its ISO_C_BINDING
    ! in such a way!!!
    ! If c_intptr_t<=0, this type is not defined and we can't do this check, of course.

    IF(c_intptr_t > 0 .AND. c_intptr_t /= MPI_ADDRESS_KIND) &
     & CALL finish(modname,'c_intptr_t /= MPI_ADDRESS_KIND, too dangerous to proceed!')

    CALL MPI_Alloc_mem(mem_bytes, MPI_INFO_NULL, c_mem_ptr, mpierr)

    ! The NEC requires a standard INTEGER array as 3rd argument for c_f_pointer,
    ! although it would make more sense to have it of size MPI_ADDRESS_KIND.

    CALL c_f_pointer(c_mem_ptr, mem_ptr, (/ INT(mem_size) /) )
#endif

    mem_ptr(:) = 0._wp

    ! Create memory window for communication

    CALL MPI_Win_create(mem_ptr,mem_bytes,nbytes_real,MPI_INFO_NULL,p_comm_work_io,mpi_win,mpierr)

    ! Send data for safety check to I/O PEs

    DO jg = 1, n_dom

      CALL p_bcast(num_output_vars(jg), bcast_root(), p_comm_work_2_io)

      ALLOCATE(ncheck(2,num_output_vars(jg)))
      DO n=1, num_output_vars(jg)
        ncheck(1,n) = outvar_desc(n,jg)%type
        ncheck(2,n) = outvar_desc(n,jg)%nlev
      ENDDO
      CALL p_bcast(ncheck, bcast_root(), p_comm_work_2_io)
      DEALLOCATE(ncheck)

    ENDDO

    ! Send first set of files to be opened.
    ! output_file_name must be set before this routine is called!

    CALL p_bcast(output_file_name(1:n_dom), bcast_root(), p_comm_work_2_io)
    new_output_files = .FALSE.

  END SUBROUTINE setup_io_procs

#ifdef USE_CRAY_POINTER
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE set_mem_ptr(arr, len)

    INTEGER len
    REAL(wp), TARGET :: arr(len)

    mem_ptr => arr

  END SUBROUTINE set_mem_ptr
#endif

  !-------------------------------------------------------------------------------------------------
  !>
  !! Shutdown of I/O procs (called from compute PEs)

  SUBROUTINE shutdown_io_procs

    INTEGER jg

    IF(my_process_is_mpi_test()) THEN
      DO jg = 1, n_dom
        CALL close_output_vlist(jg)
      ENDDO
      RETURN
    ENDIF

    ! Receive message that I/O is finished

    CALL compute_wait_for_io_ready

    ! Send message to shut down

    CALL compute_shutdown_io

  END SUBROUTINE shutdown_io_procs

  !-------------------------------------------------------------------------------------------------
  !>
  !! Sets new output file name (to be used at next output)

  SUBROUTINE set_output_file(file, jg)

    CHARACTER(LEN=*), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: jg

    output_file_name(jg) = file

    ! Please note that if this routine is called for one output file, it must be called for
    ! ALL output files since the following triggers an open/close for all files:
    new_output_files = .TRUE.

  END SUBROUTINE set_output_file

  !-------------------------------------------------------------------------------------------------
  !>
  !! Async output routine

  SUBROUTINE output_async(datetime, z_sim_time)

    TYPE(t_datetime),   INTENT(in) :: datetime
    REAL(wp), OPTIONAL, INTENT(in) :: z_sim_time(n_dom)

    INTEGER jg, jk, n, i, mpierr, n_own, n_tot, nlev_ptr, nblk_ptr
    INTEGER (KIND=MPI_ADDRESS_KIND) :: ioff ! If amount of data should exceed 2 GB
    INTEGER, POINTER :: iidx(:), iblk(:)

    LOGICAL :: reset, delete

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)
    REAL(wp)          :: p_sim_time
    CHARACTER*10 ctime


    IF(my_process_is_mpi_test()) THEN
      CALL test_pe_output(datetime, z_sim_time)
      RETURN
    ENDIF

    ! Wait for I/O PEs to get ready with previous I/O

    CALL date_and_time(TIME=ctime)
    IF(p_pe_work==0) PRINT '(a)','.................... Compute PEs ready for I/O at '//ctime
    CALL compute_wait_for_io_ready
    CALL date_and_time(TIME=ctime)
    IF(p_pe_work==0) PRINT '(a)','.................... Compute PEs starting I/O  at '//ctime

    ! Lock own window before writing to it

    CALL MPI_Win_lock(MPI_LOCK_EXCLUSIVE, p_pe_work, MPI_MODE_NOCHECK, mpi_win, mpierr)

    ! Offset into memory window

    ioff = 0_MPI_ADDRESS_KIND

    ! Go over all patches and fill memory window with data to be output

    DO jg = 1, n_dom

      IF (PRESENT(z_sim_time)) THEN
        p_sim_time = z_sim_time(jg)
      ELSE
        p_sim_time = 1.0_wp
      ENDIF
      
      DO n = 1, num_output_vars(jg)

        ! Set ptr2/ptr3 to the variable to be output

        SELECT CASE (iequations)
          CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
            CALL get_outvar_ptr_ha(outvar_desc(n,jg)%name, jg, ptr2, ptr3, reset, delete)
          CASE (inh_atmosphere)
            CALL get_outvar_ptr_nh(outvar_desc(n,jg)%name,jg,p_sim_time,ptr2, ptr3, reset, delete)
          CASE DEFAULT
            CALL finish(modname,'Unsupported value of iequations')
        END SELECT


        ! Get the number of own points and index arrays to be used

        SELECT CASE(outvar_desc(n,jg)%type)
          CASE (GATHER_C)
            n_own = patch_owner_info(jg)%n_own_cells
            iidx => patch_owner_info(jg)%own_cell_idx
            iblk => patch_owner_info(jg)%own_cell_blk
            n_tot = p_patch(jg)%n_patch_cells
          CASE (GATHER_E)
            n_own = patch_owner_info(jg)%n_own_edges
            iidx => patch_owner_info(jg)%own_edge_idx
            iblk => patch_owner_info(jg)%own_edge_blk
            n_tot = p_patch(jg)%n_patch_edges
          CASE (GATHER_V)
            n_own = patch_owner_info(jg)%n_own_verts
            iidx => patch_owner_info(jg)%own_vert_idx
            iblk => patch_owner_info(jg)%own_vert_blk
            n_tot = p_patch(jg)%n_patch_verts
          CASE DEFAULT
            CALL finish(modname, 'Illegal type in outvar_desc')
        END SELECT

        ! Make some checks if the table is set up correctly

        IF(ASSOCIATED(ptr3)) THEN
          nlev_ptr = UBOUND(ptr3,2)
          nblk_ptr = UBOUND(ptr3,3)
        ELSE
          nlev_ptr = 1
          nblk_ptr = UBOUND(ptr2,2)
        ENDIF

        IF(outvar_desc(n,jg)%nlev /= nlev_ptr) &
          & CALL finish(modname, 'Incorrect nlev in outvar_desc for '//outvar_desc(n,jg)%name)

        IF(blk_no(n_tot) /= nblk_ptr) &
          & CALL finish(modname, 'Incorrect type in outvar_desc for '//outvar_desc(n,jg)%name)

        ! Copy data into memory window (layer by layer for 3D arrays)

        IF(ASSOCIATED(ptr3)) THEN
          DO jk = 1, outvar_desc(n,jg)%nlev
            DO i = 1, n_own
              mem_ptr(i+ioff) = ptr3(iidx(i),jk,iblk(i))
            ENDDO
            ioff = ioff + n_own
          ENDDO
          IF(reset) ptr3 = 0._wp
          IF(delete) DEALLOCATE(ptr3)
        ELSE
          DO i = 1, n_own
            mem_ptr(i+ioff) = ptr2(iidx(i),iblk(i))
          ENDDO
          ioff = ioff + n_own
          IF(reset) ptr2 = 0._wp
          IF(delete) DEALLOCATE(ptr2)
        ENDIF

      ENDDO
    ENDDO

    ! Done writing to memory window, unlock it

    CALL MPI_Win_unlock(p_pe_work, mpi_win, mpierr)

    ! Send message that output is ready

    CALL compute_start_io(datetime)
    CALL date_and_time(TIME=ctime)
    IF(p_pe_work==0) PRINT '(a)','.................... Compute PEs done with I/O at '//ctime

  END SUBROUTINE output_async

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_wait_for_io_ready: Wait for a message that the I/O is ready
  !! The counterpart on the I/O side is io_send_ready_message

  SUBROUTINE compute_wait_for_io_ready

    INTEGER msg

    ! First compute PE receives message from I/O leader
    IF(p_pe_work==0) THEN
      CALL p_recv(msg, p_io_pe0, 0)
      ! Just for safety: Check if we got the correct tag
      IF(msg /= msg_io_done) CALL finish(modname,'Compute PE: Got illegal I/O tag')
    ENDIF

    ! Wait in barrier until message is here
    CALL p_barrier(comm=p_comm_work)

  END SUBROUTINE compute_wait_for_io_ready

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_start_io: Send a message to I/O PEs that they should start I/O
  !! The counterpart on the I/O side is io_wait_for_start_message

  SUBROUTINE compute_start_io(datetime)

    TYPE(t_datetime), INTENT(in) :: datetime

    INTEGER msg(4), jg

    CALL p_barrier(comm=p_comm_work) ! make sure all are here

    msg(1) = msg_io_start
    msg(2) = cdiEncodeDate(datetime%year, datetime%month, datetime%day)
    msg(3) = cdiEncodeTime(datetime%hour, datetime%minute, NINT(datetime%second))
    msg(4) = MERGE(1,0,new_output_files)

    IF(p_pe_work==0) CALL p_send(msg, p_io_pe0, 0)

    IF(new_output_files) THEN
      ! Send new output files
      DO jg = 1, n_dom
        IF(p_pe_work==0) CALL p_send(output_file_name(jg), p_io_pe0, 0)
      ENDDO
      new_output_files = .FALSE.
    ENDIF

  END SUBROUTINE compute_start_io

  !-------------------------------------------------------------------------------------------------
  !>
  !! compute_shutdown_io: Send a message to I/O PEs that they should shut down
  !! The counterpart on the I/O side is io_wait_for_start_message

  SUBROUTINE compute_shutdown_io

    INTEGER msg(4)

    CALL p_barrier(comm=p_comm_work) ! make sure all are here

    msg(1) = msg_io_shutdown
    msg(2:) = 0

    IF(p_pe_work==0) CALL p_send(msg, p_io_pe0, 0)

  END SUBROUTINE compute_shutdown_io

  !-------------------------------------------------------------------------------------------------
  !>
  !! Verification output on test PE

  SUBROUTINE test_pe_output(datetime, z_sim_time)

    TYPE(t_datetime),   INTENT(in) :: datetime
    REAL(wp), OPTIONAL, INTENT(in) :: z_sim_time(n_dom)

    INTEGER jg, jk, n, i, n_tot, idate, itime
    INTEGER (KIND=MPI_ADDRESS_KIND) :: ioff ! If amount of data should exceed 2 GB

    LOGICAL :: reset, delete

    REAL(wp), POINTER :: ptr2(:,:)
    REAL(wp), POINTER :: ptr3(:,:,:)
    REAL(wp), ALLOCATABLE :: var(:)
    REAL(wp)          :: p_sim_time

    INTEGER, SAVE :: nstep = 0

    nstep = nstep + 1

    ! Check if new output files have to be opened

    IF(new_output_files) THEN

      DO jg = 1, n_dom
        IF(nstep>1) CALL close_output_vlist(jg)
        CALL open_output_vlist('TEST-'//output_file_name(jg), jg)
      ENDDO
      new_output_files = .FALSE.

    ENDIF

    idate = cdiEncodeDate(datetime%year, datetime%month, datetime%day)
    itime = cdiEncodeTime(datetime%hour, datetime%minute, NINT(datetime%second))
    CALL vlist_start_step(idate, itime)

    ! Go over all patches and output data

    DO jg = 1, n_dom

      IF (PRESENT(z_sim_time)) THEN
        p_sim_time = z_sim_time(jg)
      ELSE
        p_sim_time = 1.0_wp
      ENDIF
      
      DO n = 1, num_output_vars(jg)

        ! Set ptr2/ptr3 to the variable to be output

        SELECT CASE (iequations)
          CASE (ishallow_water, ihs_atm_temp, ihs_atm_theta)
            CALL get_outvar_ptr_ha(outvar_desc(n,jg)%name, jg, ptr2, ptr3, reset, delete)
          CASE (inh_atmosphere)
            CALL get_outvar_ptr_nh(outvar_desc(n,jg)%name,jg,p_sim_time,ptr2, ptr3, reset, delete)
          CASE DEFAULT
            CALL finish(modname,'Unsupported value of iequations')
        END SELECT


        ! Get the number of values

        SELECT CASE(outvar_desc(n,jg)%type)
          CASE (GATHER_C)
            n_tot = p_patch(jg)%n_patch_cells
          CASE (GATHER_E)
            n_tot = p_patch(jg)%n_patch_edges
          CASE (GATHER_V)
            n_tot = p_patch(jg)%n_patch_verts
          CASE DEFAULT
            CALL finish(modname, 'Illegal type in outvar_desc')
        END SELECT

        ! Output data

        ALLOCATE(var(n_tot*outvar_desc(n,jg)%nlev))

        IF(ASSOCIATED(ptr3)) THEN
          ioff = 0
          DO jk = 1, outvar_desc(n,jg)%nlev
            DO i = 1, n_tot
              var(i+ioff) = ptr3(idx_no(i),jk,blk_no(i))
            ENDDO
            ioff = ioff + n_tot
          ENDDO
          IF(reset) ptr3 = 0._wp
          IF(delete) DEALLOCATE(ptr3)
        ELSE
          DO i = 1, n_tot
            var(i) = ptr2(idx_no(i),blk_no(i))
          ENDDO
          IF(reset) ptr2 = 0._wp
          IF(delete) DEALLOCATE(ptr2)
        ENDIF

        CALL vlist_write_var(n, jg, var)

        DEALLOCATE(var)

      ENDDO
    ENDDO

  END SUBROUTINE test_pe_output

  ! ------------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------------
  ! Remove the following when actually working with pio
  ! ------------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------------

  INTEGER FUNCTION pioInit(type, comm, my_rank, num_tasks, comm_out)
    INTEGER type, comm, my_rank, num_tasks, comm_out

    IF(p_pe_work==0) THEN
      print *
      print *
      print *,'ATTENTION:'
      print *,'If you want ro run with more I/O PEs than domains,'
      print *,'you have to remove the dummy rountines pioInit/pioFinalize'
      print *,'at the end of module mo_io_async and link with an actual'
      print *,'PIO enabled CDI version'
      CALL finish(modname,'CDI PIO not enabled')
    ENDIF

    CALL p_barrier ! prevent other PEs from going on

    ! Make the compiler happy:
    pioInit = type*comm*my_rank*num_tasks*comm_out

  END FUNCTION pioInit

  SUBROUTINE pioFinalize
  END SUBROUTINE pioFinalize

#endif

END MODULE mo_io_async

