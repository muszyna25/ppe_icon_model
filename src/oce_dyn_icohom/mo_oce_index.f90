!>
!!
!! The module <i>mo_oce_index</i> calculates indices of single cells
!! and neighbouring grids for test output of single values in the ocean core
!!
!! @par Revision History
!! Initial version by Stephan Lorenz,  MPI-M, Hamburg (2010-11)
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
MODULE mo_oce_index
!-------------------------------------------------------------------------
!
USE mo_kind,                   ONLY: wp
USE mo_mpi,                    ONLY: my_process_is_stdio
USE mo_io_units,               ONLY: nerr
USE mo_parallel_config,        ONLY: nproma, p_test_run
USE mo_grid_config,            ONLY: n_dom
USE mo_run_config,             ONLY: nsteps!, ltimer
!USE mo_timer,                  ONLY: timer_start, timer_stop, timer_print_mxmn
USE mo_sync,                   ONLY: global_max, global_min
USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
USE mo_ocean_nml,              ONLY: n_zlev, i_dbg_oce, i_dbg_inx, str_proc_tst,no_tracer, &
  &                                  i_oct_ilv, rlon_in, rlat_in
! &                                  i_oct_blk, i_oct_idx, i_oct_ilv, rlon_in, rlat_in
! &                                  i_ocv_blk, i_ocv_idx, i_ocv_ilv, t_val,  &
USE mo_dynamics_config,        ONLY: nold, nnew
USE mo_run_config,             ONLY: dtime
USE mo_loopindices,            ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_impl_constants,         ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
  &                                  min_rlcell, min_rledge, min_rlvert,                &
  &                                  max_char_length
USE mo_math_constants,         ONLY: pi
USE mo_exception,              ONLY: message, message_text, finish
USE mo_model_domain,           ONLY: t_patch
USE mo_ext_data_types,         ONLY: t_external_data
USE mo_oce_state,              ONLY: t_hydro_ocean_state, v_base


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

! indices of cells and neighbours for debug output set by namelist octst_ctl
INTEGER :: c_b, c_i, c_k, ne_b(3), ne_i(3), nc_b(3), nc_i(3), nv_b(3), nv_i(3)
INTEGER :: jkc, jkdim, ipl_src
INTEGER :: loc_nblks_c, loc_nblks_e, loc_nblks_v

! format variables for degug output
CHARACTER(len=95) :: form4ar

! control debug output
LOGICAL :: ldbg

PUBLIC :: init_index_test

! Public variables:
PUBLIC :: c_b, c_i, c_k, ne_b, ne_i, nc_b, nc_i, nv_b, nv_i
PUBLIC :: form4ar
PUBLIC :: ldbg
PUBLIC :: jkc, jkdim, ipl_src


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialization of indices
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-11)
  !!
  !! TODO: parallelize
  !! TODO: cleanup and move to shared
  !
  !
  SUBROUTINE init_index_test (ppatch, pstate_oce, p_ext_data)


  TYPE(t_patch),             TARGET, INTENT(IN)     :: ppatch(n_dom)
  TYPE(t_hydro_ocean_state), TARGET, INTENT(INOUT)  :: pstate_oce(n_dom)
  TYPE(t_external_data),     TARGET, INTENT(IN)     :: p_ext_data(n_dom)

  INTEGER :: jg, jt, i, islmval, idolic
  !INTEGER :: jb, ji
  INTEGER :: rl_start, rl_end, i_startblk, i_endblk, &
    &        i_startidxf, i_endidxf, i_startidxl, i_endidxl
  REAL(wp) :: zlon, zlat, bathy

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_oce_index:init_index_test'

  !CALL message(TRIM(routine), 'Start' )

  IF (n_dom > 1 ) THEN
    CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
  END IF
  jg = n_dom

  ! set time level to old level (nold(1)=3)
  jt = nold(jg)

  ! search for block/index of debug output cell at given lat/lon
  CALL search_latlonindex (ppatch, rlat_in, rlon_in, c_i, c_b)

  ! overwrite local block and index of test point to values of namelist
  !IF (i_oct_blk /= 0) THEN
  !  c_b = i_oct_blk
  !  c_i = i_oct_idx
  !END IF

  ! set zlevel of test point
  c_k = i_oct_ilv

  ! set level of debug output
  ldbg = .FALSE.
  ! #slo# Bugfix - more necessary
  !IF (i_dbg_oce >=1) ldbg = .TRUE.
  ldbg = .TRUE.


  IF (ldbg) THEN

    WRITE(message_text,'(a)') 'Output diverse variables:'
    CALL message (TRIM(routine), message_text)

    WRITE(message_text,'(3(a,i8))')          &
      &   '  nzlev   =', n_zlev,             &
      &   ', n_dom   =', n_dom,              &
      &   ', nproma  =', nproma
    CALL message ('', message_text)
    WRITE(message_text,'(3(a,i8))')          &
      &   '  nblks_c =', ppatch(jg)%nblks_c, &
      &   ', nblks_e =', ppatch(jg)%nblks_e, &
      &   ', nblks_v =', ppatch(jg)%nblks_v
    CALL message ('', message_text)
    WRITE(message_text,'(3(a,i8))')          &
      &   '  nold(jg)=', nold(jg),           &
      &   ', nnew(jg)=', nnew(jg),           &
      &   ', Nsteps  =', nsteps
    CALL message ('', message_text)
    WRITE(message_text,'(5(a,i8))')          &
      &   '  SEA=',     sea,                 &
      &   ', SEA_BOUNDARY=',sea_boundary,    &
      &   ', BOUNDARY=',boundary,            &
      &   ', LAND_BOUNDARY=',land_boundary,  &
      &   ', LAND=,',   land!,                &
    CALL message ('', message_text)
    WRITE(message_text,'(4(a,g18.6))')       &
      &   '  Time Step=',dtime
    CALL message ('', message_text)

    !------------------------------------------------------------------
    ! Test output of indices using rl_start = 1 and 2
    !  - differ if landpoints are not included in the grid
    !------------------------------------------------------------------

    ! values for the blocking: cells
    rl_start = 1
    rl_end = min_rlcell
    i_startblk = ppatch(jg)%cells%start_blk(rl_start,1)
    i_endblk   = ppatch(jg)%cells%end_blk(rl_end,1)
    ! first block
    CALL get_indices_c(ppatch(jg), i_startblk, i_startblk, i_endblk, i_startidxf, i_endidxf, &
      &                            rl_start, rl_end)
    ! last block
    CALL get_indices_c(ppatch(jg), i_endblk, i_startblk, i_endblk, i_startidxl, i_endidxl, &
      &                            rl_start, rl_end)
    WRITE(message_text,'(8(a,i5))')                         &
      &   ' Cells: rl_start=',rl_start,' rl_end=',rl_end,   &
      &   ' stblk=',i_startblk,' endblk=',i_endblk,         &
      &   ' fst blk: sidx=',i_startidxf,' eidx=',i_endidxf, &
      &   ' lst blk: sidx=',i_startidxl,' eidx=',i_endidxl
    CALL message ('', message_text)

    ! values for the blocking: edges
    rl_end = min_rledge
    i_startblk = ppatch(jg)%edges%start_blk(rl_start,1)
    i_endblk   = ppatch(jg)%edges%end_blk(rl_end,1)
    ! first block
    CALL get_indices_e(ppatch(jg), i_startblk, i_startblk, i_endblk, i_startidxf, i_endidxf, &
      &                            rl_start, rl_end)
    ! last block
    CALL get_indices_e(ppatch(jg), i_endblk, i_startblk, i_endblk, i_startidxl, i_endidxl, &
      &                            rl_start, rl_end)
    WRITE(message_text,'(8(a,i5))')                         &
      &   ' Edges: rl_start=',rl_start,' rl_end=',rl_end,   &
      &   ' stblk=',i_startblk,' endblk=',i_endblk,         &
      &   ' fst blk: sidx=',i_startidxf,' eidx=',i_endidxf, &
      &   ' lst blk: sidx=',i_startidxl,' eidx=',i_endidxl
    CALL message ('', message_text)

    ! values for the blocking: verts
    rl_end = min_rlvert
    i_startblk = ppatch(jg)%verts%start_blk(rl_start,1)
    i_endblk   = ppatch(jg)%verts%end_blk(rl_end,1)
    ! first block
    CALL get_indices_v(ppatch(jg), i_startblk, i_startblk, i_endblk, i_startidxf, i_endidxf, &
      &                            rl_start, rl_end)
    ! last block
    CALL get_indices_v(ppatch(jg), i_endblk, i_startblk, i_endblk, i_startidxl, i_endidxl, &
      &                            rl_start, rl_end)
    WRITE(message_text,'(8(a,i5))')                         &
      &   ' Verts: rl_start=',rl_start,' rl_end=',rl_end,   &
      &   ' stblk=',i_startblk,' endblk=',i_endblk,         &
      &   ' fst blk: sidx=',i_startidxf,' eidx=',i_endidxf, &
      &   ' lst blk: sidx=',i_startidxl,' eidx=',i_endidxl
    CALL message ('', message_text)

    ! values for the blocking using rl_start=2
    rl_start = 2
    rl_end = min_rlcell
    i_startblk = ppatch(jg)%cells%start_blk(rl_start,1)
    i_endblk   = ppatch(jg)%cells%end_blk(rl_end,1)
    WRITE(message_text,'(8(a,i5))')                     &
      &   ' Cells: rl_start=',rl_start,' rl_end=',rl_end, &
      &   ' stblk=',i_startblk,' endblk=',i_endblk
    CALL message ('', message_text)
    rl_end = min_rledge
    i_startblk = ppatch(jg)%edges%start_blk(rl_start,1)
    i_endblk   = ppatch(jg)%edges%end_blk(rl_end,1)
    WRITE(message_text,'(8(a,i5))')                     &
      &   ' Edges: rl_start=',rl_start,' rl_end=',rl_end, &
      &   ' stblk=',i_startblk,' endblk=',i_endblk
    CALL message ('', message_text)
    rl_end = min_rlvert
    i_startblk = ppatch(jg)%verts%start_blk(rl_start,1)
    i_endblk   = ppatch(jg)%verts%end_blk(rl_end,1)
    WRITE(message_text,'(8(a,i5))')                     &
      &   ' Verts: rl_start=',rl_start,' rl_end=',rl_end, &
      &   ' stblk=',i_startblk,' endblk=',i_endblk
    CALL message ('', message_text)

    !------------------------------------------------------------------
    ! Check parameters
    !------------------------------------------------------------------

    IF (c_k > n_zlev) THEN
      CALL finish(TRIM(routine), ' N_DOM > 1 is not allowed')
    END IF

    ! slm and coordinates of this point:
    islmval = v_base%lsm_oce_c(c_i,c_k,c_b)
    idolic  = v_base%dolic_c  (c_i,    c_b)
    bathy   = p_ext_data(jg)%oce%bathymetry_c (c_i,c_b)
    zlat = ppatch(jg)%cells%center(c_i,c_b)%lat * 180.0_wp / pi
    zlon = ppatch(jg)%cells%center(c_i,c_b)%lon * 180.0_wp / pi

    ! output format
    99 FORMAT(a,i4,a,i4,a,i3,a,i3,3(a,f9.2))
    97 FORMAT(a,i1,a,i4,a,i4,a,i3,a,i3,3(a,f9.2))
    form4ar = '(4(a,g20.9))'

    CALL message (TRIM(routine), 'Conditions at test cell (C), and edges/neighbouring cells:')
    WRITE(message_text,99) ' Cell C: block=',c_b,'  index=',c_i,             &
      &         '  lsm_c=', islmval,'  dolic_c=',idolic,'  bathy_c=', bathy, &
      &         '  lat=',zlat,'  lon=',zlon
    CALL message (' ', message_text)
    IF(no_tracer>=1)THEN
      WRITE(message_text,form4ar)  &
                ' Elev. h at Cell C    =', pstate_oce(jg)%p_prog(jt)%h(c_i,c_b),            &
        &                  '  Tracer 1 =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1), &
        &      '  Test level at cell C = ',REAL(c_k,wp)
      CALL message (' ', message_text)
    ELSEIF(no_tracer==0)THEN
      WRITE(message_text,form4ar)  &
                ' Elev. h at Cell C    =', pstate_oce(jg)%p_prog(jt)%h(c_i,c_b),            &
        !&                  '  Tracer 1 =', pstate_oce(jg)%p_prog(jt)%tracer(c_i,c_k,c_b,1), &
        &      '  Test level (jk) at cell C = ',REAL(c_k,wp)
      CALL message (' ', message_text)
    ENDIF

    !------------------------------------------------------------------
    ! find and print corresponding edges/verts of test cell
    !------------------------------------------------------------------

    DO i = 1, 3 ! 3 edges of cell C at (ne_i,ne_b)
      ! slm and coordinates of edges
      ne_b(i)=ppatch(jg)%cells%edge_blk(c_i,c_b,i)
      ne_i(i)=ppatch(jg)%cells%edge_idx(c_i,c_b,i)
      islmval = v_base%lsm_oce_e  (ne_i(i),c_k,ne_b(i))
      idolic  = v_base%dolic_e    (ne_i(i),ne_b(i))
      bathy   = p_ext_data(jg)%oce%bathymetry_e (ne_i(i),ne_b(i))
      zlat    = ppatch(jg)%edges%center         (ne_i(i),ne_b(i))%lat * 180.0_wp / pi
      zlon    = ppatch(jg)%edges%center         (ne_i(i),ne_b(i))%lon * 180.0_wp / pi
      ! output
      WRITE(message_text,97) ' Edge E',i,' block=',ne_b(i),'  index=',ne_i(i),              &
        &                    '  lsm_e=', islmval,'  dolic_e=',idolic,'  bathy_e=', bathy, &
        &                    '  lat=',zlat,'  lon=',zlon
      CALL message (' ', message_text)
    END DO

    DO i = 1, 3 ! 3 vertices of cell C at (nv_i,nv_b)
      ! slm and coordinates of vertices
      nv_b(i)=ppatch(jg)%cells%vertex_blk(c_i,c_b,i)
      nv_i(i)=ppatch(jg)%cells%vertex_idx(c_i,c_b,i)
      islmval = v_base%lsm_oce_c(c_i,c_k,c_b)
      idolic  = v_base%dolic_c  (c_i,    c_b)
      bathy   = p_ext_data(jg)%oce%bathymetry_c (c_i,c_b)
      zlat    = ppatch(jg)%edges%center         (nv_i(i),nv_b(i))%lat * 180.0_wp / pi
      zlon    = ppatch(jg)%edges%center         (nv_i(i),nv_b(i))%lon * 180.0_wp / pi
      ! output
      WRITE(message_text,97) ' Vert V',i,' block=',nv_b(i),'  index=',nv_i(i),              &
        &                    '  lsm_c=', islmval,'  dolic_c=',idolic,'  bathy_c=', bathy, &
        &                    '  lat=',zlat,'  lon=',zlon
      CALL message (' ', message_text)
    END DO

    DO i = 1, 3 ! 3 neighbours of cell C at (nc_i,nc_b)
      ! slm and coordinates of neighbouring cells
      ! #slo# - careful at index-boundaries: in ocean grid exist cells without neighbours !!
      nc_b(i)=ppatch(jg)%cells%neighbor_blk(c_i,c_b,i)
      nc_i(i)=ppatch(jg)%cells%neighbor_idx(c_i,c_b,i)
      IF ( nc_i(i) == 0 .OR. nc_b(i) == 0) THEN
        nc_i(i) = c_i
        nc_b(i) = c_b
        WRITE(message_text,'(a)') ' Neighbor Cell is on LAND - NOT DEFINED'
      ELSE
        islmval = v_base%lsm_oce_c  (nc_i(i),c_k,nc_b(i))
        idolic  = v_base%dolic_c    (nc_i(i),    nc_b(i))
        bathy   = p_ext_data(jg)%oce%bathymetry_c (nc_i(i),nc_b(i))
        zlat    = ppatch(jg)%cells%center         (nc_i(i),nc_b(i))%lat * 180.0_wp / pi
        zlon    = ppatch(jg)%cells%center         (nc_i(i),nc_b(i))%lon * 180.0_wp / pi
        WRITE(message_text,97) ' Neighbor  C',i,' =',nc_b(i),'  index=',nc_i(i),            &
          &                    '  lsm_c=', islmval,'  dolic_c=',idolic,'  bathy_c=', bathy, &
          &                    '  lat=',zlat,'  lon=',zlon
      END IF
      ! output
      CALL message (' ', message_text)
    END DO

  END IF

  loc_nblks_c =ppatch(jg)%nblks_c
  loc_nblks_e =ppatch(jg)%nblks_e
  loc_nblks_v =ppatch(jg)%nblks_v

  END SUBROUTINE init_index_test

  !-------------------------------------------------------------------------
  !>
  !! Search for a cell center at given longitude and latitude
  !! provided in namelist octst_nml
  !! 
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-12)
  !!
  !! TODO: parallelize
  !! TODO: cleanup and move to shared
  !
  !
  SUBROUTINE search_latlonindex (ppatch, plat_in, plon_in, iidx, iblk)

  TYPE(t_patch), TARGET, INTENT(IN)     :: ppatch(n_dom)
  REAL(wp),              INTENT(IN)     :: plat_in       ! cell latitude to search for
  REAL(wp),              INTENT(IN)     :: plon_in       ! cell longitude to search for
  INTEGER,               INTENT(OUT)    :: iidx          ! index of nearest cell
  INTEGER,               INTENT(OUT)    :: iblk          ! block of nearest cell

  INTEGER  :: jb, jc, jg, i_startidx, i_endidx, proc_id
  REAL(wp) :: zlon, zlat, zdist, zdist_cmp, ctr
  REAL(wp) :: zdst(nproma,ppatch(1)%nblks_c)
  TYPE(t_subset_range), POINTER :: all_cells
  LOGICAL  :: p_test_run_bac

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_oce_index:search_latlonindex'

  CALL message(TRIM(routine), 'Start' )

  jg = n_dom

  all_cells => ppatch(jg)%cells%all

  ! initial distance to compare
  zdist_cmp = 100000.0_wp
  zdst(:,:) = 100000.0_wp
  proc_id   = -1

  !  loop over all cells including halo
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx

      zlat    = ppatch(jg)%cells%center(jc,jb)%lat * 180.0_wp / pi
      zlon    = ppatch(jg)%cells%center(jc,jb)%lon * 180.0_wp / pi

      zdist       = sqrt((zlat-plat_in)*(zlat-plat_in) + (zlon-plon_in)*(zlon-plon_in))
      zdst(jc,jb) = zdist
      IF (zdist < zdist_cmp) THEN
        iblk = jb
        iidx = jc
        zdist_cmp = zdist
      END IF

    END DO
  END DO

  ! find PE with minimum distance, MPI-broadcast block/index from that PE - not yet

  ! disable p_test_run since global_max will be different
  p_test_run_bac = p_test_run
  p_test_run = .false.
  ctr = global_max(-zdist,proc_id)
  p_test_run = p_test_run_bac

  zlat    = ppatch(jg)%cells%center          (iidx,iblk)%lat * 180.0_wp / pi
  zlon    = ppatch(jg)%cells%center          (iidx,iblk)%lon * 180.0_wp / pi

  99 FORMAT(3a,i4,a,i4,3(a,f9.2))
  98 FORMAT(2a,3(a,f9.2))
  IF (my_process_is_stdio()) THEN
    WRITE(0,98) ' ',TRIM(routine),' Found cell nearest to          lat=', plat_in,'  lon=',plon_in
    WRITE(0,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  lat=',zlat,'  lon=',zlon
    WRITE(0,'(3a,i3)')         ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',proc_id
    WRITE(0,'(3a,2i3,a,f9.2)') ' ',TRIM(routine),' FOUND: Min dist is at idx/blk=', &
      &                  MINLOC(zdst(:,:)),' distance in deg =',MINVAL(zdst(:,:))
  END IF

  END SUBROUTINE search_latlonindex

END MODULE mo_oce_index

