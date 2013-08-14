!>
!!
!! The module <i>mo_util_dbg_prnt</i> prints out max and min as well as single
!! cell and neighbouring values of 2- and 3-dim arrays of the icon core for
!! debug purposes
!!
!! @par Revision History
!! Initial version by Stephan Lorenz,  MPI-M, Hamburg (2010-11)
!!
!! @par Revision History
!! Modified for general purpose by Stephan Lorenz, MPI-M, 2012-06
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
MODULE mo_util_dbg_prnt
!-------------------------------------------------------------------------
!
USE mo_kind,                   ONLY: wp
USE mo_mpi,                    ONLY: my_process_is_stdio, p_pe, get_my_mpi_all_id !, p_comm_work, p_bcast
USE mo_io_units,               ONLY: nerr
USE mo_parallel_config,        ONLY: nproma, p_test_run
USE mo_impl_constants,         ONLY: max_char_length
USE mo_timer,                  ONLY: timer_start, timer_stop, timer_dbg_prnt, timers_level
USE mo_sync,                   ONLY: SYNC_C, sync_patch_array, global_max
USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
USE mo_dbg_nml,                ONLY: str_mod_tst, dim_mod_tst, dbg_lon_in, dbg_lat_in, &
  &                                  idbg_mxmn, idbg_val, idbg_slev, idbg_elev, &
  &                                  idbg_idx, idbg_blk
USE mo_math_constants,         ONLY: pi
USE mo_exception,              ONLY: message, message_text
USE mo_model_domain,           ONLY: t_patch
USE mo_grid_subset,            ONLY: t_subset_range
USE mo_statistics,             ONLY: global_minmaxmean
USE mo_icon_comm_interface,    ONLY: icon_comm_barrier

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

! indices of cells and neighbours for debug output at single cell
INTEGER :: c_b, c_i, ne_b(3), ne_i(3), nc_b(3), nc_i(3), nv_b(3), nv_i(3), near_proc_id
INTEGER :: loc_nblks_c, loc_nblks_e, loc_nblks_v
LOGICAL :: p_test_run_bac
!TYPE(t_subset_range) :: v_subdom_cell !, v_suball_cell, v_subset_edge  !  part of subset to store

!
! PUBLIC INTERFACE
!

! Public subroutines:
PUBLIC :: init_dbg_index
PUBLIC :: dbg_print

! Public variables:
PUBLIC :: c_i, c_b, nc_i, nc_b
!PUBLIC :: v_subdom_cell !, v_suball_cell, v_subset_edge  !  part of subset to store

INTERFACE dbg_print
  MODULE PROCEDURE dbg_print_2d
  MODULE PROCEDURE dbg_print_3d
END INTERFACE


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialization of indices for debug output
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-11)
  !!
  !
  ! TODO: parallelize
  !
  SUBROUTINE init_dbg_index (ppatch)

    TYPE(t_patch),             TARGET, INTENT(IN)     :: ppatch

    INTEGER  :: i
    REAL(wp) :: zlon, zlat, zarea, zlength

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      &      routine = 'mo_util_dbg_prnt:init_dbg_index'

    CALL message(TRIM(routine), 'Start' )

    ! test submit via eclipse workshop / 2013-02-11

    ! fill subset for use in dbg_print without passing the patch in every call
!    v_subdom_cell = ppatch%cells%in_domain

    ! local module variables for check of cells/edges/verts
    loc_nblks_c =ppatch%nblks_c
    loc_nblks_e =ppatch%nblks_e
    loc_nblks_v =ppatch%nblks_v

 !  For a correct dicision of cells/edges/verts the number of points on domain would be better
 !    but is not available as local dimension
 !  In order to keep a difference in number of blocks please use a low number for nproma
 !    i.e nproma=1 for debugging
 !  loc_patch_c =ppatch%n_patch_cells
 !  loc_patch_e =ppatch%n_patch_edges
 !  loc_patch_v =ppatch%n_patch_verts

    ! module index/block for one cell output
    IF ((idbg_idx /= 0 ) .OR. (idbg_blk /= 0 )) THEN
      c_i = idbg_idx
      c_b = idbg_blk
    ELSE
      ! search for block/index of debug output cell at lat/lon
      ! given by namelist dbg_index_nml - not yet parallelized
      CALL find_latlonindex (ppatch, dbg_lat_in, dbg_lon_in, c_i, c_b, near_proc_id)
    END IF

    zlat = ppatch%cells%center(c_i,c_b)%lat * 180.0_wp / pi
    zlon = ppatch%cells%center(c_i,c_b)%lon * 180.0_wp / pi

    !------------------------------------------------------------------
    ! print test cell
    !------------------------------------------------------------------

    ! output format
    99 FORMAT(     2(a,i4),2(a,f9.2),a,f13.2)
    97 FORMAT(a,i1,2(a,i4),2(a,f9.2),a,f13.2)

    zarea = ppatch%cells%area(c_i,c_b)*1.0e-6_wp ! in km2
    CALL message (TRIM(routine), 'Conditions at test cell (C), and edges/verts/neighbors:')
    WRITE(message_text,99) ' Cell C: block=',c_b,'  index=',c_i,               &
      &                    '  lat=',zlat,'  lon=',zlon,                        &
      &                    '  cell-area  =', zarea
    CALL message (' ', message_text)

    !------------------------------------------------------------------
    ! find and print corresponding edges/verts/neighbors of test cell
    !------------------------------------------------------------------

    DO i = 1, 3 ! 3 edges of cell C at (ne_i,ne_b)
      ne_b(i) = ppatch%cells%edge_blk(c_i,c_b,i)
      ne_i(i) = ppatch%cells%edge_idx(c_i,c_b,i)
      zlat    = ppatch%edges%center(ne_i(i),ne_b(i))%lat * 180.0_wp / pi
      zlon    = ppatch%edges%center(ne_i(i),ne_b(i))%lon * 180.0_wp / pi
      zlength = ppatch%edges%primal_edge_length(ne_i(i),ne_b(i))*0.001_wp  ! in km
      ! output
      WRITE(message_text,97) ' Edge E',i,' block=',ne_b(i),'  index=',ne_i(i), &
        &                    '  lat=',zlat,'  lon=',zlon,                      &
        &                    '  edge-length=',zlength
      CALL message (' ', message_text)
    END DO

    DO i = 1, 3 ! 3 vertices of cell C at (nv_i,nv_b)
      nv_b(i) = ppatch%cells%vertex_blk(c_i,c_b,i)
      nv_i(i) = ppatch%cells%vertex_idx(c_i,c_b,i)
      zlat    = ppatch%edges%center(nv_i(i),nv_b(i))%lat * 180.0_wp / pi
      zlon    = ppatch%edges%center(nv_i(i),nv_b(i))%lon * 180.0_wp / pi
      ! output
      WRITE(message_text,97) ' Vert V',i,' block=',nv_b(i),'  index=',nv_i(i), &
        &                    '  lat=',zlat,'  lon=',zlon
      CALL message (' ', message_text)
    END DO

    DO i = 1, 3 ! 3 neighbours of cell C at (nc_i,nc_b)
      nc_b(i)=ppatch%cells%neighbor_blk(c_i,c_b,i)
      nc_i(i)=ppatch%cells%neighbor_idx(c_i,c_b,i)
      IF ( nc_i(i) == 0 .OR. nc_b(i) == 0) THEN
        nc_i(i) = c_i
        nc_b(i) = c_b
        WRITE(message_text,'(a)') ' Neighbor Cell is NOT DEFINED'
      ELSE
        zlat = ppatch%cells%center(nc_i(i),nc_b(i))%lat * 180.0_wp / pi
        zlon = ppatch%cells%center(nc_i(i),nc_b(i))%lon * 180.0_wp / pi
        zarea= ppatch%cells%area(nc_i(i),nc_b(i))*1.0e-6_wp  ! in km2
        WRITE(message_text,97) ' Neighbor  C',i,' =',nc_b(i),'  index=',nc_i(i),  &
          &                    '  lat=',zlat,'  lon=',zlon,                       &
          &                    '  cell-area  =', zarea
      END IF
      ! output
      CALL message (' ', message_text)
    END DO

  END SUBROUTINE init_dbg_index

  !-------------------------------------------------------------------------
  !>
  !! Search for a cell center at given longitude and latitude
  !! provided in namelist dbg_index_nml
  !!
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-12)
  !!
  !! TODO: parallelize
  !
  !
  SUBROUTINE find_latlonindex (ppatch, plat_in, plon_in, iidx, iblk, proc_id)

  TYPE(t_patch), TARGET, INTENT(IN)     :: ppatch
  REAL(wp),              INTENT(IN)     :: plat_in       ! cell latitude to search for
  REAL(wp),              INTENT(IN)     :: plon_in       ! cell longitude to search for
  INTEGER,               INTENT(OUT)    :: iidx          ! index of nearest cell
  INTEGER,               INTENT(OUT)    :: iblk          ! block of nearest cell
  INTEGER,               INTENT(OUT)    :: proc_id       ! process where nearest cell is found

  INTEGER  :: jb, jc, i_startidx, i_endidx, mproc_id   !, mpi_comm
  REAL(wp) :: zlon, zlat, zdist, zdist_cmp, xctr
  REAL(wp) :: zdst_c(nproma,ppatch%nblks_c)
  TYPE(t_subset_range), POINTER :: owned_cells !,cells_in_domain, all_cells

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_util_dbg_prnt:find_latlonindex'

  CALL message(TRIM(routine), 'Start' )

  !all_cells       => ppatch%cells%all
  !cells_in_domain => ppatch%cells%in_domain
  owned_cells      => ppatch%cells%owned

  ! initial distance to compare
  zdist_cmp   = 100000.0_wp
  zdst_c(:,:) = 100000.0_wp

  !  loop over owned cells
  DO jb = owned_cells%start_block, owned_cells%end_block
    CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)

    DO jc = i_startidx, i_endidx

      zlat    = ppatch%cells%center(jc,jb)%lat * 180.0_wp / pi
      zlon    = ppatch%cells%center(jc,jb)%lon * 180.0_wp / pi

      zdist       = sqrt((zlat-plat_in)*(zlat-plat_in) + (zlon-plon_in)*(zlon-plon_in))
      zdst_c(jc,jb) = zdist
      IF (zdist < zdist_cmp) THEN
        iblk = jb
        iidx = jc
        zdist_cmp = zdist
      END IF

    END DO
  END DO

  CALL sync_patch_array(SYNC_C, ppatch, zdst_c(:,:))

  ! find PE with minimum distance
  ! disable p_test_run since global_max will be different

  p_test_run_bac = p_test_run
  p_test_run = .false.

  ! comparing zdist_cmp over all domains
  ! local variable mproc_id must be pre-set to p_pe in all domains
  mproc_id = p_pe
  !write(20+p_pe,*) ' 0max: mproc_id=',mproc_id,'  p_pe=',p_pe,'  zdst=',zdst_c(iidx,iblk),'  idx/blk=',iidx,iblk,' cmp=',zdist_cmp
  xctr = global_max(-zdist_cmp, mproc_id)
  !write(20+p_pe,*) ' 1max: mproc_id=',mproc_id,'  p_pe=',p_pe,'  zdst=',zdst_c(iidx,iblk),'  idx/blk=',iidx,iblk,' cmp=',zdist_cmp

  ! broadcast indices and proc_id to global, mainly for output - not the best, better write info directly
  !mpi_comm = p_comm_work
  !CALL p_bcast(iidx     , mproc_id, mpi_comm)
  !CALL p_bcast(iblk     , mproc_id, mpi_comm)
  !CALL p_bcast(zdist_cmp, mproc_id, mpi_comm)
  !CALL p_bcast(zlat_min , mproc_id, mpi_comm)  !  must be set accordingly above
  !CALL p_bcast(zlon_min , mproc_id, mpi_comm)  !  must be set accordingly above

  proc_id = mproc_id
  p_test_run = p_test_run_bac

  99 FORMAT(3a,i4,a,i4,3(a,f9.3))
  98 FORMAT(2a,3(a,f9.3))

  ! write info directly by mproc_id: needs barrier to avoid merging messages from pe_io and mproc_id
  ! write(0,*) get_my_mpi_all_id(), "enter icon_comm_barrier:"
  CALL icon_comm_barrier(for_patch=ppatch)
  CALL icon_comm_barrier(for_patch=ppatch) ! twice since the pipelining of output is still sometimes out of order
  IF (p_pe .EQ. mproc_id) THEN
  !IF (my_process_is_stdio()) THEN
     zlat = ppatch%cells%center(iidx,iblk)%lat * 180.0_wp / pi
     zlon = ppatch%cells%center(iidx,iblk)%lon * 180.0_wp / pi
     WRITE(0,98) ' ',TRIM(routine),' Found  cell nearest to         latitude=', plat_in,'  longitude=',plon_in
     WRITE(0,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  latitude=',zlat,'  longitude=',zlon
     WRITE(0,'(3a,i3,a,f9.3)') ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',mproc_id, &
      &                  '; distance in degrees =', zdist_cmp
  !  WRITE(125,98) ' ',TRIM(routine),' Found  cell nearest to         latitude=', plat_in,'  longitude=',plon_in
  !  WRITE(125,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  latitude=',zlat,'  longitude=',zlon
  !  WRITE(125,'(3a,i3,a,f9.3)') ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',mproc_id, &
  !    &                  '; distance in degrees =', zdist_cmp
  ! WRITE(125,'(3a,i3,a,i3,a,f9.3)') ' ',TRIM(routine),' FOUND: at block=',iblk,' index=',iidx, &
  !   &                  ' distance in degrees =', zdist_cmp
  ! WRITE(125,'(3a,2i3,a,f9.3)') ' ',TRIM(routine),' FOUND: using MINLOC: at idx/blk=', &
  !   &                  MINLOC(zdst_c(:,:)),' distance in degrees =',MINVAL(zdst_c(:,:))
  END IF
  CALL icon_comm_barrier(for_patch=ppatch)
  CALL icon_comm_barrier(for_patch=ppatch) ! twice since the pipelining of output is still sometimes out of order
  ! write(0,*) get_my_mpi_all_id(), "leave icon_comm_barrier:"

  END SUBROUTINE find_latlonindex

  !-------------------------------------------------------------------------
  !>
  !! Print out min and max or a specific cell value and neighbors of a 3-dim array.
  !!
  !! Reduce writing effort for a simple print.
  !! The amount of prints is controlled by comparison of a fixed level of detail
  !! for output (idetail_src) with variables idbg_mxmn/idbg_val  that are
  !! given via namelist dbg_index_nml
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2012-06)
  !!
  !
  SUBROUTINE dbg_print_3d( str_prntdes, p_array, str_mod_src, idetail_src, in_subset )

  CHARACTER(len=*),      INTENT(IN) :: str_prntdes    ! description of array
  REAL(wp),              INTENT(IN) :: p_array(:,:,:) ! 3-dim array for debugging
  CHARACTER(len=*),      INTENT(IN) :: str_mod_src    ! defined string for source of current array
  INTEGER,               INTENT(IN) :: idetail_src    ! source level from module for print output
  TYPE(t_subset_range),  TARGET, OPTIONAL :: in_subset

  ! local variables
  CHARACTER(len=27) ::  strout
  CHARACTER(len=12) ::  strmod
  INTEGER           ::  slev, elev, elev_val, elev_mxmn
  INTEGER           ::  iout, icheck_str_mod, jstr, i, jk, nlev, ndimblk
  REAL(wp)          :: minmaxmean(3)

  IF (timers_level > 10) CALL timer_start(timer_dbg_prnt)


#ifdef __SX__
  ! valid g-format without offset of decimal point
  981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g12.5))
  982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g12.5))
  992 FORMAT(a,a12,':',a27,'  :',i3, 3g26.18)
#else

! ! valid g-format without offset of decimal point
! 981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3, 2g26.18)
! ! valid e-format with first digit > zero
! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pe26.18,3(a,i0,a,1pe20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pe20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2e26.18)
! ! g-format with first digit > zero, not valid for SX-compiler
  981 FORMAT(a,a12,':',a27,' C:',i3, 1pg26.18,3(a,i0,a,1pg20.12))
  982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pg20.12))
  992 FORMAT(a,a12,':',a27,'  :',i3, 1pg26.18, 1pg26.18, 1pg26.18)
#endif

  ! check print output level idetail_src (1-5) with namelist given value (idbg_val)
  ! for output at given index

  !
  ! All calculations are done inside this IF only
  !

  IF (idbg_val >= idetail_src) THEN
    !
    ! dimensions
    !                           !  index 1:nproma
    nlev    = SIZE(p_array,2)   !  vertical dimension (levels)
    ndimblk = SIZE(p_array,3)   !  blocks 1:nblks for cells/edges/verts

    ! output channel: stderr
    iout = nerr

    ! compare defined source string with namelist-given output string
    icheck_str_mod = 0
    DO jstr = 1, dim_mod_tst
      IF (str_mod_src == str_mod_tst(jstr) .OR. str_mod_tst(jstr) == 'all') &
        &  icheck_str_mod = 1
    END DO

    ! if str_mod_src not found in str_mod_tst - no output
    IF (icheck_str_mod == 0 .and. timers_level > 10) CALL timer_stop(timer_dbg_prnt)
    IF (icheck_str_mod == 0 ) RETURN

    strout=TRIM(str_prntdes)
    strmod=TRIM(str_mod_src)

    ! check start and end index for output of vertical levels via namelist
    slev = 1
    IF (idbg_slev > 1)    slev = idbg_slev
    IF (slev      > nlev) slev = nlev
    elev = nlev
    IF (idbg_elev < nlev) elev = idbg_elev

    ! idbg_val<4: one level output only (slev), apart from permanent output (init)
    elev_val = elev
    IF (idbg_val < 4 .AND. idetail_src > 0) elev_val = slev

    DO jk = slev, elev_val

      ! write value at index
      IF (ndimblk == loc_nblks_c) THEN
        IF (my_process_is_stdio()) &
          & WRITE(iout,981) '   VALUE ', strmod, strout, jk, p_array(c_i,jk,c_b), &
          &                 (' C',i,':',p_array(nc_i(i),jk,nc_b(i)),i=1,3)
      ELSE IF (ndimblk == loc_nblks_e) THEN
        IF (my_process_is_stdio()) &
          & WRITE(iout,982) '   VALUE ', strmod, strout, jk, &
          &                 (' E',i,':',p_array(ne_i(i),jk,ne_b(i)),i=1,3)
      ELSE IF (ndimblk == loc_nblks_v) THEN
        IF (my_process_is_stdio()) &
          & WRITE(iout,982) '   VALUE ', strmod, strout, jk, &
          &                 (' V',i,':',p_array(nv_i(i),jk,nv_b(i)),i=1,3)
      END IF

    END DO

  END IF

  ! check print output level idetail_src (1-5) with namelist given value (idbg_mxmn)
  ! for MIN/MAX output:

  !
  ! All calculations are done inside this IF only
  !

  IF (idbg_mxmn >= idetail_src ) THEN
    !
    ! dimensions
    !                           !  index 1:nproma
    nlev    = SIZE(p_array,2)   !  vertical dimension (levels)
    ndimblk = SIZE(p_array,3)   !  blocks 1:nblks for cells/edges/verts

    ! output channel: stderr
    iout = nerr

    ! compare defined source string with namelist-given output string
    icheck_str_mod = 0
    DO jstr = 1, dim_mod_tst
      IF (str_mod_src == str_mod_tst(jstr) .OR. str_mod_tst(jstr) == 'all') &
        &  icheck_str_mod = 1
    END DO

    ! if str_mod_src not found in str_mod_tst - no output
    IF (icheck_str_mod == 0 .and. timers_level > 10) CALL timer_stop(timer_dbg_prnt)
    IF (icheck_str_mod == 0 ) RETURN

    strout=TRIM(str_prntdes)
    strmod=TRIM(str_mod_src)

    ! check start and end index for output of vertical levels via namelist
    slev = 1
    IF (idbg_slev > 1)    slev = idbg_slev
    IF (slev      > nlev) slev = nlev
    elev = nlev
    IF (idbg_elev < nlev) elev = idbg_elev

    ! idbg_mxmn<4: one level output only (slev), independent of elev_val
    elev_mxmn = elev
    IF (idbg_mxmn < 4 .AND. idetail_src > 0) elev_mxmn = slev

    ! print out maximum and minimum value
    ! ctrn=minval(p_array(:, slev:elev_mxmn, :))
    DO jk = slev, elev_mxmn

      minmaxmean(:) = global_minmaxmean(values=p_array(:,:,:), range_subset=in_subset, start_level=jk, end_level=jk)

      IF (my_process_is_stdio()) &
        & WRITE(iout,992) ' MAX/MIN/MEAN ', strmod, strout, jk, minmaxmean(2), minmaxmean(1), minmaxmean(3)


      ! location of max/min - parallelize!
      ! WRITE(iout,983) ' LOC ',strout,jk, &
      !   &              MAXLOC(p_array(1:nproma,jk,1:ndimblk)),     &
      !   &              MINLOC(p_array(1:nproma,jk,1:ndimblk))
! 983 FORMAT(a,a12,':',a27,'  :',i3, 4i4)

    END DO

  END IF

  IF (timers_level > 10) CALL timer_stop(timer_dbg_prnt)

  END SUBROUTINE dbg_print_3d

  !-------------------------------------------------------------------------
  !>
  !! Print out min and max or a specific cell value and neighbors of a 2-dim array.
  !!

  SUBROUTINE dbg_print_2d( str_prntdes, p_array, str_mod_src, idetail_src, in_subset )

  CHARACTER(len=*),      INTENT(IN) :: str_prntdes    ! description of array
  REAL(wp),              INTENT(IN) :: p_array(:,:)   ! 2-dim array for debugging
  CHARACTER(len=*),      INTENT(IN) :: str_mod_src    ! defined string for source of current array
  INTEGER,               INTENT(IN) :: idetail_src    ! source level from module for print output
  TYPE(t_subset_range),  TARGET, OPTIONAL :: in_subset

  ! local variables
  CHARACTER(len=27) ::  strout
  CHARACTER(len=12) ::  strmod
  INTEGER           ::  iout, icheck_str_mod, jstr, i, jk, ndimblk
  REAL(wp)          ::  minmaxmean(3)

  IF (timers_level > 10) CALL timer_start(timer_dbg_prnt)

  ! dimensions - first dimension is nproma
  ndimblk = SIZE(p_array,2)

  ! output channel: stderr
  iout = nerr

  ! compare defined source string with namelist-given output string
  icheck_str_mod = 0
  DO jstr = 1, dim_mod_tst
    IF (str_mod_src == str_mod_tst(jstr) .OR. str_mod_tst(jstr) == 'all') &
      &  icheck_str_mod = 1
  END DO

  ! if str_mod_src not found in str_mod_tst - no output
  IF (icheck_str_mod == 0 .and. timers_level > 10) CALL timer_stop(timer_dbg_prnt)
  IF (icheck_str_mod == 0 ) RETURN

#ifdef __SX__
  ! valid g-format without offset of decimal point
  981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g12.5))
  982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g12.5))
  992 FORMAT(a,a12,':',a27,'  :',i3, 3g26.18)
#else

! ! valid e-format with first digit > zero
! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pe26.18,3(a,i0,a,1pe20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pe20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2e26.18)

! ! g-format with first digit > zero, not valid for SX-compiler
  981 FORMAT(a,a12,':',a27,' C:',i3, 1pg26.18,3(a,i0,a,1pg20.12))
  982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pg20.12))
  992 FORMAT(a,a12,':',a27,'  :',i3, 1pg26.18, 1pg26.18, 1pg26.18)
#endif

  strout=TRIM(str_prntdes)
  strmod=TRIM(str_mod_src)

  ! surface level output only
  jk = 1


  ! check print output level idetail_src (1-5) with namelist given value (idbg_val)
  ! for output at given index

  IF (idbg_val >= idetail_src) THEN

    !write(iout,*) ' ndimblk and loc_nblks = ',ndimblk, loc_nblks_c, loc_nblks_e, loc_nblks_v

    ! write value at index
    IF (ndimblk == loc_nblks_c) THEN
      IF (my_process_is_stdio()) &
        & WRITE(iout,981) '   VALUE ', strmod, strout, jk, p_array(c_i,c_b), &
        &                 (' C',i,':',p_array(nc_i(i),nc_b(i)),i=1,3)
    ELSE IF (ndimblk == loc_nblks_e) THEN
      IF (my_process_is_stdio()) &
        & WRITE(iout,982) '   VALUE ', strmod, strout, jk, &
        &                 (' E',i,':',p_array(ne_i(i),ne_b(i)),i=1,3)
    ELSE IF (ndimblk == loc_nblks_v) THEN
      IF (my_process_is_stdio()) &
        & WRITE(iout,982) '   VALUE ', strmod, strout, jk, &
        &                 (' V',i,':',p_array(nv_i(i),nv_b(i)),i=1,3)
    END IF

  END IF

  ! check print output level idetail_src (1-5) with namelist given value (idbg_mxmn)
  ! for MIN/MAX output:

  IF (idbg_mxmn >= idetail_src ) THEN

    minmaxmean(:) = global_minmaxmean(values=p_array(:,:), range_subset=in_subset)

    IF (my_process_is_stdio()) &
      & WRITE(iout,992) ' MAX/MIN/MEAN ', strmod, strout, jk, minmaxmean(2), minmaxmean(1), minmaxmean(3)

    ! WRITE(0,*) ' MAX/MIN/MEAN ', minmaxmean(2), minmaxmean(1), minmaxmean(3)

  END IF

  IF (timers_level > 10) CALL timer_stop(timer_dbg_prnt)

  END SUBROUTINE dbg_print_2d

END MODULE mo_util_dbg_prnt

