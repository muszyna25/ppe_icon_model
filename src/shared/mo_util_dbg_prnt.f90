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
USE mo_mpi,                    ONLY: my_process_is_stdio, p_pe
USE mo_io_units,               ONLY: nerr
USE mo_parallel_config,        ONLY: nproma, p_test_run
USE mo_impl_constants,         ONLY: max_char_length
USE mo_run_config,             ONLY: ltimer
USE mo_timer,                  ONLY: timer_start, timer_stop, timer_dbg_prnt
USE mo_sync,                   ONLY: SYNC_C, sync_patch_array, global_max, global_min
USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
USE mo_dbg_nml,                ONLY: str_mod_tst, dim_mod_tst, dbg_lon_in, dbg_lat_in, &
  &                                  idbg_mxmn, idbg_val, idbg_slev, idbg_elev, &
  &                                  idbg_idx, idbg_blk
USE mo_math_constants,         ONLY: pi
USE mo_exception,              ONLY: message, message_text
USE mo_model_domain,           ONLY: t_patch


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

! indices of cells and neighbours for debug output at single cell
INTEGER :: c_b, c_i, ne_b(3), ne_i(3), nc_b(3), nc_i(3), nv_b(3), nv_i(3)
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
      CALL find_latlonindex (ppatch, dbg_lat_in, dbg_lon_in, c_i, c_b)
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
      &                    '  cell-area = ', zarea
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
          &                    '  cell-area = ', zarea
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
  SUBROUTINE find_latlonindex (ppatch, plat_in, plon_in, iidx, iblk)

  TYPE(t_patch), TARGET, INTENT(IN)     :: ppatch
  REAL(wp),              INTENT(IN)     :: plat_in       ! cell latitude to search for
  REAL(wp),              INTENT(IN)     :: plon_in       ! cell longitude to search for
  INTEGER,               INTENT(OUT)    :: iidx          ! index of nearest cell
  INTEGER,               INTENT(OUT)    :: iblk          ! block of nearest cell

  INTEGER  :: jb, jc, i_startidx, i_endidx, proc_id
  REAL(wp) :: zlon, zlat, zdist, zdist_cmp, ctr
  REAL(wp) :: zdst_c(nproma,ppatch%nblks_c)
  TYPE(t_subset_range), POINTER :: cells_in_domain!, all_cells

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_util_dbg_prnt:find_latlonindex'

  CALL message(TRIM(routine), 'Start' )

  !all_cells       => ppatch%cells%all
  cells_in_domain => ppatch%cells%in_domain

  ! initial distance to compare
  zdist_cmp = 100000.0_wp
  zdst_c(:,:) = 100000.0_wp
  proc_id   = -1

  !  loop over all cells including halo
  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
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

  ! find PE with minimum distance, MPI-broadcast block/index from that PE - not yet
  ! disable p_test_run since global_max will be different
  p_test_run_bac = p_test_run
  p_test_run = .false.
  ctr = global_max(-zdist,proc_id)
  p_test_run = p_test_run_bac

  zlat    = ppatch%cells%center          (iidx,iblk)%lat * 180.0_wp / pi
  zlon    = ppatch%cells%center          (iidx,iblk)%lon * 180.0_wp / pi

  99 FORMAT(3a,i4,a,i4,3(a,f9.2))
  98 FORMAT(2a,3(a,f9.2))
  IF (my_process_is_stdio()) THEN
    WRITE(0,98) ' ',TRIM(routine),' Found  cell nearest to         lat=', plat_in,'  lon=',plon_in
    WRITE(0,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  lat=',zlat,'  lon=',zlon
    WRITE(0,'(3a,i3)')         ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',proc_id
    WRITE(0,'(3a,2i3,a,f9.2)') ' ',TRIM(routine),' FOUND: Min dist is at idx/blk=', &
      &                  MINLOC(zdst_c(:,:)),' distance in deg =',MINVAL(zdst_c(:,:))
  END IF

  IF (p_pe == proc_id) THEN

 !  !  loop over all cells with minimum distance
 !  DO jb = cells_in_domain%start_block, cells_in_domain%end_block
 !    CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
 !    DO jc = i_startidx, i_endidx
 !  
 !      zlat    = ppatch%cells%center(jc,jb)%lat * 180.0_wp / pi
 !      zlon    = ppatch%cells%center(jc,jb)%lon * 180.0_wp / pi
 !  
 !      zdist       = sqrt((zlat-plat_in)*(zlat-plat_in) + (zlon-plon_in)*(zlon-plon_in))
 !      zdst_c(jc,jb) = zdist
 !      IF (zdist < zdist_cmp) THEN
 !        iblk = jb
 !        iidx = jc
 !        zdist_cmp = zdist
 !      END IF
 !  
 !    END DO
 !  END DO

 !  WRITE(0,'(3a)') ' ',TRIM(routine),' Par: Found  cell now in Process proc_id'
 !  WRITE(0,99) ' ',TRIM(routine),    ' Par: Found  block=',iblk,'  index=',iidx, &
 !    &                               '  lat=',zlat,'  lon=',zlon
 !  WRITE(0,'(3a,i3)')         ' ',TRIM(routine),' Par: FOUND: proc_id for nearest cell is=',proc_id
 !  WRITE(0,'(3a,2i3,a,f9.2)') ' ',TRIM(routine),' Par: FOUND: Min dist is at idx/blk=', &
 !    &                  MINLOC(zdst_c(:,:)),' distance in deg =',MINVAL(zdst_c(:,:))

  END IF

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
  SUBROUTINE dbg_print_3d( str_prntdes, p_array, str_mod_src, idetail_src )

  CHARACTER(len=*),      INTENT(IN) :: str_prntdes    ! description of array
  REAL(wp),              INTENT(IN) :: p_array(:,:,:) ! 3-dim array for debugging
  CHARACTER(len=*),      INTENT(IN) :: str_mod_src    ! defined string for source of current array
  INTEGER,               INTENT(IN) :: idetail_src    ! source level from module for print output 

  ! local variables
  CHARACTER(len=27) ::  strout
  CHARACTER(len=12) ::  strmod
  INTEGER           ::  slev, elev, elev_val, elev_mxmn
  INTEGER           ::  iout, icheck_str_mod, jstr, i, jk, nlev, ndimblk
  REAL(wp)          ::  ctrx, ctrn, glbmx, glbmn
  !TYPE(t_subset_range), POINTER :: cells_in_domain!, all_cells
  INTEGER           ::  i_startidx, i_endidx, jb
  REAL(wp)          ::  ctrxind(nproma), ctrnind(nproma)

  IF (ltimer) CALL timer_start(timer_dbg_prnt)

  !cells_in_domain => ppatch%cells%in_domain

#ifdef __SX__
  ! valid g-format without offset of decimal point
  981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g12.5))
  982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g12.5))
  991 FORMAT(a,a12,':',a27,'  :',i3, 2g26.18)
#else

! ! valid e-format with first digit gt zero
! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pe26.18,3(a,i0,a,1pe20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pe20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2e26.18)

! ! g-format with offset for decimal point not valid for NAG compiler
! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pg26.18,3(a,i0,a,1pg20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pg20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2g26.18)

  ! valid g-format without offset of decimal point
  981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g20.12))
  982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g20.12))
  991 FORMAT(a,a12,':',a27,'  :',i3, 2g26.18)
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
    IF (icheck_str_mod == 0 .and. ltimer) CALL timer_stop(timer_dbg_prnt)
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
    IF (icheck_str_mod == 0 .and. ltimer) CALL timer_stop(timer_dbg_prnt)
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
    DO jk = slev, elev_mxmn
    
      ctrx=maxval(p_array(1:nproma,jk,1:ndimblk))
      ctrn=minval(p_array(1:nproma,jk,1:ndimblk))

      ! find max/min out of active indices only
       !DO jb = cells_in_domain%start_block, cells_in_domain%end_block
       !DO jb = v_subdom_cell%start_block, v_subdom_cell%end_block
       !  CALL get_index_range(v_subdom_cell, jb, i_startidx, i_endidx)
       !  ctrxind(jb)=maxval(p_array(i_startidx:i_endidx,jk,jb))
       !  ctrnind(jb)=minval(p_array(i_startidx:i_endidx,jk,jb))
       !END DO  
       !ctrx=maxval(ctrxind(:))
       !ctrn=minval(ctrxind(:))

      ! parallelize:
      p_test_run_bac = p_test_run
      p_test_run = .false.
      glbmx=global_max(ctrx)
      glbmn=global_min(ctrn)
      p_test_run = p_test_run_bac
    
      IF (my_process_is_stdio()) &
        & WRITE(iout,991) ' MAX/MIN ', strmod, strout, jk, glbmx, glbmn

    
      ! location of max/min - parallelize!
      ! WRITE(iout,983) ' LOC ',strout,jk, &
      !   &              MAXLOC(p_array(1:nproma,jk,1:ndimblk)),     &
      !   &              MINLOC(p_array(1:nproma,jk,1:ndimblk))
! 983 FORMAT(a,a12,':',a27,'  :',i3, 4i4)
    
    END DO

  END IF

  IF (ltimer) CALL timer_stop(timer_dbg_prnt)

  END SUBROUTINE dbg_print_3d

  !-------------------------------------------------------------------------
  !>
  !! Print out min and max or a specific cell value and neighbors of a 2-dim array.
  !!

  SUBROUTINE dbg_print_2d( str_prntdes, p_array, str_mod_src, idetail_src )

  CHARACTER(len=*),      INTENT(IN) :: str_prntdes    ! description of array
  REAL(wp),              INTENT(IN) :: p_array(:,:)   ! 2-dim array for debugging
  CHARACTER(len=*),      INTENT(IN) :: str_mod_src    ! defined string for source of current array
  INTEGER,               INTENT(IN) :: idetail_src    ! source level from module for print output 

  ! local variables
  CHARACTER(len=27) ::  strout
  CHARACTER(len=12) ::  strmod
  INTEGER           ::  iout, icheck_str_mod, jstr, i, jk, ndimblk
  REAL(wp)          ::  ctrx, ctrn, glbmx, glbmn
  !TYPE(t_subset_range), POINTER :: cells_in_domain!, all_cells
  !REAL(wp)          ::  ctrxind(nproma), ctrnind(nproma)
  !INTEGER           ::  i_startidx, i_endidx, jb

  IF (ltimer) CALL timer_start(timer_dbg_prnt)

  !cells_in_domain => ppatch%cells%in_domain

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
  IF (icheck_str_mod == 0 .and. ltimer) CALL timer_stop(timer_dbg_prnt)
  IF (icheck_str_mod == 0 ) RETURN

#ifdef __SX__
  ! valid g-format without offset of decimal point
  981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g12.5))
  982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g12.5))
  991 FORMAT(a,a12,':',a27,'  :',i3, 2g26.18)
#else

! ! valid e-format with first digit gt zero
! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pe26.18,3(a,i0,a,1pe20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pe20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2e26.18)

! ! g-format with offset for decimal point not valid for NAG compiler
! 981 FORMAT(a,a12,':',a27,' C:',i3, 1pg26.18,3(a,i0,a,1pg20.12))
! 982 FORMAT(a,a12,':',a27,'  :',i3,    26x,  3(a,i0,a,1pg20.12))
! 991 FORMAT(a,a12,':',a27,'  :',i3,1p2g26.18)

  ! valid g-format without offset of decimal point
  981 FORMAT(a,a12,':',a27,' C:',i3,  g26.18,3(a,i0,a,  g20.12))
  982 FORMAT(a,a12,':',a27,'  :',i3,   26x,  3(a,i0,a,  g20.12))
  991 FORMAT(a,a12,':',a27,'  :',i3, 2g26.18)
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

    ctrx=maxval(p_array(1:nproma,1:ndimblk))
    ctrn=minval(p_array(1:nproma,1:ndimblk))

    ! find max/min out of active indices only
 !   DO jb = cells_in_domain%start_block, cells_in_domain%end_block
 !     CALL get_index_range(cells_in_domain, jb, i_startidx, i_endidx)
 !     ctrxind(jb)=maxval(p_array(1:i_startidx,1:i_endidx))
 !     ctrnind(jb)=minval(p_array(1:i_startidx,1:i_endidx))
 !   END DO  
 !   ctrx=maxval(ctrxind(:))
 !   ctrn=minval(ctrxind(:))

    ! print out maximum and minimum value
    ! parallelize:
    p_test_run_bac = p_test_run
    p_test_run = .false.
    glbmx=global_max(ctrx)
    glbmn=global_min(ctrn)
    p_test_run = p_test_run_bac
   
    IF (my_process_is_stdio()) &
      & WRITE(iout,991) ' MAX/MIN ', strmod, strout, jk, glbmx, glbmn

  END IF

  IF (ltimer) CALL timer_stop(timer_dbg_prnt)

  END SUBROUTINE dbg_print_2d

END MODULE mo_util_dbg_prnt

