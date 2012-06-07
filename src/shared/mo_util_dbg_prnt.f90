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
USE mo_mpi,                    ONLY: my_process_is_stdio
USE mo_io_units,               ONLY: nerr
USE mo_parallel_config,        ONLY: nproma, p_test_run
USE mo_impl_constants,         ONLY: max_char_length
!USE mo_run_config,             ONLY: ltimer
!USE mo_timer,                  ONLY: timer_start, timer_stop, timer_print_mxmn
USE mo_sync,                   ONLY: SYNC_C, sync_patch_array, global_max, global_min
USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
USE mo_dbg_nml,                ONLY: str_mod_tst, dim_mod_tst, dbg_lon_in, dbg_lat_in, &
  &                                  idbg_mxmn, idbg_val
USE mo_math_constants,         ONLY: pi
USE mo_exception,              ONLY: message, message_text
USE mo_model_domain,           ONLY: t_patch


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

! indices of cells and neighbours for debug output at single cell
INTEGER :: c_b, c_i, ne_b(3), ne_i(3), nc_b(3), nc_i(3), nv_b(3), nv_i(3)
INTEGER :: loc_nblks_c, loc_nblks_e, loc_nblks_v

PUBLIC :: init_dbg_index
PUBLIC :: dbg_print

INTERFACE dbg_print
  MODULE PROCEDURE dbg_print_2d
  MODULE PROCEDURE dbg_print_3d
END INTERFACE


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Initialization of indices
  !!
  !! @par Revision History
  !! Initial release by Stephan Lorenz, MPI-M (2010-11)
  !!
  !
  ! TODO: parallelize
  !
  SUBROUTINE init_dbg_index (ppatch)


  TYPE(t_patch),             TARGET, INTENT(IN)     :: ppatch

  INTEGER :: i
  REAL(wp) :: zlon, zlat

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &      routine = 'mo_util_dbg_prnt:init_dbg_index'

  CALL message(TRIM(routine), 'Start' )

  ! module variables for check of cells/edges/verts
  loc_nblks_c =ppatch%nblks_c
  loc_nblks_e =ppatch%nblks_e
  loc_nblks_v =ppatch%nblks_v

  ! search for block/index of debug output cell at lat/lon
  ! given by namelist dbg_index_nml - not yet parallelized
  CALL find_latlonindex (ppatch, dbg_lat_in, dbg_lon_in, c_i, c_b)

  zlat = ppatch%cells%center(c_i,c_b)%lat * 180.0_wp / pi
  zlon = ppatch%cells%center(c_i,c_b)%lon * 180.0_wp / pi

  !------------------------------------------------------------------
  ! print test cell
  !------------------------------------------------------------------

  ! output format
  99 FORMAT(     2(a,i4),2(a,f9.2))
  97 FORMAT(a,i1,2(a,i4),2(a,f9.2))

  CALL message (TRIM(routine), 'Conditions at test cell (C), and edges/verts/neighbors:')
  WRITE(message_text,99) ' Cell C: block=',c_b,'  index=',c_i,               &
               &         '  lat=',zlat,'  lon=',zlon
  CALL message (' ', message_text)

  !------------------------------------------------------------------
  ! find and print corresponding edges/verts/neighbors of test cell
  !------------------------------------------------------------------

  DO i = 1, 3 ! 3 edges of cell C at (ne_i,ne_b)
    ne_b(i) = ppatch%cells%edge_blk(c_i,c_b,i)
    ne_i(i) = ppatch%cells%edge_idx(c_i,c_b,i)
    zlat    = ppatch%edges%center(ne_i(i),ne_b(i))%lat * 180.0_wp / pi
    zlon    = ppatch%edges%center(ne_i(i),ne_b(i))%lon * 180.0_wp / pi
    ! output
    WRITE(message_text,97) ' Edge E',i,' block=',ne_b(i),'  index=',ne_i(i), &
      &                    '  lat=',zlat,'  lon=',zlon
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
      WRITE(message_text,97) ' Neighbor  C',i,' =',nc_b(i),'  index=',nc_i(i),            &
        &                    '  lat=',zlat,'  lon=',zlon
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
  LOGICAL  :: p_test_run_bac

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
    WRITE(0,98) ' ',TRIM(routine),' Found cell nearest to          lat=', plat_in,'  lon=',plon_in
    WRITE(0,99) ' ',TRIM(routine),' Found  block=',iblk,'  index=',iidx,'  lat=',zlat,'  lon=',zlon
    WRITE(0,'(3a,i3)')         ' ',TRIM(routine),' FOUND: proc_id for nearest cell is=',proc_id
    WRITE(0,'(3a,2i3,a,f9.2)') ' ',TRIM(routine),' FOUND: Min dist is at idx/blk=', &
      &                  MINLOC(zdst_c(:,:)),' distance in deg =',MINVAL(zdst_c(:,:))
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
  INTEGER           ::  iout, icheck_str_mod, jstr, iper, i, jk, klev, nlev, ndimblk
  REAL(wp)          ::  ctr, glbmx, glbmn

  !IF (ltimer) CALL timer_start(timer_print_mxmn)

  ! dimensions
  !                           !  index 1:nproma
  nlev    = SIZE(p_array,2)   !  vertical dimension (levels)
  ndimblk = SIZE(p_array,3)   !  blocks 1:nblks for cells/edges/verts

  ! output channel: stderr
  iout = nerr

  ! compare defined source string with namelist-given output string ('per' for permanent output)
  icheck_str_mod = 0
  iper = 0
  DO jstr = 1, dim_mod_tst
    IF (str_mod_src == str_mod_tst(jstr) .OR. str_mod_src       == 'per'  &
      &                                  .OR. str_mod_tst(jstr) == 'all') &
      &  icheck_str_mod = 1
  END DO

  !IF (icheck_str_mod == 0 .and. ltimer) CALL timer_stop(timer_print_mxmn)

  ! if str_mod_src not found in str_mod_tst - no output
  IF (icheck_str_mod == 0 ) RETURN

! ! valid e-format with first digit gt zero
! 981 FORMAT(a,a10,a25,' C:',i3, 1pe26.18,3(a,i0,a,1pe16.8))
! 982 FORMAT(a,a10,a25,'  :',i3,    26x,  3(a,i0,a,1pe16.8))
! 991 FORMAT(a,a10,a27,  ':',i3,1p2e26.18)

! ! g-format with offset for decimal point not valid for NAG compiler
! 981 FORMAT(a,a10,a25,' C:',i3, 1pg26.18,3(a,i0,a,1pg16.8))
! 982 FORMAT(a,a10,a25,'  :',i3,    26x,  3(a,i0,a,1pg16.8))
! 991 FORMAT(a,a10,a27,  ':',i3,1p2g26.18)

  ! valid g-format without offset of decimal point
  981 FORMAT(a,a10,a25,' C:',i3,  g26.18,3(a,i0,a,  g16.8))
  982 FORMAT(a,a10,a25,'  :',i3,   26x,  3(a,i0,a,  g16.8))
  991 FORMAT(a,a10,a27,  ':',i3, 2g26.18)

  983 FORMAT(a,a10,a25,'  :',i3, 4i4)

  strout=TRIM(str_prntdes)


  ! check print output level idetail_src (1-5) with namelist given value (idbg_val)
  ! for output at given index

  IF (idbg_val >= idetail_src) THEN
    IF (my_process_is_stdio()) THEN

      ! idbg_val<4: surface level output only
      klev = nlev
      IF (idbg_val < 4) klev = 1

      DO jk = 1, klev

        ! write value at index
        IF (ndimblk == loc_nblks_c) THEN
          WRITE(iout,981) '   VALUE ', str_mod_src, strout, jk, p_array(c_i,jk,c_b), &
        &                 (' C',i,':',p_array(nc_i(i),jk,nc_b(i)),i=1,3)
        ELSE IF (ndimblk == loc_nblks_e) THEN
          WRITE(iout,982) '   VALUE ', str_mod_src, strout, jk, &
        &                 (' E',i,':',p_array(ne_i(i),jk,ne_b(i)),i=1,3)
        ELSE IF (ndimblk == loc_nblks_v) THEN
          WRITE(iout,982) '   VALUE ', str_mod_src, strout, jk, &
        &                 (' V',i,':',p_array(nv_i(i),jk,nv_b(i)),i=1,3)
        END IF

      END DO

    END IF
  END IF

  ! check print output level idetail_src (1-5) with namelist given value (idbg_mxmn)
  ! for MIN/MAX output:

  !IF (idbg_mxmn < idetail_src .and. ltimer) CALL timer_stop(timer_print_mxmn)
  IF (idbg_mxmn < idetail_src ) RETURN

  ! idbg_mxmn<4: surface level output only
  klev = nlev
  IF (idbg_mxmn < 4) klev = 1

  ! print out maximum and minimum value
  DO jk = 1, klev

    ! parallelize:
    ctr=maxval(p_array(1:nproma,jk,1:ndimblk))
    glbmx=global_max(ctr)
    ctr=minval(p_array(1:nproma,jk,1:ndimblk))
    glbmn=global_min(ctr)

    IF (my_process_is_stdio()) &
      & WRITE(iout,991) ' MAX/MIN ', str_mod_src, strout, jk, glbmx, glbmn

    ! location of max/min - parallelize!
!!$    WRITE(iout,983) ' LOC ',strout,klev, &
!!$      &              MAXLOC(p_array(1:nproma,jk,1:ndimblk)),     &
!!$      &              MINLOC(p_array(1:nproma,jk,1:ndimblk))

  END DO

  !IF (ltimer) CALL timer_stop(timer_print_mxmn)

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
  INTEGER           ::  iout, icheck_str_mod, jstr, iper, i, jk, klev, ndimblk
  REAL(wp)          ::  ctr, glbmx, glbmn

  !IF (ltimer) CALL timer_start(timer_print_mxmn)

  ! dimensions - first dimension is nproma
  ndimblk = SIZE(p_array,2)

  ! output channel: stderr
  iout = nerr

  ! compare defined source string with namelist-given output string ('per' for permanent output)
  icheck_str_mod = 0
  iper = 0
  DO jstr = 1, dim_mod_tst
    IF (str_mod_src == str_mod_tst(jstr) .OR. str_mod_src       == 'per'  &
      &                                  .OR. str_mod_tst(jstr) == 'all') &
      &  icheck_str_mod = 1
  END DO

  !IF (icheck_str_mod == 0 .and. ltimer) CALL timer_stop(timer_print_mxmn)

  ! if str_mod_src not found in str_mod_tst - no output
  IF (icheck_str_mod == 0 ) RETURN

! ! valid e-format with first digit gt zero
! 981 FORMAT(a,a10,a25,' C:',i3, 1pe26.18,3(a,i0,a,1pe16.8))
! 982 FORMAT(a,a10,a25,'  :',i3,    26x,  3(a,i0,a,1pe16.8))
! 991 FORMAT(a,a10,a27,  ':',i3,1p2e26.18)

! ! g-format with offset for decimal point not valid for NAG compiler
! 981 FORMAT(a,a10,a25,' C:',i3, 1pg26.18,3(a,i0,a,1pg16.8))
! 982 FORMAT(a,a10,a25,'  :',i3,    26x,  3(a,i0,a,1pg16.8))
! 991 FORMAT(a,a10,a27,  ':',i3,1p2g26.18)

  ! valid g-format without offset of decimal point
  981 FORMAT(a,a10,a25,' C:',i3,  g26.18,3(a,i0,a,  g16.8))
  982 FORMAT(a,a10,a25,'  :',i3,   26x,  3(a,i0,a,  g16.8))
  991 FORMAT(a,a10,a27,  ':',i3, 2g26.18)

  983 FORMAT(a,a10,a25,'  :',i3, 4i4)

  strout=TRIM(str_prntdes)


  ! check print output level idetail_src (1-5) with namelist given value (idbg_val)
  ! for output at given index

  IF (idbg_val >= idetail_src) THEN
    IF (my_process_is_stdio()) THEN

      ! surface level output only
      jk = 1

      ! write value at index
      IF (ndimblk == loc_nblks_c) THEN
        WRITE(iout,981) '   VALUE ', str_mod_src, strout, jk, p_array(c_i,c_b), &
      &                 (' C',i,':',p_array(nc_i(i),nc_b(i)),i=1,3)
      ELSE IF (ndimblk == loc_nblks_e) THEN
        WRITE(iout,982) '   VALUE ', str_mod_src, strout, jk, &
      &                 (' E',i,':',p_array(ne_i(i),ne_b(i)),i=1,3)
      ELSE IF (ndimblk == loc_nblks_v) THEN
        WRITE(iout,982) '   VALUE ', str_mod_src, strout, jk, &
      &                 (' V',i,':',p_array(nv_i(i),nv_b(i)),i=1,3)
      END IF

    END IF
  END IF

  ! check print output level idetail_src (1-5) with namelist given value (idbg_mxmn)
  ! for MIN/MAX output:

  !IF (idbg_mxmn < idetail_src .and. ltimer) CALL timer_stop(timer_print_mxmn)
  IF (idbg_mxmn < idetail_src ) RETURN

  ! surface level output only
  klev = 1

  ! print out maximum and minimum value

  ! parallelize:
  ctr=maxval(p_array(1:nproma,1:ndimblk))
  glbmx=global_max(ctr)
  ctr=minval(p_array(1:nproma,1:ndimblk))
  glbmn=global_min(ctr)

  IF (my_process_is_stdio()) &
    & WRITE(iout,991) ' MAX/MIN ', str_mod_src, strout, jk, glbmx, glbmn

  !IF (ltimer) CALL timer_stop(timer_print_mxmn)

  END SUBROUTINE dbg_print_2d

END MODULE mo_util_dbg_prnt

