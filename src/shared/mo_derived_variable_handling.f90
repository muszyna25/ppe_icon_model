!! Routines for handling proxy variables e.g. accumulation buffers
#include <omp_definitions.inc>

MODULE mo_derived_variable_handling

  USE mo_kind,                ONLY: wp, sp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_dynamics_config,     ONLY: nnow, nnew, nold
  USE mo_impl_constants,      ONLY: vname_len, REAL_T, TIMELEVEL_SUFFIX
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
                              & GRID_ZONAL, GRID_UNSTRUCTURED_VERT
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_zaxis_type,          ONLY: ZA_OCEAN_SEDIMENT
  USE mo_name_list_output_metadata, ONLY: metainfo_get_timelevel
  USE mo_var_list,            ONLY: add_var, t_var_list_ptr
  USE mo_var, ONLY: t_var, t_var_ptr, level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_var_metadata,        ONLY: get_var_name
  USE mo_var_metadata_types,  ONLY: t_var_metadata
  USE mo_var_list_register,   ONLY: vl_register
  USE mo_exception,           ONLY: finish
  USE mtime,                  ONLY: newEvent, event, isCurrentEventActive, newDatetime, datetime
  USE mo_output_event_types,  ONLY: t_sim_step_info
  USE mo_time_config,         ONLY: time_config
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_LONLAT, TSTEP_CONSTANT
  USE mo_util_texthash,       ONLY: text_hash_c
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif
#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_statistics_streams, finish_statistics_streams
  PUBLIC :: update_statistics, reset_statistics
  PUBLIC :: process_statistics_stream

  TYPE :: t_derivate_var
    TYPE(t_var_ptr) :: dst, src(3)
    INTEGER :: tls(3) = -999, counter = 0
  END TYPE t_derivate_var

  TYPE :: t_derivate_var_alloctble
    TYPE(t_derivate_var), ALLOCATABLE :: a
  END TYPE t_derivate_var_alloctble

  TYPE :: t_derivate_event
    TYPE(t_derivate_var_alloctble), ALLOCATABLE :: vars(:)
    TYPE(event), POINTER :: mtime_event
    CHARACTER(:), ALLOCATABLE :: eString
    INTEGER :: eKey = 0
  END TYPE t_derivate_event

  TYPE :: t_derivate_event_alloctble
    TYPE(t_derivate_event), ALLOCATABLE :: a
  END TYPE t_derivate_event_alloctble

  TYPE :: t_derivate_op
    TYPE(t_derivate_event_alloctble), ALLOCATABLE :: events(:)
    CHARACTER(:), ALLOCATABLE :: opname
    INTEGER :: opcode = -1
  END TYPE t_derivate_op

  TYPE(t_derivate_op) :: ops(4)
  CHARACTER(*), PARAMETER :: dlim = '|'
  CHARACTER(*), PARAMETER :: modname = 'mo_derived_variable_handling'

CONTAINS
  !! Create all needed lists and maps
  SUBROUTINE init_statistics_streams()

    ops(1)%opcode = 1
    ops(1)%opname = "mean"
    ops(2)%opcode = 2
    ops(2)%opname = "max"
    ops(3)%opcode = 3
    ops(3)%opname = "min"
    ops(4)%opcode = 4
    ops(4)%opname = "square"
  END SUBROUTINE init_statistics_streams

  SUBROUTINE finish_statistics_streams()
! just a stub
  END SUBROUTINE finish_statistics_streams

  SUBROUTINE reset_statistics()
! just a stub
  END SUBROUTINE reset_statistics

  SUBROUTINE process_statistics_stream(p_onl, i_typ, sim_step_info, patch_2d)
    TYPE(t_output_name_list), TARGET, INTENT(IN) :: p_onl
    INTEGER, INTENT(IN) :: i_typ
    TYPE(t_sim_step_info), INTENT(IN) :: sim_step_info
    TYPE(t_patch), INTENT(IN), TARGET :: patch_2d
    INTEGER :: i

    DO i = 1, 4
      CALL process_mvstream(ops(i), p_onl, i_typ, patch_2d)
    END DO
  END SUBROUTINE process_statistics_stream

  SUBROUTINE process_mvstream(bundle, p_onl, i_typ, patch_2d)
    TYPE(t_derivate_op), INTENT(INOUT), TARGET :: bundle
    TYPE(t_output_name_list), TARGET :: p_onl
    INTEGER, INTENT(IN) :: i_typ
    TYPE(t_patch), INTENT(IN), TARGET :: patch_2d
    CHARACTER(LEN=vname_len), POINTER :: in_varlist(:)
    INTEGER :: eKey, ie, iv, nv_scan, nv_new, nv_old, ne
    TYPE(t_derivate_event), POINTER :: ederiv
    TYPE(t_derivate_event_alloctble), ALLOCATABLE :: etmp(:)
    TYPE(t_derivate_var_alloctble), ALLOCATABLE :: vderiv(:), vtmp(:)
    TYPE(t_var), POINTER :: vl_elem
    CHARACTER(:), ALLOCATABLE :: dname, eString
    CHARACTER(LEN=1) :: dom_string
    TYPE(t_var_list_ptr) :: src_list
    CHARACTER(*), PARAMETER :: routine = modname//"::process_mvstream"

    IF (bundle%opname == TRIM(p_onl%operation)) THEN
      IF (ANY(1 < [p_onl%stream_partitions_ml, p_onl%stream_partitions_pl,  &
        &          p_onl%stream_partitions_hl, p_onl%stream_partitions_il])) &
        & CALL finish(routine, "only supported on global domain 1 " // &
        &                      "and without stream partitioning!")
      SELECT CASE(i_typ)
      CASE(level_type_ml)
        in_varlist => p_onl%ml_varlist
      CASE(level_type_pl)
        in_varlist => p_onl%pl_varlist
      CASE(level_type_hl)
        in_varlist => p_onl%hl_varlist
      CASE(level_type_il)
        in_varlist => p_onl%il_varlist
      END SELECT
      nv_scan = 0
      DO WHILE(.NOT.(in_varlist(nv_scan + 1) == ' '))
        nv_scan = nv_scan + 1
      END DO
      ALLOCATE(vderiv(nv_scan))
      ! uniq identifier for an event based on output start/end/interval
      eString = TRIM(p_onl%output_start(1)) // '_' // TRIM(p_onl%output_end(1)) // '_' // &
        &        TRIM(p_onl%output_interval(1))
      ! this has the advantage that we can compute a uniq id without creating the event itself
      ! fill main dictionary of variables for different event
      eKey = text_hash_c(eString)
      NULLIFY(ederiv)
      ne = 0
      IF (ALLOCATED(bundle%events)) ne = SIZE(bundle%events)
      DO ie = 2, ne
        IF (bundle%events(ie)%a%eKey .EQ. eKey) THEN
          IF (bundle%events(ie)%a%eString /= eString) THEN
            ederiv => bundle%events(ie)%a
            EXIT
          END IF
        END IF
      END DO
      nv_old = 0
      IF (.NOT.ASSOCIATED(ederiv)) THEN
        IF (ne .GT. 0) CALL MOVE_ALLOC(bundle%events, etmp)
        ALLOCATE(bundle%events(ne + 1))
        DO ie = 1, ne
          CALL MOVE_ALLOC(etmp(ie)%a, bundle%events(ie)%a)
        END DO
        ne = ne + 1
        ALLOCATE(bundle%events(ne)%a)
        ederiv => bundle%events(ne)%a
        ederiv%eString = eString
        ederiv%eKey = eKey
        ederiv%mtime_event => newEvent(eString, p_onl%output_start(1), &
          & p_onl%output_start(1), p_onl%output_end(1), p_onl%output_interval(1))
      ELSE
        IF (ALLOCATED(ederiv%vars)) THEN
          CALL MOVE_ALLOC(ederiv%vars, vtmp)
          nv_old = SIZE(vtmp)
        END IF
      END IF
      nv_new = 0
      DO iv = 1, nv_scan
        ! collect data variables only
        IF (INDEX(in_varlist(iv),':') > 0) CYCLE ! to avoid e.g. "grid:clon" stuff
        ! check for already created meanStream variable (maybe from another output_nml with the same output_interval)
        ! names consist of original spot-value names PLUS event information (start + interval of output)
        WRITE(dom_string, "(i1)") p_onl%dom
        dname = TRIM(in_varlist(iv)) // dlim // bundle%opname // dlim // &
          & TRIM(p_onl%output_interval(1)) // dlim // TRIM(p_onl%output_start(1)) &
          & // dlim // 'DOM' // dom_string
        vl_elem => vl_register%find_var_all(dname, opt_patch_id=p_onl%dom)
        IF (.NOT.ASSOCIATED(vl_elem)) THEN !not found -->> create a new one
          ALLOCATE(vderiv(nv_new+1)%a) ! staging a new var entry
          CALL find_src_element(TRIM(in_varlist(iv)), vderiv(nv_new+1)%a)
          IF (TSTEP_CONSTANT .EQ. vderiv(nv_new+1)%a%src(1)%p%info%isteptype) THEN
            DEALLOCATE(vderiv(nv_new+1)%a) ! discard staged entry
            dname = TRIM(in_varlist(iv)) ! no aggregation needed, since constant
          ELSE
            nv_new = nv_new + 1 ! new entry is valid, so keep
            CALL copy_var_to_list(vderiv(nv_new)%a) ! add_var to store accumulation
          END IF
        END IF
        in_varlist(iv) = dname
      END DO
      ALLOCATE(ederiv%vars(nv_new + nv_old))
      DO iv = 1, nv_old
        CALL MOVE_ALLOC(vtmp(iv)%a, ederiv%vars(iv)%a)
      END DO
      DO iv = 1, nv_new
        CALL MOVE_ALLOC(vderiv(iv)%a, ederiv%vars(iv+nv_old)%a)
      END DO
    END IF
  CONTAINS

  SUBROUTINE find_src_element(vname, deriv)
    CHARACTER(*), INTENT(IN) :: vname
    TYPE(t_derivate_var), INTENT(INOUT) :: deriv
    INTEGER :: k, l, tls(3)
    CHARACTER(4) :: tl_suff
    INTEGER, PARAMETER :: grids(5) = [GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
      & GRID_UNSTRUCTURED_VERT, GRID_LONLAT, GRID_ZONAL]

    DO k = 1, 5 ! scan for simple (instant) variable
      vl_elem => vl_register%find_var_all(vname, opt_patch_id=p_onl%dom, &
        & opt_hgrid=grids(k), opt_list=src_list, opt_cs=.FALSE., opt_output=.TRUE.)
      IF (ASSOCIATED(vl_elem)) THEN
        deriv%tls(1) = -1
        deriv%src(1)%p => vl_elem
        RETURN
      END IF
    END DO
    l = 0
    tls = [nold(1), nnow(1), nnew(1)]
    DO k = 1, 3 ! scan for time-levels of variable
      WRITE(tl_suff, '(a3,i1)') TIMELEVEL_SUFFIX, tls(k)
      vl_elem => vl_register%find_var_all(vname//tl_suff, &
        & opt_patch_id=p_onl%dom, opt_list=src_list, opt_output=.TRUE.)
      IF (ASSOCIATED(vl_elem)) THEN
        l = l + 1
        deriv%src(l)%p => vl_elem
        deriv%tls(l) = tls(k)
      END IF
    END DO
    IF (l .EQ. 0) &
      & CALL finish(routine, 'Could not find source variable: '//TRIM(vname))
  END SUBROUTINE find_src_element

  SUBROUTINE copy_var_to_list(deriv)
    TYPE(t_derivate_var), INTENT(INOUT) :: deriv
    TYPE(t_var_metadata), POINTER :: info

    info => deriv%src(1)%p%info
    CALL add_var(REAL_T, src_list, dname, info%hgrid, info%vgrid, info%cf, &
      & info%grib2, info%used_dimensions(1:info%ndims), vl_elem, &
      & tlev_source=info%tlev_source, isteptype=info%isteptype, &
      & post_op=info%post_op, initval_r=info%initval%rval, &
      & resetval_r=info%resetval%rval, lmiss=info%lmiss, lopenacc = .TRUE., &
      & missval_r=info%missval%rval, action_list=info%action_list, &
      & vert_interp=info%vert_interp, hor_interp=info%hor_interp, &
      & in_group=info%in_group, l_pp_scheduler_task=info%l_pp_scheduler_task, &
      & loutput=.TRUE., lrestart=.FALSE., var_class=info%var_class )
    __acc_attach(vl_elem%field%r_ptr)
    SELECT CASE(info%hgrid)
    CASE(GRID_UNSTRUCTURED_CELL)
      vl_elem%info%subset = patch_2d%cells%owned
    CASE(GRID_UNSTRUCTURED_EDGE)
      vl_elem%info%subset = patch_2d%edges%owned
    CASE(GRID_UNSTRUCTURED_VERT)
      vl_elem%info%subset = patch_2d%verts%owned
    END SELECT
    vl_elem%info%cf%datatype = &
      & MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    IF ("" == info%cf%short_name) &
      & vl_elem%info%cf%short_name = get_var_name(deriv%src(1)%p%info)
    deriv%dst%p => vl_elem
  END SUBROUTINE copy_var_to_list

  END SUBROUTINE process_mvstream

  SUBROUTINE perform_op(src, dest, opcode)
    TYPE(t_var), POINTER, INTENT(IN) :: src
    TYPE(t_var), POINTER, INTENT(INOUT) :: dest
    INTEGER, INTENT(IN) :: opcode
    INTEGER :: ic, sb, eb, br
    REAL(wp), POINTER :: sd5d(:,:,:,:,:)
    REAL(sp), POINTER :: ss5d(:,:,:,:,:)

    br = 0
    NULLIFY(sd5d, ss5d)
    sb = dest%info%subset%start_block
    eb = dest%info%subset%end_block
    IF (src%info%lcontained) THEN
      ic = src%info%ncontained
      SELECT CASE(src%info%var_ref_pos)
      CASE(1)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(ic:ic,:,:,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(ic:ic,:,:,:,:)
      CASE(2)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,ic:ic,:,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,ic:ic,:,:,:)
      CASE(3)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,ic:ic,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,ic:ic,:,:)
      CASE(4)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,:,ic:ic,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,:,ic:ic,:)
      END SELECT
    ELSE
      SELECT CASE(dest%info%ndims)
      CASE(3)
        br = 3
        ic = dest%info%used_dimensions(2)
        IF (ZA_OCEAN_SEDIMENT .EQ. dest%info%vgrid .OR. 1 .EQ. ic) THEN
          IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,1:ic,sb:eb,1:1,1:1)
          IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,1:ic,sb:eb,1:1,1:1)
        ELSE
          IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,sb:eb,1:1,1:1)
          IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,sb:eb,1:1,1:1)
        END IF
      CASE(2)
        IF (GRID_ZONAL .EQ. dest%info%hgrid) THEN
          IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,1:1,1:1,1:1)
          IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,1:1,1:1,1:1)
        ELSE
          br = 2
          IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,sb:eb,1:1,1:1,1:1)
          IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,sb:eb,1:1,1:1,1:1)
        END IF
      CASE DEFAULT
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,:,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,:,:,:)
      END SELECT
    END IF
    eb = eb + 1 - sb
    CALL perform_op_5d()
  CONTAINS

  SUBROUTINE perform_op_5d()
    INTEGER :: ni, j, k, l, m, lsi, lei, blk, si, ei
    LOGICAL :: ls, blk_is3
    REAL(wp), POINTER :: tmp1(:,:,:,:,:), tmp2(:,:,:,:,:)

    ls = br .NE. 0
    blk_is3 = br .EQ. 3
    si = dest%info%subset%start_index
    ei = dest%info%subset%end_index
    ni = SIZE(dest%r_ptr, 1)
    IF (ASSOCIATED(sd5d)) THEN
      tmp1 => sd5d
    ELSE
      ALLOCATE(tmp1(SIZE(ss5d,1), SIZE(ss5d,2), SIZE(ss5d,3), SIZE(ss5d,4), SIZE(ss5d,5)))
!ICON_OMP PARALLEL PRIVATE(j,k,l,m,blk,lsi,lei)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp1(lsi:lei,j,k,l,m) = REAL(ss5d(lsi:lei,j,k,l,m), wp)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
    END IF
    tmp2 => dest%r_ptr
!ICON_OMP PARALLEL PRIVATE(j,k,l,m,blk,lsi,lei)
    SELECT CASE(opcode)
    CASE(1)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j,k,l,m) = &
                & tmp2(lsi:lei,j,k,l,m) + tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(2)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              WHERE(tmp1(lsi:lei,j,k,l,m) .GT. tmp2(lsi:lei,j,k,l,m)) &
                & tmp2(lsi:lei,j,k,l,m) = tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(3)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              WHERE(tmp1(lsi:lei,j,k,l,m) .LT. tmp2(lsi:lei,j,k,l,m)) &
                & tmp2(lsi:lei,j,k,l,m) = tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(4)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j,k,l,m) = tmp2(lsi:lei,j,k,l,m) + &
                & tmp1(lsi:lei,j,k,l,m) * tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(5)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j,k,l,m) = tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    END SELECT
!ICON_OMP END PARALLEL
    IF (ASSOCIATED(ss5d)) DEALLOCATE(tmp1)
  END SUBROUTINE perform_op_5d

  END SUBROUTINE perform_op

  !! Execute the accumulation forall internal variables and compute mean values
  !! if the corresponding event is active
  SUBROUTINE update_statistics()
    INTEGER :: i

    DO i = 1, 4
      CALL update_mvstream(ops(i))
    END DO
  END SUBROUTINE update_statistics

  SUBROUTINE update_mvstream(bundle)
    TYPE(t_derivate_op), INTENT(INOUT), TARGET :: bundle
    INTEGER :: tl, iv, ie, it, ne, nv, j, k, l, m
    INTEGER, POINTER :: ct
    TYPE(t_var), POINTER :: src, dst
    TYPE(t_derivate_event), POINTER :: ederiv
    TYPE(datetime), POINTER :: mtime_date
    LOGICAL :: isactive
    REAL(wp) :: weight
    CHARACTER(*), PARAMETER :: routine = modname//"::update_statistics"

    ne = 0
    IF (ALLOCATED(bundle%events)) ne = SIZE(bundle%events)
    mtime_date => newDatetime(time_config%tc_current_date)
    DO ie = 1, ne
      ederiv => bundle%events(ie)%a
      isactive = LOGICAL(isCurrentEventActive(ederiv%mtime_event, mtime_date))
      nv = 0
      IF (ALLOCATED(ederiv%vars)) nv = SIZE(ederiv%vars)
      DO iv = 1, nv
        dst => ederiv%vars(iv)%a%dst%p
        src => ederiv%vars(iv)%a%src(1)%p
        IF (ederiv%vars(iv)%a%tls(1) .NE. -1) THEN
          tl = metainfo_get_timelevel(dst%info, dst%info%dom)
          DO it = 1, 3
            IF (tl .EQ. ederiv%vars(iv)%a%tls(it)) THEN
              src => ederiv%vars(iv)%a%src(it)%p
              EXIT
            END IF
          END DO
        END IF
        ct => ederiv%vars(iv)%a%counter
        IF (ct .EQ. 0) THEN
          IF (bundle%opcode .EQ. 4) THEN
!ICON_OMP PARALLEL
!ICON_OMP DO PRIVATE(j,k,l,m) COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
            DO m = 1, SIZE(dst%r_ptr,5)
              DO l = 1, SIZE(dst%r_ptr,4)
                DO k = 1, SIZE(dst%r_ptr,3)
                  DO j = 1, SIZE(dst%r_ptr,2)
                    dst%r_ptr(:,j,k,l,m) = 0.0_wp
                  END DO
                END DO
              END DO
            END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
            CALL perform_op(src, dst, 4)
          ELSE
            CALL perform_op(src, dst, 5) ! initial assignment
          END IF
        ELSE
          CALL perform_op(src, dst, bundle%opcode)
        END IF
        ct = ct + 1
        IF (isactive) THEN ! output step, so weighting applied this time
!ICON_OMP PARALLEL PRIVATE(j,k,l,m, weight)
          IF (1 .EQ. bundle%opcode .OR. 4 .EQ. bundle%opcode &
            & .AND. ct .GT. 0) THEN
            weight = 1._wp / REAL(ct, wp)
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
            DO m = 1, SIZE(dst%r_ptr,5)
              DO l = 1, SIZE(dst%r_ptr,4)
                DO k = 1, SIZE(dst%r_ptr,3)
                  DO j = 1, SIZE(dst%r_ptr,2)
                    dst%r_ptr(:,j,k,l,m) = dst%r_ptr(:,j,k,l,m) * weight
                  END DO
                END DO
              END DO
            END DO
!ICON_OMP END DO
          END IF
          IF (dst%info%lmiss) THEN ! set missval where applicable
            IF (ASSOCIATED(src%r_ptr)) THEN
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
              DO m = 1, SIZE(dst%r_ptr,5)
                DO l = 1, SIZE(dst%r_ptr,4)
                  DO k = 1, SIZE(dst%r_ptr,3)
                    DO j = 1, SIZE(dst%r_ptr,2)
                      WHERE(dst%info%missval%rval .EQ. src%r_ptr(:,j,k,l,m)) &
                        & dst%r_ptr(:,j,k,l,m) = dst%info%missval%rval
                    END DO
                  END DO
                END DO
              END DO
!ICON_OMP END DO NOWAIT
            ELSE
!ICON_OMP DO COLLAPSE(4)
!$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(4) IF(i_am_accel_node)
              DO m = 1, SIZE(dst%r_ptr,5)
                DO l = 1, SIZE(dst%r_ptr,4)
                  DO k = 1, SIZE(dst%r_ptr,3)
                    DO j = 1, SIZE(dst%r_ptr,2)
                      WHERE(src%info%missval%sval .EQ. src%s_ptr(:,j,k,l,m)) &
                        & dst%r_ptr(:,j,k,l,m) = dst%info%missval%rval
                    END DO
                  END DO
                END DO
              END DO
!ICON_OMP END DO NOWAIT
            END IF
          END IF
!ICON_OMP END PARALLEL
          ct = 0
        END IF
      END DO
    END DO
  END SUBROUTINE update_mvstream

END MODULE mo_derived_variable_handling
