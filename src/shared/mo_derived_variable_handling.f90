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
  USE mo_var,                 ONLY: t_var, t_var_ptr
  USE mo_var_metadata,        ONLY: get_var_name
  USE mo_var_metadata_types,  ONLY: t_var_metadata
  USE mo_var_list_register_utils, ONLY: vlr_find
  USE mo_exception,           ONLY: finish
  USE mtime,                  ONLY: newEvent, event, isCurrentEventActive, datetime
  USE mo_time_config,         ONLY: time_config
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_LONLAT, TSTEP_CONSTANT
  USE mo_util_texthash,       ONLY: text_hash_c
! HB: commented openACC stuff for now -- due to weird memory issues, if nproma is large
#ifdef _OPENACC
  USE mo_mpi,                 ONLY: i_am_accel_node
#endif
!#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_statistics
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
    CHARACTER(LEN=vname_len) :: eString
    INTEGER :: eKey = 0
  END TYPE t_derivate_event

  TYPE :: t_derivate_event_alloctble
    TYPE(t_derivate_event), ALLOCATABLE :: a
  END TYPE t_derivate_event_alloctble

  TYPE :: t_derivate_op
    TYPE(t_derivate_event_alloctble), ALLOCATABLE :: events(:)
  END TYPE t_derivate_op

  INTEGER, PARAMETER :: nops = 4
  TYPE(t_derivate_op), TARGET :: ops(nops)
  CHARACTER(*), PARAMETER :: dlim = '|'
  CHARACTER(*), PARAMETER :: modname = 'mo_derived_variable_handling'
  CHARACTER(*), PARAMETER :: opnames(nops) = ["mean  ", "max   ", "min   ", "square"]
  INTEGER, PARAMETER :: oplen(nops) = [4, 3, 3, 6]

CONTAINS

  SUBROUTINE process_statistics_stream(p_onl, vlist, patch_2d)
    TYPE(t_output_name_list), INTENT(IN) :: p_onl
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: vlist(:)
    TYPE(t_patch), INTENT(IN), TARGET :: patch_2d
    INTEGER :: i

    DO i = 1, nops
      IF (TRIM(p_onl%operation) == opnames(i)(1:oplen(i))) &
        & CALL process_mvstream(i, p_onl, vlist, patch_2d)
    ENd DO
  END SUBROUTINE process_statistics_stream

  SUBROUTINE process_mvstream(iop, p_onl, in_varlist, patch_2d)
    INTEGER, INTENT(IN) :: iop
    TYPE(t_output_name_list), TARGET :: p_onl
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: in_varlist(:)
    TYPE(t_patch), INTENT(IN), TARGET :: patch_2d
    INTEGER :: eKey, ie, iv, nv_scan, nv_new, nv_old, ne, dns_len, it_len, st_len
    TYPE(t_derivate_event), POINTER :: ederiv
    TYPE(t_derivate_event_alloctble), ALLOCATABLE :: etmp(:)
    TYPE(t_derivate_var_alloctble), ALLOCATABLE :: vderiv(:), vtmp(:)
    TYPE(t_var), POINTER :: vl_elem
    CHARACTER(LEN=vname_len) :: dname, eString, dname_suffix
    TYPE(t_var_list_ptr) :: src_list
    CHARACTER(*), PARAMETER :: routine = modname//"::process_mvstream"

    IF (ANY(1 < [p_onl%stream_partitions_ml, p_onl%stream_partitions_pl,  &
      &          p_onl%stream_partitions_hl, p_onl%stream_partitions_il])) &
      & CALL finish(routine, "only supported on global domain 1 " // &
      &                      "and without stream partitioning!")
    nv_scan = 0
    DO WHILE(.NOT.(in_varlist(nv_scan + 1) == ' '))
      nv_scan = nv_scan + 1
    END DO
    ALLOCATE(vderiv(nv_scan))
    ! uniq identifier for an event based on output start/end/interval
    st_len = LEN_TRIM(p_onl%output_start(1))
    it_len = LEN_TRIM(p_onl%output_interval(1))
    WRITE(eString, "(5a)") p_onl%output_start(1)(:st_len), '_', &
      & TRIM(p_onl%output_end(1)), '_', p_onl%output_interval(1)(:it_len)
    WRITE(dname_suffix, "(8a,i0)") dlim, opnames(iop)(1:oplen(iop)), dlim, &
        & p_onl%output_interval(1)(:it_len), dlim, &
        & p_onl%output_start(1)(:st_len), dlim, 'DOM', p_onl%dom
    dns_len = LEN_TRIM(dname_suffix)
    ! this has the advantage that we can compute a uniq id without creating the event itself
    ! fill main dictionary of variables for different event
    eKey = text_hash_c(eString)
    NULLIFY(ederiv)
    ne = 0
    IF (ALLOCATED(ops(iop)%events)) ne = SIZE(ops(iop)%events)
    DO ie = 2, ne
      IF (ops(iop)%events(ie)%a%eKey .EQ. eKey) THEN
        IF (ops(iop)%events(ie)%a%eString /= eString) THEN
          ederiv => ops(iop)%events(ie)%a
          EXIT
        END IF
      END IF
    END DO
    nv_old = 0
    IF (.NOT.ASSOCIATED(ederiv)) THEN
      IF (ne .GT. 0) CALL MOVE_ALLOC(ops(iop)%events, etmp)
      ALLOCATE(ops(iop)%events(ne + 1))
      DO ie = 1, ne
        CALL MOVE_ALLOC(etmp(ie)%a, ops(iop)%events(ie)%a)
      END DO
      ne = ne + 1
      ALLOCATE(ops(iop)%events(ne)%a)
      ederiv => ops(iop)%events(ne)%a
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
      WRITE(dname, "(2a)") TRIM(in_varlist(iv)), dname_suffix(:dns_len)
      vl_elem => vlr_find(dname, opt_patch_id=p_onl%dom)
      IF (.NOT.ASSOCIATED(vl_elem)) THEN !not found -->> create a new one
        ALLOCATE(vderiv(nv_new+1)%a) ! staging a new var entry
        CALL find_src_element(TRIM(in_varlist(iv)), vderiv(nv_new+1)%a)
        IF (TSTEP_CONSTANT .EQ. vderiv(nv_new+1)%a%src(1)%p%info%isteptype) THEN
          DEALLOCATE(vderiv(nv_new+1)%a) ! discard staged entry
          dname = in_varlist(iv) ! no aggregation needed, since constant
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
  CONTAINS

  SUBROUTINE find_src_element(vname, deriv)
    CHARACTER(*), INTENT(IN) :: vname
    TYPE(t_derivate_var), INTENT(INOUT) :: deriv
    INTEGER :: k, l, tls(3)
    CHARACTER(4) :: tl_suff
    INTEGER, PARAMETER :: grids(5) = [GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
      & GRID_UNSTRUCTURED_VERT, GRID_LONLAT, GRID_ZONAL]

    DO k = 1, 5 ! scan for simple (instant) variable
      vl_elem => vlr_find(vname, opt_patch_id=p_onl%dom, &
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
      vl_elem => vlr_find(vname//tl_suff, &
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
      & resetval_r=info%resetval%rval, lmiss=info%lmiss, &
      & missval_r=info%missval%rval, action_list=info%action_list, &
      & vert_interp=info%vert_interp, hor_interp=info%hor_interp, &
      & in_group=info%in_group, l_pp_scheduler_task=info%l_pp_scheduler_task, &
      & loutput=.TRUE., lrestart=.FALSE., var_class=info%var_class, &
      & lopenacc = .FALSE. )
!    __acc_attach(vl_elem%r_ptr)
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

  SUBROUTINE perform_op(src, dest, opcode, weight, miss, miss_s)
    TYPE(t_var), POINTER, INTENT(IN) :: src
    TYPE(t_var), POINTER, INTENT(INOUT) :: dest
    INTEGER, INTENT(IN) :: opcode
    REAL(wp), INTENT(IN), OPTIONAL :: weight, miss
    REAL(sp), INTENT(IN), OPTIONAL :: miss_s
    REAL(wp) :: miss_src
    INTEGER :: ic, sb, eb, br
    REAL(wp), POINTER :: sd5d(:,:,:,:,:)
    REAL(sp), POINTER :: ss5d(:,:,:,:,:)
    CHARACTER(*), PARAMETER :: routine = modname//":perform_op"

    IF (.NOT.PRESENT(weight) .AND. opcode .EQ. 7) THEN
      CALL finish(routine, "no weight factor provided")
    ELSE IF (PRESENT(weight) .AND. opcode .NE. 7) THEN
      CALL finish(routine, "no weight factor allowed for this op")
    END IF
    IF (opcode .EQ. 8 .AND. .NOT.(PRESENT(miss) .AND. PRESENT(miss_s))) THEN
      CALL finish(routine, "no missing values provided")
    ELSE IF((PRESENT(miss) .OR. PRESENT(miss_s)) .AND. opcode .NE. 8) THEN
      CALL finish(routine, "no missing values allowed for this op")
    END IF
    br = 0
    NULLIFY(sd5d, ss5d)
    IF (opcode .NE. 8) THEN
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
    ELSE
      miss_src = MERGE(miss, REAL(miss_s, wp), ASSOCIATED(src%r_ptr))
      IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,:,:,:)
      IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,:,:,:)
    END IF
    CALL perform_op_5d()
  CONTAINS

  SUBROUTINE perform_op_5d()
    INTEGER :: ni, j, k, l, m, lsi, lei, blk, si, ei, joff, koff
    LOGICAL :: ls, blk_is3
    REAL(wp), POINTER :: tmp1(:,:,:,:,:), tmp2(:,:,:,:,:)
    REAL(wp) :: miss__, weight__

    ls = br .NE. 0
    blk_is3 = br .EQ. 3
    joff = MERGE(sb - 1, 0, ls .AND. .NOT.blk_is3)
    koff = MERGE(sb - 1, 0, ls .AND. blk_is3)
    si = dest%info%subset%start_index
    ei = dest%info%subset%end_index
    ni = SIZE(dest%r_ptr, 1)
    IF (ASSOCIATED(sd5d)) THEN
      !$ACC UPDATE HOST(src%r_ptr) IF(src%info%lopenacc .AND. i_am_accel_node)
      tmp1 => sd5d
    ELSE
      !$ACC UPDATE HOST(src%s_ptr) IF(src%info%lopenacc .AND. i_am_accel_node)
      ALLOCATE(tmp1(SIZE(ss5d,1), SIZE(ss5d,2), SIZE(ss5d,3), SIZE(ss5d,4), SIZE(ss5d,5)))
!ICON_OMP PARALLEL PRIVATE(j,k,l,m,blk,lsi,lei)
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
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
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j+joff,k+koff,l,m) = &
                & tmp2(lsi:lei,j+joff,k+koff,l,m) + tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(2)
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              WHERE(tmp1(lsi:lei,j,k,l,m) .GT. tmp2(lsi:lei,j+joff,k+koff,l,m)) &
                & tmp2(lsi:lei,j+joff,k+koff,l,m) = tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(3)
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              WHERE(tmp1(lsi:lei,j,k,l,m) .LT. tmp2(lsi:lei,j+joff,k+koff,l,m)) &
                & tmp2(lsi:lei,j+joff,k+koff,l,m) = tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(4)
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j+joff,k+koff,l,m) = tmp2(lsi:lei,j+joff,k+koff,l,m) + &
                & tmp1(lsi:lei,j,k,l,m) * tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(5)
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j+joff,k+koff,l,m) = tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(6)
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j+joff,k+koff,l,m) = &
                & tmp1(lsi:lei,j,k,l,m) * tmp1(lsi:lei,j,k,l,m)
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(7)
      weight__ = weight
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp1,5)
        DO l = 1, SIZE(tmp1,4)
          DO k = 1, SIZE(tmp1,3)
            DO j = 1, SIZE(tmp1,2)
              blk = MERGE(k, j, blk_is3)
              lsi = MERGE(si, 1,  ls .AND. blk .EQ. 1)
              lei = MERGE(ei, ni, ls .AND. blk .EQ. eb)
              tmp2(lsi:lei,j+joff,k+koff,l,m) = &
                & tmp2(lsi:lei,j+joff,k+koff,l,m) * weight__
            END DO
          END DO
        END DO
      END DO
!ICON_OMP END DO NOWAIT
    CASE(8)
      miss__ = miss
!ICON_OMP DO COLLAPSE(4)
!!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) GANG VECTOR COLLAPSE(4) ASYNC(1) IF(i_am_accel_node)
      DO m = 1, SIZE(tmp2,5)
        DO l = 1, SIZE(tmp2,4)
          DO k = 1, SIZE(tmp2,3)
            DO j = 1, SIZE(tmp2,2)
              WHERE(tmp1(:,j,k,l,m) .EQ. miss_src) &
                tmp2(:,j,k,l,m) = miss__
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

    DO i = 1, nops
      CALL update_mvstream(i)
    END DO
  END SUBROUTINE update_statistics

  SUBROUTINE update_mvstream(iop)
    INTEGER, INTENT(IN) :: iop
    INTEGER :: tl, iv, ie, it, ne, nv
    INTEGER, POINTER :: ct
    TYPE(t_var), POINTER :: src, dst
    TYPE(t_derivate_event), POINTER :: ederiv
    TYPE(datetime) :: mtime_date
    LOGICAL :: isactive

    ne = 0
    IF (ALLOCATED(ops(iop)%events)) ne = SIZE(ops(iop)%events)
    mtime_date = time_config%tc_current_date
    DO ie = 1, ne
      ederiv => ops(iop)%events(ie)%a
      isactive = LOGICAL(isCurrentEventActive(ederiv%mtime_event, mtime_date))
      nv = 0
      IF (ALLOCATED(ederiv%vars)) nv = SIZE(ederiv%vars)
      DO iv = 1, nv
        dst => ederiv%vars(iv)%a%dst%p
        ct => ederiv%vars(iv)%a%counter
        src => ederiv%vars(iv)%a%src(1)%p
        IF (ederiv%vars(iv)%a%tls(1) .NE. -1) THEN
          tl = metainfo_get_timelevel(dst%info, dst%info%dom)
          it = MAXLOC(MERGE(1, 0, tl .EQ. ederiv%vars(iv)%a%tls(:)), 1)
          src => ederiv%vars(iv)%a%src(it)%p
        END IF
        IF (ct .EQ. 0) THEN ! initial assignment
          CALL perform_op(src, dst, MERGE(6, 5, iop .EQ. 4))
        ELSE ! actual update
          CALL perform_op(src, dst, iop)
        END IF
        ct = ct + 1
        IF (isactive) THEN ! output step, so weighting is applied this time
          IF ((1 .EQ. iop .OR. 4 .EQ. iop) .AND. ct .GT. 0) &
            & CALL perform_op(src, dst, 7, weight=(1._wp / REAL(ct, wp)))
          IF (dst%info%lmiss) & ! (re)set missval where applicable
            & CALL perform_op(src, dst, 8, miss=dst%info%missval%rval, &
                &             miss_s=src%info%missval%sval)
          ct = 0
        END IF
      END DO
    END DO
  END SUBROUTINE update_mvstream

END MODULE mo_derived_variable_handling
