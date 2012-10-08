
! #ifdef __xlC__
! @PROCESS HOT
! #endif
#ifdef __PGI
!pgi$g opt=1
#endif
!>
!! Contains the interpolation routines needed for grid refinement.
!!
!! These had originally been included in mo_grf_interpolation but then were
!! packed into a separate module to clean up the code
!!
!! @par Revision History
!! Created by Guenther Zaengl, DWD (2009-02-09)
!! Modification by Guenther Zaengl, DWD (2009-06-22)
!! - preparation for generalized grid refinement (affects all subroutines)
!! Modification by Constantin Junk, MPI-M (2011-05-05)
!! - moved setup_gridref to mo_gridref_nml
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
MODULE mo_grf_intp_state
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int, min_rledge_int
USE mo_model_domain,        ONLY: t_patch, p_patch_local_parent
USE mo_grid_config,         ONLY: n_dom, n_dom_start
USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e
USE mo_parallel_config,     ONLY: nproma
USE mo_mpi,                 ONLY: my_process_is_mpi_parallel, p_pe

USE mo_communication,       ONLY: t_comm_pattern, blk_no, idx_no, idx_1d, &
  &                               setup_comm_pattern, delete_comm_pattern, exchange_data


USE mo_grf_intp_data_strc
USE mo_grf_intp_coeffs
USE mo_gridref_config



IMPLICIT NONE

PRIVATE

PUBLIC ::construct_2d_gridref_state, destruct_2d_gridref_state
PUBLIC ::allocate_grf_state, deallocate_grf_state
PUBLIC ::transfer_grf_state

TYPE(t_comm_pattern) :: comm_pat_loc_to_glb_c, comm_pat_loc_to_glb_e

CONTAINS

!-------------------------------------------------------------------------
!
!
!>
!! Allocation of fields needed for grid refinement.
!!
!! @par Revision History
!! Split off from construct_2d_gridref_state, Rainer Johanni (2010-10-29)
!!
SUBROUTINE allocate_grf_state( ptr_patch, ptr_grf )
!
TYPE(t_patch), INTENT(IN) :: ptr_patch
TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: ptr_grf

TYPE(t_gridref_single_state), POINTER :: ptr_int    => NULL()

INTEGER :: jcd
INTEGER :: nblks_c, nblks_e, nblks_v, &
           nproma_grf, isb_e, ieb_e, isb_c, ieb_c, i_nchdom
INTEGER :: nlev              !< number of full levels
INTEGER :: ist
INTEGER :: grf_vec_dim_1, grf_vec_dim_2

!-----------------------------------------------------------------------

  grf_vec_dim_1 = 6
  grf_vec_dim_2 = 5

  ! determine size of arrays, i.e.
  ! values for the blocking
  !
  nblks_c  = ptr_patch%nblks_c
  nblks_e  = ptr_patch%nblks_e
  nblks_v  = ptr_patch%nblks_v

  ! number of vertical levels
  nlev   = ptr_patch%nlev

  ! determine number of child domains
  i_nchdom = MAX(1,ptr_patch%n_childdom)

  ! build data structure for grid refinement for each child domain
  ALLOCATE(ptr_grf%p_dom(i_nchdom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_intp_state:allocate_grf_state',  &
      &             'allocation for p_dom failed')
  ENDIF

  ! Now loop over all child domains of the present domain
  CD_LOOP: DO jcd = 1, i_nchdom

    ptr_int => ptr_grf%p_dom(jcd)

    ! Determine dimensions needed for grid refinement fields
    IF (ptr_patch%n_childdom > 0) THEN  ! n_childdom = 0 in the innermost nest(s)

      isb_e      = ptr_patch%edges%start_blk(grf_bdyintp_start_e,jcd)
      ieb_e      = ptr_patch%edges%end_blk(min_rledge_int,jcd)
      isb_c      = ptr_patch%cells%start_blk(grf_bdyintp_start_c,jcd)
      ieb_c      = ptr_patch%cells%end_blk(min_rlcell_int,jcd)
      nproma_grf = nproma

    ELSE

      isb_e      = 1
      ieb_e      = 1
      isb_c      = 1
      ieb_c      = 1
      nproma_grf = 1

    ENDIF

    ! index arrays (external points; global indices) for grid refinement interpolation

    ALLOCATE (ptr_int%grf_vec_ind_1a(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_ind_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_ind_1b(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_ind_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_ind_2a(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_ind_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_ind_2b(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_ind_2b failed')
    ENDIF

    ! index arrays (external points; global indices) for grid refinement interpolation

    ALLOCATE (ptr_int%grf_vec_blk_1a(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_blk_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_blk_1b(nproma_grf, grf_vec_dim_1, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_blk_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_blk_2a(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_blk_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_blk_2b(nproma_grf, grf_vec_dim_2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_blk_2b failed')
    ENDIF

    ! arrays carrying the actual number of edges available for reconstruction

    ALLOCATE (ptr_int%grf_vec_stencil_1a(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_stencil_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_stencil_1b(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_stencil_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_stencil_2a(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_stencil_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_stencil_2b(nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_stencil_2b failed')
    ENDIF

    ! arrays carrying the interpolation coefficients for reconstruction related to grid refinement

    ALLOCATE (ptr_int%grf_vec_coeff_1a(grf_vec_dim_1, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_coeff_1a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_coeff_1b(grf_vec_dim_1, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_coeff_1b failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_coeff_2a(grf_vec_dim_2, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_coeff_2a failed')
    ENDIF
    ALLOCATE (ptr_int%grf_vec_coeff_2b(grf_vec_dim_2, nproma_grf, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_vec_coeff_2b failed')
    ENDIF

    ! Distances from parent cell to child cells
    ALLOCATE (ptr_int%grf_dist_pc2cc(nproma_grf, 4, 2, isb_c:ieb_c), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_dist_pc2cc failed')
    ENDIF

    ! Normalized distances from parent edge to child edges 1 and 2
    ALLOCATE (ptr_int%grf_dist_pe2ce(nproma_grf, 2, isb_e:ieb_e), STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for grf_dist_pe2ce failed')
    ENDIF

    ptr_int%grf_vec_ind_1a   = 0
    ptr_int%grf_vec_ind_1b   = 0
    ptr_int%grf_vec_ind_2a   = 0
    ptr_int%grf_vec_ind_2b   = 0

    ptr_int%grf_vec_blk_1a   = 0
    ptr_int%grf_vec_blk_1b   = 0
    ptr_int%grf_vec_blk_2a   = 0
    ptr_int%grf_vec_blk_2b   = 0

    ptr_int%grf_vec_stencil_1a   = 0
    ptr_int%grf_vec_stencil_1b   = 0
    ptr_int%grf_vec_stencil_2a   = 0
    ptr_int%grf_vec_stencil_2b   = 0

    ptr_int%grf_vec_coeff_1a   = 0._wp
    ptr_int%grf_vec_coeff_1b   = 0._wp
    ptr_int%grf_vec_coeff_2a   = 0._wp
    ptr_int%grf_vec_coeff_2b   = 0._wp

    ptr_int%grf_dist_pc2cc   = 0._wp
    ptr_int%grf_dist_pe2ce   = 0._wp

  ENDDO CD_LOOP

  ! Feedback weights for cell-based variables
  ALLOCATE (ptr_grf%fbk_wgt_c(nproma,nblks_c,4), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for fbk_wgt_c failed')
  ENDIF
  ALLOCATE (ptr_grf%fbk_wgt_ct(nproma,nblks_c,4), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for fbk_wgt_ct failed')
  ENDIF

  ! Feedback weights for edge-based variables
  ALLOCATE (ptr_grf%fbk_wgt_e(nproma,nblks_e,6), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for fbk_wgt_e failed')
  ENDIF

  ALLOCATE (ptr_grf%fbk_dom_area(i_nchdom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for fbk_dom_area failed')
  ENDIF

  ptr_grf%fbk_wgt_c     = 0._wp
  ptr_grf%fbk_wgt_ct    = 0._wp
  ptr_grf%fbk_wgt_e     = 0._wp
  ptr_grf%fbk_dom_area  = 0._wp

END SUBROUTINE allocate_grf_state

!-------------------------------------------------------------------------
!
!
!>
!!               Allocation of fields needed for grid refinement.
!!
!!               Initialization of components.
!!
!! @par Revision History
!! Created  by  Guenther Zaengl, DWD (2009-02-11).
!! Modified by Rainer Johanni (2010-10-29): split out allocation
!!
SUBROUTINE construct_2d_gridref_state( ptr_patch, ptr_grf_state )
!
TYPE(t_patch), TARGET, INTENT(INOUT) :: ptr_patch(n_dom_start:)
TYPE(t_gridref_state), INTENT(INOUT) :: ptr_grf_state(n_dom_start:)

INTEGER :: jg

!-----------------------------------------------------------------------

CALL message('mo_grf_intp_state:construct_2d_gridref_state', &
  & 'start to construct grf_state')

DO jg = n_dom_start, n_dom
  CALL allocate_grf_state(ptr_patch(jg), ptr_grf_state(jg))
ENDDO
IF(my_process_is_mpi_parallel()) THEN
  DO jg = n_dom_start+1, n_dom
    CALL allocate_grf_state(p_patch_local_parent(jg), p_grf_state_local_parent(jg))
  ENDDO
ENDIF

CALL message ('mo_grf_intp_state:construct_2d_gridref_state',   &
  & 'memory allocation finished')

IF (n_dom_start == 0 .OR. n_dom > 1) THEN
  ! set the patch that grf_intp_coeffs will work on
  CALL grf_intp_coeffs_setpatch(ptr_patch)

  CALL gridref_info ( ptr_grf_state)

  CALL init_fbk_wgt ( ptr_grf_state)

  CALL compute_pc2cc_distances ( ptr_grf_state)
  CALL compute_pe2ce_distances ( ptr_grf_state)

  CALL grf_index( ptr_grf_state)
  IF (grf_intmethod_e == 2 .OR. grf_intmethod_e == 4) THEN
    CALL rbf_compute_coeff_grf_e ( ptr_grf_state)
  ELSE IF (grf_intmethod_e == 1 .OR. grf_intmethod_e == 3) THEN
    CALL idw_compute_coeff_grf_e ( ptr_grf_state)
  ENDIF

ENDIF

CALL message ('mo_grf_intp_state:construct_2d_gridref_state',                        &
  & 'construction of interpolation state finished')

END SUBROUTINE construct_2d_gridref_state

!-------------------------------------------------------------------------
!!
!>
!! Get local idx_no for edges
!!
ELEMENTAL INTEGER FUNCTION loc_idx_no_e(p_p, idx)
  TYPE(t_patch), INTENT(IN) :: p_p
  INTEGER, INTENT(IN) :: idx

  IF(idx<1 .OR. idx > p_p%n_patch_edges_g) THEN
    loc_idx_no_e = 0
  ELSE
    loc_idx_no_e = idx_no(p_p%edges%loc_index(idx))
  ENDIF
END FUNCTION loc_idx_no_e
!-------------------------------------------------------------------------
!!
!>
!! Get local blk_no for edges
!!
ELEMENTAL INTEGER FUNCTION loc_blk_no_e(p_p, idx)
  TYPE(t_patch), INTENT(IN) :: p_p
  INTEGER, INTENT(IN) :: idx

  IF(idx<1 .OR. idx > p_p%n_patch_edges_g) THEN
    loc_blk_no_e = 0
  ELSE
    loc_blk_no_e = blk_no(p_p%edges%loc_index(idx))
  ENDIF
END FUNCTION loc_blk_no_e
!-------------------------------------------------------------------------
!!
!>
!! Get global idx_1d for edges
!!
ELEMENTAL INTEGER FUNCTION glb_idx_1d_e(p_p, idx, blk)
  TYPE(t_patch), INTENT(IN) :: p_p
  INTEGER, INTENT(IN) :: idx, blk

  IF(idx<=0 .or. blk<=0 .or. idx_1d(idx, blk)>p_p%n_patch_edges) THEN
    glb_idx_1d_e = 0
  ELSE
    glb_idx_1d_e = p_p%edges%glb_index(idx_1d(idx, blk))
  ENDIF
END FUNCTION glb_idx_1d_e
!-------------------------------------------------------------------------
!!
!>
!!  Transfers interpolation state from local parent to global parent
!!  Some variables have to be transferred the other way round - see below
!!
!! @par Revision History
!! Created by Rainer Johanni (2011-10-26)
!!
SUBROUTINE transfer_grf_state(p_p, p_lp, p_grf, p_lgrf, jcd)
!
  TYPE(t_patch), INTENT(IN) :: p_p
  TYPE(t_patch), INTENT(IN) :: p_lp

  TYPE(t_gridref_state), INTENT(INOUT) :: p_grf
  TYPE(t_gridref_state), INTENT(INOUT) :: p_lgrf

  INTEGER, INTENT(IN) :: jcd

  INTEGER, ALLOCATABLE :: owner(:)
  INTEGER :: j, k, n
  INTEGER :: isb_e, ieb_e, isb_le, ieb_le, is_e, ie_e
  INTEGER :: isb_c, ieb_c, isb_lc, ieb_lc, is_c, ie_c

  REAL(wp), ALLOCATABLE :: z_tmp_s(:,:,:), z_tmp_r(:,:,:)

  INTEGER :: grf_vec_dim_1, grf_vec_dim_2

  grf_vec_dim_1 = 6
  grf_vec_dim_2 = 5


  isb_e      = p_p%edges%start_blk(grf_bdyintp_start_e,jcd)
  ieb_e      = p_p%edges%end_blk(min_rledge_int,jcd)
  isb_c      = p_p%cells%start_blk(grf_bdyintp_start_c,jcd)
  ieb_c      = p_p%cells%end_blk(min_rlcell_int,jcd)

  isb_le      = p_lp%edges%start_blk(grf_bdyintp_start_e,jcd)
  ieb_le      = p_lp%edges%end_blk(min_rledge_int,jcd)
  isb_lc      = p_lp%cells%start_blk(grf_bdyintp_start_c,jcd)
  ieb_lc      = p_lp%cells%end_blk(min_rlcell_int,jcd)

  is_e = idx_1d(p_p%edges%start_idx(grf_bdyintp_start_e,jcd), &
                p_p%edges%start_blk(grf_bdyintp_start_e,jcd))
  ie_e = idx_1d(p_p%edges%end_idx(min_rledge_int,jcd), &
                p_p%edges%end_blk(min_rledge_int,jcd))

  is_c = idx_1d(p_p%cells%start_idx(grf_bdyintp_start_c,jcd), &
                p_p%cells%start_blk(grf_bdyintp_start_c,jcd))
  ie_c = idx_1d(p_p%cells%end_idx(min_rlcell_int,jcd), &
                p_p%cells%end_blk(min_rlcell_int,jcd))

  ! Set up communication patterns for transferring the data to local parents.
  ! Since these communication patterns are not used elsewhere, they are
  ! stored locally and deleted at the end of the routine


  ALLOCATE(owner(p_p%n_patch_edges))
  owner = -1
  DO j = is_e, ie_e
    owner(j) = p_lp%edges%owner_g(p_p%edges%glb_index(j))
  ENDDO
  CALL setup_comm_pattern(p_p%n_patch_edges, owner, p_p%edges%glb_index, &
    & p_lp%edges%loc_index, comm_pat_loc_to_glb_e)
  DEALLOCATE(owner)

  ALLOCATE(owner(p_p%n_patch_cells))
  owner = -1
  DO j = is_c, ie_c
    owner(j) = p_lp%cells%owner_g(p_p%cells%glb_index(j))
  ENDDO
  CALL setup_comm_pattern(p_p%n_patch_cells, owner, p_p%cells%glb_index, &
    & p_lp%cells%loc_index, comm_pat_loc_to_glb_c)
  DEALLOCATE(owner)

  ALLOCATE(z_tmp_s(nproma,4,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,4,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp
  z_tmp_s(:,1,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_1a(:,:),wp)
  z_tmp_s(:,2,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_1b(:,:),wp)
  z_tmp_s(:,3,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_2a(:,:),wp)
  z_tmp_s(:,4,isb_le:ieb_le) = REAL(p_lgrf%p_dom(jcd)%grf_vec_stencil_2b(:,:),wp)

  CALL exchange_data(comm_pat_loc_to_glb_e, RECV=z_tmp_r, SEND=z_tmp_s)

  p_grf%p_dom(jcd)%grf_vec_stencil_1a(:,:) = INT(z_tmp_r(:,1,isb_e:ieb_e))
  p_grf%p_dom(jcd)%grf_vec_stencil_1b(:,:) = INT(z_tmp_r(:,2,isb_e:ieb_e))
  p_grf%p_dom(jcd)%grf_vec_stencil_2a(:,:) = INT(z_tmp_r(:,3,isb_e:ieb_e))
  p_grf%p_dom(jcd)%grf_vec_stencil_2b(:,:) = INT(z_tmp_r(:,4,isb_e:ieb_e))
  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,2*grf_vec_dim_1+2*grf_vec_dim_2,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,2*grf_vec_dim_1+2*grf_vec_dim_2,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  n = 0
  DO k = 1, grf_vec_dim_1
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_1a(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_1a(:,k,isb_le:ieb_le)),wp)
  ENDDO
  DO k = 1, grf_vec_dim_1
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_1b(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_1b(:,k,isb_le:ieb_le)),wp)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_2a(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_2a(:,k,isb_le:ieb_le)),wp)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
!CDIR IEXPAND
    z_tmp_s(:,n,isb_le:ieb_le) =REAL(glb_idx_1d_e(p_lp,p_lgrf%p_dom(jcd)%grf_vec_ind_2b(:,k,isb_le:ieb_le), &
                                                       p_lgrf%p_dom(jcd)%grf_vec_blk_2b(:,k,isb_le:ieb_le)),wp)
  ENDDO

  CALL exchange_data(comm_pat_loc_to_glb_e, RECV=z_tmp_r, SEND=z_tmp_s)

  n = 0
  DO k = 1, grf_vec_dim_1
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_ind_1a(:,k,isb_e:ieb_e) = loc_idx_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
    p_grf%p_dom(jcd)%grf_vec_blk_1a(:,k,isb_e:ieb_e) = loc_blk_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
  ENDDO
  DO k = 1, grf_vec_dim_1
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_ind_1b(:,k,isb_e:ieb_e) = loc_idx_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
    p_grf%p_dom(jcd)%grf_vec_blk_1b(:,k,isb_e:ieb_e) = loc_blk_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_ind_2a(:,k,isb_e:ieb_e) = loc_idx_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
    p_grf%p_dom(jcd)%grf_vec_blk_2a(:,k,isb_e:ieb_e) = loc_blk_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_ind_2b(:,k,isb_e:ieb_e) = loc_idx_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
    p_grf%p_dom(jcd)%grf_vec_blk_2b(:,k,isb_e:ieb_e) = loc_blk_no_e(p_p, int(z_tmp_r(:,n,isb_e:ieb_e)))
  ENDDO

  n = 0
  DO k = 1, grf_vec_dim_1
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_1a(k,:,:)
  ENDDO
  DO k = 1, grf_vec_dim_1
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_1b(k,:,:)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_2a(k,:,:)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    z_tmp_s(:,n,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_vec_coeff_2b(k,:,:)
  ENDDO

  CALL exchange_data(comm_pat_loc_to_glb_e, RECV=z_tmp_r, SEND=z_tmp_s)

  n = 0
  DO k = 1, grf_vec_dim_1
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_coeff_1a(k,:,:) = z_tmp_r(:,n,isb_e:ieb_e)
  ENDDO
  DO k = 1, grf_vec_dim_1
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_coeff_1b(k,:,:) = z_tmp_r(:,n,isb_e:ieb_e)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_coeff_2a(k,:,:) = z_tmp_r(:,n,isb_e:ieb_e)
  ENDDO
  DO k = 1, grf_vec_dim_2
    n = n+1
    p_grf%p_dom(jcd)%grf_vec_coeff_2b(k,:,:) = z_tmp_r(:,n,isb_e:ieb_e)
  ENDDO

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,2,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,2,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  z_tmp_s(:,:,isb_le:ieb_le) = p_lgrf%p_dom(jcd)%grf_dist_pe2ce(:,:,:)

  CALL exchange_data(comm_pat_loc_to_glb_e, RECV=z_tmp_r, SEND=z_tmp_s)

  p_grf%p_dom(jcd)%grf_dist_pe2ce(:,:,:) = z_tmp_r(:,:,isb_e:ieb_e)

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,8,p_lp%nblks_c))
  ALLOCATE(z_tmp_r(nproma,8,p_p%nblks_c))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  z_tmp_s(:,1,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,1,1,:)
  z_tmp_s(:,2,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,2,1,:)
  z_tmp_s(:,3,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,3,1,:)
  z_tmp_s(:,4,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,4,1,:)
  z_tmp_s(:,5,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,1,2,:)
  z_tmp_s(:,6,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,2,2,:)
  z_tmp_s(:,7,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,3,2,:)
  z_tmp_s(:,8,isb_lc:ieb_lc) = p_lgrf%p_dom(jcd)%grf_dist_pc2cc(:,4,2,:)

  CALL exchange_data(comm_pat_loc_to_glb_c, RECV=z_tmp_r, SEND=z_tmp_s)

  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,1,1,:) = z_tmp_r(:,1,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,2,1,:) = z_tmp_r(:,2,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,3,1,:) = z_tmp_r(:,3,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,4,1,:) = z_tmp_r(:,4,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,1,2,:) = z_tmp_r(:,5,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,2,2,:) = z_tmp_r(:,6,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,3,2,:) = z_tmp_r(:,7,isb_c:ieb_c)
  p_grf%p_dom(jcd)%grf_dist_pc2cc(:,4,2,:) = z_tmp_r(:,8,isb_c:ieb_c)

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_c)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_e)

  ! Now create communication patterns for the complete patch

  ALLOCATE(owner(p_p%n_patch_edges))
  DO j = 1, p_p%n_patch_edges
    owner(j) = p_lp%edges%owner_g(p_p%edges%glb_index(j))
  ENDDO
  CALL setup_comm_pattern(p_p%n_patch_edges, owner, p_p%edges%glb_index, &
    & p_lp%edges%loc_index, comm_pat_loc_to_glb_e)
  DEALLOCATE(owner)

  ALLOCATE(owner(p_p%n_patch_cells))
  DO j = 1, p_p%n_patch_cells
    owner(j) = p_lp%cells%owner_g(p_p%cells%glb_index(j))
  ENDDO
  CALL setup_comm_pattern(p_p%n_patch_cells, owner, p_p%cells%glb_index, &
    & p_lp%cells%loc_index, comm_pat_loc_to_glb_c)
  DEALLOCATE(owner)

  ALLOCATE(z_tmp_s(nproma,8,p_lp%nblks_c))
  ALLOCATE(z_tmp_r(nproma,8,p_p%nblks_c))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  n = 0
  DO k = 1, 4
    n = n+1
    z_tmp_s(:,n,:) = p_lgrf%fbk_wgt_c(:,:,k)
  ENDDO
  DO k = 1, 4
    n = n+1
    z_tmp_s(:,n,:) = p_lgrf%fbk_wgt_ct(:,:,k)
  ENDDO

  CALL exchange_data(comm_pat_loc_to_glb_c, RECV=z_tmp_r, SEND=z_tmp_s)

  n = 0
  DO k = 1, 4
    n = n+1
    p_grf%fbk_wgt_c(:,:,k) = z_tmp_r(:,n,:)
  ENDDO
  DO k = 1, 4
    n = n+1
    p_grf%fbk_wgt_ct(:,:,k) = z_tmp_r(:,n,:)
  ENDDO

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

  ALLOCATE(z_tmp_s(nproma,8,p_lp%nblks_e))
  ALLOCATE(z_tmp_r(nproma,8,p_p%nblks_e))
  z_tmp_s = 0._wp
  z_tmp_r = 0._wp

  n = 0
  DO k = 1, 6
    n = n+1
    z_tmp_s(:,n,:) = p_lgrf%fbk_wgt_e(:,:,k)
  ENDDO

  CALL exchange_data(comm_pat_loc_to_glb_e, RECV=z_tmp_r, SEND=z_tmp_s)

  n = 0
  DO k = 1, 6
    n = n+1
    p_grf%fbk_wgt_e(:,:,k) = z_tmp_r(:,n,:)
  ENDDO

  DEALLOCATE(z_tmp_s)
  DEALLOCATE(z_tmp_r)

!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_c)
!CDIR NOIEXPAND
  CALL delete_comm_pattern(comm_pat_loc_to_glb_e)

  ! Transfers from global parent to local parent

  p_lgrf%fbk_dom_area(:) = p_grf%fbk_dom_area(:)

END SUBROUTINE transfer_grf_state

!-------------------------------------------------------------------------
!
!
!>
!! Deallocation of fields needed for grid refinement for a single state.
!!
!!
!! @par Revision History
!! Split off from destruct_2d_interpol_state, Rainer Johanni (2010-10-29)
!!
SUBROUTINE deallocate_grf_state( ptr_patch, ptr_grf )
!
TYPE(t_patch), INTENT(IN) :: ptr_patch

TYPE(t_gridref_state), TARGET, INTENT(INOUT) :: ptr_grf

TYPE(t_gridref_single_state), POINTER :: ptr_int    => NULL()

INTEGER :: jcd, i_nchdom
INTEGER :: ist

!-----------------------------------------------------------------------

  i_nchdom = MAX(1,ptr_patch%n_childdom)

  DO jcd = 1, i_nchdom

    ptr_int => ptr_grf%p_dom(jcd)

    ! index arrays

    DEALLOCATE (ptr_int%grf_vec_ind_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_ind_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_ind_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_ind_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_2b failed')
    ENDIF

    ! block index arrays

    DEALLOCATE (ptr_int%grf_vec_blk_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_blk_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_blk_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_blk_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_ind_2b failed')
    ENDIF

    ! arrays holding the actual number of cells/edges available for reconstruction

    DEALLOCATE (ptr_int%grf_vec_stencil_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_stencil_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_stencil_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_stencil_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_stencil_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_stencil_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_stencil_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_stencil_2b failed')
    ENDIF

    ! arrays holding the interpolation coefficients for reconstruction related to grid refinement

    DEALLOCATE (ptr_int%grf_vec_coeff_1a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_coeff_1a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_coeff_1b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_coeff_1b failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_coeff_2a, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_coeff_2a failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_vec_coeff_2b, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_vec_coeff_2b failed')
    ENDIF

    DEALLOCATE (ptr_int%grf_dist_pc2cc, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_dist_pc2cc failed')
    ENDIF
    DEALLOCATE (ptr_int%grf_dist_pe2ce, STAT=ist )
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for grf_dist_pe2ce failed')
    ENDIF

  ENDDO

  DEALLOCATE (ptr_grf%fbk_wgt_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for fbk_wgt_c failed')
  ENDIF
  DEALLOCATE (ptr_grf%fbk_wgt_ct, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for fbk_wgt_ct failed')
  ENDIF
  DEALLOCATE (ptr_grf%fbk_wgt_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for fbk_wgt_e failed')
  ENDIF

  DEALLOCATE (ptr_grf%fbk_dom_area, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for fbk_dom_area failed')
  ENDIF

  DEALLOCATE (ptr_grf%p_dom, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for p_dom failed')
  ENDIF

END SUBROUTINE deallocate_grf_state

!-------------------------------------------------------------------------
!
!
!>
!!               Deallocation of fields needed for grid refinement.
!!
!!
!! @par Revision History
!! Created  by  Guenther Zaengl, DWD (2009-02-11).
!!
SUBROUTINE destruct_2d_gridref_state( ptr_patch, ptr_grf_state )
!
TYPE(t_patch), INTENT(IN) :: ptr_patch(n_dom_start:)
TYPE(t_gridref_state), INTENT(INOUT) :: ptr_grf_state(n_dom_start:)

INTEGER :: jg

!-----------------------------------------------------------------------

CALL message('mo_grf_intp_state:destruct_2d_gridref_state', &
  & 'start to destruct grf state')

!
! deallocate gridref_state
!
DO jg = n_dom_start, n_dom
  CALL deallocate_grf_state(ptr_patch(jg), ptr_grf_state(jg))
END DO

CALL message ('mo_grf_intp_state:destruct_2d_gridref_state', &
  & 'destruction of grf state finished')

END SUBROUTINE destruct_2d_gridref_state


!-------------------------------------------------------------------------
END MODULE mo_grf_intp_state
