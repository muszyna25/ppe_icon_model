
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
USE mo_io_units,            ONLY: nnml, nnml_output
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, min_rlcell_int, min_rledge_int
USE mo_model_domain,        ONLY: t_patch
USE mo_model_domain_import, ONLY: n_dom, n_dom_start
USE mo_namelist,            ONLY: position_nml, POSITIONED
USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e
USE mo_run_nml,             ONLY: nproma, lshallow_water
USE mo_mpi,                 ONLY: p_pe, p_io

USE mo_grf_intp_data_strc
USE mo_grf_intp_coeffs

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC ::setup_gridref, construct_2d_gridref_state, destruct_2d_gridref_state
PUBLIC ::allocate_grf_state, deallocate_grf_state

NAMELIST/gridref_ctl/  rbf_vec_kern_grf_e, rbf_scale_grf_e,             &
  &                    grf_velfbk, grf_scalfbk, grf_tracfbk,            &
  &                    grf_idw_exp_e12, grf_idw_exp_e34,                &
  &                    grf_intmethod_c, grf_intmethod_e,                &
  &                    grf_intmethod_ct, denom_diffu_v, denom_diffu_t

CONTAINS

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!

!>
!!   Set up the configuration for grid refinement.
!!
!!
!! @par Revision History
!!  Created by Guenther Zaengl, MPI-M (2009-02-09).
!!
SUBROUTINE setup_gridref
!
! !local variables
INTEGER :: i_status

CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine = 'mo_grf_interpolation: setup_gridref'

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! set default values
!-----------------------------------------------------------------------

CALL message(TRIM(routine),'default settings')

! Switch for interpolation method used for cell-based dynamical
grf_intmethod_c   = 2         ! 1: copying, 2: gradient-based interpolation
! Switch for interpolation method used for tracer variables
grf_intmethod_ct  = 2         ! 1: copying, 2: gradient-based interpolation
! Currently, grf_intmethod_c is used for temperature only; other variables are copied
grf_intmethod_e   = 4         ! 1: IDW, 2: RBF, 3: IDW/gradient-based, 4: RBF/gradient-based
! Switch for velocity feedback method.
grf_velfbk      = 1         ! 1: average over child edges 1 and 2
                            ! 2: 2nd-order method using RBF reconstruction to child vertices
! Switch for feedback method for scalar dynamical variables
grf_scalfbk     = 2         ! 1: area-weighted averaging
                            ! 2: bilinear interpolation
! Switch for feedback method for passive tracer variables
grf_tracfbk     = 2         ! 1: area-weighted averaging
                            ! 2: bilinear interpolation

! Exponents for IDW interpolation function
grf_idw_exp_e12  = 1.2_wp   ! child edges 1 and 2
grf_idw_exp_e34  = 1.7_wp   ! child edges 3 and 4

! RBF kernels for grid refinement interpolation
rbf_vec_kern_grf_e = 1        ! 1: Gaussian, 2: 1/(1+r**2), 3: inverse multiquadric

! zero whole arrays
rbf_scale_grf_e(:) = 0.0_wp

! Initialize namelist fields for scaling factors (dimension 1:n_dom); used part only
rbf_scale_grf_e(1:n_dom) = 0.5_wp  ! default setting for vector grf interpolation

! Denominator for temperature boundary diffusion
denom_diffu_t = 135._wp

! Denominator for velocity boundary diffusion
IF (lshallow_water) THEN
  denom_diffu_v = 250._wp
ELSE
  denom_diffu_v = 200._wp
ENDIF


!-----------------------------------------------------------------------
! read namelist
!-----------------------------------------------------------------------

CALL position_nml ('gridref_ctl', status=i_status)
SELECT CASE (i_status)
CASE (POSITIONED)
  READ (nnml, gridref_ctl)
END SELECT

! write the contents of the namelist to an ASCII file

   IF(p_pe == p_io) WRITE(nnml_output,nml=gridref_ctl)


END SUBROUTINE setup_gridref

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

  ! "parent-child-index", i.e. child index number a grid point has at its parent point
  ALLOCATE (ptr_grf%pc_idx_c(nproma,nblks_c), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for pc_idx_c failed')
  ENDIF
  ALLOCATE (ptr_grf%pc_idx_e(nproma,nblks_e), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for pc_idx_e failed')
  ENDIF
  ALLOCATE (ptr_grf%fbk_dom_area(i_nchdom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for fbk_dom_area failed')
  ENDIF
  ALLOCATE (ptr_grf%fbk_dom_volume(nlev,i_nchdom), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'allocation for fbk_dom_volume failed')
  ENDIF

  ptr_grf%fbk_wgt_c     = 0._wp
  ptr_grf%fbk_wgt_ct    = 0._wp
  ptr_grf%fbk_wgt_e     = 0._wp
  ptr_grf%pc_idx_c      = 0
  ptr_grf%pc_idx_e      = 0
  ptr_grf%fbk_dom_area  = 0._wp
  ptr_grf%fbk_dom_volume = 0._wp

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

  DEALLOCATE (ptr_grf%pc_idx_c, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for pc_idx_c failed')
  ENDIF
  DEALLOCATE (ptr_grf%pc_idx_e, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for pc_idx_e failed')
  ENDIF
  DEALLOCATE (ptr_grf%fbk_dom_area, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for fbk_dom_area failed')
  ENDIF
  DEALLOCATE (ptr_grf%fbk_dom_volume, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_grf_interpolation:construct_grf_state',                       &
      &             'deallocation for fbk_dom_volume failed')
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
