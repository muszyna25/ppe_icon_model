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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_grf_bdyintp
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
USE mo_model_domain,        ONLY: t_patch
USE mo_intp_data_strc,      ONLY: t_int_state
USE mo_impl_constants_grf,  ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,   &
                                  grf_bdyintp_end_c, grf_bdyintp_end_e, grf_bdywidth_e
USE mo_parallel_config,     ONLY: nproma, p_test_run
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_mpi,                 ONLY: my_process_is_mpi_parallel, my_process_is_mpi_seq
USE mo_communication,       ONLY: exchange_data, exchange_data_grf

USE mo_grf_intp_data_strc


IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: interpol_vec_grf, interpol2_vec_grf, interpol_scal_grf, interpol_scal2d_grf

CONTAINS
!-------------------------------------------------------------------------
!
!
!>
!! Performs interpolation from parent edges to child edges using inverse.
!!
!! Performs interpolation from parent edges to child edges using inverse
!! distance weighting (IDW) or RBF. The output field (p_vn_out) contains
!! the normal wind components along the lateral boundaries of nested
!! model domains.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2008-06-05)
!! Vector optimization (2009-03-20)
!!
SUBROUTINE interpol_vec_grf (ptr_pp, ptr_pc, ptr_grf, i_chidx, p_vn_in, p_vn_out)
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc

! input normal components of vectors at edge midpoints
REAL(wp), INTENT(IN) ::  &
  p_vn_in(:,:,:) ! dim: (nproma,nlev,nblks_e)

! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx


! reconstructed edge-normal wind component
REAL(wp),INTENT(INOUT) :: p_vn_out(:,:,:) ! dim: (nproma,nlev_c,nblks_e)


INTEGER :: jb, jk, je                ! loop indices
INTEGER :: js                        ! shift parameter
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

INTEGER :: nsendtot, nrecvtot        ! for MPI communication call

REAL(wp) :: vn_aux(nproma,ptr_pc%nlev,ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx):&
                   MAX(ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx),               &
                       ptr_pp%edges%end_blk(grf_bdyintp_end_e,i_chidx)),4)

! Pointers to coefficient fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx_1a, iblk_1a, iidx_1b, iblk_1b, &
                                         iidx_2a, iblk_2a, iidx_2b, iblk_2b, &
                                         icheidx, icheblk
INTEGER :: nlev_c       !< number of vertical full levels (child domain)
!-----------------------------------------------------------------------

! Set pointers to coefficient fields
iidx_1a => ptr_grf%grf_vec_ind_1a
iblk_1a => ptr_grf%grf_vec_blk_1a
iidx_1b => ptr_grf%grf_vec_ind_1b
iblk_1b => ptr_grf%grf_vec_blk_1b

iidx_2a => ptr_grf%grf_vec_ind_2a
iblk_2a => ptr_grf%grf_vec_blk_2a
iidx_2b => ptr_grf%grf_vec_ind_2b
iblk_2b => ptr_grf%grf_vec_blk_2b

icheidx => ptr_pp%edges%child_idx
icheblk => ptr_pp%edges%child_blk

! Start and end blocks for which vector interpolation is needed
i_startblk = ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx)
i_endblk   = ptr_pp%edges%end_blk(grf_bdyintp_end_e,i_chidx)

! number of vertical full levels (child domain)
nlev_c = ptr_pc%nlev
! difference between upper boundary of parent domain and 
! upper boundary of child domain (in terms of vertical levels) 
js = ptr_pc%nshift

!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx)
DO jb =  i_startblk, i_endblk

  CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_bdyintp_start_e, grf_bdyintp_end_e, i_chidx)

! Note: OMP parallelization is done over jk because for reasonable (=efficient)
! choices of nproma, the boundary interpolation zone extends over no more
! than 2 blocks.

#ifndef _OPENMP

! a) child edge 1
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
      vn_aux(je,jk,jb,1) = ptr_grf%grf_vec_coeff_1a(1,je,jb) * &
        p_vn_in(iidx_1a(je,1,jb),jk+js,iblk_1a(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(2,je,jb) *                    &
        p_vn_in(iidx_1a(je,2,jb),jk+js,iblk_1a(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(3,je,jb) *                    &
        p_vn_in(iidx_1a(je,3,jb),jk+js,iblk_1a(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(4,je,jb) *                    &
        p_vn_in(iidx_1a(je,4,jb),jk+js,iblk_1a(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(5,je,jb) *                    &
        p_vn_in(iidx_1a(je,5,jb),jk+js,iblk_1a(je,5,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(6,je,jb) *                    &
        p_vn_in(iidx_1a(je,6,jb),jk+js,iblk_1a(je,6,jb))
    ENDDO
  ENDDO

! b) child edge 2
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
      vn_aux(je,jk,jb,2) = ptr_grf%grf_vec_coeff_1b(1,je,jb) * &
        p_vn_in(iidx_1b(je,1,jb),jk+js,iblk_1b(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(2,je,jb) *                    &
        p_vn_in(iidx_1b(je,2,jb),jk+js,iblk_1b(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(3,je,jb) *                    &
        p_vn_in(iidx_1b(je,3,jb),jk+js,iblk_1b(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(4,je,jb) *                    &
        p_vn_in(iidx_1b(je,4,jb),jk+js,iblk_1b(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(5,je,jb) *                    &
        p_vn_in(iidx_1b(je,5,jb),jk+js,iblk_1b(je,5,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(6,je,jb) *                    &
        p_vn_in(iidx_1b(je,6,jb),jk+js,iblk_1b(je,6,jb))
    ENDDO
  ENDDO

! c) child edge 3
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
      vn_aux(je,jk,jb,3) = ptr_grf%grf_vec_coeff_2a(1,je,jb) * &
        p_vn_in(iidx_2a(je,1,jb),jk+js,iblk_2a(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(2,je,jb) *                    &
        p_vn_in(iidx_2a(je,2,jb),jk+js,iblk_2a(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(3,je,jb) *                    &
        p_vn_in(iidx_2a(je,3,jb),jk+js,iblk_2a(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(4,je,jb) *                    &
        p_vn_in(iidx_2a(je,4,jb),jk+js,iblk_2a(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(5,je,jb) *                    &
        p_vn_in(iidx_2a(je,5,jb),jk+js,iblk_2a(je,5,jb))
    ENDDO
  ENDDO

! d) child edge 4
!CDIR UNROLL=6
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx
      IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
        vn_aux(je,jk,jb,4) = ptr_grf%grf_vec_coeff_2b(1,je,jb) * &
          p_vn_in(iidx_2b(je,1,jb),jk+js,iblk_2b(je,1,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(2,je,jb) *                    &
          p_vn_in(iidx_2b(je,2,jb),jk+js,iblk_2b(je,2,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(3,je,jb) *                    &
          p_vn_in(iidx_2b(je,3,jb),jk+js,iblk_2b(je,3,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(4,je,jb) *                    &
          p_vn_in(iidx_2b(je,4,jb),jk+js,iblk_2b(je,4,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(5,je,jb) *                    &
          p_vn_in(iidx_2b(je,5,jb),jk+js,iblk_2b(je,5,jb))
      ENDIF
    ENDDO
  ENDDO
#else

!$OMP DO PRIVATE (jk,je) ICON_OMP_DEFAULT_SCHEDULE
  DO jk = 1, nlev_c
    DO je = i_startidx, i_endidx

      ! a) child edge 1
      vn_aux(je,jk,jb,1) = ptr_grf%grf_vec_coeff_1a(1,je,jb) * &
        p_vn_in(iidx_1a(je,1,jb),jk+js,iblk_1a(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(2,je,jb) *                    &
        p_vn_in(iidx_1a(je,2,jb),jk+js,iblk_1a(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(3,je,jb) *                    &
        p_vn_in(iidx_1a(je,3,jb),jk+js,iblk_1a(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(4,je,jb) *                    &
        p_vn_in(iidx_1a(je,4,jb),jk+js,iblk_1a(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(5,je,jb) *                    &
        p_vn_in(iidx_1a(je,5,jb),jk+js,iblk_1a(je,5,jb)) +     &
        ptr_grf%grf_vec_coeff_1a(6,je,jb) *                    &
        p_vn_in(iidx_1a(je,6,jb),jk+js,iblk_1a(je,6,jb))

      ! b) child edge 2
      vn_aux(je,jk,jb,2) = ptr_grf%grf_vec_coeff_1b(1,je,jb) * &
        p_vn_in(iidx_1b(je,1,jb),jk+js,iblk_1b(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(2,je,jb) *                    &
        p_vn_in(iidx_1b(je,2,jb),jk+js,iblk_1b(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(3,je,jb) *                    &
        p_vn_in(iidx_1b(je,3,jb),jk+js,iblk_1b(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(4,je,jb) *                    &
        p_vn_in(iidx_1b(je,4,jb),jk+js,iblk_1b(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(5,je,jb) *                    &
        p_vn_in(iidx_1b(je,5,jb),jk+js,iblk_1b(je,5,jb)) +     &
        ptr_grf%grf_vec_coeff_1b(6,je,jb) *                    &
        p_vn_in(iidx_1b(je,6,jb),jk+js,iblk_1b(je,6,jb))

      ! c) child edge 3
      vn_aux(je,jk,jb,3) = ptr_grf%grf_vec_coeff_2a(1,je,jb) * &
        p_vn_in(iidx_2a(je,1,jb),jk+js,iblk_2a(je,1,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(2,je,jb) *                    &
        p_vn_in(iidx_2a(je,2,jb),jk+js,iblk_2a(je,2,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(3,je,jb) *                    &
        p_vn_in(iidx_2a(je,3,jb),jk+js,iblk_2a(je,3,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(4,je,jb) *                    &
        p_vn_in(iidx_2a(je,4,jb),jk+js,iblk_2a(je,4,jb)) +     &
        ptr_grf%grf_vec_coeff_2a(5,je,jb) *                    &
        p_vn_in(iidx_2a(je,5,jb),jk+js,iblk_2a(je,5,jb))

      ! d) child edge 4
      IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
        vn_aux(je,jk,jb,4) = ptr_grf%grf_vec_coeff_2b(1,je,jb) * &
          p_vn_in(iidx_2b(je,1,jb),jk+js,iblk_2b(je,1,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(2,je,jb) *                    &
          p_vn_in(iidx_2b(je,2,jb),jk+js,iblk_2b(je,2,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(3,je,jb) *                    &
          p_vn_in(iidx_2b(je,3,jb),jk+js,iblk_2b(je,3,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(4,je,jb) *                    &
          p_vn_in(iidx_2b(je,4,jb),jk+js,iblk_2b(je,4,jb)) +     &
          ptr_grf%grf_vec_coeff_2b(5,je,jb) *                    &
          p_vn_in(iidx_2b(je,5,jb),jk+js,iblk_2b(je,5,jb))
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
#endif

ENDDO ! blocks


! Store results in p_vn_out

IF (my_process_is_mpi_seq()) THEN

  DO jb =  i_startblk, i_endblk

    CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
      i_startidx, i_endidx, grf_bdyintp_start_e, grf_bdyintp_end_e, i_chidx)

!$OMP DO PRIVATE (jk,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jk = 1, nlev_c
!CDIR NODEP,VOVERTAKE,VOB
      DO je = i_startidx, i_endidx

        p_vn_out(icheidx(je,jb,1),jk,icheblk(je,jb,1))   = vn_aux(je,jk,jb,1)
        p_vn_out(icheidx(je,jb,2),jk,icheblk(je,jb,2))   = vn_aux(je,jk,jb,2)

        ! Avoid overwriting prognostic grid points in case of applying this to a prognostic variable
        IF (ptr_pc%edges%refin_ctrl(icheidx(je,jb,3),icheblk(je,jb,3)) <= grf_bdywidth_e) &
          p_vn_out(icheidx(je,jb,3),jk,icheblk(je,jb,3))   = vn_aux(je,jk,jb,3)

        IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1 .AND.                                    &
            ptr_pc%edges%refin_ctrl(icheidx(je,jb,4),icheblk(je,jb,4)) <= grf_bdywidth_e) &
          p_vn_out(icheidx(je,jb,4),jk,icheblk(je,jb,4)) = vn_aux(je,jk,jb,4)

      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

ENDIF
!$OMP END PARALLEL

IF (my_process_is_mpi_parallel()) THEN

  nsendtot = SUM(ptr_pc%comm_pat_interpol_vec_grf(1:4)%n_send)
  nrecvtot = SUM(ptr_pc%comm_pat_interpol_vec_grf(1:4)%n_recv)
  CALL exchange_data_grf(ptr_pc%comm_pat_interpol_vec_grf,1,nlev_c,nsendtot,nrecvtot, &
                         RECV1=p_vn_out,SEND1=vn_aux,SEND_LBOUND3=LBOUND(vn_aux,3))

ENDIF

END SUBROUTINE interpol_vec_grf

!
!>
!! Performs gradient-based interpolation from parent edges to child edges 1 and 2
!! and RBF/IDW-based interpolation to child edges 3 and 4
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2010-03-12)
!!
SUBROUTINE interpol2_vec_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, i_chidx, nfields,&
                              f3din1, f3dout1, f3din2, f3dout2, f3din3, f3dout3, &
                              f3din4, f3dout4)

TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc

! input normal components of vectors at edge midpoints; dim: (nproma,nlev,nblks_e)
REAL(wp), INTENT(IN), DIMENSION(:,:,:), OPTIONAL, TARGET :: f3din1, f3din2, f3din3, f3din4

! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf
TYPE(t_int_state),              TARGET, INTENT(IN)    ::  ptr_int

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! number of fields provided on input (needed for aux fields and pointer allocation)
INTEGER, INTENT(IN) :: nfields

! reconstructed edge-normal vector component; dim: (nproma,nlev_c,nblks_e)
REAL(wp), INTENT(INOUT), DIMENSION(:,:,:), OPTIONAL, TARGET :: f3dout1, f3dout2, f3dout3, f3dout4


INTEGER :: jb, jk, je, jv, jn        ! loop indices
INTEGER :: js                        ! shift parameter
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

INTEGER :: nsendtot, nrecvtot, nlevtot ! for MPI communication call

REAL(wp) :: dvn_tang(nproma)

! Auxiliary fields
REAL(wp) :: vn_aux(nproma,ptr_pc%nlev,ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx):&
                   MAX(ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx),               &
                       ptr_pp%edges%end_blk(grf_bdyintp_end_e,i_chidx)),4,nfields)

REAL(wp), DIMENSION(nproma,ptr_pc%nlev,ptr_pp%verts%start_blk(grf_bdyintp_start_c,i_chidx):&
                    MAX(ptr_pp%verts%start_blk(grf_bdyintp_start_c,i_chidx),               &
                        ptr_pp%verts%end_blk(grf_bdyintp_end_c-1,i_chidx))) :: u_vert,v_vert

! Pointers to index fields
INTEGER,  DIMENSION(:,:,:), POINTER :: iidx_2a, iblk_2a, iidx_2b, iblk_2b, ividx, ivblk, &
                                       irvidx, irvblk, icheidx, icheblk

REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_rvcoeff
REAL(wp), DIMENSION(:,:,:), POINTER   :: ptr_dist

! Allocatable pointer to input and output fields
TYPE t_fieldptr
  REAL(wp), POINTER :: fld(:,:,:)
END TYPE t_fieldptr
TYPE(t_fieldptr) :: p_in(nfields), p_out(nfields)

INTEGER :: nlev_c       !< number of vertical full levels (child domain)
!-----------------------------------------------------------------------
IF (p_test_run) vn_aux=0._wp

IF (PRESENT(f3din1)) THEN
  p_in(1)%fld  => f3din1
  p_out(1)%fld => f3dout1
ENDIF
IF (PRESENT(f3din2)) THEN
  p_in(2)%fld  => f3din2
  p_out(2)%fld => f3dout2
ENDIF
IF (PRESENT(f3din3)) THEN
  p_in(3)%fld  => f3din3
  p_out(3)%fld => f3dout3
ENDIF
IF (PRESENT(f3din4)) THEN
  p_in(4)%fld  => f3din4
  p_out(4)%fld => f3dout4
ENDIF

! Set pointers to index and coefficient fields

iidx_2a => ptr_grf%grf_vec_ind_2a
iblk_2a => ptr_grf%grf_vec_blk_2a
iidx_2b => ptr_grf%grf_vec_ind_2b
iblk_2b => ptr_grf%grf_vec_blk_2b

! vertices of edges
ividx => ptr_pp%edges%vertex_idx
ivblk => ptr_pp%edges%vertex_blk

! for RBF interpolation to vertices
irvidx => ptr_int%rbf_vec_idx_v
irvblk => ptr_int%rbf_vec_blk_v
ptr_rvcoeff => ptr_int%rbf_vec_coeff_v

! normalized distances from parent edges to child edges
ptr_dist  => ptr_grf%grf_dist_pe2ce

! child edge indices
icheidx => ptr_pp%edges%child_idx
icheblk => ptr_pp%edges%child_blk

! number of vertical full levels (child domain)
nlev_c = ptr_pc%nlev
! difference between upper boundary of parent domain and 
! upper boundary of child domain (in terms of vertical levels
js = ptr_pc%nshift

!$OMP PARALLEL PRIVATE (jb,jn,i_startblk,i_endblk,i_startidx,i_endidx)

DO jn = 1, nfields

  ! Start and end blocks for which vector reconstruction at vertices is needed
  ! Note: the use of grf_bdyintp_start_c (=-1) and grf_bdyintp_end_c-1 (=-3) is intended
  i_startblk = ptr_pp%verts%start_blk(grf_bdyintp_start_c,i_chidx)
  i_endblk   = ptr_pp%verts%end_blk(grf_bdyintp_end_c-1,i_chidx)

!$OMP DO PRIVATE(jk,jv) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_pp, jb, i_startblk, i_endblk, &
     i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c-1, i_chidx)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = 1, nlev_c
#else
!CDIR UNROLL=6
    DO jk = 1, nlev_c
      DO jv = i_startidx, i_endidx
#endif

        u_vert(jv,jk,jb) =  &
          ptr_rvcoeff(1,1,jv,jb)*p_in(jn)%fld(irvidx(1,jv,jb),jk+js,irvblk(1,jv,jb)) + &
          ptr_rvcoeff(2,1,jv,jb)*p_in(jn)%fld(irvidx(2,jv,jb),jk+js,irvblk(2,jv,jb)) + &
          ptr_rvcoeff(3,1,jv,jb)*p_in(jn)%fld(irvidx(3,jv,jb),jk+js,irvblk(3,jv,jb)) + &
          ptr_rvcoeff(4,1,jv,jb)*p_in(jn)%fld(irvidx(4,jv,jb),jk+js,irvblk(4,jv,jb)) + &
          ptr_rvcoeff(5,1,jv,jb)*p_in(jn)%fld(irvidx(5,jv,jb),jk+js,irvblk(5,jv,jb)) + &
          ptr_rvcoeff(6,1,jv,jb)*p_in(jn)%fld(irvidx(6,jv,jb),jk+js,irvblk(6,jv,jb))
        v_vert(jv,jk,jb) =  &
          ptr_rvcoeff(1,2,jv,jb)*p_in(jn)%fld(irvidx(1,jv,jb),jk+js,irvblk(1,jv,jb)) + &
          ptr_rvcoeff(2,2,jv,jb)*p_in(jn)%fld(irvidx(2,jv,jb),jk+js,irvblk(2,jv,jb)) + &
          ptr_rvcoeff(3,2,jv,jb)*p_in(jn)%fld(irvidx(3,jv,jb),jk+js,irvblk(3,jv,jb)) + &
          ptr_rvcoeff(4,2,jv,jb)*p_in(jn)%fld(irvidx(4,jv,jb),jk+js,irvblk(4,jv,jb)) + &
          ptr_rvcoeff(5,2,jv,jb)*p_in(jn)%fld(irvidx(5,jv,jb),jk+js,irvblk(5,jv,jb)) + &
          ptr_rvcoeff(6,2,jv,jb)*p_in(jn)%fld(irvidx(6,jv,jb),jk+js,irvblk(6,jv,jb))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

  ! Start and end blocks for which interpolation to boundary edge points is needed
  i_startblk = ptr_pp%edges%start_blk(grf_bdyintp_start_e,i_chidx)
  i_endblk   = ptr_pp%edges%end_blk(grf_bdyintp_end_e,i_chidx)

!$OMP DO PRIVATE (jk,je,dvn_tang) ICON_OMP_DEFAULT_SCHEDULE
  DO jb =  i_startblk, i_endblk

    CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
     i_startidx, i_endidx, grf_bdyintp_start_e, grf_bdyintp_end_e, i_chidx)

#ifndef __LOOP_EXCHANGE
    ! child edges 1 and 2
!CDIR UNROLL=6
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx

        dvn_tang(je) = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,2)%v2 - &
                      (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,1)%v2 )

        vn_aux(je,jk,jb,1,jn) = p_in(jn)%fld(je,jk+js,jb) + dvn_tang(je)*ptr_dist(je,1,jb)
        vn_aux(je,jk,jb,2,jn) = p_in(jn)%fld(je,jk+js,jb) + dvn_tang(je)*ptr_dist(je,2,jb)
      ENDDO
    ENDDO

    ! child edge 3
!CDIR UNROLL=6
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx
        vn_aux(je,jk,jb,3,jn) = ptr_grf%grf_vec_coeff_2a(1,je,jb) * &
          p_in(jn)%fld(iidx_2a(je,1,jb),jk+js,iblk_2a(je,1,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(2,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,2,jb),jk+js,iblk_2a(je,2,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(3,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,3,jb),jk+js,iblk_2a(je,3,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(4,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,4,jb),jk+js,iblk_2a(je,4,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(5,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,5,jb),jk+js,iblk_2a(je,5,jb))
      ENDDO
    ENDDO

    ! child edge 4
!CDIR UNROLL=6
    DO jk = 1, nlev_c
      DO je = i_startidx, i_endidx
        IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
          vn_aux(je,jk,jb,4,jn) = ptr_grf%grf_vec_coeff_2b(1,je,jb) * &
            p_in(jn)%fld(iidx_2b(je,1,jb),jk+js,iblk_2b(je,1,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(2,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,2,jb),jk+js,iblk_2b(je,2,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(3,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,3,jb),jk+js,iblk_2b(je,3,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(4,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,4,jb),jk+js,iblk_2b(je,4,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(5,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,5,jb),jk+js,iblk_2b(je,5,jb))
        ENDIF
      ENDDO
    ENDDO

#else

    DO je = i_startidx, i_endidx
      DO jk = 1, nlev_c

        ! child edges 1 and 2
        dvn_tang(je) = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,2)%v1 + &
                       v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,2)%v2 - &
                      (u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,1)%v1 + &
                       v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) *    &
                       ptr_pp%edges%primal_normal_vert(je,jb,1)%v2 )

        vn_aux(je,jk,jb,1,jn) = p_in(jn)%fld(je,jk+js,jb) + dvn_tang(je)*ptr_dist(je,1,jb)
        vn_aux(je,jk,jb,2,jn) = p_in(jn)%fld(je,jk+js,jb) + dvn_tang(je)*ptr_dist(je,2,jb)

        ! child edge 3
        vn_aux(je,jk,jb,3,jn) = ptr_grf%grf_vec_coeff_2a(1,je,jb) * &
          p_in(jn)%fld(iidx_2a(je,1,jb),jk+js,iblk_2a(je,1,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(2,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,2,jb),jk+js,iblk_2a(je,2,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(3,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,3,jb),jk+js,iblk_2a(je,3,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(4,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,4,jb),jk+js,iblk_2a(je,4,jb)) +   &
          ptr_grf%grf_vec_coeff_2a(5,je,jb) *                       &
          p_in(jn)%fld(iidx_2a(je,5,jb),jk+js,iblk_2a(je,5,jb))

        ! child edge 4
        IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
          vn_aux(je,jk,jb,4,jn) = ptr_grf%grf_vec_coeff_2b(1,je,jb) * &
            p_in(jn)%fld(iidx_2b(je,1,jb),jk+js,iblk_2b(je,1,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(2,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,2,jb),jk+js,iblk_2b(je,2,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(3,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,3,jb),jk+js,iblk_2b(je,3,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(4,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,4,jb),jk+js,iblk_2b(je,4,jb)) +   &
            ptr_grf%grf_vec_coeff_2b(5,je,jb) *                       &
            p_in(jn)%fld(iidx_2b(je,5,jb),jk+js,iblk_2b(je,5,jb))
        ENDIF

      ENDDO
    ENDDO

#endif

  ENDDO ! blocks
!$OMP END DO

ENDDO ! fields


! Store results in p_vn_out

IF (my_process_is_mpi_seq()) THEN

  DO jn = 1, nfields

!$OMP DO PRIVATE (jk,je) ICON_OMP_DEFAULT_SCHEDULE
      DO jb =  i_startblk, i_endblk

      CALL get_indices_e(ptr_pp, jb, i_startblk, i_endblk, &
        i_startidx, i_endidx, grf_bdyintp_start_e, grf_bdyintp_end_e, i_chidx)

#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev_c
#else
      DO jk = 1, nlev_c
!CDIR NODEP,VOVERTAKE,VOB
        DO je = i_startidx, i_endidx
#endif

          p_out(jn)%fld(icheidx(je,jb,1),jk,icheblk(je,jb,1))   = vn_aux(je,jk,jb,1,jn)
          p_out(jn)%fld(icheidx(je,jb,2),jk,icheblk(je,jb,2))   = vn_aux(je,jk,jb,2,jn)

          ! Avoid overwriting prognostic grid points in case of applying this to a prognostic variable
          IF (ptr_pc%edges%refin_ctrl(icheidx(je,jb,3),icheblk(je,jb,3)) <= grf_bdywidth_e) &
            p_out(jn)%fld(icheidx(je,jb,3),jk,icheblk(je,jb,3))   = vn_aux(je,jk,jb,3,jn)

          IF (ptr_pp%edges%refin_ctrl(je,jb) /= -1) THEN
            IF (ptr_pc%edges%refin_ctrl(icheidx(je,jb,4),icheblk(je,jb,4)) <= grf_bdywidth_e) &
              & p_out(jn)%fld(icheidx(je,jb,4),jk,icheblk(je,jb,4)) = vn_aux(je,jk,jb,4,jn)
          ENDIF

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

  ENDDO
ENDIF ! not MPI-parallel

!$OMP END PARALLEL

IF (my_process_is_mpi_parallel()) THEN

  nsendtot = SUM(ptr_pc%comm_pat_interpol_vec_grf(1:4)%n_send)
  nrecvtot = SUM(ptr_pc%comm_pat_interpol_vec_grf(1:4)%n_recv)

  IF (nfields == 1) THEN
    nlevtot = nlev_c
    CALL exchange_data_grf(ptr_pc%comm_pat_interpol_vec_grf,nfields,nlevtot,nsendtot, &
                           nrecvtot,RECV1=f3dout1,SEND1=vn_aux(:,:,:,:,1),            &
                           SEND_LBOUND3=LBOUND(vn_aux,3))

  ELSE IF (nfields == 2) THEN
    nlevtot = 2*nlev_c
    CALL exchange_data_grf(ptr_pc%comm_pat_interpol_vec_grf,nfields,nlevtot,nsendtot,   &
                           nrecvtot,RECV1=f3dout1,SEND1=vn_aux(:,:,:,:,1),RECV2=f3dout2,&
                           SEND2=vn_aux(:,:,:,:,2),SEND_LBOUND3=LBOUND(vn_aux,3))

  ELSE IF (nfields == 3) THEN
    nlevtot = 3*nlev_c
    CALL exchange_data_grf(ptr_pc%comm_pat_interpol_vec_grf,nfields,nlevtot,nsendtot,    &
                           nrecvtot,RECV1=f3dout1,SEND1=vn_aux(:,:,:,:,1),RECV2=f3dout2, &
                           SEND2=vn_aux(:,:,:,:,2),RECV3=f3dout3,SEND3=vn_aux(:,:,:,:,3),&
                           SEND_LBOUND3=LBOUND(vn_aux,3))

  ELSE IF (nfields == 4) THEN
    nlevtot = 4*nlev_c
    CALL exchange_data_grf(ptr_pc%comm_pat_interpol_vec_grf,nfields,nlevtot,nsendtot,         &
                           nrecvtot,RECV1=f3dout1,SEND1=vn_aux(:,:,:,:,1),RECV2=f3dout2,      &
                           SEND2=vn_aux(:,:,:,:,2),RECV3=f3dout3,SEND3=vn_aux(:,:,:,:,3),     &
                           RECV4=f3dout4,SEND4=vn_aux(:,:,:,:,4),SEND_LBOUND3=LBOUND(vn_aux,3))
  ENDIF
ENDIF

END SUBROUTINE interpol2_vec_grf


!-------------------------------------------------------------------------
!
!>
!! Performs interpolation of scalar fields from parent cells to child
!! cells using the 2D gradient at the cell center.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2009-12-16)
!!
SUBROUTINE interpol_scal_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, i_chidx, nfields,&
                              f3din1, f3dout1, f3din2, f3dout2, f3din3, f3dout3, &
                              f3din4, f3dout4, f3din5, f3dout5, f3din6, f3dout6, &
                              f4din1, f4dout1, f4din2, f4dout2,                  &
                              lpar_fields, llimit_nneg, lnoshift )
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc


! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf
TYPE(t_int_state),              TARGET, INTENT(IN)    ::  ptr_int

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! number of fields provided on input (needed for aux fields and pointer allocation)
INTEGER, INTENT(IN) :: nfields

! logical switch: if present and true, do OpenMP parallelization over nfields rather than nlev
LOGICAL, INTENT(IN), OPTIONAL :: lpar_fields

! logical switch: if present and true, limit horizontal gradient so as to avoid negative values
LOGICAL, INTENT(IN), OPTIONAL :: llimit_nneg(nfields)

! logical switch: if present and true, turn off accounting for shifts related to vertical nesting
! (needed for boundary interpolation of 2D diagnostic variables combined into a 3D field)
LOGICAL, INTENT(IN), OPTIONAL :: lnoshift

! input scalar fields at cell points (up to six 3D fields or two 4D fields)
REAL(wp), INTENT(IN), OPTIONAL, TARGET ::  & ! dim: (nproma,nlev,nblks_c)
  f3din1(:,:,:), f3din2(:,:,:), f3din3(:,:,:), f3din4(:,:,:), f3din5(:,:,:), f3din6(:,:,:), &
  f4din1(:,:,:,:), f4din2(:,:,:,:)

! reconstructed scalar output fields (up to six 3D fields or two 4D fields)
REAL(wp), INTENT(INOUT), OPTIONAL, TARGET ::  & ! dim: (nproma,nlev_c,nblks_c)
  f3dout1(:,:,:), f3dout2(:,:,:), f3dout3(:,:,:), f3dout4(:,:,:), f3dout5(:,:,:), f3dout6(:,:,:), &
  f4dout1(:,:,:,:), f4dout2(:,:,:,:)

INTEGER :: jb, jk, jc, jn, n         ! loop indices
INTEGER :: js                        ! shift parameter
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index
INTEGER :: elev                      ! end index of vertical loop

INTEGER :: nsendtot, nrecvtot, nlevtot  ! for MPI communication call
INTEGER :: n4d

LOGICAL :: l_par_fields              ! local variable corresponding to lpar_fields
LOGICAL :: l_limit_nneg(nfields)     ! local variable corresponding to llimit_nneg
LOGICAL :: l_noshift                 ! local variable corresponding to lnoshift
LOGICAL :: l4d                       ! 4D field is provided as input

! Auxiliary fields
REAL(wp), DIMENSION(nproma,MAX(32,ptr_pc%nlevp1)) :: grad_x, grad_y, maxval_neighb, minval_neighb
REAL(wp) :: h_aux(nproma,MAX(32,ptr_pc%nlevp1), &
                      ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx):  &
                  MAX(ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx),  &
                      ptr_pp%cells%end_blk(grf_bdyintp_end_c,i_chidx)),4,nfields)
REAL(wp) :: limfac1, limfac2, limfac, min_expval, max_expval, epsi, ovsht_fac, r_ovsht_fac, &
            relaxed_minval, relaxed_maxval

! Pointers to index and coefficient fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk, ichcidx, ichcblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff, ptr_dist

! Allocatable pointer to input and output fields
TYPE t_fieldptr
  REAL(wp), POINTER :: fld(:,:,:)
END TYPE t_fieldptr
TYPE(t_fieldptr) :: p_in(nfields), p_out(nfields)


!-----------------------------------------------------------------------

IF (PRESENT(f4din1) .AND. .NOT. PRESENT(f4din2)) THEN
  DO n = 1, nfields
    p_in(n)%fld  => f4din1(:,:,:,n)
    p_out(n)%fld => f4dout1(:,:,:,n)
  ENDDO
  l4d = .TRUE.
ELSE IF (PRESENT(f4din1) .AND. PRESENT(f4din2)) THEN
  n4d = nfields/2
  DO n = 1, n4d
    p_in(n)%fld  => f4din1(:,:,:,n)
    p_out(n)%fld => f4dout1(:,:,:,n)
  ENDDO
  DO n = 1, n4d
    p_in(n4d+n)%fld  => f4din2(:,:,:,n)
    p_out(n4d+n)%fld => f4dout2(:,:,:,n)
  ENDDO
  l4d = .TRUE.
ELSE
  IF (PRESENT(f3din1)) THEN
    p_in(1)%fld  => f3din1
    p_out(1)%fld => f3dout1
  ENDIF
  IF (PRESENT(f3din2)) THEN
    p_in(2)%fld  => f3din2
    p_out(2)%fld => f3dout2
  ENDIF
  IF (PRESENT(f3din3)) THEN
    p_in(3)%fld  => f3din3
    p_out(3)%fld => f3dout3
  ENDIF
  IF (PRESENT(f3din4)) THEN
    p_in(4)%fld  => f3din4
    p_out(4)%fld => f3dout4
  ENDIF
  IF (PRESENT(f3din5)) THEN
    p_in(5)%fld  => f3din5
    p_out(5)%fld => f3dout5
  ENDIF
  IF (PRESENT(f3din6)) THEN
    p_in(6)%fld  => f3din6
    p_out(6)%fld => f3dout6
  ENDIF
  l4d = .FALSE.
ENDIF

! Check if jk loop or jn loop shall be OpenMP parallelized
IF (PRESENT(lpar_fields)) THEN
  l_par_fields = lpar_fields
ELSE
  l_par_fields = .FALSE.
ENDIF

! Check if gradient limiting is required
IF (PRESENT(llimit_nneg)) THEN
  l_limit_nneg(:) = llimit_nneg(:)
ELSE
  l_limit_nneg(:) = .FALSE.
ENDIF

! Check if vertical shifting is to be turned off
IF (PRESENT(lnoshift)) THEN
  l_noshift = lnoshift
ELSE
  l_noshift = .FALSE.
ENDIF

epsi = 1.e-75_wp
ovsht_fac = 1.05_wp ! factor of allowed overshooting
r_ovsht_fac = 1._wp/ovsht_fac

! Start and end blocks for which scalar interpolation is needed
i_startblk = ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx)
i_endblk   = ptr_pp%cells%end_blk(grf_bdyintp_end_c,i_chidx)

! Pointers to index lists and interpolation coefficients for gradient computation
iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff
ptr_dist  => ptr_grf%grf_dist_pc2cc

! child cell indices and blocks for non-MPI parent-to-child communication
ichcidx => ptr_pp%cells%child_idx
ichcblk => ptr_pp%cells%child_blk

IF (l_noshift) THEN
  js = 0
ELSE
  js = ptr_pc%nshift
ENDIF

IF (p_test_run) h_aux = 0._wp

!$OMP PARALLEL PRIVATE(jb,jn,elev,i_startidx,i_endidx)

IF (l_par_fields) THEN ! parallelization over fields
!$OMP DO PRIVATE (jk,jc,grad_x,grad_y,min_expval,max_expval,limfac1, &
!$OMP   limfac2,limfac,maxval_neighb,minval_neighb,relaxed_minval,   &
!$OMP   relaxed_maxval) ICON_OMP_DEFAULT_SCHEDULE

  DO jn = 1, nfields

    elev   = UBOUND(p_out(jn)%fld,2)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, elev
#else
      DO jk = 1, elev
        DO jc = i_startidx, i_endidx
#endif

          grad_x(jc,jk) =  &
            ptr_coeff(1,1,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,1,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,1,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,1,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,1,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,1,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,1,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,1,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,1,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,1,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          grad_y(jc,jk) =  &
            ptr_coeff(1,2,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,2,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,2,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,2,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,2,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,2,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,2,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,2,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,2,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,2,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          maxval_neighb(jc,jk) =                               &
            MAX(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))
          minval_neighb(jc,jk) =                               &
            MIN(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))

        ENDDO
      ENDDO

      DO jk = 1, elev
        DO jc = i_startidx, i_endidx
          min_expval = MIN(grad_x(jc,jk)*ptr_dist(jc,1,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,1,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,2,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,2,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,3,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,3,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,4,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,4,2,jb),  &
                           -1.e-80_wp )
          max_expval = MAX(grad_x(jc,jk)*ptr_dist(jc,1,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,1,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,2,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,2,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,3,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,3,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,4,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,4,2,jb),  &
                           1.e-80_wp )

          limfac1 = 1._wp
          limfac2 = 1._wp
          ! Allow a limited amount of over-/undershooting in the downscaled fields
          IF (minval_neighb(jc,jk) > 0._wp) THEN
            relaxed_minval = r_ovsht_fac*minval_neighb(jc,jk)
          ELSE
            relaxed_minval = ovsht_fac*minval_neighb(jc,jk)
          ENDIF
          IF (maxval_neighb(jc,jk) > 0._wp) THEN
            relaxed_maxval = ovsht_fac*maxval_neighb(jc,jk)
          ELSE
            relaxed_maxval = r_ovsht_fac*maxval_neighb(jc,jk)
          ENDIF

          IF (p_in(jn)%fld(jc,jk+js,jb) + min_expval < relaxed_minval-epsi) THEN
            limfac1 = ABS((relaxed_minval-p_in(jn)%fld(jc,jk+js,jb))/min_expval)
          ENDIF
          IF (p_in(jn)%fld(jc,jk+js,jb) + max_expval > relaxed_maxval+epsi) THEN
            limfac2 = ABS((relaxed_maxval-p_in(jn)%fld(jc,jk+js,jb))/max_expval)
          ENDIF
          limfac = MIN(limfac1,limfac2)

          grad_x(jc,jk) = grad_x(jc,jk)*limfac
          grad_y(jc,jk) = grad_y(jc,jk)*limfac

        ENDDO
      ENDDO

      IF (l_limit_nneg(jn)) THEN
        DO jk = 1, elev
          DO jc = i_startidx, i_endidx

            h_aux(jc,jk,jb,1,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,1,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,1,2,jb))
            h_aux(jc,jk,jb,2,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,2,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,2,2,jb))
            h_aux(jc,jk,jb,3,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,3,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,3,2,jb))
            h_aux(jc,jk,jb,4,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,4,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,4,2,jb))

          ENDDO
        ENDDO
      ELSE
        DO jk = 1, elev
          DO jc = i_startidx, i_endidx

            h_aux(jc,jk,jb,1,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,1,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,1,2,jb)
            h_aux(jc,jk,jb,2,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,2,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,2,2,jb)
            h_aux(jc,jk,jb,3,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,3,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,3,2,jb)
            h_aux(jc,jk,jb,4,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,4,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,4,2,jb)

          ENDDO
        ENDDO
      ENDIF

    ENDDO ! blocks
  ENDDO ! fields
!$OMP END DO

ELSE ! parallelization over jb loop

  DO jn = 1, nfields

    elev = UBOUND(p_out(jn)%fld,2)

!$OMP DO PRIVATE (jk,jc,grad_x,grad_y,min_expval,max_expval,limfac1,  &
!$OMP   limfac2,limfac,maxval_neighb,minval_neighb,relaxed_minval,    &
!$OMP   relaxed_maxval) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, elev
#else
      DO jk = 1, elev
        DO jc = i_startidx, i_endidx
#endif

          grad_x(jc,jk) =  &
            ptr_coeff(1,1,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,1,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,1,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,1,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,1,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,1,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,1,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,1,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,1,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,1,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          grad_y(jc,jk) =  &
            ptr_coeff(1,2,jc,jb)*p_in(jn)%fld(jc,jk+js,jb) + &
            ptr_coeff(2,2,jc,jb)*p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)) + &
            ptr_coeff(3,2,jc,jb)*p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)) + &
            ptr_coeff(4,2,jc,jb)*p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)) + &
            ptr_coeff(5,2,jc,jb)*p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)) + &
            ptr_coeff(6,2,jc,jb)*p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)) + &
            ptr_coeff(7,2,jc,jb)*p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)) + &
            ptr_coeff(8,2,jc,jb)*p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)) + &
            ptr_coeff(9,2,jc,jb)*p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)) + &
            ptr_coeff(10,2,jc,jb)*p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb))
          maxval_neighb(jc,jk) =                               &
            MAX(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))
          minval_neighb(jc,jk) =                              &
            MIN(p_in(jn)%fld(jc,jk+js,jb),                       &
                p_in(jn)%fld(iidx(2,jc,jb),jk+js,iblk(2,jc,jb)), &
                p_in(jn)%fld(iidx(3,jc,jb),jk+js,iblk(3,jc,jb)), &
                p_in(jn)%fld(iidx(4,jc,jb),jk+js,iblk(4,jc,jb)), &
                p_in(jn)%fld(iidx(5,jc,jb),jk+js,iblk(5,jc,jb)), &
                p_in(jn)%fld(iidx(6,jc,jb),jk+js,iblk(6,jc,jb)), &
                p_in(jn)%fld(iidx(7,jc,jb),jk+js,iblk(7,jc,jb)), &
                p_in(jn)%fld(iidx(8,jc,jb),jk+js,iblk(8,jc,jb)), &
                p_in(jn)%fld(iidx(9,jc,jb),jk+js,iblk(9,jc,jb)), &
                p_in(jn)%fld(iidx(10,jc,jb),jk+js,iblk(10,jc,jb)))

        ENDDO
      ENDDO

      DO jk = 1, elev
        DO jc = i_startidx, i_endidx
          min_expval = MIN(grad_x(jc,jk)*ptr_dist(jc,1,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,1,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,2,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,2,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,3,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,3,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,4,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,4,2,jb),  &
                           -1.e-80_wp )
          max_expval = MAX(grad_x(jc,jk)*ptr_dist(jc,1,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,1,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,2,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,2,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,3,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,3,2,jb),  &
                           grad_x(jc,jk)*ptr_dist(jc,4,1,jb) + &
                           grad_y(jc,jk)*ptr_dist(jc,4,2,jb),  &
                           1.e-80_wp )

          limfac1 = 1._wp
          limfac2 = 1._wp
          ! Allow a limited amount of over-/undershooting in the downscaled fields
          IF (minval_neighb(jc,jk) > 0._wp) THEN
            relaxed_minval = r_ovsht_fac*minval_neighb(jc,jk)
          ELSE
            relaxed_minval = ovsht_fac*minval_neighb(jc,jk)
          ENDIF
          IF (maxval_neighb(jc,jk) > 0._wp) THEN
            relaxed_maxval = ovsht_fac*maxval_neighb(jc,jk)
          ELSE
            relaxed_maxval = r_ovsht_fac*maxval_neighb(jc,jk)
          ENDIF

          IF (p_in(jn)%fld(jc,jk+js,jb) + min_expval < relaxed_minval-epsi) THEN
            limfac1 = ABS((relaxed_minval-p_in(jn)%fld(jc,jk+js,jb))/min_expval)
          ENDIF
          IF (p_in(jn)%fld(jc,jk+js,jb) + max_expval > relaxed_maxval+epsi) THEN
            limfac2 = ABS((relaxed_maxval-p_in(jn)%fld(jc,jk+js,jb))/max_expval)
          ENDIF
          limfac = MIN(limfac1,limfac2)

          grad_x(jc,jk) = grad_x(jc,jk)*limfac
          grad_y(jc,jk) = grad_y(jc,jk)*limfac

        ENDDO
      ENDDO

      IF (l_limit_nneg(jn)) THEN
        DO jk = 1, elev
          DO jc = i_startidx, i_endidx

            h_aux(jc,jk,jb,1,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,1,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,1,2,jb))
            h_aux(jc,jk,jb,2,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,2,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,2,2,jb))
            h_aux(jc,jk,jb,3,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,3,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,3,2,jb))
            h_aux(jc,jk,jb,4,jn) = MAX(0._wp, p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,4,1,jb)                         + &
              grad_y(jc,jk)*ptr_dist(jc,4,2,jb))

          ENDDO
        ENDDO

      ELSE

        DO jk = 1, elev
          DO jc = i_startidx, i_endidx

            h_aux(jc,jk,jb,1,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,1,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,1,2,jb)
            h_aux(jc,jk,jb,2,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,2,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,2,2,jb)
            h_aux(jc,jk,jb,3,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,3,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,3,2,jb)
            h_aux(jc,jk,jb,4,jn) = p_in(jn)%fld(jc,jk+js,jb) + &
              grad_x(jc,jk)*ptr_dist(jc,4,1,jb)              + &
              grad_y(jc,jk)*ptr_dist(jc,4,2,jb)

          ENDDO
        ENDDO
      ENDIF

    ENDDO ! blocks
!$OMP END DO
  ENDDO ! fields

ENDIF


! Store results in p_out

IF (my_process_is_mpi_seq()) THEN

  DO jn = 1, nfields

    elev = UBOUND(p_out(jn)%fld,2)

!$OMP DO PRIVATE (jk,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb =  i_startblk, i_endblk

      CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
           i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)

#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 1, elev
#else
      DO jk = 1, elev
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = i_startidx, i_endidx
#endif

          p_out(jn)%fld(ichcidx(jc,jb,1),jk,ichcblk(jc,jb,1)) = h_aux(jc,jk,jb,1,jn)
          p_out(jn)%fld(ichcidx(jc,jb,2),jk,ichcblk(jc,jb,2)) = h_aux(jc,jk,jb,2,jn)
          p_out(jn)%fld(ichcidx(jc,jb,3),jk,ichcblk(jc,jb,3)) = h_aux(jc,jk,jb,3,jn)
          p_out(jn)%fld(ichcidx(jc,jb,4),jk,ichcblk(jc,jb,4)) = h_aux(jc,jk,jb,4,jn)

        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDDO

ENDIF
!$OMP END PARALLEL

IF (my_process_is_mpi_parallel()) THEN

  nsendtot = SUM(ptr_pc%comm_pat_interpol_scal_grf(1:4)%n_send)
  nrecvtot = SUM(ptr_pc%comm_pat_interpol_scal_grf(1:4)%n_recv)

  IF (l4d .AND. .NOT. PRESENT(f4din2)) THEN

    nlevtot = SIZE(f4dout1,2)*nfields
    CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot,       &
                           nrecvtot,RECV4D1=f4dout1,SEND4D1=h_aux,SEND_LBOUND3=LBOUND(h_aux,3))

  ELSE IF (l4d) THEN

    nlevtot = SIZE(f4dout1,2)*n4d+SIZE(f4dout2,2)*n4d
    CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot,  &
                           nrecvtot,RECV4D1=f4dout1,SEND4D1=h_aux(:,:,:,:,1:n4d),       &
                           RECV4D2=f4dout2,SEND4D2=h_aux(:,:,:,:,n4d+1:nfields),        &
                           SEND_LBOUND3=LBOUND(h_aux,3))

  ELSE
    IF (nfields == 1) THEN
      nlevtot = SIZE(f3dout1,2)
      CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot, &
                             nrecvtot,RECV1=f3dout1,SEND1=h_aux(:,:,:,:,1),              &
                             SEND_LBOUND3=LBOUND(h_aux,3))

    ELSE IF (nfields == 2) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)
      CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot, &
                             nrecvtot,RECV1=f3dout1,SEND1=h_aux(:,:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,:,2),SEND_LBOUND3=LBOUND(h_aux,3))

    ELSE IF (nfields == 3) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2)
      CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot, &
                             nrecvtot,RECV1=f3dout1,SEND1=h_aux(:,:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,:,3),&
                             SEND_LBOUND3=LBOUND(h_aux,3))

    ELSE IF (nfields == 4) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2)+SIZE(f3dout4,2)
      CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot,      &
                             nrecvtot,RECV1=f3dout1,SEND1=h_aux(:,:,:,:,1),RECV2=f3dout2,     &
                             SEND2=h_aux(:,:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,:,3),     &
                             RECV4=f3dout4,SEND4=h_aux(:,:,:,:,4),SEND_LBOUND3=LBOUND(h_aux,3))

    ELSE IF (nfields == 5) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2) + &
                SIZE(f3dout4,2)+SIZE(f3dout5,2)
      CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot, &
                             nrecvtot,RECV1=f3dout1,SEND1=h_aux(:,:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,:,3),&
                             RECV4=f3dout4,SEND4=h_aux(:,:,:,:,4),RECV5=f3dout5,         &
                             SEND5=h_aux(:,:,:,:,5),SEND_LBOUND3=LBOUND(h_aux,3))

    ELSE IF (nfields == 6) THEN
      nlevtot = SIZE(f3dout1,2)+SIZE(f3dout2,2)+SIZE(f3dout3,2) + &
                SIZE(f3dout4,2)+SIZE(f3dout5,2)+SIZE(f3dout6,2)
      CALL exchange_data_grf(ptr_pc%comm_pat_interpol_scal_grf,nfields,nlevtot,nsendtot, &
                             nrecvtot,RECV1=f3dout1,SEND1=h_aux(:,:,:,:,1),RECV2=f3dout2,&
                             SEND2=h_aux(:,:,:,:,2),RECV3=f3dout3,SEND3=h_aux(:,:,:,:,3),&
                             RECV4=f3dout4,SEND4=h_aux(:,:,:,:,4),RECV5=f3dout5,         &
                             SEND5=h_aux(:,:,:,:,5),RECV6=f3dout6,SEND6=h_aux(:,:,:,:,6),&
                             SEND_LBOUND3=LBOUND(h_aux,3))

    ENDIF
  ENDIF
ENDIF

END SUBROUTINE interpol_scal_grf


!-------------------------------------------------------------------------
!
!>
!! Performs interpolation of scalar fields from parent cells to child
!! cells using the 2D gradient at the cell center.
!! Version for 2D input and output fields (i.e. psfc in hydrostatic model)
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD (2009-12-16)
!!
SUBROUTINE interpol_scal2d_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, &
                                i_chidx, p_h_in, p_h_out, llimit_nneg )
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pp
TYPE(t_patch), TARGET, INTENT(in) :: ptr_pc

! input scalar field at cell points
REAL(wp), INTENT(IN) ::  &
  p_h_in(:,:) ! dim: (nproma,nblks_c)

! Indices of source points and interpolation coefficients
TYPE(t_gridref_single_state),   TARGET, INTENT(IN)    ::  ptr_grf
TYPE(t_int_state),              TARGET, INTENT(IN)    ::  ptr_int

! child domain index as seen from parent domain
INTEGER, INTENT(IN) :: i_chidx

! logical switch: if present and true, limit horizontal gradient so as to avoid negative values
LOGICAL, INTENT(IN), OPTIONAL :: llimit_nneg

! reconstructed scalar output field
REAL(wp),INTENT(INOUT) :: p_h_out(:,:) ! dim: (nproma,nblks_c)

INTEGER :: jb, jc, n                 ! loop indices
INTEGER :: i_startblk                ! start block
INTEGER :: i_endblk                  ! end index
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

LOGICAL :: l_limit_nneg              ! local variable corresponding to llimit_nneg

! Auxiliary fields
REAL(wp), DIMENSION(nproma,ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx):&
                   MAX(ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx),    &
                       ptr_pp%cells%end_blk(grf_bdyintp_end_c,i_chidx))) :: grad_x, grad_y
REAL(wp) :: h_aux(nproma,ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx):&
                  MAX(ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx),   &
                      ptr_pp%cells%end_blk(grf_bdyintp_end_c,i_chidx)),4)

REAL(wp) :: limfac, min_expval

! Pointers to index and coefficient fields
INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk, ichcidx, ichcblk
REAL(wp), DIMENSION(:,:,:,:), POINTER :: ptr_coeff, ptr_dist

!-----------------------------------------------------------------------

! Check if gradient limiting is required
IF (PRESENT(llimit_nneg)) THEN
  l_limit_nneg = llimit_nneg
ELSE
  l_limit_nneg = .FALSE.
ENDIF


! Start and end blocks for which scalar interpolation is needed
i_startblk = ptr_pp%cells%start_blk(grf_bdyintp_start_c,i_chidx)
i_endblk   = ptr_pp%cells%end_blk(grf_bdyintp_end_c,i_chidx)

! Pointers to index lists and interpolation coefficients for gradient computation
iidx => ptr_int%rbf_c2grad_idx
iblk => ptr_int%rbf_c2grad_blk

ptr_coeff => ptr_int%rbf_c2grad_coeff
ptr_dist  => ptr_grf%grf_dist_pc2cc

! child cell indices and blocks for non-MPI parent-to-child communication
ichcidx => ptr_pp%cells%child_idx
ichcblk => ptr_pp%cells%child_blk


DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)


  DO jc = i_startidx, i_endidx

    grad_x(jc,jb) =  &
      ptr_coeff(1,1,jc,jb)*p_h_in(jc,jb) + &
      ptr_coeff(2,1,jc,jb)*p_h_in(iidx(2,jc,jb),iblk(2,jc,jb)) + &
      ptr_coeff(3,1,jc,jb)*p_h_in(iidx(3,jc,jb),iblk(3,jc,jb)) + &
      ptr_coeff(4,1,jc,jb)*p_h_in(iidx(4,jc,jb),iblk(4,jc,jb)) + &
      ptr_coeff(5,1,jc,jb)*p_h_in(iidx(5,jc,jb),iblk(5,jc,jb)) + &
      ptr_coeff(6,1,jc,jb)*p_h_in(iidx(6,jc,jb),iblk(6,jc,jb)) + &
      ptr_coeff(7,1,jc,jb)*p_h_in(iidx(7,jc,jb),iblk(7,jc,jb)) + &
      ptr_coeff(8,1,jc,jb)*p_h_in(iidx(8,jc,jb),iblk(8,jc,jb)) + &
      ptr_coeff(9,1,jc,jb)*p_h_in(iidx(9,jc,jb),iblk(9,jc,jb)) + &
      ptr_coeff(10,1,jc,jb)*p_h_in(iidx(10,jc,jb),iblk(10,jc,jb))
    grad_y(jc,jb) =  &
      ptr_coeff(1,2,jc,jb)*p_h_in(jc,jb) + &
      ptr_coeff(2,2,jc,jb)*p_h_in(iidx(2,jc,jb),iblk(2,jc,jb)) + &
      ptr_coeff(3,2,jc,jb)*p_h_in(iidx(3,jc,jb),iblk(3,jc,jb)) + &
      ptr_coeff(4,2,jc,jb)*p_h_in(iidx(4,jc,jb),iblk(4,jc,jb)) + &
      ptr_coeff(5,2,jc,jb)*p_h_in(iidx(5,jc,jb),iblk(5,jc,jb)) + &
      ptr_coeff(6,2,jc,jb)*p_h_in(iidx(6,jc,jb),iblk(6,jc,jb)) + &
      ptr_coeff(7,2,jc,jb)*p_h_in(iidx(7,jc,jb),iblk(7,jc,jb)) + &
      ptr_coeff(8,2,jc,jb)*p_h_in(iidx(8,jc,jb),iblk(8,jc,jb)) + &
      ptr_coeff(9,2,jc,jb)*p_h_in(iidx(9,jc,jb),iblk(9,jc,jb)) + &
      ptr_coeff(10,2,jc,jb)*p_h_in(iidx(10,jc,jb),iblk(10,jc,jb))

  ENDDO

  ! Limit horizontal gradients to avoid generation of negative values
  IF (l_limit_nneg) THEN

    DO jc = i_startidx, i_endidx
      min_expval = MIN(grad_x(jc,jb)*ptr_dist(jc,1,1,jb) + &
                       grad_y(jc,jb)*ptr_dist(jc,1,2,jb),  &
                       grad_x(jc,jb)*ptr_dist(jc,2,1,jb) + &
                       grad_y(jc,jb)*ptr_dist(jc,2,2,jb),  &
                       grad_x(jc,jb)*ptr_dist(jc,3,1,jb) + &
                       grad_y(jc,jb)*ptr_dist(jc,3,2,jb),  &
                       grad_x(jc,jb)*ptr_dist(jc,4,1,jb) + &
                       grad_y(jc,jb)*ptr_dist(jc,4,2,jb),  &
                       -1.e-100_wp )

      IF (p_h_in(jc,jb) + min_expval < 0._wp) THEN
        limfac = p_h_in(jc,jb)/ABS(min_expval)
        grad_x(jc,jb) = grad_x(jc,jb)*limfac
        grad_y(jc,jb) = grad_y(jc,jb)*limfac
      ENDIF

    ENDDO

  ENDIF !l_limit_nneg

ENDDO

DO jb =  i_startblk, i_endblk

  CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)

  DO jc = i_startidx, i_endidx

    h_aux(jc,jb,1) = p_h_in(jc,jb) + &
      grad_x(jc,jb)*ptr_dist(jc,1,1,jb) + &
      grad_y(jc,jb)*ptr_dist(jc,1,2,jb)
    h_aux(jc,jb,2) = p_h_in(jc,jb) + &
      grad_x(jc,jb)*ptr_dist(jc,2,1,jb) + &
      grad_y(jc,jb)*ptr_dist(jc,2,2,jb)
    h_aux(jc,jb,3) = p_h_in(jc,jb) + &
      grad_x(jc,jb)*ptr_dist(jc,3,1,jb) + &
      grad_y(jc,jb)*ptr_dist(jc,3,2,jb)
    h_aux(jc,jb,4) = p_h_in(jc,jb) + &
      grad_x(jc,jb)*ptr_dist(jc,4,1,jb) + &
      grad_y(jc,jb)*ptr_dist(jc,4,2,jb)

  ENDDO

ENDDO ! blocks

! Store results in p_h_out

IF (my_process_is_mpi_seq()) THEN

  DO jb =  i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdyintp_start_c, grf_bdyintp_end_c, i_chidx)

!CDIR NODEP,VOVERTAKE,VOB
    DO jc = i_startidx, i_endidx

      p_h_out(ichcidx(jc,jb,1),ichcblk(jc,jb,1)) = h_aux(jc,jb,1)
      p_h_out(ichcidx(jc,jb,2),ichcblk(jc,jb,2)) = h_aux(jc,jb,2)
      p_h_out(ichcidx(jc,jb,3),ichcblk(jc,jb,3)) = h_aux(jc,jb,3)
      p_h_out(ichcidx(jc,jb,4),ichcblk(jc,jb,4)) = h_aux(jc,jb,4)

    ENDDO
  ENDDO

ENDIF

IF (my_process_is_mpi_parallel()) THEN

  DO n = 1, 4
    CALL exchange_data(ptr_pc%comm_pat_interpol_scal_grf(n),RECV=p_h_out,SEND=h_aux(:,:,n), &
                       SEND_LBOUND2=LBOUND(h_aux,2))
  ENDDO

ENDIF


END SUBROUTINE interpol_scal2d_grf


!-------------------------------------------------------------------------
END MODULE mo_grf_bdyintp
