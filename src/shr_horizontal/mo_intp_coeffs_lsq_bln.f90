
#ifdef __xlC__
! @PROCESS nosmp
! @PROCESS NOOPTimize
! @PROCESS smp=noopt
@PROCESS noopt
#endif
#ifdef __PGI
!pgi$g opt=1
#endif
!>
!! Contains the implementation of interpolation and reconstruction.
!!
!! Contains the implementation of interpolation and reconstruction
!! routines used by the shallow water model, including the RBF
!! reconstruction routines.
!!
!! @par Revision History
!! Developed  by Luca Bonaventura and Will Sawyer (2002-4).
!! Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!! Adapted to new data structure by Thomas Heinze,
!! Peter Korn and Luca Bonaventura (2005).
!! Modification by Thomas Heinze (2006-02-21):
!! - renamed m_modules to mo_modules
!! Modification by Thomas Heinze (2006-07-05):
!! - modified cell2edge_lin_int_coeff
!! - created cc_dot_product
!! Modification by Peter Korn and Luca Bonaventura(2006-07-28):
!! - moved several auxiliary functions to mo_math_utilities
!! - introduced recoded rbf interpolation for vector fields
!! - added lraviart switch to force RT interpolation to be used
!! Modification by Thomas Heinze  and Luca Bonaventura(2006-10-05):
!! - merged with 'Milano' version by P. Korn
!! Modification by Pilar Ripodas (2006-11):
!! - new subroutine rbf_vec_interpol_car with the cartesian
!!   coordinates as output
!! Modification by Peter Korn, MPI-M, (2006-11-23):
!! - replacements in TYPE patch: ic by l2g_c, ie by l2g_e, iv by l2g_v,
!!   iic by g2l_c, iie by g2l_e, iiv by g2l_v
!! - replaced edge_index by edge_idx
!! - replaced vertex_index by vertex_idx
!! - replaced cell_index by cell_idx
!! - replaced neighbor_index by neighbor_idx
!! Modification by Pilar Ripodas (2006-12):
!! - dt_tan_vec and dt_tan_rt_vec are wrong. They are renamed to
!!   dt_tan_vec_old and dt_tan_rt_vec_old and should not be used
!! - New subroutines dt_tan_vec_h and dt_tan_vec_kin and
!!   dt_tan_vec_gen are produced and
!!   moved to mo_sw_state.f90
!!  Modification by Peter Korn, MPI-M (2007-02)
!!  Modification by Hui Wan, MPI-M (2007-02-22)
!!  - changes in the USE section because
!!    the coordinate types had been move from mo_model_domain
!!    to mo_math_utilities;
!!  Modification by Almut Gassmann, MPI-M (2007-04)
!!  - removed reference to unused halo_verts
!!  - summing over all halos of the various parallel patches (Quick and Dirty!)
!!  Modification by Almut Gassmann, MPI-M (2007-04)
!!  - abandon grid for the sake of patch
!!  Modification by Thomas Heinze, DWD (2007-07-26)
!!  - including all the improvements of Tobias Ruppert's diploma thesis
!!  - several changes according to the programming guide
!!  Modification by Pilar Ripodas, DWD (2007-07):
!!  - substruct the outgoing component of the reconstructed
!!    vector in subroutine "rbf_vec_interpol_car"
!!  Modification by Thomas Heinze, DWD (2007-08-02)
!!  - replaced rbf_kern_dim by rbf_kern_dim_c
!!  - replaced rbf_vec_dim by rbf_vec_dim_c
!!  - replaced rbf_mat_dim by rbf_mat_dim_c
!!  - replaced rbf_vec_scale by rbf_vec_scale_c
!!  - replaced rbf_vec_pdeg_c by rbf_vec_rbf_vec_pdeg_c_c
!!  Modification by Hui Wan, MPI-M (2007-08-02; 2007-11-30)
!!  - added interpolation coefficients c_aw_e and e_aw_c
!!    and the initialization subroutine aw_int_coeff.
!!  - added subroutine edges2cells_scalar
!!  Modification by Jochen Foerstner, DWD (2008-05-05)
!!  - four new subroutines
!!      rbf_vec_index_vertex
!!      rbf_vec_compute_coeff_vertex
!!      rbf_vec_interpol_car_vertex
!!      prepare_simpson
!!    to reconstruct a Cartesian vector at the vertices using
!!    RBF interpolation and to prepare quadrature via the
!!    Simpson's rule.
!!  Modification by Marco Restelli, MPI (2008-07-17)
!!  - included the subroutines
!!      cells2vertex_scalar, cells2vertex_coeff, ravtom_normgrad2,
!!      ls_normgrad2, ls_normgrad2_ii, edges2points_vector
!!    to compute polynomial fitting with sufficient accuracy as
!!    required in SW-alpha model.
!!  Modification by Jochen Foerstner, DWD (2008-09-12)
!!  - moved SUBROUTINE ravtom_normgrad2 to mo_math_operators
!!    because of conflicting use statements.
!!  Modification by Almut Gassmann, MPI-M (2008-10-09)
!!  - added features for helicity bracket reconstruction
!!  Modification by Guenther Zaengl, DWD (2008-10-23)
!!  - added interpolation routines needed for mesh refinement
!!  Modification by Almut Gassmann, MPI-M (2009-01-29)
!!  - conforming scalar interpolation routines and adjusting coefficients
!!  Modification by Guenther Zaengl, DWD (2009-02-11)
!!  - all routines needed for grid refinement are moved into the new
!!    module mo_grf_interpolation
!!  Modification by Guenther Zaengl, DWD (2009-02-13)
!!  - RBFs are changed to direct reconstruction of velocity components on
!!    the sphere, avoiding the detour over the 3D Cartesian space
!!  Modification by Almut Gassmann, DWD (2009-03-17)
!!  - remove lraviart
!!  Modification by Almut Gassmann, MPI-M (2009-04-23)
!!  - remove all Raviart Thomas stuff, add edge to verts averaging
!!  Modification by Daniel Reinert, DWD (2009-07-20)
!!  - added subroutine grad_lsq_compute_coeff_cell to prepare
!!    (2D) gradient reconstruction at circumcenter via the least squares
!!    method.
!!  Modification by Almut Gassmann, MPI-M (2009-10-05)
!!  - set RBF vec dimensions to predefined values (edges:4,vertices:6,cells:9);
!!    All other switches and belongings are deleted. The reason is that
!!    the Hollingsworth instability requires 4 edges, cell reconstruction
!!    is only needed for output and vertices are only used in the bracket
!!    version, where the dimension at the vertices should be 6
!!  Modification by Daniel Reinert, DWD (2009-12-10)
!!  - replaced grad_lsq_compute_coeff_cell by lsq_compute_coeff_cell
!!    which initializes either a second order or a third order least squares
!!    reconstruction.
!!  Modification by Almut Gassmann, MPI-M (2010-01-12)
!!  - generalize p_int%primal_normal_ec and p_int%edge_cell_length to hexagons
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_intp_coeffs_lsq_bln
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
USE mo_math_constants,      ONLY: pi2
USE mo_exception,           ONLY: message, finish
USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell
USE mo_model_domain,        ONLY: t_patch
USE mo_math_types,          ONLY: t_cartesian_coordinates
USE mo_math_utilities,      ONLY: gnomonic_proj, rotate_latlon, &
                                  plane_torus_closest_coordinates
USE mo_math_utility_solvers, ONLY: qrdec
USE mo_parallel_config,     ONLY: nproma
USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
USE mo_advection_config,    ONLY: advection_config
USE mo_sync,                ONLY: SYNC_C, SYNC_E, SYNC_V, sync_patch_array, sync_idx
USE mo_grid_config,         ONLY: grid_sphere_radius
USE mo_grid_geometry_info,  ONLY: planar_torus_geometry, sphere_geometry
USE mo_intp_data_strc,      ONLY: t_lsq, t_int_state
USE mo_fortran_tools,       ONLY: copy

IMPLICIT NONE

PRIVATE

PUBLIC :: lsq_stencil_create
PUBLIC :: lsq_compute_coeff_cell
PUBLIC :: scalar_int_coeff
PUBLIC :: bln_int_coeff_e2c

CONTAINS

!-------------------------------------------------------------------------
!
!
!>
!! This routine initializes the indices used to define the stencil.
!!
!! This routine initializes the indices used to define the stencil
!! of the lsq reconstruction. The stencil is cell based and includes
!! a variable number of cells (lsq_dim_c) around each control volume
!! (currently 3 or 9)
!!
!! @par Revision History
!! Developed and tested by Daniel Reinert (2009-11-11)
!! Modification by Daniel Reinert, DWD (2010-09-06)
!! - added 12-point stencil
!! Modification by Rainer Johanni 2010-10-26: used for one patch at a time only
!!
!!
SUBROUTINE lsq_stencil_create( ptr_patch, ptr_int_lsq, lsq_dim_c)
!
TYPE(t_patch), INTENT(IN) :: ptr_patch

TYPE(t_lsq), INTENT(INOUT) :: ptr_int_lsq

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER :: ilc, ibc                 ! line and block index
INTEGER :: ilc_n(3), ibc_n(3)       ! line and block index for neighbors of
                                    ! direct neighbors
INTEGER :: ilv(3), ibv(3)           ! vertex line and block indices
INTEGER :: ilc_v(3,6), ibc_v(3,6)   ! cell line and block indices
                                    ! around each of the three vertices
INTEGER :: jb                       ! loop index blocks
INTEGER :: jc                       ! loop index cells
INTEGER :: jj                       ! loop index
INTEGER :: jec                      ! loop index
INTEGER :: jtri                     ! loop index
INTEGER :: cnt                      ! counter
INTEGER :: nblks_c
INTEGER :: i_startblk               ! start block
INTEGER :: i_startidx               ! start index
INTEGER :: i_endidx                 ! end index
INTEGER :: i_rlstart                ! refinement control start level

REAL(wp) :: z_stencil(UBOUND(ptr_int_lsq%lsq_dim_stencil,1),UBOUND(ptr_int_lsq%lsq_dim_stencil,2))

!--------------------------------------------------------------------

  CALL message('mo_interpolation:lsq_stencil_create', '')

  i_rlstart = 2

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rlstart,1)

!$OMP PARALLEL
  IF ( lsq_dim_c == ptr_patch%geometry_info%cell_type ) THEN
    ! The stencil consists of 3 cells surrounding the control volume
    ! i.e. the direct neighbors are taken.

    ! the cell and block indices are just copied from ptr_patch%cells%neighbor_idx
    ! and ptr_patch%cells%neighbor_blk
    CALL copy(ptr_patch%cells%neighbor_idx(:,:,:), &
         ptr_int_lsq%lsq_idx_c(:,:,:))
    CALL copy(ptr_patch%cells%neighbor_blk(:,:,:), &
         ptr_int_lsq%lsq_blk_c(:,:,:))
    CALL copy(ptr_patch%cells%num_edges(:,:), &
         ptr_int_lsq%lsq_dim_stencil(:,:))
!$OMP BARRIER

  ELSE IF (lsq_dim_c == 9) THEN

    !
    ! The stencil consists of 9 cells surrounding the control volume. The 3 direct
    ! neighbors and the neighbors of the direct neighbors are taken.
    !
!$OMP DO PRIVATE(jb,jc,jec,jj,i_startidx,i_endidx,cnt,ilc,ibc,ilc_n,&
!$OMP ibc_n) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,     &
                         i_startidx, i_endidx, i_rlstart)

      DO jc = i_startidx, i_endidx

        IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        cnt = 1

        DO jec = 1, 3
          ilc = ptr_patch%cells%neighbor_idx(jc,jb,jec)
          ibc = ptr_patch%cells%neighbor_blk(jc,jb,jec)

          ! direct neighbors
          ptr_int_lsq%lsq_idx_c(jc,jb,cnt) = ilc
          ptr_int_lsq%lsq_blk_c(jc,jb,cnt) = ibc

          cnt = cnt + 1

          ! neighbors of direct neighbors
          DO jj = 1,3
            ilc_n(jj) = ptr_patch%cells%neighbor_idx(ilc,ibc,jj)
            ibc_n(jj) = ptr_patch%cells%neighbor_blk(ilc,ibc,jj)

            IF (ilc_n(jj) /= jc .or. ibc_n(jj) /= jb) THEN
              ptr_int_lsq%lsq_idx_c(jc,jb,cnt) = ilc_n(jj)
              ptr_int_lsq%lsq_blk_c(jc,jb,cnt) = ibc_n(jj)
              cnt = cnt + 1
            ENDIF
          ENDDO

        ENDDO ! jec loop

        ptr_int_lsq%lsq_dim_stencil(jc,jb) = cnt - 1

      ENDDO ! loop over cells

    ENDDO ! loop over blocks
!$OMP END DO

  ELSE

    !
    ! The stencil consists of 12 cells surrounding the control volume.
    ! This one is similar to the 9-point stencil except that it
    ! is more isotropic. The stencil includes the 3 direct neighbors,
    ! the neighbors of the direct neighbors and 3 additional cells
    ! which share 1 vertex with the CV under consideration.
    !
    ! Note: At pentagon points the size of the stencil reduces to 11.
    !
!$OMP DO PRIVATE(jb,jc,jec,jj,jtri,i_startidx,i_endidx,cnt,ilv,ibv, &
!$OMP            ilc_v,ibc_v,ilc_n,ibc_n) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,     &
                         i_startidx, i_endidx, i_rlstart)

      DO jc = i_startidx, i_endidx

        IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

        cnt = 1

        ! get get line and block indices of cell vertices
        ilv(1:3) = ptr_patch%cells%vertex_idx(jc,jb,1:3)
        ibv(1:3) = ptr_patch%cells%vertex_blk(jc,jb,1:3)

        ! for each vertex: get all the cells which share this vertex
        DO jj = 1,3
          ilc_v(jj,:)=ptr_patch%verts%cell_idx(ilv(jj),ibv(jj),:)
          ibc_v(jj,:)=ptr_patch%verts%cell_blk(ilv(jj),ibv(jj),:)
        ENDDO

        !
        ! 1. add the 3 direct neighbors to the stencil
        !
        DO jec = 1, 3
          ! get line and block indices of direct neighbors
          ilc_n(jec) = ptr_patch%cells%neighbor_idx(jc,jb,jec)
          ibc_n(jec) = ptr_patch%cells%neighbor_blk(jc,jb,jec)

          ptr_int_lsq%lsq_idx_c(jc,jb,cnt) = ilc_n(jec)
          ptr_int_lsq%lsq_blk_c(jc,jb,cnt) = ibc_n(jec)

          cnt = cnt + 1
        ENDDO

        !
        ! 2. loop over the vertices and add all the cells
        !    that are no direct neighbors and not our CV.
        !
        DO jj = 1,3   ! loop over vertices
          DO jtri=1,6 ! loop over cells around each vertex

            IF (.NOT.( (ilc_v(jj,jtri) == ilc_n(1) .AND. ibc_v(jj,jtri) == ibc_n(1))  &
              &  .OR.  (ilc_v(jj,jtri) == ilc_n(2) .AND. ibc_v(jj,jtri) == ibc_n(2))  &
              &  .OR.  (ilc_v(jj,jtri) == ilc_n(3) .AND. ibc_v(jj,jtri) == ibc_n(3))  &
              &  .OR.  (ilc_v(jj,jtri) == jc       .AND. ibc_v(jj,jtri) == jb)        &
              &  .OR.  (ilc_v(jj,jtri) == 0        .AND. ibc_v(jj,jtri) == 0 ) ) ) THEN

              ptr_int_lsq%lsq_idx_c(jc,jb,cnt) = ilc_v(jj,jtri)
              ptr_int_lsq%lsq_blk_c(jc,jb,cnt) = ibc_v(jj,jtri)

              cnt = cnt + 1
            ENDIF
          ENDDO
        ENDDO

        ptr_int_lsq%lsq_dim_stencil(jc,jb) = cnt - 1

      ENDDO ! loop over cells

    ENDDO ! loop over blocks
!$OMP END DO NOWAIT

  ENDIF
!$OMP END PARALLEL

  DO cnt = 1, lsq_dim_c
    CALL sync_idx(SYNC_C, SYNC_C, ptr_patch, ptr_int_lsq%lsq_idx_c(:,:,cnt), &
                                           & ptr_int_lsq%lsq_blk_c(:,:,cnt))
  ENDDO

  z_stencil(:,:) = REAL(ptr_int_lsq%lsq_dim_stencil(:,:),wp)
  CALL sync_patch_array(SYNC_C,ptr_patch,z_stencil)
  ptr_int_lsq%lsq_dim_stencil(:,:) = NINT(z_stencil(:,:))


END SUBROUTINE lsq_stencil_create

!-------------------------------------------------------------------------
!
!
!>
!! This routine bifurcates into lsq_compute_coeff_cell based on geometry type
!!
!! Anurag Dipankar, MPIM (2013-04)
!!
SUBROUTINE lsq_compute_coeff_cell( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
  &                                lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )
!

!
TYPE(t_patch), INTENT(IN) ::  ptr_patch

TYPE(t_lsq), TARGET, INTENT(INOUT) ::  ptr_int_lsq

LOGICAL, INTENT(IN) ::   &  ! flag determining whether the least
  &  llsq_rec_consv         ! squares reconstruction should be conservative

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER, INTENT(IN)  ::  &  ! parameter determining the dimension of the solution
  &  lsq_dim_unk

INTEGER, INTENT(IN)  ::  &  ! least squares weighting exponent
  &  lsq_wgt_exp

 !Local variable
 CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs_lsq_bln:lsq_compute_coeff_cell'

    !
    SELECT CASE(ptr_patch%geometry_info%geometry_type)

    CASE (planar_torus_geometry)

     CALL lsq_compute_coeff_cell_torus( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
                                        lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )

    CASE (sphere_geometry)

     CALL lsq_compute_coeff_cell_sphere( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
                                         lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )

    CASE DEFAULT

      CALL finish(method_name, "Undefined geometry type")

    END SELECT


END SUBROUTINE lsq_compute_coeff_cell


!-------------------------------------------------------------------------
!
!
!>
!! AD: This routine has been just renamed with affix "_sphere"
!!
!! This routine computes the coefficients needed for a weighted least-squares.
!!
!! This routine computes the coefficients needed for a weighted least-squares
!! reconstruction at cell centers. Optionally, the reconstruction can be
!! enforced to be conservative in the sense that, when integrated over the
!! control volume, it recovers the area average stored at the mass point.
!! Works for triangular and hexagonal control volumes.
!!
!! @par Revision History
!! Developed and tested by Daniel Reinert (2009-09-29)
!! Modification by Daniel Reinert, DWD (2009-11-02)
!! - application of gnomonic projection for calculation of distance
!!   vectors between points. Replaces call of rotate_latlon
!! Modification by Daniel Reinert, DWD (2009-11-11)
!! - generalization to arbitrary order of reconstruction (yet 2nd or 3rd order)
!!
SUBROUTINE lsq_compute_coeff_cell_sphere( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
  &                                       lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )
!

!
TYPE(t_patch), INTENT(IN) ::  ptr_patch

TYPE(t_lsq), TARGET, INTENT(INOUT) ::  ptr_int_lsq

LOGICAL, INTENT(IN) ::   &  ! flag determining whether the least
  &  llsq_rec_consv         ! squares reconstruction should be conservative

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER, INTENT(IN)  ::  &  ! parameter determining the dimension of the solution
  &  lsq_dim_unk

INTEGER, INTENT(IN)  ::  &  ! least squares weighting exponent
  &  lsq_wgt_exp

REAL(wp), DIMENSION(lsq_dim_c,2) ::  &      ! geographical coordinates of all cell centers
  & xytemp_c                                ! in the stencil

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type,2)   ::  &  ! geogr. coordinates of vertices of the
  & xytemp_v                                ! control volume

REAL(wp), ALLOCATABLE,DIMENSION(:,:,:,:) ::  &
  & z_dist_g                                ! distance vectors to neighbouring cell
                                            ! centers stored for each cell

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type,2)   ::  &  ! lat/lon distance vector edge midpoint -> cvertex
  & distxy_v

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) :: dely, delx ! difference in latitude and longitude between
                                               ! vertices

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk) ::  & ! lsq matrix
  & z_lsq_mat_c

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk)   ::  &  ! Q matrix of QR-factorization
  & z_qmat

REAL(wp), DIMENSION(nproma,lsq_dim_unk,lsq_dim_unk) ::  &  ! R matrix of QR-factorization
  & z_rmat

REAL(wp) :: z_rcarea                   ! reciprocal of cell area

REAL(wp) :: xloc, yloc                 ! geographical coordinates of
                                       ! point under consideration

REAL(wp) :: z_norm                     ! vector length (distance between control volume
                                       ! center and cell centers in the stencil on tangent
                                       ! plane) (also used for normalization)

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) ::  & ! integrand for each edge
  & fx, fy, fxx, fyy, fxy,           & ! for analytical calculation of moments
  & fxxx, fyyy, fxxy, fxyy

INTEGER, POINTER  ::       &           ! pointer to stencil size (cell dependent)
  & ptr_ncells(:,:)

INTEGER, DIMENSION(lsq_dim_c) ::  &    ! line and block indices of cells in the stencil
  & ilc_s, ibc_s
INTEGER, DIMENSION(ptr_patch%geometry_info%cell_type) :: jlv, jbv      ! line and block indices of vertex
INTEGER :: cnt                         ! counter
INTEGER :: jrow                        ! matrix row-identifier
INTEGER :: nel                         ! number of matrix elements
INTEGER :: nblks_c
INTEGER :: pid                         ! patch ID
INTEGER :: jb                          ! index of current block
INTEGER :: jc                          ! index of current cell
INTEGER :: js                          ! index of current control volume in the stencil
INTEGER :: ju                          ! loop index for column of lsq matrix
INTEGER :: jec                         ! loop index for cell's edge
INTEGER :: i_startblk                  ! start block
INTEGER :: i_startidx                  ! start index
INTEGER :: i_endidx                    ! end index
INTEGER :: i_rcstartlev                ! refinement control start level
INTEGER :: ist, icheck                 ! status
INTEGER :: nverts
INTEGER :: jecp
INTEGER :: jja, jjb, jjk               ! loop indices for Moore-Penrose inverse

REAL(wp) ::   &                        ! singular values of lsq design matrix A
  &  zs(lsq_dim_unk,nproma)            ! min(lsq_dim_c,lsq_dim_unk)

REAL(wp) ::   &                        ! U matrix of SVD. Columns of U are the left
  &  zu  (lsq_dim_c,lsq_dim_c,nproma)  ! singular vectors of A

REAL(wp) ::   &                        ! TRANSPOSE of V matrix of SVD. Columns of V are
  &  zv_t(lsq_dim_unk,lsq_dim_unk,nproma) ! the right singular vectors of A.


INTEGER, PARAMETER  :: &     ! size of work array for SVD lapack routine
  &  lwork=10000
REAL(wp) ::   &              ! work array for SVD lapack routine
  &  zwork(lwork)
INTEGER  ::   &              ! work array for SVD lapack routine
  & ziwork(8*min(lsq_dim_c,lsq_dim_unk))


!DR for DEBUG purposes
! #ifdef DEBUG_COEFF LL it's used in openmp directives,
REAL(wp) :: za_debug(nproma,lsq_dim_c,lsq_dim_unk)
! #endif
!--------------------------------------------------------------------


  CALL message('mo_interpolation:lsq_compute_coeff_cell_sphere', '')

  i_rcstartlev = 2

  ! get patch id
  pid = ptr_patch%id

  ! stencil size
  ptr_ncells => ptr_int_lsq%lsq_dim_stencil(:,:)

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

  ! allocate array in which the distance vectors between the
  ! cell center of the control volume and the cell centers of the
  ! neighboring control volumes are stored.
  ALLOCATE (z_dist_g(nproma,nblks_c,lsq_dim_c,2), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:lsq_compute_coeff_cell_sphere',   &
      &             'allocation for z_dist_g failed')
  ENDIF


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jc,js,jec,i_startidx,i_endidx,jlv,jbv,ilc_s,ibc_s, &
!!$OMP            xloc,yloc,xytemp_c,xytemp_v,z_norm,distxy_v,z_rcarea, &
!!$OMP            delx,dely,fx,fy,fxx,fyy,fxy,jecp,nverts) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,     &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each cell, calculate weights, moments, matrix coefficients
    ! and QR decomposition of normal equation matrix
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      IF (ptr_patch%geometry_info%cell_type == 3 )THEN
        nverts  = 3
      ELSE
        nverts = ptr_patch%cells%num_edges(jc,jb)
      ENDIF
    !
    ! Gather some information about the control volume and
    ! all the cells in the stencil.
    !
    ! get line and block indices of edge vertices
      jlv(1:nverts) = ptr_patch%cells%vertex_idx(jc,jb,1:nverts)
      jbv(1:nverts) = ptr_patch%cells%vertex_blk(jc,jb,1:nverts)

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    !
    ! 1. Get geographical coordinates of the control volume center
    !    and all cell centers in the stencil. In addition, get geographical
    !    coordinates of the vertices of the control volume.
    !
    ! get geographical coordinates of control volume center
      xloc = ptr_patch%cells%center(jc,jb)%lon
      yloc = ptr_patch%cells%center(jc,jb)%lat

    ! get geogr. coordinates of all cell centers in the stencil
      DO js = 1, ptr_ncells(jc,jb)
        xytemp_c(js,1) = ptr_patch%cells%center(ilc_s(js),ibc_s(js))%lon
        xytemp_c(js,2) = ptr_patch%cells%center(ilc_s(js),ibc_s(js))%lat
      ENDDO

    ! get geogr. coordinates of edge-vertices (only for control volume)
      DO jec = 1, nverts
        xytemp_v(jec,1) = ptr_patch%verts%vertex(jlv(jec),jbv(jec))%lon
        xytemp_v(jec,2) = ptr_patch%verts%vertex(jlv(jec),jbv(jec))%lat
      ENDDO


    !
    ! 2. Now, all points (centers, vertices) are projected onto a tangent
    !    plane having its center at the cell center of the control volume.
    !    At the same time distance vectors (r_x,r_y) are calculated between
    !    the control volume center and all other mass-points in the stencil.
    !    From these the lsq weights can be deduced.

    !
    ! a: Project centers and calculate distance vectors between control volume
    !    center and all cell centers in the stencil in geographical coordinates.
    !    Since the control volume center has the local coordinates (x,y)=(0,0),
    !    the coordinates of the other cell centers equal the distance
    !    vectors.
      DO js = 1, ptr_ncells(jc,jb)
        CALL gnomonic_proj( xloc, yloc, xytemp_c(js,1), xytemp_c(js,2),  &! in
         &                  z_dist_g(jc,jb,js,1), z_dist_g(jc,jb,js,2) )  ! out

      ENDDO
      ! multiply with earth radius and store
      z_dist_g(jc,jb,1:ptr_ncells(jc,jb),:) = grid_sphere_radius * &
        & z_dist_g(jc,jb,1:ptr_ncells(jc,jb),:)



    !
    ! b: compute normalized weights for weighted least squares system
    !    The closest cell circumcenter is assigned weight of w=1.
    !
      DO js = 1, ptr_ncells(jc,jb)
        z_norm = SQRT(DOT_PRODUCT(z_dist_g(jc,jb,js,1:2),z_dist_g(jc,jb,js,1:2)))

        !
        ! weights for weighted least squares system
        !
        ptr_int_lsq%lsq_weights_c(jc,js,jb)= 1._wp/(z_norm**lsq_wgt_exp)

      ENDDO
      !
      ! Normalization
      !
      ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)=                &
        &     ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)         &
        &   / MAXVAL(ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb))



    ! 3. the following part (including calculation of moments) will only
    !    be called, if a conservative least squares reconstruction
    !    is chosen. Otherwise all moments will be equal to zero and the
    !    reconstruction simplifies to the standard non-conservative
    !    reconstruction.
      IF (llsq_rec_consv) THEN
      !
      ! a: Project control volume vertices and calculate distance vectors
      !    between cell center and vertices
      !
        DO jec=1,nverts
          CALL gnomonic_proj( xloc, yloc, xytemp_v(jec,1), xytemp_v(jec,2),  &
           &                 distxy_v(jec,1), distxy_v(jec,2) )

        ENDDO
      ! multiply with earth radius
        distxy_v(1:nverts,1:2) = grid_sphere_radius * distxy_v(1:nverts,1:2)


      !
      ! b: calculate moments for given cell
      !    (calculated analytically; see Lauritzen CSLAM 09 for first 5 moments)
      !
      ! !DR: Those moments have been re-rechecked, using an alternative,
      !      quadrature-based formulation. Results have been identical up to
      !      roundoff-errors. Similarly the hat-moments have been checked. The
      !      inconsistency caused by the different projections involved do not
      !      seem to negatively effect the results.

      ! Storage docu for x^ny^m:
      ! lsq_moments(:,:,1) : x^1y^0
      ! lsq_moments(:,:,2) : x^0y^1
      ! lsq_moments(:,:,3) : x^2y^0
      ! lsq_moments(:,:,4) : x^0y^2
      ! lsq_moments(:,:,5) : x^1y^1
      ! lsq_moments(:,:,6) : x^3y^0
      ! lsq_moments(:,:,7) : x^0y^3
      ! lsq_moments(:,:,8) : x^2y^1
      ! lsq_moments(:,:,9) : x^1y^2
      !

        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! note that the distance vector distxy_v between each vertex and
          ! the center of the tangent plane are identical to the coordinates of
          ! each vertex on the tangent plane. Thus the distances between the
          ! vertices in x and y direction can be derived as follows:
          !
          ! longitudinal-distance between vertices on tangent plane
          delx(jec) = distxy_v(jecp,1) - distxy_v(jec,1)

          ! latitudinal-distance between vertices on tangent plane
          dely(jec) = distxy_v(jecp,2) - distxy_v(jec,2)


          !
          ! analytic moment calculation
          !
          ! 0: control volume area (reciprocal value)
          fx(jec) = distxy_v(jecp,1) + distxy_v(jec,1)

        ENDDO

        z_rcarea = 2._wp/DOT_PRODUCT(fx(1:nverts),dely(1:nverts))


        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! I. x^1y^0
          fx(jec) = distxy_v(jecp,1)**2              &
          &       + distxy_v(jecp,1)*distxy_v(jec,1) &
          &       + distxy_v(jec ,1)**2

          ! II. x^0y^1
          fy(jec) = distxy_v(jecp,2)**2              &
          &       + distxy_v(jecp,2)*distxy_v(jec,2) &
          &       + distxy_v(jec,2)**2

          IF ( lsq_dim_unk > 2 ) THEN

            ! III. x^2y^0
            fxx(jec) = (distxy_v(jecp,1)    + distxy_v(jec,1)       ) &
            &        * (distxy_v(jecp,1)**2 + distxy_v(jec,1)**2)

            ! IV. x^0y^2
            fyy(jec) = (distxy_v(jecp,2)    + distxy_v(jec,2)       ) &
            &        * (distxy_v(jecp,2)**2 + distxy_v(jec,2)**2)

            ! V. x^1y^1
            fxy(jec) = distxy_v(jecp,2) * (3._wp*distxy_v(jecp,1)**2                         &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + distxy_v(jec,1)**2 )     &
            &        + distxy_v(jec,2) * ( distxy_v(jecp,1)**2                               &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + 3._wp*distxy_v(jec,1)**2 )

          ENDIF ! lsq_dim_unk > 2

          IF ( lsq_dim_unk > 5 ) THEN

            ! VI.  x^3y^0
            fxxx(jec) = 5._wp*distxy_v(jec,1)**4                   &
              &       + 10._wp*distxy_v(jec,1)**3 * delx(jec)      &
              &       + 10._wp*distxy_v(jec,1)**2 * delx(jec)**2   &
              &       + 5._wp *distxy_v(jec,1)    * delx(jec)**3   &
              &       + delx(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fxxx(jec) = distxy_v(jecp,1)**4                       &
            !DR    &       + distxy_v(jec,1) * distxy_v(jecp,1)**3     &
            !DR    &       + distxy_v(jec,1)**2 * distxy_v(jecp,1)**2  &
            !DR    &       + distxy_v(jec,1)**3 * distxy_v(jecp,1)     &
            !DR    &       + distxy_v(jec,1)**4


            ! VII. x^0y^3
            fyyy(jec) = 5._wp*distxy_v(jec,2)**4                   &
              &       + 10._wp*distxy_v(jec,2)**3 * dely(jec)      &
              &       + 10._wp*distxy_v(jec,2)**2 * dely(jec)**2   &
              &       + 5._wp *distxy_v(jec,2)    * dely(jec)**3   &
              &       + dely(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fyyy(jec) = distxy_v(jecp,2)**4                       &
            !DR    &       + distxy_v(jec,2)    * distxy_v(jecp,2)**3  &
            !DR    &       + distxy_v(jec,2)**2 * distxy_v(jecp,2)**2  &
            !DR    &       + distxy_v(jec,2)**3 * distxy_v(jecp,2)     &
            !DR    &       + distxy_v(jec,2)**4

          ENDIF ! lsq_dim_unk > 5

          IF ( lsq_dim_unk > 7 ) THEN

            ! VIII. x^2y^1
            fxxy(jec) = 4._wp*distxy_v(jecp,1)**3 *distxy_v(jecp,2)                    &
              &       + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jecp,2)   &
              &       + 2._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jecp,2) &
              &       +       distxy_v(jec,1)**3 * distxy_v(jecp,2)                    &
              &       +       distxy_v(jecp,1)**3 * distxy_v(jec,2)                    &
              &       + 2._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jec,2)    &
              &       + 3._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jec,2)  &
              &       + 4._wp*distxy_v(jec,1)**3 * distxy_v(jec,2)

            ! IX. x^1y^2
            fxyy(jec) = 6._wp*distxy_v(jecp,1)**2 * distxy_v(jecp,2)**2                     &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jecp,2)**2              &
              &   +       distxy_v(jec,1)**2 * distxy_v(jecp,2)**2                          &
              &   + 3._wp*distxy_v(jecp,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)            &
              &   + 4._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)*distxy_v(jecp,2) &
              &   + 3._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)             &
              &   +       distxy_v(jecp,1)**2 * distxy_v(jec,2)**2                          &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)**2               &
              &   + 6._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)**2

          ENDIF ! lsq_dim_unk > 7


        ENDDO ! loop over nverts

        ptr_int_lsq%lsq_moments(jc,jb,1)= z_rcarea/6._wp*DOT_PRODUCT(fx(1:nverts),dely(1:nverts))
        ptr_int_lsq%lsq_moments(jc,jb,2)=-z_rcarea/6._wp*DOT_PRODUCT(fy(1:nverts),delx(1:nverts))

        IF ( lsq_dim_unk > 2 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,3) =  z_rcarea/12._wp * &
            & DOT_PRODUCT(fxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,4) = -z_rcarea/12._wp * &
            & DOT_PRODUCT(fyy(1:nverts),delx(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,5) =  z_rcarea/24._wp * &
            & DOT_PRODUCT(fxy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 2


        IF ( lsq_dim_unk > 5 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,6) =  z_rcarea/20._wp * &
            & DOT_PRODUCT(fxxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,7) = -z_rcarea/20._wp * &
            & DOT_PRODUCT(fyyy(1:nverts),delx(1:nverts))

        END IF  ! lsq_dim_unk > 5


        IF ( lsq_dim_unk > 7 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,8) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxxy(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,9) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxyy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 7

      END IF  ! llsq_rec_consv

    END DO  ! loop over cells

  END DO  ! loop over blocks
!!$OMP END DO NOWAIT
!! For unknown reasons, closing the parallel section here is needed to get the above
!! loop parallelized.
!!$OMP END PARALLEL

  DO jb = 1, lsq_dim_c
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_weights_c(:,jb,:))
  ENDDO
  DO jb = 1, lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments(:,:,jb))
  ENDDO


!$OMP PARALLEL PRIVATE(jb,jc,js,ju,jja,jjb,jjk,i_startidx,i_endidx,ilc_s,ibc_s, &
!$OMP            z_lsq_mat_c,zs,zu,zv_t,zwork,ziwork,ist,icheck,za_debug, &
!$OMP            z_qmat,z_rmat,cnt,jrow,nel)
!$OMP DO ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! 4. for each cell, calculate LSQ design matrix A
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
        ! Take care that z_lsq_mat_c isn't singular
        z_lsq_mat_c(jc,:,:) = 0.0_wp
        DO js = 1, MIN(lsq_dim_unk, lsq_dim_c)
          z_lsq_mat_c(jc,js,js) = 1.0_wp
        ENDDO
        CYCLE
      ENDIF

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    ! Calculate full moments lsq_moments_hat(ilc_s(js),ibc_s(js),ju)
    !
    ! Storage docu for x^ny^m:
    ! lsq_moments_hat(:,:,:,1) : \hat{x^1y^0}
    ! lsq_moments_hat(:,:,:,2) : \hat{x^0y^1}
    ! lsq_moments_hat(:,:,:,3) : \hat{x^2y^0}
    ! lsq_moments_hat(:,:,:,4) : \hat{x^0y^2}
    ! lsq_moments_hat(:,:,:,5) : \hat{x^1y^1}
    ! lsq_moments_hat(:,:,:,6) : \hat{x^3y^0}
    ! lsq_moments_hat(:,:,:,7) : \hat{x^0y^3}
    ! lsq_moments_hat(:,:,:,8) : \hat{x^2y^1}
    ! lsq_moments_hat(:,:,:,9) : \hat{x^1y^2}
    !
      DO js = 1, ptr_ncells(jc,jb)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,1) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1) + z_dist_g(jc,jb,js,1)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,2) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2) + z_dist_g(jc,jb,js,2)

        IF (lsq_dim_unk > 2) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,3) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
           &    + z_dist_g(jc,jb,js,1)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,4) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2) &
           &    + z_dist_g(jc,jb,js,2)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,5) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)                              &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)        &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)        &
           &    + z_dist_g(jc,jb,js,1) * z_dist_g(jc,jb,js,2)
        ENDIF

        IF ( lsq_dim_unk > 5 ) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,6) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),6)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,1) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)                       &
            &    * z_dist_g(jc,jb,js,1)**2                                                    &
            &    + z_dist_g(jc,jb,js,1)**3

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,7) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),7)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,2) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)                       &
            &    * z_dist_g(jc,jb,js,2)**2                                                    &
            &    + z_dist_g(jc,jb,js,2)**3
        ENDIF

        IF ( lsq_dim_unk > 7 ) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,8) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),8)                             &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,2)       &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,1) &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
            &    * z_dist_g(jc,jb,js,2)                                                       &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)**2    &
            &    + z_dist_g(jc,jb,js,1)**2 * z_dist_g(jc,jb,js,2)


          ptr_int_lsq%lsq_moments_hat(jc,jb,js,9) =                                          &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),9)                            &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,1)      &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,2)&
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2)&
            &    * z_dist_g(jc,jb,js,1)                                                      &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)**2   &
            &    + z_dist_g(jc,jb,js,2)**2 * z_dist_g(jc,jb,js,1)
        ENDIF
      ENDDO


      ! loop over rows of lsq design matrix (all cells in the stencil)
      DO js = 1, ptr_ncells(jc,jb)
        ! loop over columns of lsq design matrix (number of unknowns)
        DO ju = 1,lsq_dim_unk

          z_lsq_mat_c(jc,js,ju) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                         &
           &    * (ptr_int_lsq%lsq_moments_hat(jc,jb,js,ju)-ptr_int_lsq%lsq_moments(jc,jb,ju))

        END DO
       END DO
       IF(ptr_ncells(jc,jb) < lsq_dim_c) THEN
         z_lsq_mat_c(jc,lsq_dim_c,:) = 0.0_wp
       ENDIF

    ENDDO   ! loop over cells



    !
    ! compute QR decomposition and Singular Value Decomposition (SVD)
    ! of least squares design matrix A. For the time being both methods are
    ! retained.

    !
    ! 5a. QR-factorization of design matrix A
    !
    IF (.NOT. advection_config(pid)%llsq_svd) THEN
!CDIR NOIEXPAND
    CALL qrdec(lsq_dim_c, lsq_dim_unk, i_startidx, & ! in
     &         i_endidx, z_lsq_mat_c,              & ! in
     &         z_qmat, z_rmat)                       ! out


    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! 7. Save transposed Q-Matrix
      ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,1:lsq_dim_c,jb)  =  &
       &      TRANSPOSE(z_qmat(jc,1:lsq_dim_c,1:lsq_dim_unk))

      ! 8. Save R-Matrix
      !
      ! a. Save reciprocal values of the diagonal elements
      !
      DO ju = 1,lsq_dim_unk
        ptr_int_lsq%lsq_rmat_rdiag_c(jc,ju,jb) = 1._wp/z_rmat(jc,ju,ju)
      ENDDO

      !
      ! b. Save upper triangular elements without diagonal elements in a 1D-array
      !    (starting from the bottom right)
      !
      cnt = 1
      DO jrow = lsq_dim_unk-1,1,-1
        ! number of elements to store
        nel = lsq_dim_unk - jrow
        ptr_int_lsq%lsq_rmat_utri_c(jc,cnt:cnt+nel-1,jb) = z_rmat(jc,jrow,jrow+1:lsq_dim_unk)
        cnt = cnt + nel
      ENDDO


      ! Multiply ith column of the transposed Q-matrix (corresponds to the
      ! different members of the stencil) with the ith weight. This avoids
      ! multiplication of the RHS of the LSQ-System with this weight during
      ! runtime.
      DO js = 1,lsq_dim_c
        ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  =                  &
          &                ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  &
          &              * ptr_int_lsq%lsq_weights_c(jc,js,jb)
      ENDDO

    END DO  ! loop over cells

    ELSE   ! llsq_svd=.TRUE.

    !
    ! 5b. Singular value decomposition of lsq design matrix A
    ! !!! does not vectorize, unfortunately !!!
    ist = 0
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! A = U * SIGMA * transpose(V)
      !
      ! z_lsq_mat_c : M x N least squares design matrix A            (IN)
      ! zu          : M x M orthogonal matrix U                      (OUT)
      ! zv_t        : N x N orthogonal matrix V                      (OUT)
      ! zs          : min(M,N) Singular values of A                  (OUT)
      ! zwork       : workspace(1,LWORK)                             (OUT)
      ! lwork       : 3*min(M,N)                                     (IN)
      !              + max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N))  (IN)
      ! iwork       : workspace(8*min(M,N))                          (IN)

      CALL DGESDD('A',                 & !in
        &         lsq_dim_c,           & !in
        &         lsq_dim_unk,         & !in
        &         z_lsq_mat_c(jc,:,:), & !inout Note: destroyed on output
        &         lsq_dim_c,           & !in
        &         zs(:,jc),            & !out
        &         zu(:,:,jc),          & !out
        &         lsq_dim_c,           & !in
        &         zv_t(:,:,jc),        & !out
        &         lsq_dim_unk,         & !in
        &         zwork,               & !out
        &         lwork,               & !in
        &         ziwork,              & !inout
        &         icheck               ) !out
      ist = ist + icheck
    ENDDO
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:lsq_compute_coeff_cell_sphere',   &
        &             'singular value decomposition failed')
    ENDIF

    ! compute Moore-Penrose inverse
    ! INVERSE(A):: V * INVERSE(SIGMA) * TRANSPOSE(U) and store
    ! note that the ith column is multiplied with the ith weight
    ! in order to avoid the weighting of the r.h.s. during runtime.
    DO jja = 1, lsq_dim_unk
      DO jjb = 1, lsq_dim_c
        DO jjk = 1, lsq_dim_unk
          DO jc = i_startidx, i_endidx
            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
            ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb) =            &
              &  ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb)         &
              &  + zv_t(jjk,jja,jc) /zs(jjk,jc) * zu(jjb,jjk,jc)  &
              &  * ptr_int_lsq%lsq_weights_c(jc,jjb,jb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!!!!!DEBUG !!!
#ifdef DEBUG_COEFF
    za_debug(:,:,:)  = 0._wp
    ! re-COMPUTE A:: U * SIGMA * TRANSPOSE(V)  !!! Funktioniert
    DO jja = 1, lsq_dim_c
      DO jjb = 1, lsq_dim_unk
        DO jjk = 1, lsq_dim_unk  !lsq_dim_c
          DO jc = i_startidx, i_endidx
            za_debug(jc,jja,jjb) = za_debug(jc,jja,jjb)  &
              &  + zu(jja,jjk,jc) * zs(jjk,jc) * zv_t(jjk,jjb,jc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    write(0,*) "za_debug(10,:,:)- z_lsq_mat_c(10,:,:)",za_debug(10,:,:)- z_lsq_mat_c(10,:,:)
    write(0,*) "za_debug(55,:,:)- z_lsq_mat_c(55,:,:)",za_debug(55,:,:)- z_lsq_mat_c(55,:,:)
    write(0,*) "zs(:,10) ",zs(:,10)
    write(0,*) "zs(:,55) ",zs(:,55)
    write(0,*) "ptr_int_lsq%lsq_weights_c(10,:,jb)",ptr_int_lsq%lsq_weights_c(10,:,jb)
    write(0,*) "ptr_int_lsq%lsq_weights_c(55,:,jb)",ptr_int_lsq%lsq_weights_c(55,:,jb)
#endif
!!!! END DEBUG !!!

    ENDIF  ! llsq_svd

  END DO  ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  DO ju = 1, lsq_dim_unk
    DO jc = 1, lsq_dim_c
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments_hat(:,:,jc,ju))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_pseudoinv(:,ju,jc,:))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_qtmat_c(:,ju,jc,:))
    ENDDO
  ENDDO

  DO jc = 1, UBOUND(ptr_int_lsq%lsq_rmat_utri_c, 2)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_utri_c(:,jc,:))
  ENDDO
  DO ju = 1,lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_rdiag_c(:,ju,:))
  ENDDO

  DEALLOCATE (z_dist_g, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:lsq_compute_coeff_cell_sphere',   &
      &             'deallocation for z_dist_g failed')
  ENDIF


END SUBROUTINE lsq_compute_coeff_cell_sphere


!-------------------------------------------------------------------------
!
!! This is same routine as lsq_compute_coeff_cell_sphere just modified for
!! flat geometry
!>
!!
!! @par Revision History
!! Initial version by Daniel Reinert for sphere geometry modified
!! for torus geometry by Anurag Dipankar, MPIM (2013-04)
!!
SUBROUTINE lsq_compute_coeff_cell_torus( ptr_patch, ptr_int_lsq, llsq_rec_consv, &
  &                                      lsq_dim_c, lsq_dim_unk, lsq_wgt_exp )
!

!
TYPE(t_patch), INTENT(IN) ::  ptr_patch

TYPE(t_lsq), TARGET, INTENT(INOUT) ::  ptr_int_lsq

LOGICAL, INTENT(IN) ::   &  ! flag determining whether the least
  &  llsq_rec_consv         ! squares reconstruction should be conservative

INTEGER, INTENT(IN)  ::  &  ! parameter determining the size of the lsq stencil
  &  lsq_dim_c

INTEGER, INTENT(IN)  ::  &  ! parameter determining the dimension of the solution
  &  lsq_dim_unk

INTEGER, INTENT(IN)  ::  &  ! least squares weighting exponent
  &  lsq_wgt_exp


!CC of points in the stencil
TYPE(t_cartesian_coordinates) :: cc_cv, cc_cell(lsq_dim_c), cc_vert(ptr_patch%geometry_info%cell_type)

REAL(wp), ALLOCATABLE,DIMENSION(:,:,:,:) ::  &
  & z_dist_g                                ! distance vectors to neighbouring cell
                                            ! centers stored for each cell

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type,2)   ::  &  ! lat/lon distance vector edge midpoint -> cvertex
  & distxy_v

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) :: dely, delx ! difference in latitude and longitude between
                                               ! vertices

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk) ::  & ! lsq matrix
  & z_lsq_mat_c

REAL(wp), DIMENSION(nproma,lsq_dim_c,lsq_dim_unk)   ::  &  ! Q matrix of QR-factorization
  & z_qmat

REAL(wp), DIMENSION(nproma,lsq_dim_unk,lsq_dim_unk) ::  &  ! R matrix of QR-factorization
  & z_rmat

REAL(wp) :: z_rcarea                   ! reciprocal of cell area

REAL(wp) :: z_norm                     ! vector length (distance between control volume
                                       ! center and cell centers in the stencil on tangent
                                       ! plane) (also used for normalization)

REAL(wp), DIMENSION(ptr_patch%geometry_info%cell_type) ::  & ! integrand for each edge
  & fx, fy, fxx, fyy, fxy,           & ! for analytical calculation of moments
  & fxxx, fyyy, fxxy, fxyy

INTEGER, POINTER  ::       &           ! pointer to stencil size (cell dependent)
  & ptr_ncells(:,:)

INTEGER, DIMENSION(lsq_dim_c) ::  &    ! line and block indices of cells in the stencil
  & ilc_s, ibc_s
INTEGER, DIMENSION(ptr_patch%geometry_info%cell_type) :: jlv, jbv      ! line and block indices of vertex
INTEGER :: cnt                         ! counter
INTEGER :: jrow                        ! matrix row-identifier
INTEGER :: nel                         ! number of matrix elements
INTEGER :: nblks_c
INTEGER :: pid                         ! patch ID
INTEGER :: jb                          ! index of current block
INTEGER :: jc                          ! index of current cell
INTEGER :: js                          ! index of current control volume in the stencil
INTEGER :: ju                          ! loop index for column of lsq matrix
INTEGER :: jec                         ! loop index for cell's edge
INTEGER :: i_startblk                  ! start block
INTEGER :: i_startidx                  ! start index
INTEGER :: i_endidx                    ! end index
INTEGER :: i_rcstartlev                ! refinement control start level
INTEGER :: ist, icheck                 ! status
INTEGER :: nverts
INTEGER :: jecp
INTEGER :: jja, jjb, jjk               ! loop indices for Moore-Penrose inverse

REAL(wp) ::   &                        ! singular values of lsq design matrix A
  &  zs(lsq_dim_unk,nproma)            ! min(lsq_dim_c,lsq_dim_unk)

REAL(wp) ::   &                        ! U matrix of SVD. Columns of U are the left
  &  zu  (lsq_dim_c,lsq_dim_c,nproma)  ! singular vectors of A

REAL(wp) ::   &                        ! TRANSPOSE of V matrix of SVD. Columns of V are
  &  zv_t(lsq_dim_unk,lsq_dim_unk,nproma) ! the right singular vectors of A.


INTEGER, PARAMETER  :: &     ! size of work array for SVD lapack routine
  &  lwork=10000
REAL(wp) ::   &              ! work array for SVD lapack routine
  &  zwork(lwork)
INTEGER  ::   &              ! work array for SVD lapack routine
  & ziwork(8*min(lsq_dim_c,lsq_dim_unk))


!DR for DEBUG purposes
! #ifdef DEBUG_COEFF LL it's used in openmp directives,
REAL(wp) :: za_debug(nproma,lsq_dim_c,lsq_dim_unk)
! #endif
!--------------------------------------------------------------------


  CALL message('mo_interpolation:lsq_compute_coeff_cell_torus', '')

  i_rcstartlev = 2

  ! get patch id
  pid = ptr_patch%id

  ! stencil size
  ptr_ncells => ptr_int_lsq%lsq_dim_stencil(:,:)

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c

  ! The start block depends on the width of the stencil
  i_startblk = ptr_patch%cells%start_blk(i_rcstartlev,1)

  ! allocate array in which the distance vectors between the
  ! cell center of the control volume and the cell centers of the
  ! neighboring control volumes are stored.
  ALLOCATE (z_dist_g(nproma,nblks_c,lsq_dim_c,2), STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:lsq_compute_coeff_cell_torus',   &
      &             'allocation for z_dist_g failed')
  ENDIF


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jc,js,jec,i_startidx,i_endidx,jlv,jbv,ilc_s,ibc_s, &
!!$OMP            cc_cv,cc_cell,cc_vert,z_norm,distxy_v,z_rcarea, &
!!$OMP            delx,dely,fx,fy,fxx,fyy,fxy,jecp,nverts) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c,     &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! for each cell, calculate weights, moments, matrix coefficients
    ! and QR decomposition of normal equation matrix
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      IF (ptr_patch%geometry_info%cell_type == 3 )THEN
        nverts  = 3
      ELSE
        nverts = ptr_patch%cells%num_edges(jc,jb)
      ENDIF
    !
    ! Gather some information about the control volume and
    ! all the cells in the stencil.
    !
    ! get line and block indices of edge vertices
      jlv(1:nverts) = ptr_patch%cells%vertex_idx(jc,jb,1:nverts)
      jbv(1:nverts) = ptr_patch%cells%vertex_blk(jc,jb,1:nverts)

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    !
    ! 1. Get CC of the control volume center and all cell centers in the stencil.
    !    In addition, get cartesian coordinates of the vertices of the control volume.
    !

    ! get CC of control volume center
      cc_cv = ptr_patch%cells%cartesian_center(jc,jb)

    ! get CC of all cell centers in the stencil
      DO js = 1, ptr_ncells(jc,jb)
        cc_cell(js) = ptr_patch%cells%cartesian_center(ilc_s(js),ibc_s(js))
      ENDDO

    ! get cartesian coordinates of edge-vertices (only for control volume)
      DO jec = 1, nverts
        cc_vert(jec) = ptr_patch%verts%cartesian(jlv(jec),jbv(jec))
      ENDDO


    !
    ! 2. Now, all points (centers, vertices) are projected onto a tangent
    !    plane having its center at the cell center of the control volume.
    !    At the same time distance vectors (r_x,r_y) are calculated between
    !    the control volume center and all other mass-points in the stencil.
    !    From these the lsq weights can be deduced.

    !
    ! a: Calculate distance vectors between control volume center and all cell
    !    centers in the stencil.
      DO js = 1, ptr_ncells(jc,jb)
        !Get the actual location of the cell w.r.t the cc_cv
        cc_cell(js) = plane_torus_closest_coordinates(cc_cv%x, cc_cell(js)%x, &
                                                      ptr_patch%geometry_info)

        !the distance vector: z coord is 0
        z_dist_g(jc,jb,js,1) = cc_cell(js)%x(1) - cc_cv%x(1)
        z_dist_g(jc,jb,js,2) = cc_cell(js)%x(2) - cc_cv%x(2)
      ENDDO


    !
    ! b: compute normalized weights for weighted least squares system
    !    The closest cell circumcenter is assigned weight of w=1.
    !
      DO js = 1, ptr_ncells(jc,jb)
        z_norm = SQRT(DOT_PRODUCT(z_dist_g(jc,jb,js,1:2),z_dist_g(jc,jb,js,1:2)))

        !
        ! weights for weighted least squares system
        !
        ptr_int_lsq%lsq_weights_c(jc,js,jb)= 1._wp/(z_norm**lsq_wgt_exp)

      ENDDO
      !
      ! Normalization
      !
      ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)=                &
        &     ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb)         &
        &   / MAXVAL(ptr_int_lsq%lsq_weights_c(jc,1:ptr_ncells(jc,jb),jb))



    ! 3. the following part (including calculation of moments) will only
    !    be called, if a conservative least squares reconstruction
    !    is chosen. Otherwise all moments will be equal to zero and the
    !    reconstruction simplifies to the standard non-conservative
    !    reconstruction.
      IF (llsq_rec_consv) THEN
      !
      ! a: Project control volume vertices and calculate distance vectors
      !    between cell center and vertices
      !
        DO jec=1,nverts

          !Get the actual location of the cell w.r.t the cc_cv
          cc_vert(jec) = plane_torus_closest_coordinates(cc_cv%x, cc_vert(jec)%x, &
                                                         ptr_patch%geometry_info)

          !the distance vector: z coord is 0
          distxy_v(jec,1) = cc_vert(jec)%x(1) - cc_cv%x(1)
          distxy_v(jec,2) = cc_vert(jec)%x(2) - cc_cv%x(2)

        ENDDO

      !AD: Remaining part of the code is same as the spherical part

      !
      ! b: calculate moments for given cell
      !    (calculated analytically; see Lauritzen CSLAM 09 for first 5 moments)
      !
      ! !DR: Those moments have been re-rechecked, using an alternative,
      !      quadrature-based formulation. Results have been identical up to
      !      roundoff-errors. Similarly the hat-moments have been checked. The
      !      inconsistency caused by the different projections involved do not
      !      seem to negatively effect the results.

      ! Storage docu for x^ny^m:
      ! lsq_moments(:,:,1) : x^1y^0
      ! lsq_moments(:,:,2) : x^0y^1
      ! lsq_moments(:,:,3) : x^2y^0
      ! lsq_moments(:,:,4) : x^0y^2
      ! lsq_moments(:,:,5) : x^1y^1
      ! lsq_moments(:,:,6) : x^3y^0
      ! lsq_moments(:,:,7) : x^0y^3
      ! lsq_moments(:,:,8) : x^2y^1
      ! lsq_moments(:,:,9) : x^1y^2
      !

        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! note that the distance vector distxy_v between each vertex and
          ! the center of the tangent plane are identical to the coordinates of
          ! each vertex on the tangent plane. Thus the distances between the
          ! vertices in x and y direction can be derived as follows:
          !
          ! longitudinal-distance between vertices on tangent plane
          delx(jec) = distxy_v(jecp,1) - distxy_v(jec,1)

          ! latitudinal-distance between vertices on tangent plane
          dely(jec) = distxy_v(jecp,2) - distxy_v(jec,2)


          !
          ! analytic moment calculation
          !
          ! 0: control volume area (reciprocal value)
          fx(jec) = distxy_v(jecp,1) + distxy_v(jec,1)

        ENDDO

        z_rcarea = 2._wp/DOT_PRODUCT(fx(1:nverts),dely(1:nverts))


        DO jec=1,nverts

          jecp = jec + 1
          IF(jec==nverts) jecp=1

          ! I. x^1y^0
          fx(jec) = distxy_v(jecp,1)**2              &
          &       + distxy_v(jecp,1)*distxy_v(jec,1) &
          &       + distxy_v(jec ,1)**2

          ! II. x^0y^1
          fy(jec) = distxy_v(jecp,2)**2              &
          &       + distxy_v(jecp,2)*distxy_v(jec,2) &
          &       + distxy_v(jec,2)**2

          IF ( lsq_dim_unk > 2 ) THEN

            ! III. x^2y^0
            fxx(jec) = (distxy_v(jecp,1)    + distxy_v(jec,1)       ) &
            &        * (distxy_v(jecp,1)**2 + distxy_v(jec,1)**2)

            ! IV. x^0y^2
            fyy(jec) = (distxy_v(jecp,2)    + distxy_v(jec,2)       ) &
            &        * (distxy_v(jecp,2)**2 + distxy_v(jec,2)**2)

            ! V. x^1y^1
            fxy(jec) = distxy_v(jecp,2) * (3._wp*distxy_v(jecp,1)**2                         &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + distxy_v(jec,1)**2 )     &
            &        + distxy_v(jec,2) * ( distxy_v(jecp,1)**2                               &
            &        + 2._wp * distxy_v(jecp,1) * distxy_v(jec,1) + 3._wp*distxy_v(jec,1)**2 )

          ENDIF ! lsq_dim_unk > 2

          IF ( lsq_dim_unk > 5 ) THEN

            ! VI.  x^3y^0
            fxxx(jec) = 5._wp*distxy_v(jec,1)**4                   &
              &       + 10._wp*distxy_v(jec,1)**3 * delx(jec)      &
              &       + 10._wp*distxy_v(jec,1)**2 * delx(jec)**2   &
              &       + 5._wp *distxy_v(jec,1)    * delx(jec)**3   &
              &       + delx(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fxxx(jec) = distxy_v(jecp,1)**4                       &
            !DR    &       + distxy_v(jec,1) * distxy_v(jecp,1)**3     &
            !DR    &       + distxy_v(jec,1)**2 * distxy_v(jecp,1)**2  &
            !DR    &       + distxy_v(jec,1)**3 * distxy_v(jecp,1)     &
            !DR    &       + distxy_v(jec,1)**4


            ! VII. x^0y^3
            fyyy(jec) = 5._wp*distxy_v(jec,2)**4                   &
              &       + 10._wp*distxy_v(jec,2)**3 * dely(jec)      &
              &       + 10._wp*distxy_v(jec,2)**2 * dely(jec)**2   &
              &       + 5._wp *distxy_v(jec,2)    * dely(jec)**3   &
              &       + dely(jec)**4

            !DR equivalent to the following MAPLE result
            !DR marginally more accurate, when compared against quadrature-based
            !DR reference solution.
            !DR  fyyy(jec) = distxy_v(jecp,2)**4                       &
            !DR    &       + distxy_v(jec,2)    * distxy_v(jecp,2)**3  &
            !DR    &       + distxy_v(jec,2)**2 * distxy_v(jecp,2)**2  &
            !DR    &       + distxy_v(jec,2)**3 * distxy_v(jecp,2)     &
            !DR    &       + distxy_v(jec,2)**4

          ENDIF ! lsq_dim_unk > 5

          IF ( lsq_dim_unk > 7 ) THEN

            ! VIII. x^2y^1
            fxxy(jec) = 4._wp*distxy_v(jecp,1)**3 *distxy_v(jecp,2)                    &
              &       + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jecp,2)   &
              &       + 2._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jecp,2) &
              &       +       distxy_v(jec,1)**3 * distxy_v(jecp,2)                    &
              &       +       distxy_v(jecp,1)**3 * distxy_v(jec,2)                    &
              &       + 2._wp*distxy_v(jec,1)*distxy_v(jecp,1)**2 * distxy_v(jec,2)    &
              &       + 3._wp*distxy_v(jec,1)**2 * distxy_v(jecp,1) * distxy_v(jec,2)  &
              &       + 4._wp*distxy_v(jec,1)**3 * distxy_v(jec,2)

            ! IX. x^1y^2
            fxyy(jec) = 6._wp*distxy_v(jecp,1)**2 * distxy_v(jecp,2)**2                     &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jecp,2)**2              &
              &   +       distxy_v(jec,1)**2 * distxy_v(jecp,2)**2                          &
              &   + 3._wp*distxy_v(jecp,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)            &
              &   + 4._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)*distxy_v(jecp,2) &
              &   + 3._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)*distxy_v(jecp,2)             &
              &   +       distxy_v(jecp,1)**2 * distxy_v(jec,2)**2                          &
              &   + 3._wp*distxy_v(jec,1)*distxy_v(jecp,1)*distxy_v(jec,2)**2               &
              &   + 6._wp*distxy_v(jec,1)**2 * distxy_v(jec,2)**2

          ENDIF ! lsq_dim_unk > 7


        ENDDO ! loop over nverts

        ptr_int_lsq%lsq_moments(jc,jb,1)= z_rcarea/6._wp*DOT_PRODUCT(fx(1:nverts),dely(1:nverts))
        ptr_int_lsq%lsq_moments(jc,jb,2)=-z_rcarea/6._wp*DOT_PRODUCT(fy(1:nverts),delx(1:nverts))

        IF ( lsq_dim_unk > 2 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,3) =  z_rcarea/12._wp * &
            & DOT_PRODUCT(fxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,4) = -z_rcarea/12._wp * &
            & DOT_PRODUCT(fyy(1:nverts),delx(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,5) =  z_rcarea/24._wp * &
            & DOT_PRODUCT(fxy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 2


        IF ( lsq_dim_unk > 5 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,6) =  z_rcarea/20._wp * &
            & DOT_PRODUCT(fxxx(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,7) = -z_rcarea/20._wp * &
            & DOT_PRODUCT(fyyy(1:nverts),delx(1:nverts))

        END IF  ! lsq_dim_unk > 5


        IF ( lsq_dim_unk > 7 ) THEN

          ptr_int_lsq%lsq_moments(jc,jb,8) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxxy(1:nverts),dely(1:nverts))
          ptr_int_lsq%lsq_moments(jc,jb,9) = z_rcarea/60._wp * &
            & DOT_PRODUCT(fxyy(1:nverts),dely(1:nverts))

        END IF  ! lsq_dim_unk > 7

      END IF  ! llsq_rec_consv

    END DO  ! loop over cells

  END DO  ! loop over blocks
!!$OMP END DO NOWAIT
!! For unknown reasons, closing the parallel section here is needed to get the above
!! loop parallelized.
!!$OMP END PARALLEL

  DO jb = 1, lsq_dim_c
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_weights_c(:,jb,:))
  ENDDO
  DO jb = 1, lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments(:,:,jb))
  ENDDO


!$OMP PARALLEL PRIVATE(jb,jc,js,ju,jja,jjb,jjk,i_startidx,i_endidx,ilc_s,ibc_s, &
!$OMP            z_lsq_mat_c,zs,zu,zv_t,zwork,ziwork,ist,icheck,za_debug, &
!$OMP            z_qmat,z_rmat,cnt,jrow,nel)
!$OMP DO ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_c

    CALL get_indices_c(ptr_patch, jb, i_startblk, nblks_c, &
                       i_startidx, i_endidx, i_rcstartlev)

    !
    ! 4. for each cell, calculate LSQ design matrix A
    !
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) THEN
        ! Take care that z_lsq_mat_c isn't singular
        z_lsq_mat_c(jc,:,:) = 0.0_wp
        DO js = 1, MIN(lsq_dim_unk, lsq_dim_c)
          z_lsq_mat_c(jc,js,js) = 1.0_wp
        ENDDO
        CYCLE
      ENDIF

    ! line and block indices of cells in the stencil
      ilc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_idx_c(jc,jb,1:ptr_ncells(jc,jb))
      ibc_s(1:ptr_ncells(jc,jb)) = ptr_int_lsq%lsq_blk_c(jc,jb,1:ptr_ncells(jc,jb))


    ! Calculate full moments lsq_moments_hat(ilc_s(js),ibc_s(js),ju)
    !
    ! Storage docu for x^ny^m:
    ! lsq_moments_hat(:,:,:,1) : \hat{x^1y^0}
    ! lsq_moments_hat(:,:,:,2) : \hat{x^0y^1}
    ! lsq_moments_hat(:,:,:,3) : \hat{x^2y^0}
    ! lsq_moments_hat(:,:,:,4) : \hat{x^0y^2}
    ! lsq_moments_hat(:,:,:,5) : \hat{x^1y^1}
    ! lsq_moments_hat(:,:,:,6) : \hat{x^3y^0}
    ! lsq_moments_hat(:,:,:,7) : \hat{x^0y^3}
    ! lsq_moments_hat(:,:,:,8) : \hat{x^2y^1}
    ! lsq_moments_hat(:,:,:,9) : \hat{x^1y^2}
    !
      DO js = 1, ptr_ncells(jc,jb)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,1) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1) + z_dist_g(jc,jb,js,1)

        ptr_int_lsq%lsq_moments_hat(jc,jb,js,2) =                                             &
         &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2) + z_dist_g(jc,jb,js,2)

        IF (lsq_dim_unk > 2) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,3) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
           &    + z_dist_g(jc,jb,js,1)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,4) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)                              &
           &    + 2._wp* ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2) &
           &    + z_dist_g(jc,jb,js,2)**2

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,5) =                                           &
           &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)                              &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)        &
           &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)        &
           &    + z_dist_g(jc,jb,js,1) * z_dist_g(jc,jb,js,2)
        ENDIF

        IF ( lsq_dim_unk > 5 ) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,6) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),6)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,1) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)                       &
            &    * z_dist_g(jc,jb,js,1)**2                                                    &
            &    + z_dist_g(jc,jb,js,1)**3

          ptr_int_lsq%lsq_moments_hat(jc,jb,js,7) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),7)                             &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,2) &
            &    + 3._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)                       &
            &    * z_dist_g(jc,jb,js,2)**2                                                    &
            &    + z_dist_g(jc,jb,js,2)**3
        ENDIF

        IF ( lsq_dim_unk > 7 ) THEN
          ptr_int_lsq%lsq_moments_hat(jc,jb,js,8) =                                           &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),8)                             &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),3)* z_dist_g(jc,jb,js,2)       &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,1) &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,1) &
            &    * z_dist_g(jc,jb,js,2)                                                       &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,1)**2    &
            &    + z_dist_g(jc,jb,js,1)**2 * z_dist_g(jc,jb,js,2)


          ptr_int_lsq%lsq_moments_hat(jc,jb,js,9) =                                          &
            &      ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),9)                            &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),4)* z_dist_g(jc,jb,js,1)      &
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),5)* z_dist_g(jc,jb,js,2)&
            &    + 2._wp*ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),2)* z_dist_g(jc,jb,js,2)&
            &    * z_dist_g(jc,jb,js,1)                                                      &
            &    + ptr_int_lsq%lsq_moments(ilc_s(js),ibc_s(js),1)* z_dist_g(jc,jb,js,2)**2   &
            &    + z_dist_g(jc,jb,js,2)**2 * z_dist_g(jc,jb,js,1)
        ENDIF
      ENDDO


      ! loop over rows of lsq design matrix (all cells in the stencil)
      DO js = 1, ptr_ncells(jc,jb)
        ! loop over columns of lsq design matrix (number of unknowns)
        DO ju = 1,lsq_dim_unk

          z_lsq_mat_c(jc,js,ju) = ptr_int_lsq%lsq_weights_c(jc,js,jb)                         &
           &    * (ptr_int_lsq%lsq_moments_hat(jc,jb,js,ju)-ptr_int_lsq%lsq_moments(jc,jb,ju))

        END DO
       END DO
       IF(ptr_ncells(jc,jb) < lsq_dim_c) THEN
         z_lsq_mat_c(jc,lsq_dim_c,:) = 0.0_wp
       ENDIF

    ENDDO   ! loop over cells



    !
    ! compute QR decomposition and Singular Value Decomposition (SVD)
    ! of least squares design matrix A. For the time being both methods are
    ! retained.

    !
    ! 5a. QR-factorization of design matrix A
    !
    IF (.NOT. advection_config(pid)%llsq_svd) THEN
!CDIR NOIEXPAND
    CALL qrdec(lsq_dim_c, lsq_dim_unk, i_startidx, & ! in
     &         i_endidx, z_lsq_mat_c,              & ! in
     &         z_qmat, z_rmat)                       ! out


    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! 7. Save transposed Q-Matrix
      ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,1:lsq_dim_c,jb)  =  &
       &      TRANSPOSE(z_qmat(jc,1:lsq_dim_c,1:lsq_dim_unk))

      ! 8. Save R-Matrix
      !
      ! a. Save reciprocal values of the diagonal elements
      !
      DO ju = 1,lsq_dim_unk
        ptr_int_lsq%lsq_rmat_rdiag_c(jc,ju,jb) = 1._wp/z_rmat(jc,ju,ju)
      ENDDO

      !
      ! b. Save upper triangular elements without diagonal elements in a 1D-array
      !    (starting from the bottom right)
      !
      cnt = 1
      DO jrow = lsq_dim_unk-1,1,-1
        ! number of elements to store
        nel = lsq_dim_unk - jrow
        ptr_int_lsq%lsq_rmat_utri_c(jc,cnt:cnt+nel-1,jb) = z_rmat(jc,jrow,jrow+1:lsq_dim_unk)
        cnt = cnt + nel
      ENDDO


      ! Multiply ith column of the transposed Q-matrix (corresponds to the
      ! different members of the stencil) with the ith weight. This avoids
      ! multiplication of the RHS of the LSQ-System with this weight during
      ! runtime.
      DO js = 1,lsq_dim_c
        ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  =                  &
          &                ptr_int_lsq%lsq_qtmat_c(jc,1:lsq_dim_unk,js,jb)  &
          &              * ptr_int_lsq%lsq_weights_c(jc,js,jb)
      ENDDO

    END DO  ! loop over cells

    ELSE   ! llsq_svd=.TRUE.

    !
    ! 5b. Singular value decomposition of lsq design matrix A
    ! !!! does not vectorize, unfortunately !!!
    ist = 0
    DO jc = i_startidx, i_endidx

      IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

      ! A = U * SIGMA * transpose(V)
      !
      ! z_lsq_mat_c : M x N least squares design matrix A            (IN)
      ! zu          : M x M orthogonal matrix U                      (OUT)
      ! zv_t        : N x N orthogonal matrix V                      (OUT)
      ! zs          : min(M,N) Singular values of A                  (OUT)
      ! zwork       : workspace(1,LWORK)                             (OUT)
      ! lwork       : 3*min(M,N)                                     (IN)
      !              + max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N))  (IN)
      ! iwork       : workspace(8*min(M,N))                          (IN)

      CALL DGESDD('A',                 & !in
        &         lsq_dim_c,           & !in
        &         lsq_dim_unk,         & !in
        &         z_lsq_mat_c(jc,:,:), & !inout Note: destroyed on output
        &         lsq_dim_c,           & !in
        &         zs(:,jc),            & !out
        &         zu(:,:,jc),          & !out
        &         lsq_dim_c,           & !in
        &         zv_t(:,:,jc),        & !out
        &         lsq_dim_unk,         & !in
        &         zwork,               & !out
        &         lwork,               & !in
        &         ziwork,              & !inout
        &         icheck               ) !out
      ist = ist + icheck
    ENDDO
    IF (ist /= SUCCESS) THEN
      CALL finish ('mo_interpolation:lsq_compute_coeff_cell_torus',   &
        &             'singular value decomposition failed')
    ENDIF

    ! compute Moore-Penrose inverse
    ! INVERSE(A):: V * INVERSE(SIGMA) * TRANSPOSE(U) and store
    ! note that the ith column is multiplied with the ith weight
    ! in order to avoid the weighting of the r.h.s. during runtime.
    DO jja = 1, lsq_dim_unk
      DO jjb = 1, lsq_dim_c
        DO jjk = 1, lsq_dim_unk
          DO jc = i_startidx, i_endidx
            IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE
            ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb) =            &
              &  ptr_int_lsq%lsq_pseudoinv(jc,jja,jjb,jb)         &
              &  + zv_t(jjk,jja,jc) /zs(jjk,jc) * zu(jjb,jjk,jc)  &
              &  * ptr_int_lsq%lsq_weights_c(jc,jjb,jb)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!!!!!DEBUG !!!
#ifdef DEBUG_COEFF
    za_debug(:,:,:)  = 0._wp
    ! re-COMPUTE A:: U * SIGMA * TRANSPOSE(V)  !!! Funktioniert
    DO jja = 1, lsq_dim_c
      DO jjb = 1, lsq_dim_unk
        DO jjk = 1, lsq_dim_unk  !lsq_dim_c
          DO jc = i_startidx, i_endidx
            za_debug(jc,jja,jjb) = za_debug(jc,jja,jjb)  &
              &  + zu(jja,jjk,jc) * zs(jjk,jc) * zv_t(jjk,jjb,jc)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    write(0,*) "za_debug(10,:,:)- z_lsq_mat_c(10,:,:)",za_debug(10,:,:)- z_lsq_mat_c(10,:,:)
    write(0,*) "za_debug(55,:,:)- z_lsq_mat_c(55,:,:)",za_debug(55,:,:)- z_lsq_mat_c(55,:,:)
    write(0,*) "zs(:,10) ",zs(:,10)
    write(0,*) "zs(:,55) ",zs(:,55)
    write(0,*) "ptr_int_lsq%lsq_weights_c(10,:,jb)",ptr_int_lsq%lsq_weights_c(10,:,jb)
    write(0,*) "ptr_int_lsq%lsq_weights_c(55,:,jb)",ptr_int_lsq%lsq_weights_c(55,:,jb)
#endif
!!!! END DEBUG !!!

    ENDIF  ! llsq_svd

  END DO  ! loop over blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  DO ju = 1, lsq_dim_unk
    DO jc = 1, lsq_dim_c
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_moments_hat(:,:,jc,ju))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_pseudoinv(:,ju,jc,:))
      CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_qtmat_c(:,ju,jc,:))
    ENDDO
  ENDDO

  DO jc = 1, UBOUND(ptr_int_lsq%lsq_rmat_utri_c, 2)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_utri_c(:,jc,:))
  ENDDO
  DO ju = 1,lsq_dim_unk
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_lsq%lsq_rmat_rdiag_c(:,ju,:))
  ENDDO

  DEALLOCATE (z_dist_g, STAT=ist )
  IF (ist /= SUCCESS) THEN
    CALL finish ('mo_interpolation:lsq_compute_coeff_cell_torus',   &
      &             'deallocation for z_dist_g failed')
  ENDIF


END SUBROUTINE lsq_compute_coeff_cell_torus



!-------------------------------------------------------------------------
!>
!! This routine initializes the coefficients used.
!!
!! This routine initializes the coefficients used
!! for interpolations needed for scalars. The original routines were aw_int_coeff
!! and cell2edge_lin_int_coeff
!!
!! @par Revision History
!!  Original version by Hui Wan, MPI-M (2007-08-02)
!!  Modified by Almut Gassmann, MPI-M (2009-01-05)
!!  - added interpolation weights for edges to verts
!!  Renamed by Almut Gassmann, MPI-M (2009-01-27)
!!  - joining aw_int_coeff and cell2edge_lin_int_coeff to this routine
!!  - use of different edge-midpoint to cell distances
!!  - implement area weightings
!!  Modification by Daniel Reinert, DWD (2013-07-19)
!!  - added coefficients for pseudo-Laplacian weighted averaging (PLWA)
!!    from cells to verts. This type of averaging is second order accurate
!!    on general grids.
!!
!!  Literature for PLWA:
!!  - Holmes and Connel (1993): Solution of the 2D Navier-Stokes equations on
!!    unstructured adaptive grids. AIAA Paper 89-1932, AIAA 9th Computational
!!    Fluid Dynamics Conference
!!  - Kim et al (2003) : A multi-dimensional linear reconstruction scheme for
!!    arbitrary unstructured grids. AIAA Paper 2003-3990, 16th AIAA Computational
!!    Fluid Dynamics Conference.
!!  - Miura, H. (2013) : An Upwind-Biased Conservative Transport Scheme for Multi-Stage
!!                       Temporal Integrations on Spherical Icosahedral Grids.
!!                       Monthly Weather Review, in press
!!
SUBROUTINE scalar_int_coeff( ptr_patch, ptr_int_state )
!

TYPE(t_patch), INTENT(inout) :: ptr_patch

TYPE(t_int_state), INTENT(inout) :: ptr_int_state

INTEGER :: nlen, nblks_c, npromz_c, nblks_e, npromz_e, nblks_v, npromz_v
INTEGER :: jc, je, jb, jv  ! integer over edges, blocks and levels
INTEGER :: ile, ibe, &
           ilc, ibc, ilc1, ilc2, ibc1, ibc2, idx_ve,&
           ilv, ibv, ilv1, ilv2, ibv1, ibv2, idx_ce
INTEGER :: i_startblk                ! start block
INTEGER :: i_startidx                ! start index
INTEGER :: i_endidx                  ! end index

REAL(wp) :: z_sum

REAL(wp) :: vlon, vlat                 ! vertex coordinates
REAL(wp) :: clon, clat                 ! cell center coordinates
REAL(wp) :: r_x, r_y, i_xx, i_yy, i_xy ! moments for PLWA
REAL(wp) :: lambda_x, lambda_y         ! Lagrangesch multipliers for PLWA
REAL(wp) :: delx(6), dely(6)           ! distance cell center - vertex in x and y
REAL(wp) :: dist(6)                    ! sqrt(delx**2 + dely**2)
REAL(wp) :: wgt_sum                    ! sum of weights

!DR Will be removed after a first testing phase
!DRREAL(wp) :: test_l_x, test_l_y


!--------------------------------------------------------------------

  ! values for the blocking
  nblks_c  = ptr_patch%nblks_c
  npromz_c = ptr_patch%npromz_c
  nblks_e  = ptr_patch%nblks_e
  npromz_e = ptr_patch%npromz_e
  nblks_v  = ptr_patch%nblks_v
  npromz_v = ptr_patch%npromz_v

  ! a) the control volume associated to each edge is defined as the
  ! quadrilateral whose edges are the primal edge and the associated dual edge
  !----------------------------------------------------------------------------
  ! loop over all blocks and edges

!$OMP PARALLEL PRIVATE(i_startblk)
!$OMP DO PRIVATE(jb,je,nlen) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_e

    IF (jb /= nblks_e) THEN
      nlen = nproma
    ELSE
      nlen = npromz_e
    ENDIF

    DO je = 1, nlen
      ptr_patch%edges%area_edge(je,jb) =  &
        &    ptr_patch%edges%primal_edge_length(je,jb)  &
        &  * ptr_patch%edges%dual_edge_length(je,jb)
    END DO

  END DO
!$OMP END DO

  ! b1) cell to edge averages
  !-------------------------
  ! The calculation cannot be done for boundary edges
  i_startblk = ptr_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,ilc1,ilc2,ibc1,ibc2,ilv1,ilv2,&
!$OMP            ibv1,ibv2) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_e

    CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                       i_startidx, i_endidx, 2)

    DO je = i_startidx, i_endidx

      IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

      ! inverse distance averaging (former subroutine cell2edge_lin_int_coeff)
      ! For the hexagonal grid, this is also the direct distance averaging
      ! because the edge is exactly half way between the cells
      ! (needed for proper bracket formalism)
      ptr_int_state%c_lin_e(je,1,jb) = ptr_patch%edges%edge_cell_length(je,jb,2)/&
                                           ptr_patch%edges%dual_edge_length(je,jb)
      ptr_int_state%c_lin_e(je,2,jb) = 1._wp - ptr_int_state%c_lin_e(je,1,jb)

      IF (ptr_patch%geometry_info%cell_type == 6) THEN
        ilv1 = ptr_patch%edges%vertex_idx(je,jb,1)
        ilv2 = ptr_patch%edges%vertex_idx(je,jb,2)
        ibv1 = ptr_patch%edges%vertex_blk(je,jb,1)
        ibv2 = ptr_patch%edges%vertex_blk(je,jb,2)
        ptr_int_state%tria_aw_rhom(je,1,jb)=ptr_patch%verts%dual_area(ilv1,ibv1)/&
              (ptr_patch%verts%dual_area(ilv1,ibv1)+ptr_patch%verts%dual_area(ilv2,ibv2))
        ptr_int_state%tria_aw_rhom(je,2,jb)=ptr_patch%verts%dual_area(ilv2,ibv2)/&
              (ptr_patch%verts%dual_area(ilv1,ibv1)+ptr_patch%verts%dual_area(ilv2,ibv2))
      ENDIF

    ENDDO

  ENDDO
!$OMP END DO



  ! b2) vert to edge averages
  !-------------------------
  IF (ptr_patch%geometry_info%cell_type == 6) THEN
    ! The calculation cannot be done for boundary edges
    i_startblk = ptr_patch%edges%start_blk(2,1)
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_e

      CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

      DO je = i_startidx, i_endidx

        IF(.NOT. ptr_patch%edges%decomp_info%owner_mask(je,jb)) CYCLE

        ! distance averaging
        ptr_int_state%v_1o2_e(je,1,jb) = ptr_patch%edges%edge_vert_length(je,jb,1)/&
                                             ptr_patch%edges%primal_edge_length(je,jb)
        ptr_int_state%v_1o2_e(je,2,jb) = ptr_patch%edges%edge_vert_length(je,jb,2)/&
                                             ptr_patch%edges%primal_edge_length(je,jb)

      ENDDO

    ENDDO
!$OMP END DO
  ENDIF

  ! c) vert to cell averagings, edge to cell inner product
  !-------------------------------------------------------
  ! loop over all blocks and cells

!$OMP DO PRIVATE(jb,jc,je,jv,nlen,ile,ibe,idx_ce,ilv1,ilv2,ibv1,ibv2,&
!$OMP ilv,ibv,z_sum) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = 1, nblks_c
    IF (jb /= nblks_c) THEN
      nlen = nproma
    ELSE
      nlen = npromz_c
    ENDIF

    DO jc = 1, nlen

       IF(.NOT. ptr_patch%cells%decomp_info%owner_mask(jc,jb)) CYCLE

       ptr_int_state%verts_aw_cells(jc,:,jb) = 0.0_wp

       IF (ptr_patch%geometry_info%cell_type == 6) z_sum = 0.0_wp

       DO je = 1, ptr_patch%cells%num_edges(jc,jb)

          ile = ptr_patch%cells%edge_idx(jc,jb,je)
          ibe = ptr_patch%cells%edge_blk(jc,jb,je)
          IF ( ptr_patch%edges%cell_idx(ile,ibe,1) == jc .AND. &
               ptr_patch%edges%cell_blk(ile,ibe,1) == jb ) THEN
               idx_ce = 1
          ELSE
               idx_ce = 2
          ENDIF

          ptr_int_state%e_inn_c(jc,je,jb) = &
              ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
              ptr_patch%edges%primal_edge_length(ile,ibe)/&
              ptr_patch%cells%area(jc,jb)

          IF (ptr_patch%geometry_info%cell_type == 6) THEN
            ptr_int_state%e_aw_c(jc,je,jb) = 0.5_wp*&
              ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
              ptr_patch%edges%primal_edge_length(ile,ibe)/&
              ptr_patch%cells%area(jc,jb)
            ptr_int_state%r_aw_c(jc,je,jb) = &
              ptr_patch%edges%quad_area(ile,ibe)
            z_sum = z_sum+ptr_patch%edges%quad_area(ile,ibe)
          ENDIF

          ilv1 = ptr_patch%edges%vertex_idx(ile,ibe,1)
          ibv1 = ptr_patch%edges%vertex_blk(ile,ibe,1)
          ilv2 = ptr_patch%edges%vertex_idx(ile,ibe,2)
          ibv2 = ptr_patch%edges%vertex_blk(ile,ibe,2)

          DO jv = 1, ptr_patch%cells%num_edges(jc,jb)
            ilv = ptr_patch%cells%vertex_idx(jc,jb,jv)
            ibv = ptr_patch%cells%vertex_blk(jc,jb,jv)

            IF (ilv == ilv1 .AND. ibv == ibv1) THEN
              ptr_int_state%verts_aw_cells(jc,jv,jb) =   &
                ptr_int_state%verts_aw_cells(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,1)
            ENDIF
            IF (ilv == ilv2 .AND. ibv == ibv2) THEN
              ptr_int_state%verts_aw_cells(jc,jv,jb)  =  &
                ptr_int_state%verts_aw_cells(jc,jv,jb) + &
                0.5_wp/ptr_patch%cells%area(jc,jb) *             &
                ptr_patch%edges%edge_cell_length(ile,ibe,idx_ce)*&
                ptr_patch%edges%edge_vert_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO
       IF (ptr_patch%geometry_info%cell_type == 6) THEN
         ptr_int_state%r_aw_c(jc,:,jb)=ptr_int_state%r_aw_c(jc,:,jb)/z_sum
       ENDIF

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!$OMP END DO



  ! d) cells to verts averagings, edge to verts averagings
  !-------------------------------------------------------
  ! loop over all blocks and verts

  i_startblk = ptr_patch%verts%start_blk(2,1)
!$OMP DO PRIVATE(jb,jc,je,jv,i_startidx,i_endidx,ile,ibe,idx_ve,ilc,ibc, &
!$OMP            ilc1,ilc2,ibc1,ibc2 ) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, 2)

    DO jv = i_startidx, i_endidx

       IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE

       ptr_int_state%cells_aw_verts(jv,:,jb) = 0.0_wp

       DO je = 1, ptr_patch%verts%num_edges(jv,jb)

          ile = ptr_patch%verts%edge_idx(jv,jb,je)
          ibe = ptr_patch%verts%edge_blk(jv,jb,je)
          IF ( ptr_patch%edges%vertex_idx(ile,ibe,1) == jv .AND. &
               ptr_patch%edges%vertex_blk(ile,ibe,1) == jb ) THEN
               idx_ve = 1
          ELSE
               idx_ve = 2
          ENDIF

          ptr_int_state%e_aw_v(jv,je,jb) = 0.5_wp*&
            & ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve) &
            &*ptr_patch%edges%dual_edge_length(ile,ibe) &
            &/ptr_patch%verts%dual_area(jv,jb)

          IF (ptr_patch%geometry_info%cell_type == 6 ) THEN
            ptr_int_state%e_inn_v(jv,je,jb) = &
            & ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve) &
            &*ptr_patch%edges%dual_edge_length(ile,ibe) &
            &/ptr_patch%verts%dual_area(jv,jb)
          ENDIF

          ilc1 = ptr_patch%edges%cell_idx(ile,ibe,1)
          ibc1 = ptr_patch%edges%cell_blk(ile,ibe,1)
          ilc2 = ptr_patch%edges%cell_idx(ile,ibe,2)
          ibc2 = ptr_patch%edges%cell_blk(ile,ibe,2)

          DO jc = 1, ptr_patch%verts%num_edges(jv,jb)
            ilc = ptr_patch%verts%cell_idx(jv,jb,jc)
            ibc = ptr_patch%verts%cell_blk(jv,jb,jc)

            IF (ilc == ilc1 .AND. ibc == ibc1) THEN
              ptr_int_state%cells_aw_verts(jv,jc,jb) =   &
                ptr_int_state%cells_aw_verts(jv,jc,jb) + &
                0.5_wp/ptr_patch%verts%dual_area(jv,jb) *             &
                ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve)*&
                ptr_patch%edges%edge_cell_length(ile,ibe,1)
            ENDIF
            IF (ilc == ilc2 .AND. ibc == ibc2) THEN
              ptr_int_state%cells_aw_verts(jv,jc,jb)  =  &
                ptr_int_state%cells_aw_verts(jv,jc,jb) + &
                0.5_wp/ptr_patch%verts%dual_area(jv,jb) *             &
                ptr_patch%edges%edge_vert_length(ile,ibe,idx_ve)*&
                ptr_patch%edges%edge_cell_length(ile,ibe,2)
            ENDIF

          ENDDO

       ENDDO

    ENDDO !loop over all cells

  ENDDO   !loop over all blocks
!$OMP END DO NOWAIT


  ! e) cells to verts pseudo-Laplacian weighted averaging (PLWA)
  !-------------------------------------------------------------

  i_startblk = ptr_patch%verts%start_blk(2,1)
!$OMP DO PRIVATE(jb,jv,jc,i_startidx,i_endidx,vlon,vlat,ilc,ibc,clon,clat, &
!$OMP            delx,dely,dist,r_x,r_y,i_xx,i_yy,i_xy,lambda_x,lambda_y,  &
!$OMP            wgt_sum) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, nblks_v

    CALL get_indices_v(ptr_patch, jb, i_startblk, nblks_v, &
                       i_startidx, i_endidx, 2)

    DO jv = i_startidx, i_endidx

      IF(.NOT. ptr_patch%verts%decomp_info%owner_mask(jv,jb)) CYCLE


      ! initialize moments
      r_x     = 0._wp
      r_y     = 0._wp
      i_xx    = 0._wp
      i_yy    = 0._wp
      i_xy    = 0._wp
      wgt_sum = 0._wp
!DR      test_l_x = 0._wp
!DR      test_l_y = 0._wp

      ptr_int_state%cells_plwa_verts(jv,:,jb) = 0.0_wp


      ! get geographical coordinates of vertex
      vlon = ptr_patch%verts%vertex(jv,jb)%lon
      vlat = ptr_patch%verts%vertex(jv,jb)%lat

      ! loop over all cells adjacent to this vertex
      DO jc = 1, ptr_patch%verts%num_edges(jv,jb)
        ilc = ptr_patch%verts%cell_idx(jv,jb,jc)
        ibc = ptr_patch%verts%cell_blk(jv,jb,jc)

        ! get geographical coordinates of cell center
        clon = ptr_patch%cells%center(ilc,ibc)%lon
        clat = ptr_patch%cells%center(ilc,ibc)%lat


!DR could be replaced by rotation to equator as done at various other places
!DR (not tested yet)
        ! project cell center onto local 2D cartesian system
        ! tangent to the vertex
        CALL gnomonic_proj( vlon, vlat, clon, clat,  &! in
          &                 delx(jc), dely(jc))       ! out


        ! compute distance to vertex
        dist(jc) = SQRT(delx(jc)**2 + dely(jc)**2)


        ! sum up moments
        r_x  = r_x  + delx(jc)
        r_y  = r_y  + dely(jc)
        i_xx = i_xx + delx(jc)*delx(jc)/dist(jc)**2
        i_yy = i_yy + dely(jc)*dely(jc)/dist(jc)**2
        i_xy = i_xy + delx(jc)*dely(jc)/dist(jc)**2
      ENDDO  ! jc


      ! compute Lagrangesch multipliers
      lambda_x = (r_y*i_xy - r_x*i_yy)/(i_xx*i_yy - i_xy**2)
      lambda_y = (r_x*i_xy - r_y*i_xx)/(i_xx*i_yy - i_xy**2)



      ! compute weights
      ! loop over all cells adjacent to this vertex
      DO jc = 1, ptr_patch%verts%num_edges(jv,jb)

        ! weights
        ptr_int_state%cells_plwa_verts(jv,jc,jb) = 1._wp                                    &
          &                                      + ((lambda_x*delx(jc) + lambda_y*dely(jc)) &
          &                                      / dist(jc)**2 )

        ! sum up weights for normalization
        wgt_sum = wgt_sum + ptr_int_state%cells_plwa_verts(jv,jc,jb)

!DR        test_l_x = test_l_x + ptr_int_state%cells_plwa_verts(jv,jc,jb) * delx(jc)
!DR        test_l_y = test_l_y + ptr_int_state%cells_plwa_verts(jv,jc,jb) * dely(jc)
      ENDDO  ! jc

      ! Normalization
      DO jc = 1, ptr_patch%verts%num_edges(jv,jb)
        ptr_int_state%cells_plwa_verts(jv,jc,jb) = ptr_int_state%cells_plwa_verts(jv,jc,jb)  &
          &                                      / wgt_sum
      ENDDO
!!$write(0,*) "test_l_x: ", test_l_x
!!$write(0,*) "test_l_y: ", test_l_y
    ENDDO  ! jv

!!$write(0,*) "MAX(cells_plwa_verts), jb: ", MAXVAL(ptr_int_state%cells_plwa_verts(:,:,jb)), jb
!!$write(0,*) "MIN(cells_plwa_verts), jb: ", MINVAL(ptr_int_state%cells_plwa_verts(:,:,jb)), jb
  ENDDO  ! loop over all blocks
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_patch%edges%area_edge)
  CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int_state%c_lin_e)
  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%verts_aw_cells)
  CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_inn_c)
  CALL sync_patch_array(SYNC_V,ptr_patch,ptr_int_state%cells_aw_verts)
  CALL sync_patch_array(SYNC_V,ptr_patch,ptr_int_state%cells_plwa_verts)
  CALL sync_patch_array(SYNC_V,ptr_patch,ptr_int_state%e_aw_v)
  IF (ptr_patch%geometry_info%cell_type == 6) THEN
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int_state%tria_aw_rhom)
    CALL sync_patch_array(SYNC_E,ptr_patch,ptr_int_state%v_1o2_e)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_aw_c)
    CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%r_aw_c)
    CALL sync_patch_array(SYNC_V,ptr_patch,ptr_int_state%e_inn_v)
  ENDIF

END SUBROUTINE scalar_int_coeff
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!>
!! Calls routines to calculate the weighting coefficients for bilinear
!! edge-to-cell interpolation depending on grid geometry.
!!
SUBROUTINE bln_int_coeff_e2c( ptr_patch, ptr_int_state )
  TYPE(t_patch),     INTENT(inout) :: ptr_patch
  TYPE(t_int_state), INTENT(inout) :: ptr_int_state

  CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_intp_coeffs_lsq_bln:bln_int_coeff_e2c'

  !
  SELECT CASE(ptr_patch%geometry_info%geometry_type)

  CASE (planar_torus_geometry)
    CALL flat_scalar_coeffs( ptr_patch, ptr_int_state )
    CALL vector_coeffs(ptr_patch, ptr_int_state)

  CASE (sphere_geometry)
    CALL spherical_scalar_coeffs( ptr_patch, ptr_int_state )
    CALL vector_coeffs(ptr_patch, ptr_int_state)

  CASE DEFAULT
    CALL finish(method_name, "Undefined geometry type")

  END SELECT

END SUBROUTINE bln_int_coeff_e2c
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!
!
!>
!! Computes the weighting coefficients for bilinear edge-to-cell interpolation.
!!
!! Results are stored in ptr_int_state\\%e_bln_c_s
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2009-01-06
!!
SUBROUTINE spherical_scalar_coeffs ( ptr_patch, ptr_int_state )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3

REAL(wp) :: xtemp,ytemp,wgt(3),xloc,yloc,x(3),y(3), &
            pollat,pollon

!-----------------------------------------------------------------------

rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)
!
! loop through all patch cells (and blocks)
!
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,yloc,xloc,pollat,pollon,ile1,ibe1,&
!$OMP            ile2,ibe2,ile3,ibe3,xtemp,ytemp,wgt,x,y) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    yloc = ptr_patch%cells%center(jc,jb)%lat
    xloc = ptr_patch%cells%center(jc,jb)%lon

    ! Rotate local point into the equator for better accuracy of bilinear weights
    IF (yloc >= 0._wp) THEN
      pollat = yloc - pi2/4._wp
    ELSE
      pollat = yloc + pi2/4._wp
    ENDIF
    pollon = xloc

    CALL rotate_latlon( yloc, xloc, pollat, pollon )

    !  get the line and block indices of the cell edges

    ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
    ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
    ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
    ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
    ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
    ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

    ! x and y are the zonal and meridional distances from the local
    ! cell point (ignoring the earth's radius, which drops out anyway)

    xtemp = ptr_patch%edges%center(ile1,ibe1)%lon
    ytemp = ptr_patch%edges%center(ile1,ibe1)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(1)  = ytemp-yloc
    x(1)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(1) >  3.5_wp) x(1) = x(1) - pi2
    IF (x(1) < -3.5_wp) x(1) = x(1) + pi2

    xtemp = ptr_patch%edges%center(ile2,ibe2)%lon
    ytemp = ptr_patch%edges%center(ile2,ibe2)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(2)  = ytemp-yloc
    x(2)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(2) >  3.5_wp) x(2) = x(2) - pi2
    IF (x(2) < -3.5_wp) x(2) = x(2) + pi2

    xtemp = ptr_patch%edges%center(ile3,ibe3)%lon
    ytemp = ptr_patch%edges%center(ile3,ibe3)%lat
    CALL rotate_latlon( ytemp, xtemp, pollat, pollon )

    y(3)  = ytemp-yloc
    x(3)  = xtemp-xloc
    ! This is needed when the date line is crossed
    IF (x(3) >  3.5_wp) x(3) = x(3) - pi2
    IF (x(3) < -3.5_wp) x(3) = x(3) + pi2

    ! The weighting factors are based on the requirement that sum(w(i)*x(i)) = 0
    ! and sum(w(i)*y(i)) = 0, which ensures that linear horizontal gradients
    ! are not aliased into a checkerboard pattern between upward- and downward
    ! directed cells. The third condition is sum(w(i)) = 1. Analytical elimination yields...

    IF (ABS(x(2)-x(1)) > 1.e-11_wp .AND. ABS(y(3)-y(1)) > 1.e-11_wp ) THEN
      wgt(3) = 1.0_wp/( (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1))/(x(2)-x(1)) ) * &
                  ( -y(1) + x(1)*(y(2)-y(1))/(x(2)-x(1)) )
      wgt(2) = (-x(1) - wgt(3)*(x(3)-x(1)))/(x(2)-x(1))
      wgt(1) = 1.0_wp - wgt(2) - wgt(3)
    ELSE
      wgt(2) = 1.0_wp/( (y(2)-y(1)) - (x(2)-x(1))*(y(3)-y(1))/(x(3)-x(1)) ) * &
                  ( -y(1) + x(1)*(y(3)-y(1))/(x(3)-x(1)) )
      wgt(3) = (-x(1) - wgt(2)*(x(2)-x(1)))/(x(3)-x(1))
      wgt(1) = 1.0_wp - wgt(2) - wgt(3)
    ENDIF

    ! Store results in ptr_int_state%e_bln_c_s
    ptr_int_state%e_bln_c_s(jc,1,jb) = wgt(1)
    ptr_int_state%e_bln_c_s(jc,2,jb) = wgt(2)
    ptr_int_state%e_bln_c_s(jc,3,jb) = wgt(3)

  ENDDO !cell loop

END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_s)


END SUBROUTINE spherical_scalar_coeffs
!-------------------------------------------------------------------------


!-------------------------------------------------------------------------
!
!
!>
!! Computes the vector weighting coefficients using the scalar part computed earlier
!!
!! @par Revision History
!!  developed by Guenther Zaengl, 2009-01-06
!!
SUBROUTINE vector_coeffs ( ptr_patch, ptr_int_state )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER :: jc, jb
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom

INTEGER :: ile1, ibe1, ile2, ibe2, ile3, ibe3

!-----------------------------------------------------------------------

! Now compute vector interpolation weights: These take the normal and tangential
! wind components at the edges and reconstruct u and v at the cell midpoints

rl_start = 2
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ile1,ibe1,ile2,ibe2,&
!$OMP ile3,ibe3) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)

  DO jc = i_startidx, i_endidx

    !  get the line and block indices of the cell edges

    ile1 = ptr_patch%cells%edge_idx(jc,jb,1)
    ibe1 = ptr_patch%cells%edge_blk(jc,jb,1)
    ile2 = ptr_patch%cells%edge_idx(jc,jb,2)
    ibe2 = ptr_patch%cells%edge_blk(jc,jb,2)
    ile3 = ptr_patch%cells%edge_idx(jc,jb,3)
    ibe3 = ptr_patch%cells%edge_blk(jc,jb,3)

    IF (ptr_patch%edges%cell_idx(ile1,ibe1,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile1,ibe1,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v1
      ptr_int_state%e_bln_c_u(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,1)%v1

      ptr_int_state%e_bln_c_v(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,1)%v2
      ptr_int_state%e_bln_c_v(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile1,ibe1,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile1,ibe1,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v1
      ptr_int_state%e_bln_c_u(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,2)%v1

      ptr_int_state%e_bln_c_v(jc,1,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%primal_normal_cell(ile1,ibe1,2)%v2
      ptr_int_state%e_bln_c_v(jc,2,jb) = ptr_int_state%e_bln_c_s(jc,1,jb)* &
        ptr_patch%edges%dual_normal_cell(ile1,ibe1,2)%v2

    ENDIF

    IF (ptr_patch%edges%cell_idx(ile2,ibe2,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile2,ibe2,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,1)%v1
      ptr_int_state%e_bln_c_u(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,1)%v1

      ptr_int_state%e_bln_c_v(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,1)%v2
      ptr_int_state%e_bln_c_v(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile2,ibe2,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile2,ibe2,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,2)%v1
      ptr_int_state%e_bln_c_u(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,2)%v1

      ptr_int_state%e_bln_c_v(jc,3,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%primal_normal_cell(ile2,ibe2,2)%v2
      ptr_int_state%e_bln_c_v(jc,4,jb) = ptr_int_state%e_bln_c_s(jc,2,jb)* &
        ptr_patch%edges%dual_normal_cell(ile2,ibe2,2)%v2

    ENDIF

    IF (ptr_patch%edges%cell_idx(ile3,ibe3,1) == jc .AND. &
        ptr_patch%edges%cell_blk(ile3,ibe3,1) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,1)%v1
      ptr_int_state%e_bln_c_u(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,1)%v1

      ptr_int_state%e_bln_c_v(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,1)%v2
      ptr_int_state%e_bln_c_v(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,1)%v2

    ELSE IF (ptr_patch%edges%cell_idx(ile3,ibe3,2) == jc .AND. &
        ptr_patch%edges%cell_blk(ile3,ibe3,2) == jb) THEN

      ptr_int_state%e_bln_c_u(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,2)%v1
      ptr_int_state%e_bln_c_u(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,2)%v1

      ptr_int_state%e_bln_c_v(jc,5,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%primal_normal_cell(ile3,ibe3,2)%v2
      ptr_int_state%e_bln_c_v(jc,6,jb) = ptr_int_state%e_bln_c_s(jc,3,jb)* &
        ptr_patch%edges%dual_normal_cell(ile3,ibe3,2)%v2

    ENDIF
  ENDDO !cell loop

END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_u)
CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_v)

END SUBROUTINE vector_coeffs
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!>
!! Computes the weighting coefficients for bilinear edge-to-cell interpolation for
!! flat geometry with equilateral triangular cells
!!
!! @par Revision History
!!  developed by Anurag Dipankar, 2012-28-12 (taken from calculate_spherical_scalar_intp_coeffs)
!!
SUBROUTINE flat_scalar_coeffs ( ptr_patch, ptr_int_state )
!
!
!  patch on which computation is performed
!
TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! Interpolation coefficients
TYPE(t_int_state), TARGET, INTENT(inout) :: ptr_int_state
!

INTEGER  :: jc, jb
INTEGER  :: rl_start, rl_end
INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
REAL(wp) :: wgt

!-----------------------------------------------------------------------

rl_start = 1
rl_end = min_rlcell

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)

wgt = 1._wp/3._wp

!
! loop through all patch cells (and blocks)
!
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)
    DO jc = i_startidx, i_endidx

      ! Simple for torus
      ptr_int_state%e_bln_c_s(jc,1:3,jb) = wgt

    ENDDO !cell loop
  END DO !block loop
!$OMP END DO NOWAIT
!$OMP END PARALLEL

CALL sync_patch_array(SYNC_C,ptr_patch,ptr_int_state%e_bln_c_s)


END SUBROUTINE flat_scalar_coeffs
!-------------------------------------------------------------------------


END MODULE mo_intp_coeffs_lsq_bln
