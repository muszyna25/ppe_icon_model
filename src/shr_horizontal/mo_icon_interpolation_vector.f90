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
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_icon_interpolation_vector
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert, min_rlcell_int
  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: edges2cells_vector

CONTAINS



!------------------------------------------------------------------------
!
!>
!!  Bilinear interpolation of normal and tangential velocity components.
!!
!!  Bilinear interpolation of normal and tangential velocity components
!!  at the edges to u and v at the cells
!!  Works only for triangles (bilinear interpolation weights are not implemented
!!  for hexagons)
!!
!! @par Revision History
!!  Developed by Guenther Zaengl, DWD (2009-05-08)
!!
SUBROUTINE edges2cells_vector( p_vn_in, p_vt_in, ptr_patch, p_int, &
  &  p_u_out, p_v_out, opt_slev, opt_elev, opt_rlstart, opt_rlend  )
!


TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! normal velocity component at edges
REAL(wp), INTENT(in) ::  p_vn_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! (reconstructed) tangential velocity component at edges
REAL(wp), INTENT(in) ::  p_vt_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! Interpolation state
TYPE(t_int_state), INTENT(IN) :: p_int

! cell based output fields: u and v
REAL(wp), INTENT(inout) :: p_u_out(:,:,:), p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: rl_start, rl_end

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-------------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_vn_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell_int
END IF

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


IF (ltimer) CALL timer_start(timer_intp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
  DO jc = i_startidx, i_endidx
    DO jk = slev, elev
#else
!CDIR UNROLL=6
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx
#endif

      p_u_out(jc,jk,jb) =  &
        p_int%e_bln_c_u(jc,1,jb)*p_vn_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_u(jc,2,jb)*p_vt_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_u(jc,3,jb)*p_vn_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_u(jc,4,jb)*p_vt_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_u(jc,5,jb)*p_vn_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
        p_int%e_bln_c_u(jc,6,jb)*p_vt_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      p_v_out(jc,jk,jb) =  &
        p_int%e_bln_c_v(jc,1,jb)*p_vn_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_v(jc,2,jb)*p_vt_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_v(jc,3,jb)*p_vn_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_v(jc,4,jb)*p_vt_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_v(jc,5,jb)*p_vn_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
        p_int%e_bln_c_v(jc,6,jb)*p_vt_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

     ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE edges2cells_vector


END MODULE mo_icon_interpolation_vector
