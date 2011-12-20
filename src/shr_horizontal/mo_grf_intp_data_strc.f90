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
!! - moved gridref_ctl namelist variables to mo_gridref_nml
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
MODULE mo_grf_intp_data_strc

USE mo_kind, ONLY: wp

PUBLIC

!>
!----------------------------------------------------------------------------
!                                    Types
!----------------------------------------------------------------------------
TYPE t_gridref_single_state

  INTEGER, ALLOCATABLE  :: grf_vec_ind_1a (:,:,:), & ! index arrays defining the stencil
                           grf_vec_ind_1b (:,:,:), & ! of surrounding edges for RBF/IDW
                           grf_vec_ind_2a (:,:,:), & ! interpolation to lateral boundaries of
                           grf_vec_ind_2b (:,:,:)    ! refined grids

  INTEGER, ALLOCATABLE  :: grf_vec_blk_1a (:,:,:), & ! index arrays defining the stencil
                           grf_vec_blk_1b (:,:,:), & ! of surrounding edges for RBF/IDW
                           grf_vec_blk_2a (:,:,:), & ! interpolation to lateral boundaries of
                           grf_vec_blk_2b (:,:,:)    ! refined grids

  INTEGER, ALLOCATABLE  :: grf_vec_stencil_1a (:,:), & ! arrays defining the number of surrounding
                           grf_vec_stencil_1b (:,:), & ! edges/cells in the stencil for RBF/IDW
                           grf_vec_stencil_2a (:,:), & ! interpolation to lateral boundaries of
                           grf_vec_stencil_2b (:,:)   ! refined grids

  REAL(wp), ALLOCATABLE :: grf_vec_coeff_1a(:,:,:), & ! arrays containing the coefficients
                           grf_vec_coeff_1b(:,:,:), & ! used for RBF/IDW vector interpolation
                           grf_vec_coeff_2a(:,:,:), & ! to lateral boundaries of refined grids
                           grf_vec_coeff_2b(:,:,:)

  ! distances between parent cells and child cells needed for gradient-based interpolation
  REAL(wp), ALLOCATABLE :: grf_dist_pc2cc(:,:,:,:), grf_dist_pe2ce(:,:,:)

END TYPE t_gridref_single_state

TYPE t_gridref_state

  ! p_dom holds fields that are allocated for the boundary interpolation zone only
  ! They therefore need to be defined separately for each child domain
  TYPE(t_gridref_single_state), ALLOCATABLE :: p_dom(:)

  ! These fields are allocated for the full parent domain and thus do not need
  ! to be held separately for each child domain
  REAL(wp), ALLOCATABLE :: fbk_wgt_c(:,:,:)     ! Feedback weights for cell-based variables
                                                !  index1=1,nproma, index2=nblks_c, index3=1,4
  REAL(wp), ALLOCATABLE :: fbk_wgt_ct(:,:,:)    ! Feedback weights for cell-based tracer variables
                                                !  index1=1,nproma, index2=nblks_c, index3=1,4
  REAL(wp), ALLOCATABLE :: fbk_wgt_e(:,:,:)     ! The same for edge-based variables

  REAL(wp), ALLOCATABLE :: fbk_dom_area(:) ! Area of subdomain for which feedback is performed
                                           ! dimension: n_childdom

  REAL(wp), ALLOCATABLE :: fbk_dom_volume(:,:)! Area of subdomain for which feedback is performed
                                              ! dimensions: (nlev,n_childdom)

END TYPE t_gridref_state

  !----------------------------------------------------------------------------
  !                                Variables
  !----------------------------------------------------------------------------
  TYPE(t_gridref_state), TARGET, ALLOCATABLE :: p_grf_state(:), p_grf_state_local_parent(:)


END MODULE mo_grf_intp_data_strc
