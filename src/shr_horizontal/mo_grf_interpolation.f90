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
MODULE mo_grf_interpolation
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------

  USE mo_grf_intp_data_strc, ONLY: rbf_vec_kern_grf_e, rbf_scale_grf_e, &
       &                           grf_intmethod_c, grf_intmethod_e,    &
       &                           grf_intmethod_ct, grf_velfbk,        &
       &                           grf_scalfbk, grf_tracfbk,            &
       &                           grf_idw_exp_e12, grf_idw_exp_e34,    &
       &                           denom_diffu_v, denom_diffu_t,        &
       &                           t_gridref_single_state, t_gridref_state

  USE mo_grf_intp_state,     ONLY: setup_gridref,              &
       &                           construct_2d_gridref_state, &
       &                           destruct_2d_gridref_state

  PUBLIC  ! all that is USEd above

!-------------------------------------------------------------------------
END MODULE mo_grf_interpolation
