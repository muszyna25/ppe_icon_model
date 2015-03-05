!>
!!  Module contains some constants needed for grid refinement implementation.
!!
!!
!! @par Revision History
!!  Created by Guenther Zaengl, DWD (2009-12-14)
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
MODULE mo_impl_constants_grf
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

  IMPLICIT NONE

  PUBLIC

  ! Width of lateral boundary zones (as seen from the child domain) for which
  ! tendencies are interpolated from the parent domain
  INTEGER, PARAMETER :: grf_bdywidth_c = 4
  INTEGER, PARAMETER :: grf_bdywidth_e = 9

  ! Start and end refin_ctrl levels for boundary interpolation as seen from
  ! the parent domain
  INTEGER, PARAMETER :: grf_bdyintp_start_c = -1
  INTEGER, PARAMETER :: grf_bdyintp_start_e = -1
  INTEGER, PARAMETER :: grf_bdyintp_end_c   = - grf_bdywidth_c/2
  INTEGER, PARAMETER :: grf_bdyintp_end_e   = -(grf_bdywidth_e+1)/2

  ! Start refin_ctrl level of overlap area with nested domain (prognostic part plus one halo row)
  INTEGER, PARAMETER :: grf_ovlparea_start_c = -2

  ! Start refin_ctrl levels for feedback (as seen from the parent domain)
  INTEGER, PARAMETER :: grf_fbk_start_c = -3
  INTEGER, PARAMETER :: grf_fbk_start_e = -5

  ! Parameters for boundary nudging (needed for 1-way nesting and limited-area version)

  ! Width of nudging zone in units of child cell rows; MUST be a multiple of 2,
  ! AND grf_nudge_start_c + grf_width_nudgezone MUST NOT EXCEED the bdy_indexing_depth
  ! namelist variable in prepare_gridref
  INTEGER, PARAMETER :: grf_nudgezone_width = 8 ! to be determined later

  ! Start and end refin_ctrl levels for boundary nudging (as seen from the child
  ! domain)
  INTEGER, PARAMETER :: grf_nudge_start_c = grf_bdywidth_c + 1
  INTEGER, PARAMETER :: grf_nudge_start_e = grf_bdywidth_e + 1
  INTEGER, PARAMETER :: grf_nudge_end_c   = grf_nudge_start_c + grf_nudgezone_width
  INTEGER, PARAMETER :: grf_nudge_end_e   = grf_nudge_start_e + grf_nudgezone_width*2

  ! Start and end refin_ctrl levels for boundary nudging as seen from the parent
  ! domain, denoting cells/edges that have to do parent-to-child interpolation
  INTEGER, PARAMETER :: grf_nudgintp_start_c = grf_bdyintp_end_c - 1
  INTEGER, PARAMETER :: grf_nudgintp_start_e = grf_bdyintp_end_e
  INTEGER, PARAMETER :: grf_nudgintp_end_c   = - grf_nudge_end_c/2
  INTEGER, PARAMETER :: grf_nudgintp_end_e   = -(grf_nudge_end_e+1)/2

!--------------------------------------------------------------------
END MODULE mo_impl_constants_grf
