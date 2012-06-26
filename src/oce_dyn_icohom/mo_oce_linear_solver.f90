!>
!!  General comments: provide an implementation of the GMRES solver for.
!!
!!  General comments: provide an implementation of the GMRES solver for
!!  linear systems with the modified Gram-Schmidt orthogonalization.
!!  Implementation as in "Templates for the Solution of Linear Systems:
!!  Building Blocks for Iterative Methods
!!
!! @par Revision History
!!  Original version from http://www.netlib.org/templates/templates.pdf
!!  F90 rewriting by Marco Restelli.
!!  Adapted for use in ICOHDC by Hui Wan, MPI-M (2007-12-15)
!!  Cosmetic changes following the ICON programming guide by Hui Wan,
!!  MPI-M (2007-12-18)
!!  Included blocking by Marco Restelli (2008-09-23)
!!  Dummy argument for preconditioner shifted to last position
!!  and made optional, Marco Giorgetta (2009-02-14)
!!  Inlining of former functions active_dot and norm2 to simplify OpenMP
!!  parallelization, reduction of parallel sections, Guenther Zaengl (2009-06-16)
!! - For use within ocean model input parameter of type "interpolation-state" removed
!!   and adapted to 2D arrays (2010-04)
!!
!! mpi parallelized LL
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
!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_oce_linear_solver
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
!
!
!
USE mo_kind,                ONLY: wp
USE mo_parallel_config,     ONLY: nproma, p_test_run
USE mo_run_config,          ONLY: ltimer
!USE mo_impl_constants,      ONLY: sea_boundary
USE mo_model_domain,        ONLY: t_patch
#ifndef __SX__
USE mo_timer,               ONLY: timer_start, timer_stop, timer_gmres
#endif
USE mo_sync,                ONLY: omp_global_sum_array
USE mo_sync,                ONLY: sync_e, sync_c, sync_v, sync_patch_array
USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

! !PUBLIC MEMBER FUNCTIONS/SUBROUTINES
CHARACTER(len=*), PARAMETER :: this_mod_name = 'mo_oce_linear_solver'




!-------------------------------------------------------------------------
END MODULE mo_oce_linear_solver
