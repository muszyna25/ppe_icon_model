!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_gridref_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_gridref_config, gridref_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Basic configuration setup for grid refinement
  !--------------------------------------------------------------------------
  TYPE t_gridref_config

    INTEGER  :: rbf_vec_kern_grf_e ! rbf kernel for vector interpolation

    ! scale factors for rbf grid refinement interpolation
    REAL(wp) :: rbf_scale_grf_e(max_dom)

    INTEGER  :: grf_intmethod_c,  &  ! switch for type of grid refinement interpolation
      &         grf_intmethod_ct, &  ! (see below for explanation of options)
      &         grf_intmethod_e

    INTEGER  :: grf_velfbk     ! switch for velocity feedback method
                               ! 1 = averaging over child edges 1 and 2;
                               ! 2 = 2nd-order method using RBF reconstruction to child vertices
  
    INTEGER  :: grf_scalfbk    ! switch for feedback method of scalar dynamical variables
                               ! 1 = area-weighted averaging
                               ! 2 = bilinear interpolation

    INTEGER  :: grf_tracfbk    ! switch for feedback method of passive tracer variables
                               ! 1 = area-weighted averaging
                               ! 2 = bilinear interpolation

    ! Exponents for IDW interpolation in idw_compute_coeff_grf
    REAL(wp) :: grf_idw_exp_e12, grf_idw_exp_e34

    ! Denominators of normalized diffusion coefficients for boundary diffusion
    REAL(wp) :: denom_diffu_v, denom_diffu_t

  END TYPE t_gridref_config
  !>
  !!
  TYPE(t_gridref_config) :: gridref_config(max_dom)

END MODULE mo_gridref_config
