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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gridref_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: t_gridref_config, gridref_config

  PUBLIC :: rbf_vec_kern_grf_e, rbf_scale_grf_e,                                &
    &                    grf_velfbk, grf_scalfbk, grf_tracfbk,                  &
    &                    grf_idw_exp_e12, grf_idw_exp_e34,                      &
    &                    grf_intmethod_c, grf_intmethod_e,                      &
    &                    grf_intmethod_ct, denom_diffu_v, denom_diffu_t,        &
    &                    l_mass_consvcorr, l_density_nudging, fbk_relax_timescale

  !--------------------------------------------------------------------------
  ! Basic configuration setup for grid refinement
  !--------------------------------------------------------------------------
!  TYPE t_gridref_config

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

    LOGICAL  :: l_mass_consvcorr  ! .true.: apply mass conservation correction
    LOGICAL  :: l_density_nudging ! .true.: apply density nudging near lateral nest boundaries if feedback is turned on
                                  ! (in case of one-way nesting, all prognostic variables are nudged irrespective of this switch)

    ! Relaxation time scale for feedback in case of ifeedback_type = 2
    REAL(wp) :: fbk_relax_timescale
!  END TYPE t_gridref_config
  !>
  !!
!  TYPE(t_gridref_config) :: gridref_config(max_dom)


END MODULE mo_gridref_config
