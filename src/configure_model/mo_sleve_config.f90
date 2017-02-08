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
MODULE mo_sleve_config

  USE mo_kind,                ONLY: wp
 !USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: itype_laydistr, min_lay_thckn, max_lay_thckn, htop_thcknlimit, stretch_fac, top_height
  PUBLIC :: decay_scale_1, decay_scale_2, decay_exp, flat_height
  PUBLIC :: lread_smt
  !>
  !!--------------------------------------------------------------------------
  !! Type definition 
  !!--------------------------------------------------------------------------
  !TYPE :: t_sleve_config

    ! a) Parameters specifying the distrubution of the coordinate surfaces

    INTEGER :: itype_laydistr ! Type of analytical function used for computing the coordinate surface distribution
    REAL(wp):: min_lay_thckn  ! Layer thickness of lowermost level
    REAL(wp):: max_lay_thckn  ! Maximum layer thickness below htop_thcknlimit
    REAL(wp):: htop_thcknlimit! Height below which the layer thickness must not exceed max_lay_thckn
    REAL(wp):: stretch_fac    ! Factor for stretching/squeezing the model layer distribution
    REAL(wp):: top_height     ! Height of model top

    ! b) Parameters for SLEVE definition

    REAL(wp):: decay_scale_1  ! Decay scale for large-scale topography component
    REAL(wp):: decay_scale_2  ! Decay scale for small-scale topography component
    REAL(wp):: decay_exp      ! Exponent for decay function
    REAL(wp):: flat_height    ! Height above which the coordinate surfaces are exactly flat
                              ! additional feature not available in the standard
                              ! SLEVE definition

    ! c) Parameter for reading in smoothed topography
    LOGICAL :: lread_smt

  !END TYPE t_sleve_config
  !>
  !!
  !TYPE(t_sleve_config) :: sleve_config(max_dom)

END MODULE mo_sleve_config
