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
MODULE mo_sleve_config

  USE mo_kind,                ONLY: wp
 !USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: min_lay_thckn, stretch_fac, top_height
  PUBLIC :: decay_scale_1, decay_scale_2, decay_exp, flat_height

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'
  !>
  !!--------------------------------------------------------------------------
  !! Type definition 
  !!--------------------------------------------------------------------------
  !TYPE :: t_sleve_config

    ! a) Parameters specifying the distrubution of the coordinate surfaces
    !     (the initializations are a workaround for a NEC compiler bug)

    REAL(wp):: min_lay_thckn  ! Layer thickness of lowermost level
    REAL(wp):: stretch_fac    ! Factor for stretching/squeezing the model layer distribution
    REAL(wp):: top_height     ! Height of model top

    ! b) Parameters for SLEVE definition

    REAL(wp):: decay_scale_1  ! Decay scale for large-scale topography component
    REAL(wp):: decay_scale_2  ! Decay scale for small-scale topography component
    REAL(wp):: decay_exp      ! Exponent for decay function
    REAL(wp):: flat_height    ! Height above which the coordinate surfaces are exactly flat
                              ! additional feature not available in the standard
                              ! SLEVE definition
  !END TYPE t_sleve_config
  !>
  !!
  !TYPE(t_sleve_config) :: sleve_config(max_dom)

END MODULE mo_sleve_config
