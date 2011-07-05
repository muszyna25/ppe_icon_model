!>
!! Contains the setup of the sleve coordinate
!!
!!        
!! @par Revision History
!!   Revision History in mo_global_variables.f90 (r3919)
!!   Modification by Constantin Junk (2011-03-29)
!!   - added new module mo_sleve_nml.f90 which includes
!!     setting up the SLEVE namelist.
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_sleve_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
!
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish
  USE mo_impl_constants,     ONLY: max_char_length
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_namelist,           ONLY: position_nml, positioned
  USE mo_mpi,                ONLY: p_pe, p_io

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC


  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the SLEVE coordinate
  ! ----------------------------------------------------------------------------
  !
  ! a) Parameters specifying the distrubution of the coordinate surfaces
  !     (the initializations are a workaround for a NEC compiler bug)
  REAL(wp):: min_lay_thckn = 1._wp  ! Layer thickness of lowermost level
  REAL(wp):: stretch_fac   = 1._wp  ! Factor for stretching/squeezing the model layer distribution
  REAL(wp):: top_height    = 1._wp  ! Height of model top

  ! b) Parameters for SLEVE definition
  REAL(wp):: decay_scale_1 = 1._wp  ! Decay scale for large-scale topography component
  REAL(wp):: decay_scale_2 = 1._wp  ! Decay scale for small-scale topography component
  REAL(wp):: decay_exp     = 1._wp  ! Exponent for decay function
  REAL(wp):: flat_height   = 1._wp  ! Height above which the coordinate surfaces are exactly flat
                            ! additional feature not available in the standard SLEVE definition

  NAMELIST /sleve_ctl/ min_lay_thckn, top_height, decay_scale_1,     &
                       decay_scale_2, decay_exp, flat_height, stretch_fac
  !
  !
  !

CONTAINS

!-------------------------------------------------------------------------
!
!
 !>
 !!  Initialization of the SLEVE coordinate namelist
 !!
 !!
 !! @par Revision History
 !!  Initial version by Guenther Zaengl (2010-07-21)

 SUBROUTINE sleve_nml_setup

!   CHARACTER(len=max_char_length), PARAMETER :: &
!             routine = 'mo_sleve_nml/sleve_nml_setup:'


  !local variable
  INTEGER :: i_status

  !------------------------------------------------------------
  ! 2.0 set up the default values for dynamics_ctl
  !------------------------------------------------------------
  !
  !
  ! a) Parameters determining the distribution of model layers
  !    (if not read in from a table)
  min_lay_thckn   = 50._wp      ! Layer thickness of lowermost layer
  top_height      = 23500._wp   ! Height of model top
  stretch_fac     = 1._wp       ! Scaling factor for stretching/squeezing the model layer distribution

  ! b) Parameters setting up the decay function of the topographic signal
  decay_scale_1   = 4000._wp    ! Decay scale of large-scale topography component
  decay_scale_2   = 2500._wp    ! Decay scale of small-scale topography component
  decay_exp       = 1.2_wp      ! Exponent for decay function
  flat_height     = 16000._wp   ! Height above which the coordinate surfaces are flat
  !
  !
  !------------------------------------------------------------
  ! 3.0 Read the nonhydrostatic namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL position_nml ('sleve_ctl', status=i_status)
  SELECT CASE (i_status)
  CASE (positioned)
     READ (nnml, sleve_ctl)
  END SELECT
  !
  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  !currently no consistency check...

  ! write the contents of the namelist to an ASCII file

  IF(p_pe == p_io) WRITE(nnml_output,nml=sleve_ctl)

END SUBROUTINE sleve_nml_setup

END MODULE mo_sleve_nml
