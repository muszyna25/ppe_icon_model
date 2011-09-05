!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-09-05)
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
MODULE mo_nh_pzlev_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for p/p-level output
  !!--------------------------------------------------------------------------
  TYPE :: t_nh_pzlev_config

    ! namelist variables
    !
    LOGICAL :: lwrite_zlev           !< if .TRUE.: write output on z-levels

    LOGICAL :: lwrite_plev           !< if .TRUE.: write output on z-levels

    INTEGER :: nzlev                 !< number of z-levels

    INTEGER :: nplev                 !< number of p-levels

    INTEGER :: zlevels(100)          !< zlevel heights [m] 

    INTEGER :: plevels(100)          !< plevel heights [m] 

    ! derived variables
    !
                   
  END TYPE t_nh_pzlev_config

  !>
  !!
  TYPE(t_nh_pzlev_config), TARGET :: nh_pzlev_config(0:max_dom)



!!$CONTAINS


END MODULE mo_nh_pzlev_config
