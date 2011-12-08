!>
!! @brief configuration setup for ART-package
!!
!! configuration setup for ART-package
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-12-08)
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
MODULE mo_art_config

  USE mo_impl_constants     ,ONLY: max_dom

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-ART
  !!--------------------------------------------------------------------------
  TYPE :: t_art_config

    ! namelist variables
    !    
    LOGICAL :: lart                 !< main switch for using the ART-package
                                    !< .TRUE.: switch ON
                                    !<.FALSE.: switch OFF

    ! derived variables

  END TYPE t_art_config

  !>
  !!
  TYPE(t_art_config) :: art_config(0:max_dom)


!!$CONTAINS

!!$  !>
!!$  !! setup components of ICON-ART depending on this namelist
!!$  !!
!!$  !! Setup of additional ICON-ART control variables depending on the 
!!$  !! art-NAMELIST and potentially other namelists. This routine is 
!!$  !! called, after all namelists have been read and a synoptic consistency 
!!$  !! check has been done.
!!$  !!
!!$  !! @par Revision History
!!$  !! Initial revision by Daniel Reinert, DWD (2011-12-08)
!!$  !!
!!$  SUBROUTINE configure_art( )
!!$  !
!!$    INTEGER, INTENT(IN) :: jg           !< patch 
!!$
!!$    !-----------------------------------------------------------------------
!!$
!!$
!!$  END SUBROUTINE configure_art


END MODULE mo_art_config
