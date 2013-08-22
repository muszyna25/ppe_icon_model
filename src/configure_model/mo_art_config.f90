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
!! Modifications by Kristina Lundgren, KIT (2012-07-03)
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
  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: max_dom
  USE mo_math_utilities,       ONLY: t_geographical_coordinates  
  IMPLICIT NONE


  PUBLIC 


  CHARACTER(len=*),PARAMETER :: version = '$Id$'


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-ART
  !!--------------------------------------------------------------------------
   INTEGER, PARAMETER  :: max_volc_input  = 20 !Maximum number of volcanoes in input namelist art_volclist_tot
   INTEGER             :: nart_tendphy  = 0    !Maximum number of tracers that are effected by deep convective transport 

  TYPE t_volc_list
    CHARACTER(len=20)                :: zname    ! < name of volcanoe or location
    TYPE(t_geographical_coordinates) :: location !< geographical position
  END TYPE t_volc_list

  TYPE t_art_config

    ! namelist variables
    ! General
    LOGICAL :: lart                        !< main switch for using the ART-package
                                           !< .TRUE.: switch ON
                                           !< .FALSE.: switch OFF
    ! Sea Salt Aerosol
    LOGICAL :: lart_seasalt        !< Treatment of sea salt aerosol (TRUE/FALSE)
    
    ! Mineral Dust Aerosol
    LOGICAL :: lart_dust           !< Treatment of mineral dust aerosol (TRUE/FALSE)
    
    ! Processes
    LOGICAL :: lart_emiss          !< Emission of aerosol (TRUE/FALSE)

    LOGICAL :: lart_conv           !< Convection of aerosol (TRUE/FALSE)

    LOGICAL :: lart_wash           !< Washout of aerosol (TRUE/FALSE)

    LOGICAL :: lart_rad            !< Radiative impact of aerosol (TRUE/FALSE)

    LOGICAL :: lart_cloud          !< Cloud aerosol interaction (TRUE/FALSE)

    ! Volcanic Ash Aerosol
    LOGICAL :: lart_volcano        !< Treatment of volcanic ash (TRUE/FALSE)

    INTEGER :: nart_emis_volcano_update    !< Time interval for reading volcano emission file

    LOGICAL :: lart_volclist    !< Use list of volcanoes from namelist (TRUE/FALSE)
    
   INTEGER                     :: nvolc             !< Number ov volcanoes
   
   TYPE(t_volc_list), POINTER  :: volclist(:,:)    !< (idx,blk)
   
   CHARACTER (LEN=120) :: volcanofile_path  !< Path of volcano file
   
   ! for radioactive tracers

   LOGICAL :: lart_radioact                !< Treatment of radioactive nuclides (TRUE/FALSE)

   LOGICAL :: lart_decay_radioact          !< Treatment of radioactive decay (TRUE/FALSE)   
   
   CHARACTER (LEN=120) :: radioactfile_path !< Path of emission file for radiactive nuclides

   !For specification of locations.

   INTEGER                     :: nblks,npromz              
   INTEGER                     :: nconv_tracer     ! number of tracers in convection 

  END TYPE t_art_config

  !>
  !!
  TYPE(t_art_config), TARGET :: art_config(0:max_dom)


CONTAINS

  !>
  !! setup components of ICON-ART depending on this namelist
  !!
  !! Setup of additional ICON-ART control variables depending on the 
  !! art-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-12-08)
  !! Modification by kristina Lundgren, KIT (2012-11-27)
  !
  SUBROUTINE configure_art(jg)
  !
    INTEGER, INTENT(IN) :: jg           !< patch 

    !-----------------------------------------------------------------------
  art_config(jg)%nconv_tracer=0
 
  END SUBROUTINE configure_art


END MODULE mo_art_config
