!>
!! @brief configuration setup for Grib output
!!
!! configuration setup for tracer transport
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2013-01-29)
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
MODULE mo_gribout_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom


  IMPLICIT NONE
  PRIVATE


  PUBLIC :: t_gribout_config
  PUBLIC :: gribout_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for grib output
  !!--------------------------------------------------------------------------
  TYPE :: t_gribout_config

    ! namelist variables

    INTEGER :: &                          ! Table 1.2
      & significanceOfReferenceTime       ! 0: Analysis
                                          ! 1: Start of forecast
                                          ! 2: Verifying time of forecast
                                          ! 4: ... 

    INTEGER :: &                          ! Table 1.3
      & productionStatusOfProcessedData   ! 0: Oper. products
                                          ! 1: Oper. test products
                                          ! 2: Research products
                                          ! 3: ...

    INTEGER :: &                          ! Table 1.4
      & typeOfProcessedData               ! 0: Analysis products
                                          ! 1: Forecast products
                                          ! 2: Analysis and forecast products
                                          ! 3: ...

    INTEGER :: &                          ! Table 4.3
      & typeOfGeneratingProcess           ! 0: Analysis
                                          ! 1: Initialization
                                          ! 2: Forecast
                                          ! 3: ...

    INTEGER :: &                          ! Table: backgroundProcess
      & backgroundProcess                 ! 0: main run
                                          ! 1: pre-assimilation
                                          ! 2: assimilation
                                          ! 3: ... 

    INTEGER :: &                          ! Table: generatingProcessIdentifier
      & generatingProcessIdentifier       ! 1: icogl
                                          ! 2: icrgl
                                          ! 3: icoeu
                                          ! 4: ...

    LOGICAL :: ldate_grib_act             ! add Creation date to GRIB file
                                          ! .TRUE. : activated
                                          ! .FALSE.: deactivated (use dummy date/time) 


    ! derived variables
 
  END TYPE t_gribout_config

  !>
  !!
  TYPE(t_gribout_config), TARGET :: gribout_config(0:max_dom)

END MODULE mo_gribout_config
