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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gribout_config

  USE mo_impl_constants,     ONLY: max_phys_dom
  USE mo_grib2_tile,         ONLY: t_grib2_template_tile

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: t_gribout_config
  PUBLIC :: gribout_config
  PUBLIC :: configure_gribout

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for grib output
  !!--------------------------------------------------------------------------
  TYPE :: t_gribout_config

    ! namelist variables

    INTEGER :: tablesVersion              ! Main switch for table version

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
                                          ! 196: invariant data 

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

    INTEGER :: &                          ! Table: local.78.254.def
      & localDefinitionNumber             ! 252: Ensemble system incl. postprocessing 
                                          ! 253: Ensemble system 
                                          ! 254: Deterministic system

    INTEGER :: &                          ! Table: local.78.254.def
      & localNumberOfExperiment           !  


    INTEGER :: &                          ! Output generating center
      & generatingCenter                  ! DWD  : 78
                                          ! ECMWF: 98 


    INTEGER :: &                          ! Output generating subcenter
      & generatingSubcenter               ! DWD  : 255
                                          ! ECMWF: 0 

    LOGICAL :: lspecialdate_invar         ! .TRUE.: use special date 10101 for encoding
                                          ! invariant and climatological fields
                                          ! .FALSE.: no special treatment of invariant
                                          ! and climatological fields.

    LOGICAL :: ldate_grib_act             ! add Creation date to GRIB file
                                          ! .TRUE. : activated
                                          ! .FALSE.: deactivated (use dummy date/time) 

    LOGICAL :: lgribout_24bit             ! write thermodynamic fields rho, theta_v, T, p
                                          ! with 24bit precision


    INTEGER :: typeOfEnsembleForecast,        &
      &        localTypeOfEnsembleForecast,   &
      &        numberOfForecastsInEnsemble,   &
      &        perturbationNumber

    CHARACTER(len=32) ::  &               !< type of GRIB2 templates used for surface tile fields
      &  typeOfGrib2TileTemplate          !  'wmo': official WMO templates 55, 59, ...
                                          !  'dwd': local 'DWD' templates 40455, 40456, ...

    ! derived variables
    !
    TYPE(t_grib2_template_tile)::  &      !< defines set of employed GRIB2 tile templates for writing
      &  grib2_template_tile              !< contains tile template numbers as well as the corresponding 
                                          !< set GRIB2 tile keys.
                                          !< variant 1: templates 40455, 40456, etc
                                          !<            DWD's local tile templates
                                          !< variant 2: 55, 59, etc
                                          !< official WMO tile templates
                                          !< the variant in use is determined by the namelist parameter itype_tiletemplate  


  END TYPE t_gribout_config

  !>
  !!
  TYPE(t_gribout_config), TARGET :: gribout_config(1:max_phys_dom)


CONTAINS


  !>
  !! potentially modify generatingCenter and generatingSubcenter 
  !!
  !! If generatingCenter and generatingSubcenter are not set via namelist, 
  !! they are filled with values read from the grid file. 
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-04-26)
  !!
  SUBROUTINE configure_gribout(grid_generatingCenter, grid_generatingSubcenter, &
    &                          n_dom)
  !
    INTEGER,       INTENT(IN)  :: grid_generatingCenter(0:)
    INTEGER,       INTENT(IN)  :: grid_generatingSubcenter(0:)
    INTEGER,       INTENT(IN)  :: n_dom

    ! local fields
    INTEGER  :: jg

    !-----------------------------------------------------------------------

    DO jg = 1, n_dom
      !
      ! check, whether generatingCenter was set in gribout_nml
      !
      IF ( gribout_config(jg)%generatingCenter == -1 ) THEN
        ! If not, then fill with grid generating center
        gribout_config(jg)%generatingCenter = grid_generatingCenter(jg)
      ENDIF

      ! check, whether generatingSubcenter was set in gribout_nml
      !
      IF ( gribout_config(jg)%generatingSubcenter == -1 ) THEN
        ! If not, then fill with grid generating subcenter
        gribout_config(jg)%generatingSubcenter = grid_generatingSubcenter(jg)
      ENDIF

    ENDDO  ! jg

  END SUBROUTINE configure_gribout

END MODULE mo_gribout_config
