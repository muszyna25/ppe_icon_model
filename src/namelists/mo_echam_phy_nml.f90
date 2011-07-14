!>
!! Main switches of the ECHAM physics package, for turning on/off
!! the main parameterized processes.
!!
!! All the switches have default values, but can be changed by model
!! user via namelist "echam_phy_nml". This module contains a subroutine
!! that reads the namelist and makes necessary modifications of the
!! default or user-specified values. Note that output of the namelist values
!! to the stdout and the ASCII file is not done here, but after
!! all namelists have been read in and the cross-check has been finished.
!!
!! @author Hui Wan, MPI-M
!! @author Marco, Giorgetta, MPI-M
!!
!! @par Revision History
!! -Variables taken from mo_control and mo_param_switches of ECHAM6
!!  by Hui Wan, MPI (2010-07)
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
MODULE mo_echam_phy_nml

  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_namelist,           ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml
  USE mo_master_nml,         ONLY: lrestart
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_echam_phy_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  LOGICAL :: nml_lrad       !< .true. for radiation.
  LOGICAL :: nml_lvdiff     !< .true. for vertical diffusion.
  LOGICAL :: nml_lconv      !< .true. for moist convection
  LOGICAL :: nml_lcond      !< .true. for large scale condensation
  LOGICAL :: nml_lcover     !< .true. for prognostic cloud cover scheme
  LOGICAL :: nml_llandsurf  !< .true. for surface exchanges. (lsurf in ECHAM6)
  LOGICAL :: nml_lssodrag   !< .true. for subgrid scale orographic drag,
                            !< by blocking and gravity waves (lgwdrag in ECHAM6)
  LOGICAL :: nml_lgw_hines  !< .true. for atmospheric gravity wave drag
  LOGICAL :: nml_lice       !< .true. for sea-ice temperature calculation
  LOGICAL :: nml_lmeltpond  !< .true. for calculation of meltponds
  LOGICAL :: nml_lmlo       !< .true. for mixed layer ocean
  LOGICAL :: nml_lhd        !< .true. for hydrologic discharge model
  LOGICAL :: nml_lmidatm    !< .true. for middle atmosphere model version

  NAMELIST /echam_phy_nml/ nml_lrad, nml_lvdiff, nml_lconv, nml_lcond, &
                         & nml_lcover, nml_lssodrag, nml_lgw_hines,    &
                         & nml_llandsurf, nml_lice, nml_lmeltpond,     &
                         & nml_lmlo, nml_lhd, nml_lmidatm

CONTAINS
  !>
  !!
  SUBROUTINE read_echam_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    !----------------------------------------------------------------
    ! Set default values
    !----------------------------------------------------------------
    nml_lrad      = .TRUE.
    nml_lvdiff    = .TRUE.
    nml_lconv     = .TRUE.
    nml_lcond     = .TRUE.
    nml_lcover    = .FALSE.
    nml_llandsurf = .FALSE.
    nml_lssodrag  = .FALSE.
    nml_lgw_hines = .FALSE.
    nml_lice      = .FALSE.
    nml_lmeltpond = .FALSE.
    nml_lmlo      = .FALSE.
    nml_lhd       = .FALSE.
    nml_lmidatm   = .FALSE.

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !----------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('echam_phy_nml')
      READ(funit,NML=echam_phy_nml)
      CALL close_tmpfile(funit)
    END IF
                                                                                          
    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processors)
    !---------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('echam_phy_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, echam_phy_nml)
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=echam_phy_nml)
    CALL store_and_close_namelist(funit, 'echam_phy_nml')

    !-----------------------------------------------------
    ! Fill the configuration state
    !-----------------------------------------------------
    echam_phy_config% lrad      = nml_lrad                                                
    echam_phy_config% lvdiff    = nml_lvdiff                                              
    echam_phy_config% lconv     = nml_lconv                                               
    echam_phy_config% lcond     = nml_lcond                                               
    echam_phy_config% lcover    = nml_lcover                                              
    echam_phy_config% llandsurf = nml_llandsurf                                           
    echam_phy_config% lssodrag  = nml_lssodrag                                            
    echam_phy_config% lgw_hines = nml_lgw_hines                                           
    echam_phy_config% lice      = nml_lice                                                
    echam_phy_config% lmeltpond = nml_lmeltpond                                           
    echam_phy_config% lmlo      = nml_lmlo                                                
    echam_phy_config% lhd       = nml_lhd                                                 
    echam_phy_config% lmidatm   = nml_lmidatm  

  END SUBROUTINE read_echam_phy_namelist

END MODULE mo_echam_phy_nml
