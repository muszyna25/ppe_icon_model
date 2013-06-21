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

  USE mo_kind,               ONLY: wp
  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_master_control,     ONLY: is_restart_run
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings, log_nml_settings 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_echam_phy_namelist

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  LOGICAL :: lrad       !< .true. for radiation.
  LOGICAL :: lvdiff     !< .true. for vertical diffusion.
  LOGICAL :: lconv      !< .true. for moist convection
  LOGICAL :: lcond      !< .true. for large scale condensation
  LOGICAL :: lcover     !< .true. for prognostic cloud cover scheme
  LOGICAL :: llandsurf  !< .true. for surface exchanges. (lsurf in ECHAM6)
  LOGICAL :: lssodrag   !< .true. for subgrid scale orographic drag,
                        !< by blocking and gravity waves (lgwdrag in ECHAM6)
  LOGICAL :: lgw_hines  !< .true. for atmospheric gravity wave drag
  LOGICAL :: lice       !< .true. for sea-ice temperature calculation
  LOGICAL :: lmeltpond  !< .true. for calculation of meltponds
  LOGICAL :: lamip      !< .true. for AMIP simulations
  LOGICAL :: lmlo       !< .true. for mixed layer ocean
  LOGICAL :: ljsbach    !< .true. for calculating the JSBACH land surface
  LOGICAL :: lhd        !< .true. for hydrologic discharge model
!!$  LOGICAL :: lmidatm    !< .true. for middle atmosphere model version
  REAL(wp) :: dt_rad   !! "-"                     radiation


  NAMELIST /echam_phy_nml/ lrad, lvdiff, lconv, lcond,  &
                         & lcover, lssodrag, lgw_hines, &
                         & llandsurf, lice, lmeltpond,  &
                         & lamip, lmlo, ljsbach, lhd, dt_rad!!$, lmidatm

CONTAINS
  !>
  !!
  SUBROUTINE read_echam_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    !------------------------------------------------------------------
    ! 1. Set default values
    !------------------------------------------------------------------
    lrad      = .TRUE.
    lvdiff    = .TRUE.
    lconv     = .TRUE.
    lcond     = .TRUE.
    lcover    = .FALSE.
    llandsurf = .FALSE.
    lssodrag  = .FALSE.
    lgw_hines = .FALSE.
    lice      = .FALSE.
    lmeltpond = .FALSE.
    lamip     = .FALSE.
    lmlo      = .FALSE.
    ljsbach   = .FALSE.
    lhd       = .FALSE.
    dt_rad    = 3600.0_wp
!!$    lmidatm   = .FALSE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('echam_phy_nml')
      READ(funit,NML=echam_phy_nml)
      CALL close_tmpfile(funit)
    END IF
                                                                                          
    !------------------------------------------------------------------
    ! 3. Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('echam_phy_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (positioned)
      WRITE(temp_defaults(), echam_phy_nml)                     ! write defaults to temporary text file
      READ (nnml, echam_phy_nml, iostat=istat)                  ! overwrite default settings
      WRITE(temp_settings(), echam_phy_nml)                     ! write settings to temporary text file
      CALL log_nml_settings("nml.log")
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! 4. Sanity Check
    !------------------------------------------------------------------
    ! nothing to be done

    !------------------------------------------------------------------
    ! 5. Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=echam_phy_nml)
      CALL store_and_close_namelist(funit, 'echam_phy_nml')
    ENDIF

    !------------------------------------------------------------------
    ! 6. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=echam_phy_nml)

    !------------------------------------------------------------------
    ! 7. Fill the configuration state
    !------------------------------------------------------------------
    echam_phy_config% lrad      = lrad                                                
    echam_phy_config% lvdiff    = lvdiff                                              
    echam_phy_config% lconv     = lconv                                               
    echam_phy_config% lcond     = lcond                                               
    echam_phy_config% lcover    = lcover                                              
    echam_phy_config% llandsurf = llandsurf                                           
    echam_phy_config% lssodrag  = lssodrag                                            
    echam_phy_config% lgw_hines = lgw_hines                                           
    echam_phy_config% lice      = lice                                                
    echam_phy_config% lmeltpond = lmeltpond                                           
    echam_phy_config% lamip     = lamip                                                
    echam_phy_config% lmlo      = lmlo                                                
    echam_phy_config% ljsbach   = ljsbach
    echam_phy_config% lhd       = lhd                                                 
    echam_phy_config% dt_rad    = dt_rad
!!$    echam_phy_config% lmidatm   = lmidatm  

  END SUBROUTINE read_echam_phy_namelist

END MODULE mo_echam_phy_nml
