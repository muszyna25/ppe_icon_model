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
!! 2002-2010 by DWD and MPI-M
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

  USE mo_namelist,       ONLY: position_nml, POSITIONED
  USE mo_io_units,       ONLY: nnml
! USE mo_master_nml,     ONLY: lrestart
  USE mo_io_restart_namelist,ONLY: open_tmpfile, store_and_close_namelist !,   &                                             
!                                & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  LOGICAL :: lrad      = .TRUE.   !<  .true. for radiation.
  LOGICAL :: lvdiff    = .TRUE.   !<  .true. for vertical diffusion.
  LOGICAL :: lconv     = .TRUE.   !<  .true. for moist convection
  LOGICAL :: lcond     = .TRUE.   !<  .true. for large scale condensation
  LOGICAL :: lcover    = .FALSE.  !<  .true. for prognostic cloud cover scheme
  LOGICAL :: llandsurf = .FALSE.  !<  .true. for surface exchanges. (lsurf in ECHAM6)
  LOGICAL :: lssodrag  = .FALSE.  !<  .true. for subgrid scale orographic drag,
  !                                    by blocking and gravity waves (lgwdrag in ECHAM6)
  LOGICAL :: lgw_hines  = .FALSE. !<  .true. for atmospheric gravity wave drag
  LOGICAL :: lice      = .FALSE.  !<  .true. for sea-ice temperature calculation
  LOGICAL :: lmeltpond = .FALSE.  !<  .true. for calculation of meltponds
  LOGICAL :: lmlo      = .FALSE.  !<  .true. for mixed layer ocean
  LOGICAL :: lhd       = .FALSE.  !<  .true. for hydrologic discharge model
  LOGICAL :: lmidatm   = .FALSE.  !<  .true. for middle atmosphere model version

  NAMELIST /echam_phy_nml/ lrad, lvdiff, lconv, lcond, lcover, lssodrag, lgw_hines, &
    &                      llandsurf, lice, lmeltpond, lmlo, lhd, lmidatm

CONTAINS
  !>
  !!
  SUBROUTINE read_echam_phy_nml

    INTEGER :: istat, funit

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values used in the previous integration.
    !----------------------------------------------------------------
!   IF (lrestart) THEN
!     funit = open_and_restore_namelist('echam_phy_nml')
!     READ(funit,NML=echam_phy_nml)
!     CALL close_tmpfile(funit)
!    ! for testing
!     WRITE (0,*) 'contents of namelist ...'
!     WRITE (0,NML=echam_phy_nml)
!   END IF
                                                                                          
    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !---------------------------------------------------------------------
    CALL position_nml ('echam_phy_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, echam_phy_nml)
    END SELECT

    !-----------------------------------------------------                                
    ! Store the namelist for restart                                                      
    !-----------------------------------------------------                                
    funit = open_tmpfile()                                                                
    WRITE(funit,NML=echam_phy_nml)                                                             
    CALL store_and_close_namelist(funit, 'echam_phy_nml')                                      
    write(0,*) 'stored echam_phy_nml'

  END SUBROUTINE read_echam_phy_nml

END MODULE mo_echam_phy_nml
