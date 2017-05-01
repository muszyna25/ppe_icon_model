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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_phy_nml

  USE mo_kind,               ONLY: wp
  USE mo_echam_phy_config,   ONLY: echam_phy_config
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &
                                 & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_echam_phy_namelist

  INTEGER  :: idcphycpl  !< determines the coupling between the dynamical core and the
                         !  phyiscs package
                         !  1: dynamics and physics update sequentially
                         !  2: dynamics uses physics forcing for updating
  LOGICAL  :: ldrymoist  !   .true.  physics assumes moist air dynamics
  LOGICAL  :: lrad       !< .true. for radiation.
  REAL(wp) :: dt_rad   !! "-"                     radiation
  LOGICAL  :: lvdiff     !< .true. for vertical diffusion.
  LOGICAL  :: lconv      !< .true. for moist convection
  LOGICAL  :: lcond      !< .true. for large scale condensation
  LOGICAL  :: lssodrag   !< .true. for subgrid scale orographic drag,
                         !< by blocking and gravity waves (lgwdrag in ECHAM6)
  LOGICAL  :: lgw_hines  !< .true. for atmospheric gravity wave drag
  LOGICAL  :: lcariolle  !< .true. for Cariolle ozone scheme (you need a
                         !< transported ozone tracer then!)
  LOGICAL  :: lmlo       !< .true. for mixed layer ocean
  LOGICAL  :: lice       !< .true. for sea-ice temperature calculation
  LOGICAL  :: ljsbach    !< .true. for calculating the JSBACH land surface
  LOGICAL  :: llake      !< .true. for using lakes in JSBACH
  LOGICAL  :: lamip      !< .true. for AMIP simulations
  LOGICAL  :: lebudget   !< .true. for echam physcics energy budget calculation


  NAMELIST /echam_phy_nml/ idcphycpl, ldrymoist,    &
    &                      lrad, dt_rad, lvdiff,    &
    &                      lconv, lcond,            &
    &                      lssodrag, lgw_hines,     &
    &                      lcariolle,               &
    &                      lmlo, lice, ljsbach,     &
    &                      llake, lamip, lebudget

CONTAINS
  !>
  !!
  SUBROUTINE read_echam_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: istat
    INTEGER :: funit
    INTEGER :: iunit

    !------------------------------------------------------------------
    ! 1. Set default values
    !------------------------------------------------------------------
    idcphycpl = 1         ! 1: dynamics and physics update sequentially
    ldrymoist = .TRUE.    ! TRUE.: dry air mass conservation; .FALSE.: moist air mass cons.
    lrad      = .TRUE.
    dt_rad    = 3600.0_wp ! [s]
    lvdiff    = .TRUE.
    lconv     = .TRUE.
    lcond     = .TRUE.
    lssodrag  = .TRUE.
    lgw_hines = .TRUE.
    lcariolle = .FALSE.
    lmlo      = .FALSE.
    lice      = .FALSE.
    ljsbach   = .FALSE.
    llake     = .FALSE.
    lamip     = .FALSE.
    lebudget  = .FALSE.

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('echam_phy_nml')
      READ(funit,NML=echam_phy_nml)
      CALL close_tmpfile(funit)
    END IF
                                                                                          
    !------------------------------------------------------------------
    ! 3. Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('echam_phy_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, echam_phy_nml)      ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, echam_phy_nml)       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, echam_phy_nml)    ! write settings to temporary text file
      END IF
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
    echam_phy_config% idcphycpl = idcphycpl
    echam_phy_config% ldrymoist = ldrymoist
    echam_phy_config% lrad      = lrad                                                
    echam_phy_config% dt_rad    = dt_rad
    echam_phy_config% lvdiff    = lvdiff                                              
    echam_phy_config% lconv     = lconv                                               
    echam_phy_config% lcond     = lcond                                               
    echam_phy_config% lssodrag  = lssodrag                                            
    echam_phy_config% lgw_hines = lgw_hines                                           
    echam_phy_config% lcariolle = lcariolle
    echam_phy_config% lmlo      = lmlo                                                
    echam_phy_config% lice      = lice                                                
    echam_phy_config% ljsbach   = ljsbach
    echam_phy_config% llake     = llake
    echam_phy_config% lamip     = lamip                                                
    echam_phy_config% lebudget  = lebudget

  END SUBROUTINE read_echam_phy_namelist

END MODULE mo_echam_phy_nml
