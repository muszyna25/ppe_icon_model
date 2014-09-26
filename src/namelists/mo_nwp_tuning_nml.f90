!>
!! @brief Namelist for tuning and/or perturbing nwp physics
!!
!! These subroutines are called by read_atmo_namelists and do some 
!! nwp physics tuning 
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2014-09-25)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_tuning_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_nwp_tuning_config,   ONLY: config_tune_gkwake    => tune_gkwake,    &
    &                               config_tune_gkdrag    => tune_gkdrag,    &
    &                               config_tune_gfluxlaun => tune_gfluxlaun, &
    &                               config_tune_zceff_min => tune_zceff_min, &
    &                               config_tune_v0snow    => tune_v0snow,    &
    &                               config_tune_zvz0i     => tune_zvz0i !!$,     &  
!!$    &                               config_itune_albedo   => itune_albedo  
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_nwp_tuning_namelist


  !-----------------------------------!
  ! nwp_tuning_nml namelist variables !
  !-----------------------------------!

  REAL(wp) :: &                    !< low level wake drag constant
    &  tune_gkwake

  REAL(wp) :: &                    !< gravity wave drag constant
    &  tune_gkdrag

  REAL(wp) :: &                    !< total launch momentum flux in each azimuth (rho_o x F_o)
    &  tune_gfluxlaun

  REAL(wp) :: &                    !< Minimum value for sticking efficiency
    &  tune_zceff_min

  REAL(wp) :: &                    !< factor in the terminal velocity for snow
    &  tune_v0snow

  REAL(wp) :: &                    !< Terminal fall velocity of ice 
    &  tune_zvz0i

!!$  REAL(wp) :: &                    !< (MODIS) albedo tuning
!!$    &  itune_albedo                ! 1: dimmed sahara

  NAMELIST/nwp_tuning_nml/ tune_gkwake, tune_gkdrag, tune_gfluxlaun, &
    &                      tune_zceff_min, tune_v0snow, tune_zvz0i !!$, &
!!$    &                      itune_albedo


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP physics tuning. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics tuning
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !! Initial Revision by Daniel Reinert, DWD (2014-09-25)
  !!
  SUBROUTINE read_nwp_tuning_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nwp_tuning_nml: read_tuning_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    ! Comment: In case we want to draw from a normal distribution, the namelist 
    ! parameters could be extended to arrays of size 2. The first value is the mean, 
    ! while the second one is the satndard deviation. 

    ! SSO tuning
    tune_gkwake     = 1.333_wp     ! original COSMO value 0.5
    tune_gkdrag     = 0.1_wp       ! original COSMO value 0.075
    !
    ! GWD tuning
    tune_gfluxlaun  = 2.50e-3_wp   ! original IFS value 3.75e-3
    !
    ! grid scale microphysics
    tune_zceff_min  = 0.075_wp
    tune_v0snow     = 25.0_wp      ! previous ICON value was 20
    tune_zvz0i      = 1.25_wp      ! original value of Heymsfield+Donner 1990: 3.29 

!!$    itune_albedo    = 0            ! original (measured) albedo


    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('nwp_tuning_nml')
      READ(funit,NML=nwp_tuning_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_tuning_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nwp_tuning_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nwp_tuning_nml)             ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nwp_tuning_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------



    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      config_tune_gkwake    = tune_gkwake
      config_tune_gkdrag    = tune_gkdrag
      config_tune_gfluxlaun = tune_gfluxlaun
      config_tune_zceff_min = tune_zceff_min 
      config_tune_v0snow    = tune_v0snow
      config_tune_zvz0i     = tune_zvz0i
!!$      config_itune_albedo   = itune_albedo
    ENDDO



    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nwp_tuning_nml)                    
      CALL store_and_close_namelist(funit, 'nwp_tuning_nml')             
    ENDIF

    !--------------------------------------------------------
    ! 7. write the contents of the namelist to an ASCII file
    !--------------------------------------------------------
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_tuning_nml)


  END SUBROUTINE read_nwp_tuning_namelist

END MODULE mo_nwp_tuning_nml
