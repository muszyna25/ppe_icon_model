!>
!! @brief namelist setup for HAMOCC
!!
!! Namelist setup for HAMOCC, the HAMburg Ocean Carbon Cycle model
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_hamocc_nml

  USE mo_kind,                ONLY: wp
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output, find_next_free_unit
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist, &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_exception,           ONLY: finish, message
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_hamocc_namelist

  INTEGER,PUBLIC :: i_settling                  !< Switch for POC settling parameterisation
                                                !< 0 constant sinking speed
                                                !< 1 variable sinking speed following the 'Martin curve'
                                                !< 2 sinking via aggregation
  INTEGER, PUBLIC :: isac 
  REAL(wp), PUBLIC :: sinkspeed_opal 
  REAL(wp), PUBLIC :: sinkspeed_calc
  REAL(wp), PUBLIC :: sinkspeed_poc
  REAL(wp), PUBLIC :: sinkspeed_martin_ez
  REAL(wp), PUBLIC :: mc_fac
  REAL(wp), PUBLIC :: mc_depth
  REAL(wp), PUBLIC, ALLOCATABLE :: dzs(:)
  REAL(wp), PUBLIC :: dzsed(100),inpw(100)
  REAL(wp), PUBLIC, ALLOCATABLE :: porwat(:)
  REAL(wp), PUBLIC :: deltacalc
  REAL(wp), PUBLIC :: deltaorg
  REAL(wp), PUBLIC :: deltasil
 
  INTEGER, PUBLIC  :: io_stdo_bgc        !<  io unit for HAMOCC LOG file
  INTEGER, PUBLIC  :: ks,ksp

  LOGICAL, PUBLIC :: l_cyadyn         = .TRUE.   !  prognostic cyanobacteria
  LOGICAL, PUBLIC :: l_cpl_co2        = .FALSE.   !  co2 coupling to atm
  LOGICAL, PUBLIC :: l_diffat         = .FALSE.   !  diffusive atm
  LOGICAL, PUBLIC :: l_bgc_check      = .FALSE.   ! MASS check at every time step?
  LOGICAL, PUBLIC :: l_up_sedshi      = .FALSE.   ! Upward sediment shifting
  LOGICAL, PUBLIC :: l_implsed        = .FALSE.   ! Implicit sediment formulation
  LOGICAL, PUBLIC :: l_dynamic_pi      = .FALSE.    ! Depth dependent pi_alpha 
  LOGICAL, PUBLIC :: l_PDM_settling   = .FALSE.   ! PDM scheme for particle settling

  REAL(wp), PUBLIC :: denit_sed, disso_po

  !LOGICAL, PUBLIC :: l_avflux         = .TRUE.   ! flux redistribution
  
  REAL(wp), PUBLIC :: atm_co2, atm_o2, atm_n2
  INTEGER         :: iunit
 

  NAMELIST /hamocc_nml/ &
    &  i_settling, &
    &  ks, dzsed,inpw,&
    &  l_cyadyn, &
    &  sinkspeed_opal, &
    &  sinkspeed_calc, &
    &  sinkspeed_poc, &
    &  sinkspeed_martin_ez, &
    &  atm_co2, &
    &  atm_o2, &
    &  atm_n2, &
    &  mc_fac, &
    &  mc_depth, &
    &  deltacalc, &
    &  deltaorg, &
    &  deltasil, &
    &  l_cpl_co2, &
    &  l_bgc_check, &
    &  l_up_sedshi, &
    &  l_implsed,  &
    &  l_dynamic_pi, &
    &  l_PDM_settling

CONTAINS
  !>
  !!
  SUBROUTINE read_hamocc_namelist(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_hamocc_nml:read_hamocc_namelist'

    !------------------------------------------------------------------
    ! Set default values
    !------------------------------------------------------------------
    i_settling        = 0             ! constant sinking
   
    isac = 1       ! no sediment acceleration
    l_cyadyn = .TRUE.

    sinkspeed_opal =30._wp             ! m/d
    sinkspeed_calc = 30._wp            ! m/d
    sinkspeed_poc = 5._wp            ! m/d
    sinkspeed_martin_ez = 3.5_wp       ! m/d 
  
  ! sediment nml parameters
    denit_sed = 0.01_wp  
    disso_po = 0.01_wp


  ! Martin curve sinking
    mc_fac = 2.0_wp       !0.858_wp default value from Martin ea 1987 
    mc_depth = 100._wp

  ! Weathering fluxes
    deltacalc = 0._wp
    deltaorg = 0._wp
    deltasil = 0._wp

  ! Atmospheri concencentrations
   IF (l_diffat) THEN
       ! all concentrations will be calculated in carchm
    ELSE
       atm_co2 = 278._wp
       atm_o2  = 196800._wp
       atm_n2  = 802000._wp
    ENDIF

   ks=12
   dzsed(:) =-1._wp
   dzsed(1) = 0.001_wp
   dzsed(2) = 0.003_wp
   dzsed(3) = 0.005_wp
   dzsed(4) = 0.007_wp
   dzsed(5) = 0.009_wp
   dzsed(6) = 0.011_wp
   dzsed(7) = 0.013_wp
   dzsed(8) = 0.015_wp
   dzsed(9) = 0.017_wp
   dzsed(10) = 0.019_wp
   dzsed(11) = 0.021_wp
   dzsed(12) = 0.023_wp
   dzsed(13) = 0.025_wp

   inpw(:)=1._wp
   inpw(1) = 0.85_wp
   inpw(2) = 0.83_wp
   inpw(3) = 0.8_wp
   inpw(4) = 0.79_wp
   inpw(5) = 0.77_wp
   inpw(6) = 0.75_wp
   inpw(7) = 0.73_wp
   inpw(8) = 0.7_wp
   inpw(9) = 0.68_wp
   inpw(10) = 0.66_wp
   inpw(11) = 0.64_wp
   inpw(12) = 0.62_wp




    !------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('hamocc_nml')
      READ(funit,NML=hamocc_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------
    ! Read user's (new) specifications (done by all MPI processes)
    !------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('hamocc_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, hamocc_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (positioned)
      READ (nnml, hamocc_nml)                                        ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, hamocc_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------------
    ! Sanity Check
    !------------------------------------------------------------------
    IF (i_settling > 1 ) THEN
      CALL finish(TRIM(routine), 'Aggregation not yet implemented, i_settling must be 0 or 1')
    END IF

    ksp=ks+1
    ALLOCATE(porwat(ks))
    ALLOCATE(dzs(ksp))
    porwat(1:ks)=inpw(1:ks)
    dzs(1:ksp)=dzsed(1:ksp)
    dzsed(1:ksp)=dzs(1:ksp)*1000._wp ! for zaxis in output

    !------------------------------------------------------------------
    ! Store the namelist for restart
    !------------------------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=hamocc_nml)
      CALL store_and_close_namelist(funit, 'hamocc_nml')
    ENDIF

    !------------------------------------------------------------------
    ! Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=hamocc_nml)

   ! OPEN bgcout
   IF ( my_process_is_stdio() ) CALL open_bgcout 
    !------------------------------------------------------------------
    ! Fill the configuration state
    !------------------------------------------------------------------

  END SUBROUTINE read_hamocc_namelist

  !>
  !!  opens an ASCII file with HAMOCC debug messages

  SUBROUTINE open_bgcout 



    INTEGER :: istat

!-------------------------------------------------------------------------

    io_stdo_bgc = find_next_free_unit(10,20)

    OPEN (io_stdo_bgc, FILE='bgcout', IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish ('open_bgcout','Could not open bgcout')
    END IF

  END SUBROUTINE 



  
END MODULE mo_hamocc_nml
