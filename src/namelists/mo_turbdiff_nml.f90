!>
!! Namelist for turbulent diffusion (turbdiff)
!!
!! These subroutines are called by read_atmo_namelists and do the turbulent
!! diffusion setup (for turbdiff).
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-09-19)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_turbdiff_nml

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_config,       ONLY: isRestart
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_namelist,    ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_turbdiff_config,     ONLY: turbdiff_config 

  USE mo_data_turbdiff,       ONLY: &
    & itype_sher, imode_shshear, itype_wcld, &
    & imode_tran, imode_turb, icldm_tran, icldm_turb, &
    & ltkesso, ltkecon, lexpcor, ltmpcor, lprfcor, lnonloc, lcpfluc, lsflcnd, &
    & tur_len, pat_len, a_stab, tkhmin, tkmmin, c_diff, ltkeshs, &
    & rlam_heat, rlam_mom, rat_sea, tkesmot, frcsmot, impl_s, impl_t, &
    & a_hshr, imode_frcsmot, lfreeslip, alpha0, alpha0_max, alpha0_pert, tkhmin_strat, tkmmin_strat
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_turbdiff_namelist

  !----------------------------------!
  ! additional turbdiff_nml namelist variables  !
  !----------------------------------!

  LOGICAL :: lconst_z0   ! TRUE: horizontally homogeneous roughness length 
                         ! (for idealized testcases)

  REAL(wp) :: const_z0   ! horizontally homogeneous roughness length 
                         ! (for idealized testcases)

  LOGICAL :: ldiff_qi    ! turbulent diffusion of cloud ice QI
                         ! .TRUE.: ON

  NAMELIST/turbdiff_nml/ &
    & itype_sher, imode_shshear, itype_wcld, &
    & imode_tran, imode_turb, icldm_tran, icldm_turb, &
    & ltkesso, ltkecon, lexpcor, ltmpcor, lprfcor, lnonloc, lcpfluc, lsflcnd, &
    & tur_len, pat_len, a_stab, tkhmin, tkmmin, c_diff, ltkeshs, &
    & rlam_heat, rlam_mom, rat_sea, tkesmot, frcsmot, impl_s, impl_t, &
    & a_hshr, imode_frcsmot, alpha0, alpha0_max, tkhmin_strat, tkmmin_strat, &
!   additional namelist parameters:
    & lconst_z0, const_z0, lfreeslip, ldiff_qi

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for turbulent diffusion. 
  !!
  !! This subroutine 
  !! - reads the Namelist for turbulent diffusion
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-09-19)
  !!
  SUBROUTINE read_turbdiff_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_turbdiff_nml: read_turbdiff_nml'

    ! 0. default settings of internal turbdiff namelist variables
    !    is done by initialization in MODULE 'mo_data_turbdiff'

    !-----------------------
    ! 1. default settings of additional namelist variables:
    !-----------------------

    lconst_z0    =.FALSE. ! horizontally homogeneous roughness length 
                          ! (for idealized testcases)

    const_z0     = 0.001_wp ! horizontally homogeneous roughness length
                            ! (for idealized testcases)

    ldiff_qi     = .FALSE.  ! no turbulent diffusion of QI  

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('turbdiff_nml')
      READ(funit,NML=turbdiff_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('turbdiff_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults() 
      WRITE(iunit, turbdiff_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, turbdiff_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, turbdiff_nml)  ! write settings to temporary text file
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
      turbdiff_config(jg)%imode_tran   = imode_tran
      turbdiff_config(jg)%icldm_tran   = icldm_tran
      turbdiff_config(jg)%imode_turb   = imode_turb
      turbdiff_config(jg)%icldm_turb   = icldm_turb
      turbdiff_config(jg)%itype_sher   = itype_sher
      turbdiff_config(jg)%imode_shshear= imode_shshear
      turbdiff_config(jg)%imode_frcsmot= imode_frcsmot
      turbdiff_config(jg)%ltkesso      = ltkesso
      turbdiff_config(jg)%ltkeshs      = ltkeshs
      turbdiff_config(jg)%ltkecon      = ltkecon
      turbdiff_config(jg)%lexpcor      = lexpcor
      turbdiff_config(jg)%ltmpcor      = ltmpcor
      turbdiff_config(jg)%lprfcor      = lprfcor
      turbdiff_config(jg)%lnonloc      = lnonloc
      turbdiff_config(jg)%lfreeslip    = lfreeslip
      turbdiff_config(jg)%lcpfluc      = lcpfluc
      turbdiff_config(jg)%lsflcnd      = lsflcnd
      turbdiff_config(jg)%itype_wcld   = itype_wcld
      turbdiff_config(jg)%tur_len      = tur_len
      turbdiff_config(jg)%pat_len      = pat_len
      turbdiff_config(jg)%a_stab       = a_stab
      turbdiff_config(jg)%alpha0       = alpha0
      turbdiff_config(jg)%alpha0_max   = alpha0_max
      turbdiff_config(jg)%alpha0_pert  = alpha0_pert
      turbdiff_config(jg)%tkhmin       = tkhmin
      turbdiff_config(jg)%tkmmin       = tkmmin
      turbdiff_config(jg)%tkhmin_strat = tkhmin_strat
      turbdiff_config(jg)%tkmmin_strat = tkmmin_strat
      turbdiff_config(jg)%c_diff       = c_diff
      turbdiff_config(jg)%rlam_heat    = rlam_heat
      turbdiff_config(jg)%rlam_mom     = rlam_mom
      turbdiff_config(jg)%rat_sea      = rat_sea
      turbdiff_config(jg)%tkesmot      = tkesmot
      turbdiff_config(jg)%frcsmot      = frcsmot
      turbdiff_config(jg)%impl_s       = impl_s
      turbdiff_config(jg)%impl_t       = impl_t
      turbdiff_config(jg)%a_hshr       = a_hshr

      turbdiff_config(jg)%lconst_z0    = lconst_z0
      turbdiff_config(jg)%const_z0     = const_z0
      turbdiff_config(jg)%ldiff_qi     = ldiff_qi
    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=turbdiff_nml)                    
      CALL store_and_close_namelist(funit, 'turbdiff_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=turbdiff_nml)


  END SUBROUTINE read_turbdiff_namelist

END MODULE mo_turbdiff_nml
